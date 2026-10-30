"""CoSimulation solver wrapper deploying a PhysicsNeMo surrogate as a solver.

Lets a trained model participate as a *solver* in Kratos co-simulation
(CoSimulationApplication): each coupling iteration, SolveSolutionStep gathers
the configured input fields from the wrapper's own model part, runs one
forward pass and writes the output fields back - the coupled solver moves
the interface data to/from the other solvers through its data-transfer
operators and convergence accelerators, exactly as with any other wrapper.

Reference it in the co-simulation ProjectParameters by full module path
(resolved through the solver-wrapper factory's PYTHONPATH fallback):

```json
"surrogate": {
    "type" : "KratosMultiphysics.PhysicsNeMoApplication.deployment.cosim_surrogate_solver_wrapper",
    "solver_wrapper_settings" : {
        "mdpa_file"       : "surrogate_interface",
        "model_settings"  : { "checkpoint_file" : "surrogate.pt" },
        "input_fields"    : [ { "variable_name" : "POINT_LOAD",   "data_location" : "node_historical" } ],
        "output_fields"   : [ { "variable_name" : "DISPLACEMENT", "data_location" : "node_historical" } ]
    },
    "data" : {
        "load" : { "model_part_name" : "Surrogate", "variable_name" : "POINT_LOAD",   "dimension" : 3 },
        "disp" : { "model_part_name" : "Surrogate", "variable_name" : "DISPLACEMENT", "dimension" : 3 }
    }
}
```

Model interfaces: "flat" (default) is InferenceProcess's contract - the
input fields concatenate to one (n_entities, total_input_width) tensor and
the model returns (n_entities, total_output_width). The point-cloud
interfaces ("generic", "transolver", "flare", "geotransolver",
"figconvnet") dispatch through point_cloud_inference_process's
RunPointCloudForward with the node coordinates, so any point-cloud
surrogate deploys unchanged.

A "time_step" > 0 makes the wrapper own time (AdvanceInTime returns
current_time + time_step); the default 0.0 means another solver drives the
time (the coupled solver takes the max over all wrappers).

MPI. By default the wrapper is serial (rank-zero data communicator, like the
sdof wrapper) and behaves exactly as before. Setting "distributed" : true
hands the data communicator back to the base class, so the wrapper runs on
all ranks, or - when "mpi_settings" names a "num_processes"/
"data_communicator_name" - on the first N of them. The mesh is then read
through DistributedImportModelPartUtility (Metis), or, for meshes without
elements that Metis cannot partition, through the "partition_mdpa" : false
path which reads on every rank and keeps the nodes whose round-robin
PARTITION_INDEX matches. Inference runs on the *owned* nodes only (the
communicator's LocalMesh, which is the layout CouplingInterfaceData uses)
and ghosts are refreshed with SynchronizeVariable afterwards.

Only the "flat" interface is partition-safe by construction: it maps each
entity row independently, so per-rank inference reproduces the serial result
node for node. The point-cloud interfaces mix information across the whole
cloud (attention in transolver/flare/geotransolver, and the per-rank
coordinate normalization differs too), so a partitioned run does NOT
reproduce the serial answer; they are rejected under "distributed" unless
"assume_partition_safe" : true says the model really is pointwise.

Known upstream hazard: if *both* coupled wrappers use the rank-zero data
communicator in an MPI run, CoSimulation's "kratos_mapping" data-transfer
operator deadlocks - it picks the serial mapper on rank 0 (which sees two
plain model parts) and the MPI mapper on the other ranks (which see the
rank-zero dummy parts, and those report IsDistributed() == True). Keep at
least one side of such a coupling distributed.

This module is only ever imported by the CoSimulation factory, so importing
CoSimulationApplication eagerly here does not affect the application's
bare-import contract; torch stays lazy via model_registry.
"""

import KratosMultiphysics as Kratos
from KratosMultiphysics.CoSimulationApplication.base_classes.co_simulation_solver_wrapper import CoSimulationSolverWrapper
from KratosMultiphysics.CoSimulationApplication.utilities.data_communicator_utilities import GetRankZeroDataCommunicator

from KratosMultiphysics.PhysicsNeMoApplication.deployment import model_registry
from KratosMultiphysics.PhysicsNeMoApplication.bridges import torch_bridge
from KratosMultiphysics.PhysicsNeMoApplication.processes.inference.inference_process import (
    GatherInputFields, SynchronizeOutputFields, WriteOutputFields)
from KratosMultiphysics.PhysicsNeMoApplication.processes.inference.point_cloud_inference_process import (
    GatherPointCloudCoordinates, RunPointCloudForward, _MODEL_INTERFACES)
from KratosMultiphysics.PhysicsNeMoApplication.utilities.nvtx_utils import NvtxRange

_WRAPPER_INTERFACES = ("flat",) + _MODEL_INTERFACES

# Interfaces whose output for an entity depends only on that entity's own
# row, hence reproduce the serial result when the mesh is partitioned. The
# point-cloud interfaces mix across the cloud and are not on this list.
_PARTITION_SAFE_INTERFACES = ("flat",)


def Create(settings, model, solver_name):
    return CoSimSurrogateSolverWrapper(settings, model, solver_name)


class CoSimSurrogateSolverWrapper(CoSimulationSolverWrapper):
    """Deploys a trained PhysicsNeMo/TorchScript model as a co-simulation solver."""

    def __init__(self, settings, model, solver_name):
        super().__init__(settings, model, solver_name)

        wrapper_defaults = Kratos.Parameters("""{
            "mdpa_file"             : "PLEASE_SPECIFY_MDPA_FILE",
            "main_model_part_name"  : "Surrogate",
            "domain_size"           : 3,
            "buffer_size"           : 1,
            "time_step"             : 0.0,
            "model_settings"        : {},
            "model_interface"       : "flat",
            "normalize_coordinates" : true,
            "pass_geometry"         : true,
            "distributed"           : false,
            "partition_mdpa"        : true,
            "assume_partition_safe" : false,
            "input_fields"          : [
                {
                    "variable_name" : "PLEASE_SPECIFY_VARIABLE_NAME",
                    "data_location" : "node_historical"
                }
            ],
            "output_fields"         : [
                {
                    "variable_name" : "PLEASE_SPECIFY_VARIABLE_NAME",
                    "data_location" : "node_historical"
                }
            ]
        }""")
        wrapper_settings = self.settings["solver_wrapper_settings"]
        wrapper_settings.ValidateAndAssignDefaults(wrapper_defaults)
        for key in ("input_fields", "output_fields"):
            for i in range(wrapper_settings[key].size()):
                wrapper_settings[key][i].ValidateAndAssignDefaults(wrapper_defaults[key][0])

        self.model_settings = wrapper_settings["model_settings"].Clone()
        self.model_interface = wrapper_settings["model_interface"].GetString()
        self.normalize_coordinates = wrapper_settings["normalize_coordinates"].GetBool()
        self.pass_geometry = wrapper_settings["pass_geometry"].GetBool()
        self.time_step = wrapper_settings["time_step"].GetDouble()
        self.input_specs = self._ReadFieldSpecs(wrapper_settings["input_fields"])
        self.output_specs = self._ReadFieldSpecs(wrapper_settings["output_fields"])
        self.last_scalar_prediction = None  # figconvnet-style global output

        if self.model_interface not in _WRAPPER_INTERFACES:
            raise ValueError(
                f"Unsupported model interface \"{self.model_interface}\". "
                f"Supported: {', '.join(_WRAPPER_INTERFACES)}.")

        if self.distributed and self.model_interface not in _PARTITION_SAFE_INTERFACES \
                and not wrapper_settings["assume_partition_safe"].GetBool():
            raise ValueError(
                f"Model interface \"{self.model_interface}\" is not partition-safe: it mixes "
                "information across the whole point cloud (attention, and a per-rank "
                "coordinate normalization), so a distributed run would not reproduce the "
                "serial prediction. Partition-safe interfaces: "
                f"{', '.join(_PARTITION_SAFE_INTERFACES)}. Set "
                "\"assume_partition_safe\" : true if this model really is pointwise.")

        # Off-group ranks own no mesh at all: the coupled solver replaces the
        # wrapper with an UndefinedSolver, and building here would both waste
        # the read and risk joining collectives the group does not make.
        if not self.data_communicator.IsDefinedOnThisRank():
            self.model_part = None
            self._model = self._device = None
            return

        self.model_part = self.model.CreateModelPart(
            wrapper_settings["main_model_part_name"].GetString())
        self.model_part.ProcessInfo[Kratos.DOMAIN_SIZE] = wrapper_settings["domain_size"].GetInt()

        # Historical variables must be allocated BEFORE the mesh is read
        # (Kratos memory layout); collect them from the field specs and from
        # the coupling "data" blocks (components resolve to their source).
        for variable_name, data_location in self.input_specs + self.output_specs:
            if data_location == "node_historical":
                self._AddHistoricalVariable(variable_name)
        for data_name in self.settings["data"].keys():
            self._AddHistoricalVariable(
                self.settings["data"][data_name]["variable_name"].GetString())

        mdpa_file = wrapper_settings["mdpa_file"].GetString()
        if self.distributed:
            self.model_part.AddNodalSolutionStepVariable(Kratos.PARTITION_INDEX)
            self._ReadDistributedModelPart(
                mdpa_file, wrapper_settings["partition_mdpa"].GetBool())
        else:
            Kratos.ModelPartIO(mdpa_file).ReadModelPart(self.model_part)
        self.model_part.SetBufferSize(wrapper_settings["buffer_size"].GetInt())

        # Model loading is deferred to the first solve so that merely
        # constructing the wrapper never requires torch.
        self._model = None
        self._device = None

    @staticmethod
    def _ReadFieldSpecs(fields: Kratos.Parameters):
        return [
            (fields[i]["variable_name"].GetString(), fields[i]["data_location"].GetString())
            for i in range(fields.size())
        ]

    def _AddHistoricalVariable(self, variable_name: str) -> None:
        variable = Kratos.KratosGlobals.GetVariable(variable_name)
        if hasattr(variable, "GetSourceVariable"):  # component -> full variable
            variable = variable.GetSourceVariable()
        self.model_part.AddNodalSolutionStepVariable(variable)

    def _ReadDistributedModelPart(self, mdpa_file: str, partition_mdpa: bool) -> None:
        """Fills self.model_part across self.data_communicator's ranks.

        With partition_mdpa the mdpa is Metis-partitioned in memory. Metis
        needs a connectivity graph, so a mesh with no elements or conditions
        is rejected outright ("number of connected nodes = 0"); such meshes -
        the pure point clouds a surrogate is often deployed on - must use the
        partition-free path, which reads the whole mdpa on every rank and
        then assigns ownership round-robin. That path trades memory (every
        rank holds every node) for the ability to partition an edgeless
        cloud, and is exact for pointwise models either way.
        """
        import KratosMultiphysics.mpi as KratosMPI
        from KratosMultiphysics.mpi.distributed_import_model_part_utility import (
            DistributedImportModelPartUtility)

        communicator = self.data_communicator
        if partition_mdpa:
            import_settings = Kratos.Parameters("""{
                "echo_level" : 0,
                "model_import_settings" : {
                    "input_type"          : "mdpa",
                    "input_filename"      : "",
                    "partition_in_memory" : true
                }
            }""")
            import_settings["model_import_settings"]["input_filename"].SetString(mdpa_file)
            # Partition over exactly this wrapper's ranks, not the world -
            # the same alignment kratos_base_wrapper performs for the Kratos
            # solver wrappers. DataCommunicator exposes no name, so take it
            # from the mpi_settings that created the sub-communicator.
            mpi_settings = self.settings["mpi_settings"]
            if mpi_settings.Has("data_communicator_name"):
                import_settings["model_import_settings"].AddString(
                    "data_communicator_name",
                    mpi_settings["data_communicator_name"].GetString())
            utility = DistributedImportModelPartUtility(self.model_part, import_settings)
            utility.ImportModelPart()
            utility.CreateCommunicators()
            return

        Kratos.ModelPartIO(mdpa_file).ReadModelPart(self.model_part)
        size = communicator.Size()
        for offset, node in enumerate(self.model_part.Nodes):
            node.SetSolutionStepValue(Kratos.PARTITION_INDEX, offset % size)
        KratosMPI.ModelPartCommunicatorUtilities.SetMPICommunicator(self.model_part, communicator)
        KratosMPI.ParallelFillCommunicator(self.model_part, communicator).Execute()

    def _GetModel(self):
        if self._model is None:
            self._model, self._device = model_registry.LoadModelWithCardCheck(
                self.model_settings, self.input_specs, self.output_specs, type(self).__name__)
        return self._model

    def _GetNormalization(self):
        """The card's output de-normalization, or None - the same contract
        as InferenceProcess._GetNormalization. This wrapper writes through
        WriteOutputFields like every other deployment path and was
        documented as covered, but passed no normalization, so a card was
        silently ignored here."""
        if not hasattr(self, "_normalization"):
            self._normalization = model_registry.LoadOutputNormalization(self.model_settings)
        return self._normalization

    def _GetInputNormalization(self):
        """The card's input normalization, or None (the symmetric half)."""
        if not hasattr(self, "_input_normalization"):
            self._input_normalization = model_registry.LoadInputNormalization(self.model_settings)
        return self._input_normalization

    def AdvanceInTime(self, current_time):
        if self.time_step <= 0.0:
            return 0.0  # another solver owns the time
        new_time = current_time + self.time_step
        self.model_part.ProcessInfo[Kratos.STEP] += 1
        self.model_part.CloneTimeStep(new_time)
        return new_time

    def SolveSolutionStep(self):
        model = self._GetModel()
        torch = torch_bridge._TryImportTorch()

        # On a distributed model part, read and write the OWNED entities only
        # (the communicator's LocalMesh). model_part.Nodes would also carry
        # ghosts, which both wastes forward passes and puts the rows out of
        # step with CouplingInterfaceData, whose layout is ghost-free.
        local_only = self.model_part.IsDistributed()

        with NvtxRange("PhysicsNeMo::GatherInputs"):
            inputs, n_entities = GatherInputFields(
                self.model_part, self.input_specs, local_only=local_only)
            features = model_registry.ApplyInputNormalization(
                torch.cat(inputs, dim=-1), self._GetInputNormalization())  # (N, C_in)

        if self.model_interface == "flat":
            with torch.no_grad(), NvtxRange("PhysicsNeMo::Forward"):
                prediction = model(features.to(self._device)).cpu()
        else:
            coordinates = torch.from_numpy(GatherPointCloudCoordinates(
                self.model_part, self.normalize_coordinates, local_only=local_only))  # (N, 3)
            prediction, self.last_scalar_prediction = RunPointCloudForward(
                model, self._device, self.model_interface, features, coordinates,
                self.pass_geometry)

        with NvtxRange("PhysicsNeMo::WriteOutputs"):
            WriteOutputFields(self.model_part, self.output_specs, prediction, n_entities,
                              local_only=local_only, normalization=self._GetNormalization())

        if local_only:
            # owned values are authoritative; refresh every rank's ghosts
            SynchronizeOutputFields(self.model_part, self.output_specs)

    def Check(self):
        field_variables = {name for name, _ in self.input_specs + self.output_specs}
        for data in self.data_dict.values():
            variable = data.variable
            if hasattr(variable, "GetSourceVariable"):
                variable = variable.GetSourceVariable()
            if variable.Name() not in field_variables:
                Kratos.Logger.PrintWarning(
                    type(self).__name__,
                    f"Interface data \"{data.name}\" uses variable \"{variable.Name()}\" "
                    "which is neither an input nor an output field of the surrogate; "
                    "the model will not read or write it.")

    def _GetDataCommunicator(self):
        # Called from the base __init__, i.e. before solver_wrapper_settings
        # have been validated against the defaults - read the raw value.
        wrapper_settings = self.settings["solver_wrapper_settings"]
        self.distributed = (wrapper_settings.Has("distributed")
                            and wrapper_settings["distributed"].GetBool())

        if not self.distributed:
            # the default: this wrapper is serial (like the sdof wrapper) and
            # the surrogate runs on rank zero only
            return GetRankZeroDataCommunicator()

        if not Kratos.IsDistributedRun():
            # asking for a distributed surrogate in a non-MPI run degrades to
            # serial rather than failing, mirroring how the Kratos wrappers
            # treat a solver configured "OpenMP"
            self.distributed = False
            return Kratos.ParallelEnvironment.GetDefaultDataCommunicator()

        # honours "mpi_settings": the world comm, or the first N ranks under
        # a named sub-communicator
        return super()._GetDataCommunicator()
