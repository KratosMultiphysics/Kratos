"""Process putting an adjoint sensitivity field onto the model part.

The gradient counterpart of the export processes, and deliberately narrower
than one: it computes dJ/dX at every node and writes it into an ordinary
Kratos nodal variable (``SHAPE_SENSITIVITY`` by default). It writes no files
of its own, because it does not need to - once the field is a normal
variable, ``DatasetExportProcess``, ``StreamingDatasetExportProcess``,
``MeshExportProcess`` and the vtk/vtu outputs all pick it up unchanged, and

    CreateNpzDataset(directory, input_keys,
                     output_keys=[..., "SHAPE_SENSITIVITY__node_non_historical"])

already yields gradient-carrying training targets with no new dataset code.
That is what makes derivative-informed (Sobolev) training a configuration
rather than a new pipeline; see ``training.sobolev_training``.

Two sources of the same field:

``"shipped"``
    ``sensitivity_utils.ComputeShapeSensitivityField`` - the discretely exact
    dJ/dX from one element-local pass over the mesh, central-differenced.
    Needs nothing but the solved model part.
``"response_function"``
    Kratos's own adjoint stack through ``adjoint_bridge``: an application's
    response function runs its own primal and adjoint analyses and
    ``SensitivityBuilder`` produces the field. Needs that application, and
    the response's own case files.

Attach it *after* the solve of a step and *before* any export process, and
never after ``analysis.Finalize()`` - see the adjoint_bridge docstring.
"""

import numpy

import KratosMultiphysics as Kratos

from KratosMultiphysics.PhysicsNeMoApplication.bridges import adjoint_bridge
from KratosMultiphysics.PhysicsNeMoApplication.physics import differentiable_residual
from KratosMultiphysics.PhysicsNeMoApplication.physics import sensitivity_utils

_EXECUTION_POINTS = ("initialize_solution_step", "finalize_solution_step", "before_output_step")
_SOURCES = ("shipped", "response_function")


def Factory(settings: Kratos.Parameters, model: Kratos.Model) -> "AdjointSensitivityProcess":
    if not isinstance(settings, Kratos.Parameters):
        raise TypeError("Expected input shall be a Parameters object, encapsulating a json string")
    return AdjointSensitivityProcess(model, settings["Parameters"])


class AdjointSensitivityProcess(Kratos.Process):
    """Writes dJ/dX into a nodal variable at a step interval.

    Settings:
        {
            "model_part_name"    : "",
            "sensitivity_source" : "shipped",
            "objective"          : { },
            "dof_fields"         : [ { "variable_name" : "", "data_location" : "node_historical" } ],
            "fd_step"            : 1e-6,
            "design_sub_model_part_name" : "",
            "output_variable"    : "SHAPE_SENSITIVITY",
            "output_location"    : "node_non_historical",
            "response_settings"  : { },
            "execution_point"    : "finalize_solution_step",
            "output_interval"    : 1,
            "echo_level"         : 0
        }

    ``objective`` is an ``adjoint_bridge.MakeObjectiveWeights`` block and is
    what defines J for the ``"shipped"`` source; ``dof_fields`` names the
    solved unknowns (the DOF fields, historical by definition).
    ``response_settings`` is an ``adjoint_bridge.CreateResponseFunction``
    block and is what defines J for the ``"response_function"`` source - the
    two sources take the objective from different places because Kratos's
    response functions carry their own.
    """

    def __init__(self, model: Kratos.Model, settings: Kratos.Parameters) -> None:
        super().__init__()

        defaults = Kratos.Parameters("""{
            "model_part_name"            : "PLEASE_SPECIFY_MODEL_PART_NAME",
            "sensitivity_source"         : "shipped",
            "objective"                  : {},
            "dof_fields"                 : [],
            "fd_step"                    : 1e-6,
            "design_sub_model_part_name" : "",
            "output_variable"            : "SHAPE_SENSITIVITY",
            "output_location"            : "node_non_historical",
            "response_settings"          : {},
            "execution_point"            : "finalize_solution_step",
            "output_interval"            : 1,
            "echo_level"                 : 0
        }""")
        settings.ValidateAndAssignDefaults(defaults)

        self.model = model
        self.model_part = model[settings["model_part_name"].GetString()]

        self.source = settings["sensitivity_source"].GetString()
        if self.source not in _SOURCES:
            raise ValueError(
                f"Unsupported \"sensitivity_source\" \"{self.source}\". Use "
                f"{' or '.join(chr(34) + s + chr(34) for s in _SOURCES)}.")

        self.execution_point = settings["execution_point"].GetString()
        if self.execution_point not in _EXECUTION_POINTS:
            raise ValueError(
                f"Unsupported \"execution_point\" \"{self.execution_point}\". Supported: "
                f"{', '.join(_EXECUTION_POINTS)}.")

        self.output_interval = settings["output_interval"].GetInt()
        if self.output_interval < 1:
            raise ValueError(f"\"output_interval\" must be >= 1, got {self.output_interval}.")

        self.output_variable = Kratos.KratosGlobals.GetVariable(
            settings["output_variable"].GetString())
        self.output_location = settings["output_location"].GetString()
        self.fd_step = settings["fd_step"].GetDouble()
        self.echo_level = settings["echo_level"].GetInt()

        design_name = settings["design_sub_model_part_name"].GetString()
        self.design_node_ids = None
        if design_name:
            self.design_node_ids = [node.Id
                                    for node in self.model_part.GetSubModelPart(design_name).Nodes]

        self.objective_settings = settings["objective"].Clone()
        self.response_settings = settings["response_settings"].Clone()
        self.dof_field_specs = [
            (settings["dof_fields"][i]["variable_name"].GetString(),
             settings["dof_fields"][i]["data_location"].GetString())
            for i in range(settings["dof_fields"].size())]

        if self.source == "shipped" and not self.dof_field_specs:
            raise ValueError(
                "\"sensitivity_source\" : \"shipped\" needs \"dof_fields\" - the solved "
                "unknowns whose DOFs the adjoint system is posed on.")

        # built on first use: the DOF set does not exist until the solver has
        # initialized, which has not happened when processes are constructed
        self._assembler = None
        self._dof_map = None
        self._response = None

        #: The most recent (n_nodes, 3) field and its objective value.
        self.last_field = None
        self.last_value = None

    def ExecuteInitializeSolutionStep(self) -> None:
        if self.execution_point == "initialize_solution_step":
            self._MaybeExecute()

    def ExecuteFinalizeSolutionStep(self) -> None:
        if self.execution_point == "finalize_solution_step":
            self._MaybeExecute()

    def ExecuteBeforeOutputStep(self) -> None:
        if self.execution_point == "before_output_step":
            self._MaybeExecute()

    def _MaybeExecute(self) -> None:
        if self.model_part.ProcessInfo[Kratos.STEP] % self.output_interval == 0:
            self.Execute()

    def Execute(self) -> None:
        if self.source == "shipped":
            self.last_field = self._ComputeShipped()
        else:
            self.last_field = self._ComputeFromResponseFunction()

        if self.echo_level > 0:
            Kratos.Logger.PrintInfo(
                "AdjointSensitivityProcess",
                f"J = {self.last_value}, max |dJ/dX| = "
                f"{float(numpy.abs(self.last_field).max()) if self.last_field.size else 0.0}")

    def _EnsureAssembler(self):
        if self._assembler is None:
            self._assembler = differentiable_residual.TangentAssembler(self.model_part)
            self._dof_map = differentiable_residual.DofFieldMap(
                self._assembler, self.dof_field_specs)
        return self._assembler, self._dof_map

    def _ObjectiveGradientVector(self, dof_map):
        """dJ/du in equation order, from the objective's nodal weights."""
        variable_name = self.objective_settings["variable_name"].GetString()
        variable = Kratos.KratosGlobals.GetVariable(variable_name)
        weights = adjoint_bridge.MakeObjectiveWeights(self.objective_settings, self.model_part)

        selector = numpy.zeros((dof_map.n_nodes, dof_map.total_width), dtype=numpy.float64)
        if weights.shape[1] == 1:
            columns = [variable_name]
        else:
            columns = [f"{variable_name}_{axis}" for axis in "XYZ"]
        for column_index, name in enumerate(columns):
            if name not in dof_map.column_of:
                raise ValueError(
                    f"The objective reads \"{name}\", which is not among the \"dof_fields\" "
                    f"{[spec[0] for spec in self.dof_field_specs]} - an objective must be a "
                    "function of the solved unknowns.")
            selector[:, dof_map.column_of[name]] = weights[:, column_index]

        self.last_value = adjoint_bridge.EvaluateObjective(
            self.model_part, weights, variable, "node_historical")
        return dof_map.FieldsToDofVector(selector)

    def _ComputeShipped(self):
        assembler, dof_map = self._EnsureAssembler()
        dJ_du = self._ObjectiveGradientVector(dof_map)
        return sensitivity_utils.ComputeShapeSensitivityField(
            assembler, dof_map, dJ_du, fd_step=self.fd_step,
            design_node_ids=self.design_node_ids,
            output_variable=self.output_variable, output_location=self.output_location)

    def _ComputeFromResponseFunction(self):
        # The response owns its analyses, so it is built and Initialize()d
        # once and re-driven per call - rebuilding it would re-read the mdpa
        # and re-run the primal from scratch every step.
        first_call = self._response is None
        if first_call:
            self._response = adjoint_bridge.CreateResponseFunction(
                self.response_settings, self.model)

        fields = adjoint_bridge.EvaluateResponse(
            self._response, self.model_part,
            nodal_variables=(self.output_variable,), initialize=first_call)
        self.last_value = fields.value
        field = fields.nodal[self.output_variable.Name()]
        adjoint_bridge.WriteSensitivityField(
            self.model_part, field, self.output_variable, self.output_location)
        return field
