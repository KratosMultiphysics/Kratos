"""Process deploying DoMINO surrogates on the CAE datapipe exports.

DoMINO (physicsnemo.models.domino) consumes the preprocessed sample dict of
DoMINODataPipe - geometry encodings, SDF grids, kNN surface neighborhoods -
not raw tensors, so deployment reuses the exact training-time pipeline: per
execution, the current model-part state is exported as a single-case .npz
into a scratch directory through the CaeDatasetExportProcess machinery, a
DoMINODataPipe (sampling=False: node/cell order preserved) preprocesses it,
and the model's per-node volume / per-triangle surface predictions are
written back onto the Kratos entities:

- volume outputs land on the volume model part's NODES (volume mesh centers
  are exported in tensor-adaptor node order);
- surface outputs are predicted per tessellated TRIANGLE and collapsed to
  their parent conditions/elements by mean (the mesh-bridge provenance maps
  each triangle to the entity it came from), so surface output fields must
  live on "condition"/"element" data locations.

Pretrained checkpoints (e.g. DoMINO-Automotive-Aero) are external - this
process runs any .mdlus DoMINO through model_registry, tiny randomly
initialized configs included.

torch/physicsnemo/warp are imported lazily on first execution.
"""

import tempfile
from pathlib import Path

import numpy

import KratosMultiphysics as Kratos
from KratosMultiphysics.PhysicsNeMoApplication.deployment import model_registry
from KratosMultiphysics.PhysicsNeMoApplication.bridges import torch_bridge
from KratosMultiphysics.PhysicsNeMoApplication.processes.export.cae_dataset_export_process import CaeDatasetExportProcess
from KratosMultiphysics.PhysicsNeMoApplication.bridges.mesh_bridge import domain_mesh_builder
from KratosMultiphysics.PhysicsNeMoApplication.processes.inference.inference_process import WriteOutputFields
from KratosMultiphysics.PhysicsNeMoApplication.utilities.tensor_adaptor_dataset_utils import GetTensorAdaptor
from KratosMultiphysics.PhysicsNeMoApplication.utilities.nvtx_utils import NvtxRange

_EXECUTION_POINTS = ("initialize_solution_step", "finalize_solution_step")
_NORMALIZATIONS = ("none", "min_max_scaling", "mean_std_scaling")
_SURFACE_LOCATIONS = ("condition", "element")


def Factory(settings: Kratos.Parameters, model: Kratos.Model) -> "DominoInferenceProcess":
    if not isinstance(settings, Kratos.Parameters):
        raise TypeError("Expected input shall be a Parameters object, encapsulating a json string")
    return DominoInferenceProcess(model, settings["Parameters"])


def _ParametersToKwargs(settings: Kratos.Parameters) -> dict:
    """Free-form Parameters block -> python kwargs (numbers, bools, strings,
    flat arrays)."""
    kwargs = {}
    for key in settings.keys():
        value = settings[key]
        if value.IsBool():
            kwargs[key] = value.GetBool()
        elif value.IsInt():
            kwargs[key] = value.GetInt()
        elif value.IsDouble():
            kwargs[key] = value.GetDouble()
        elif value.IsString():
            kwargs[key] = value.GetString()
        elif value.IsVector() or value.IsArray():
            entries = [value[i].GetDouble() if value[i].IsDouble() else value[i].GetInt()
                       for i in range(value.size())]
            kwargs[key] = entries
        else:
            raise ValueError(f"Unsupported override type for \"{key}\".")
    return kwargs


class DominoInferenceProcess(Kratos.Process):
    """Runs DoMINO inference each output_interval steps."""

    def __init__(self, model: Kratos.Model, settings: Kratos.Parameters) -> None:
        super().__init__()

        # free-form blocks leave the settings before the schema validation
        overrides = Kratos.Parameters("{}")
        if settings.Has("datapipe_overrides"):
            overrides = settings["datapipe_overrides"].Clone()
            settings.RemoveValue("datapipe_overrides")
        self.datapipe_overrides = _ParametersToKwargs(overrides)
        global_params = Kratos.Parameters("{}")
        if settings.Has("global_params"):
            global_params = settings["global_params"].Clone()
            settings.RemoveValue("global_params")
        global_params_reference = Kratos.Parameters("{}")
        if settings.Has("global_params_reference"):
            global_params_reference = settings["global_params_reference"].Clone()
            settings.RemoveValue("global_params_reference")

        default_settings = Kratos.Parameters("""{
            "volume_model_part_name"   : "",
            "surface_model_part_name"  : "PLEASE_SPECIFY_SURFACE_MODEL_PART_NAME",
            "surface_source_container" : "Conditions",
            "model_settings"           : {},
            "model_type"               : "surface",
            "bounding_box"             : [],
            "bounding_box_surface"     : [],
            "global_params_order"      : [],
            "tessellation_mode"        : "smallest_id_diagonal",
            "higher_order_mode"        : "reduce",
            "curved_refinement_levels" : 2,
            "output_fields_surface"    : [],
            "output_fields_volume"     : [],
            "scaling_factors_file"     : "",
            "normalization"            : "none",
            "redimensionalize"         : false,
            "scratch_directory"        : "",
            "execution_point"          : "finalize_solution_step",
            "output_interval"          : 1
        }""")
        settings.ValidateAndAssignDefaults(default_settings)
        field_defaults = Kratos.Parameters("""{
            "variable_name" : "PLEASE_SPECIFY_VARIABLE_NAME",
            "data_location" : "condition"
        }""")
        for key in ("output_fields_surface", "output_fields_volume"):
            for i in range(settings[key].size()):
                settings[key][i].ValidateAndAssignDefaults(field_defaults)

        self.model_type = settings["model_type"].GetString()
        if self.model_type not in ("surface", "volume", "combined"):
            raise ValueError(
                f"Unsupported model type \"{self.model_type}\". "
                "Use \"surface\", \"volume\" or \"combined\".")

        box = settings["bounding_box"].GetVector()
        if len(box) != 6:
            raise ValueError(
                "\"bounding_box\" [x0,y0,z0,x1,y1,z1] is required - DoMINO's preprocessing "
                "needs the volume bounding box even for surface-only models.")
        self.bounding_box = (tuple(box[:3]), tuple(box[3:]))
        surface_box = settings["bounding_box_surface"].GetVector()
        self.bounding_box_surface_was_given = len(surface_box) == 6
        if len(surface_box) == 6:
            self.bounding_box_surface = (tuple(surface_box[:3]), tuple(surface_box[3:]))
        elif len(surface_box) == 0:
            self.bounding_box_surface = None
        else:
            raise ValueError("\"bounding_box_surface\" must be empty or [x0,y0,z0,x1,y1,z1].")

        self.surface_model_part = model[settings["surface_model_part_name"].GetString()]
        volume_name = settings["volume_model_part_name"].GetString()
        self.volume_model_part = model[volume_name] if volume_name else None
        self.surface_source_container = settings["surface_source_container"].GetString()
        self.tessellation_mode = settings["tessellation_mode"].GetString()
        self.higher_order_mode = settings["higher_order_mode"].GetString()
        self.curved_refinement_levels = settings["curved_refinement_levels"].GetInt()
        self.model_settings = settings["model_settings"].Clone()

        # De-normalization. A pretrained DoMINO emits DIMENSIONLESS, min-max
        # normalized values; upstream's own wrapper unnormalizes against
        # scaling_factors.pkl and then multiplies by U^2 * rho. Writing the
        # raw output onto Kratos entities is wrong by ~3 orders of magnitude,
        # so this is opt-in but loudly documented.
        self.scaling_factors_file = settings["scaling_factors_file"].GetString()
        self.normalization = settings["normalization"].GetString()
        if self.normalization not in _NORMALIZATIONS:
            raise ValueError(
                f"Unsupported normalization \"{self.normalization}\". "
                f"Supported: {', '.join(_NORMALIZATIONS)}.")
        if self.normalization != "none" and not self.scaling_factors_file:
            raise ValueError(
                f"normalization \"{self.normalization}\" needs \"scaling_factors_file\" "
                "(the scaling_factors.pkl shipped beside the checkpoint).")
        self.redimensionalize = settings["redimensionalize"].GetBool()
        self._scaling = None
        self._global_params_values = _ParametersToKwargs(global_params)

        def read_specs(key):
            return [(settings[key][i]["variable_name"].GetString(),
                     settings[key][i]["data_location"].GetString())
                    for i in range(settings[key].size())]

        self.surface_output_specs = read_specs("output_fields_surface")
        self.volume_output_specs = read_specs("output_fields_volume")
        for variable_name, data_location in self.surface_output_specs:
            if data_location not in _SURFACE_LOCATIONS:
                raise ValueError(
                    f"Surface output \"{variable_name}\" has location \"{data_location}\"; "
                    "DoMINO predicts per tessellated triangle, collapsed onto the parent "
                    f"entities, so only {_SURFACE_LOCATIONS} data locations are supported.")
        if self.model_type in ("surface", "combined") and not self.surface_output_specs:
            raise ValueError("\"output_fields_surface\" must be configured for surface models.")
        if self.model_type in ("volume", "combined"):
            if self.volume_model_part is None:
                raise ValueError("\"volume_model_part_name\" must be configured for volume models.")
            if not self.volume_output_specs:
                raise ValueError("\"output_fields_volume\" must be configured for volume models.")

        scratch = settings["scratch_directory"].GetString()
        self.scratch_directory = Path(scratch) if scratch else Path(
            tempfile.mkdtemp(prefix="domino_inference_"))

        # the exporter that writes the single-case .npz the datapipe reads
        exporter_settings = Kratos.Parameters("{}")
        exporter_settings.AddString(
            "model_part_name",
            volume_name if volume_name else settings["surface_model_part_name"].GetString())
        exporter_settings.AddString(
            "surface_model_part_name", settings["surface_model_part_name"].GetString())
        exporter_settings.AddString("surface_source_container", self.surface_source_container)
        exporter_settings.AddString("tessellation_mode", self.tessellation_mode)
        exporter_settings.AddString("higher_order_mode", self.higher_order_mode)
        exporter_settings.AddInt("curved_refinement_levels", self.curved_refinement_levels)
        exporter_settings.AddString("output_path", str(self.scratch_directory))
        exporter_settings.AddString("file_prefix", "inference_case")
        exporter_settings.AddInt("case_id", 0)
        exporter_settings.AddValue("global_params", global_params)
        exporter_settings.AddValue("global_params_reference", global_params_reference)
        exporter_settings.AddValue("global_params_order", settings["global_params_order"])
        wrapped = Kratos.Parameters("{}")
        wrapped.AddValue("Parameters", exporter_settings)
        from KratosMultiphysics.PhysicsNeMoApplication.processes.export import cae_dataset_export_process
        self._exporter = cae_dataset_export_process.Factory(wrapped, model)

        self.execution_point = settings["execution_point"].GetString()
        if self.execution_point not in _EXECUTION_POINTS:
            raise ValueError(
                f"Unsupported execution point \"{self.execution_point}\". "
                f"Supported: {', '.join(_EXECUTION_POINTS)}.")
        self.output_interval = settings["output_interval"].GetInt()
        if self.output_interval < 1:
            raise ValueError(f"\"output_interval\" must be >= 1 [ output_interval = {self.output_interval} ].")

        self._model = None
        self._device = None
        # provenance is purely topological, so it is built once and reused;
        # the entity count it was built for is kept to detect a remesh
        self._surface_rows = None
        self._surface_n_entities = None

    def ExecuteInitialize(self) -> None:
        self.scratch_directory.mkdir(parents=True, exist_ok=True)
        self._exporter.ExecuteInitialize()

    def ExecuteInitializeSolutionStep(self) -> None:
        if self.execution_point == "initialize_solution_step":
            self._ExecuteIfDue()

    def ExecuteFinalizeSolutionStep(self) -> None:
        if self.execution_point == "finalize_solution_step":
            self._ExecuteIfDue()

    def _ExecuteIfDue(self) -> None:
        if self.surface_model_part.ProcessInfo[Kratos.STEP] % self.output_interval == 0:
            self.RunInference()

    def _PreprocessCase(self):
        """Exports the current state and preprocesses it through DoMINODataPipe."""
        from KratosMultiphysics.PhysicsNeMoApplication.training.torch_dataset import (
            CreateDoMINODataPipe, DOMINO_SURFACE_KEYS, DOMINO_VOLUME_KEYS)

        with NvtxRange("PhysicsNeMo::CaeExport"):
            self._exporter.Export()

        # inference cases carry no label fields - trim them from the reader keys
        keys = {
            "surface": DOMINO_SURFACE_KEYS,
            "volume": DOMINO_VOLUME_KEYS,
            "combined": tuple(dict.fromkeys(DOMINO_SURFACE_KEYS + DOMINO_VOLUME_KEYS)),
        }[self.model_type]
        keys = tuple(k for k in keys if k not in ("surface_fields", "volume_fields"))

        pipe = CreateDoMINODataPipe(
            self.scratch_directory, self.model_type, self.bounding_box,
            bounding_box_surface=self.bounding_box_surface,
            keys_to_read=keys, phase="test",
            device=str(self._device) if self._device is not None else None,
            **self.datapipe_overrides)
        with NvtxRange("PhysicsNeMo::DoMINOPreprocessing"):
            return pipe[0]

    def RunInference(self) -> None:
        torch = torch_bridge._TryImportTorch()
        if self._model is None:
            self._model, self._device = model_registry.LoadModelWithCardCheck(
                self.model_settings, [],
                self.surface_output_specs + self.volume_output_specs, type(self).__name__)

        self._CheckCheckpointAgreesWithSettings()
        sample = self._PreprocessCase()
        parameter = next(self._model.parameters(), None)
        dtype = parameter.dtype if parameter is not None else torch.float32
        batch = {}
        for key, value in sample.items():
            # DoMINODataPipe samples come pre-batched (leading dim 1)
            tensor = torch.as_tensor(numpy.asarray(value))
            if tensor.is_floating_point():
                tensor = tensor.to(dtype)
            batch[key] = tensor.to(self._device)

        with torch.no_grad(), NvtxRange("PhysicsNeMo::Forward"):
            output_vol, output_surf = self._model(batch)

        with NvtxRange("PhysicsNeMo::WriteOutputs"):
            if output_surf is not None and self.surface_output_specs:
                surface = self._Denormalize(
                    output_surf[0].cpu().double().numpy(), "surface_fields")
                self._WriteSurfaceOutputs(surface)
            if output_vol is not None and self.volume_output_specs:
                volume = self._Denormalize(
                    output_vol[0].cpu().double().numpy(), "volume_fields")
                prediction = torch.as_tensor(volume, dtype=torch.float64)
                WriteOutputFields(self.volume_model_part, self.volume_output_specs,
                                  prediction, prediction.shape[0])

    def _CheckCheckpointAgreesWithSettings(self) -> None:
        """Catches the three ways a pretrained checkpoint is silently misfed.

        Run once, after the model is loaded. Two of these are warnings rather
        than errors because a deliberately-trained-here model may legitimately
        differ; the grid resolution is an error because the alternative is an
        opaque reshape failure deep inside the geometry encoder.
        """
        if getattr(self, "_checked_checkpoint", False):
            return
        self._checked_checkpoint = True

        expected = getattr(self._model, "grid_resolution", None)
        if expected is not None:
            configured = self.datapipe_overrides.get("grid_resolution")
            if configured is not None and list(configured) != list(expected):
                raise ValueError(
                    f"datapipe_overrides.grid_resolution {list(configured)} does not match "
                    f"the checkpoint's own interp_res {list(expected)}. DoMINO would fail "
                    "with an opaque reshape error inside its geometry encoder; set them "
                    "equal.")
            if configured is None:
                Kratos.Logger.PrintWarning(
                    type(self).__name__,
                    f"the checkpoint was trained at grid resolution {list(expected)} but "
                    "\"grid_resolution\" is not set in datapipe_overrides, so the datapipe "
                    "default applies. Set it to the checkpoint's value.")

        if "normalize_coordinates" not in self.datapipe_overrides:
            Kratos.Logger.PrintWarning(
                type(self).__name__,
                "\"normalize_coordinates\" is not set in datapipe_overrides and the "
                "datapipe default is false, while pretrained DoMINO checkpoints are "
                "trained with it true. It is what emits surface_min_max, and getting it "
                "wrong produces no error - only wrong numbers.")

        if not self.bounding_box_surface_was_given:
            Kratos.Logger.PrintWarning(
                type(self).__name__,
                "\"bounding_box_surface\" was not given, so it defaults to \"bounding_box\". "
                "Pretrained checkpoints generally declare a tighter surface box than the "
                "volume one; check the checkpoint's config.yaml.")

    def _GetScalingFactors(self):
        """Loads scaling_factors.pkl through physicsnemo-cfd's restricted
        unpickler (a plain pickle.load cannot read it)."""
        if self._scaling is None:
            from KratosMultiphysics.PhysicsNeMoApplication.bridges import cfd_bridge
            _, ScalingUnpickler = cfd_bridge._TryImportDominoScaling()
            with open(self.scaling_factors_file, "rb") as source:
                self._scaling = ScalingUnpickler(source).load()
        return self._scaling

    def _Denormalize(self, prediction, field_key: str):
        """Turns a raw DoMINO prediction into physical values.

        A pretrained checkpoint predicts dimensionless, normalized fields.
        Upstream applies two steps, and BOTH matter: the inverse of the
        training normalization (from scaling_factors.pkl), then a
        redimensionalization by stream_velocity^2 * air_density. Skipping
        them leaves the values wrong by roughly three orders of magnitude
        and shifted, with nothing to indicate it.

        With the default "none"/false this is the identity, so existing
        configurations - which were written against raw output - are
        unaffected.
        """
        if self.normalization != "none":
            scaling = self._GetScalingFactors()
            if self.normalization == "min_max_scaling":
                maximum = numpy.asarray(scaling.max_val[field_key], dtype=numpy.float64)
                minimum = numpy.asarray(scaling.min_val[field_key], dtype=numpy.float64)
                if maximum.shape[-1] != prediction.shape[-1]:
                    raise ValueError(
                        f"scaling factors for \"{field_key}\" have "
                        f"{maximum.shape[-1]} channels but the model returned "
                        f"{prediction.shape[-1]}; the scaling_factors.pkl does not "
                        "belong to this checkpoint.")
                # upstream's unnormalize: x * (max - min) / 2 + (max + min) / 2
                prediction = (prediction * (maximum - minimum) / 2.0
                              + (maximum + minimum) / 2.0)
            else:
                # mean_std_scaling: these attributes are "mean"/"std", not the
                # "*_val" naming the min/max pair uses
                mean = numpy.asarray(scaling.mean[field_key], dtype=numpy.float64)
                std = numpy.asarray(scaling.std[field_key], dtype=numpy.float64)
                prediction = prediction * std + mean

        if self.redimensionalize:
            prediction = prediction * self._DynamicPressure()
        return prediction

    def _DynamicPressure(self) -> float:
        """U^2 * rho from the configured global parameters.

        The names follow the exporter's global_params; DoMINO's own config
        calls them inlet_velocity and air_density.
        """
        values = self._global_params_values
        velocity = None
        density = None
        for name, value in values.items():
            lowered = name.lower()
            if "velocity" in lowered:
                velocity = float(numpy.asarray(value).reshape(-1)[0])
            elif "density" in lowered:
                density = float(value)
        if velocity is None or density is None:
            raise ValueError(
                "redimensionalize needs a velocity and a density among \"global_params\"; "
                f"got {sorted(values)}. DoMINO's own names are inlet_velocity/air_density.")
        return velocity ** 2 * density

    def _SurfaceRows(self):
        """Tessellated-cell index -> parent entity row, built once.

        `BuildProvenance` is the most expensive per-entity path in
        `benchmarks/benchmark_bridges.py` (4 us/entity on tets, ~39 on
        hexes), and it maps cells to parent entities: topology only. Nothing
        about it changes between solution steps while the mesh does not, so
        rebuilding it every step - as this process used to - was 232 ms per
        step of repeated work on a 28k-triangle surface.

        The application ships `AdaptiveRemeshProcess`, so the mesh genuinely
        can change mid-run; the entity count is kept and the map rebuilt when
        it moves. A remesh that preserved the exact entity count would not be
        detected, which covers remeshing as it behaves rather than every
        conceivable mutation.
        """
        container = (self.surface_model_part.Conditions
                     if self.surface_source_container == "Conditions"
                     else self.surface_model_part.Elements)
        n_entities = len(container)
        if self._surface_rows is not None and self._surface_n_entities == n_entities:
            return self._surface_rows, n_entities

        provenance = domain_mesh_builder.BuildProvenance(
            self.surface_model_part, self.surface_source_container,
            self.tessellation_mode, self.higher_order_mode, self.curved_refinement_levels)
        parent_ids = provenance.cell_provenance[:, 0]
        row_of_entity = {entity.Id: row for row, entity in enumerate(container)}
        self._surface_rows = numpy.fromiter(
            (row_of_entity[int(pid)] for pid in parent_ids),
            dtype=numpy.int64, count=len(parent_ids))
        self._surface_n_entities = n_entities
        return self._surface_rows, n_entities

    def _WriteSurfaceOutputs(self, per_triangle) -> None:
        """Collapses (T, C_surf) per-triangle predictions onto the parent
        conditions/elements by mean and writes the entity fields."""
        rows, n_entities = self._SurfaceRows()
        sums = numpy.zeros((n_entities, per_triangle.shape[1]))
        counts = numpy.zeros(n_entities)
        numpy.add.at(sums, rows, per_triangle)
        numpy.add.at(counts, rows, 1.0)
        entity_values = sums / counts[:, None]

        offset = 0
        for variable_name, data_location in self.surface_output_specs:
            variable = Kratos.KratosGlobals.GetVariable(variable_name)
            tensor_adaptor = GetTensorAdaptor(self.surface_model_part, data_location, variable)
            width = int(numpy.prod(tensor_adaptor.data.shape[1:], dtype=int))
            chunk = entity_values[:, offset:offset + width]
            tensor_adaptor.data[:] = chunk.reshape(tensor_adaptor.data.shape)
            tensor_adaptor.StoreData()
            offset += width
        if offset != per_triangle.shape[1]:
            raise ValueError(
                f"DoMINO returned {per_triangle.shape[1]} surface channels but the configured "
                f"surface output fields consume {offset}.")
