"""Process calling a PhysicsNeMo NIM microservice from the solution loop.

The serving counterpart of DominoInferenceProcess: instead of loading a
local checkpoint (weights, GPU, torch), the model part's triangle surface
is packed into a binary STL and POSTed to a running NIM
(deployment.nim_client, documented contract of
nvcr.io/nim/nvidia/domino-automotive-aero), and the returned point-sampled
fields are mapped back onto the Kratos nodes by nearest neighbour - the
NIM samples its own point cloud, so its points are NOT the mesh nodes and
an exact provenance scatter does not exist on this path.

Surface outputs (e.g. pressure_surface, wall_shear_stress) are mapped onto
the nodes of the surface that was sent; volume outputs (e.g. pressure,
velocity) onto all nodes of the model part. Scalar outputs (drag_force,
lift_force, ...) are collected into ``last_scalars`` and logged.

The solver host needs numpy and scipy only. Like TritonInferenceProcess,
``SetClient`` is the test seam - the payload contract is pinned against a
stub transport, since running a real NIM needs an NGC API key and docker.
"""

import os

import numpy

import KratosMultiphysics as Kratos
from KratosMultiphysics.PhysicsNeMoApplication.deployment import nim_client

_EXECUTION_POINTS = ("initialize_solution_step", "finalize_solution_step")


def Factory(settings: Kratos.Parameters, model: Kratos.Model) -> "NimInferenceProcess":
    if not isinstance(settings, Kratos.Parameters):
        raise TypeError("Expected input shall be a Parameters object, encapsulating a json string")
    return NimInferenceProcess(model, settings["Parameters"])


class NimInferenceProcess(Kratos.Process):
    """Runs one NIM inference every output_interval steps and scatters the answer."""

    def __init__(self, model: Kratos.Model, settings: Kratos.Parameters) -> None:
        super().__init__()

        # free-form string->string block, validated by hand (a schema template
        # default would be inherited as a real entry - the cfd_metrics lesson)
        extra_form_fields = {}
        if settings.Has("extra_form_fields"):
            block = settings["extra_form_fields"]
            for key in block.keys():
                extra_form_fields[key] = block[key].GetString()
            settings.RemoveValue("extra_form_fields")

        default_settings = Kratos.Parameters("""{
            "model_part_name"       : "PLEASE_SPECIFY_MODEL_PART_NAME",
            "base_url"              : "http://localhost:8000",
            "endpoint"              : "infer",
            "api_key"               : "",
            "api_key_env"           : "",
            "timeout"               : 120.0,
            "stream_velocity"       : 30.0,
            "stencil_size"          : 1,
            "point_cloud_size"      : 500000,
            "source_container"      : "Conditions",
            "surface_output_fields" : [],
            "volume_output_fields"  : [],
            "execution_point"       : "finalize_solution_step",
            "output_interval"       : 1
        }""")
        field_defaults = Kratos.Parameters("""{
            "nim_key"       : "PLEASE_SPECIFY_NIM_KEY",
            "variable_name" : "PLEASE_SPECIFY_VARIABLE_NAME",
            "data_location" : "node_historical"
        }""")
        settings.ValidateAndAssignDefaults(default_settings)
        for key in ("surface_output_fields", "volume_output_fields"):
            for i in range(settings[key].size()):
                settings[key][i].ValidateAndAssignDefaults(field_defaults)

        self.model_part = model[settings["model_part_name"].GetString()]
        self.endpoint = settings["endpoint"].GetString()
        api_key = settings["api_key"].GetString()
        api_key_env = settings["api_key_env"].GetString()
        if not api_key and api_key_env:
            api_key = os.environ.get(api_key_env, "")
        self._client = nim_client.NimClient(
            settings["base_url"].GetString(),
            timeout=settings["timeout"].GetDouble(),
            api_key=api_key)

        self.form_fields = {
            "stream_velocity": str(settings["stream_velocity"].GetDouble()),
            "stencil_size": str(settings["stencil_size"].GetInt()),
            "point_cloud_size": str(settings["point_cloud_size"].GetInt()),
        }
        self.form_fields.update(extra_form_fields)

        self.source_container = settings["source_container"].GetString()
        if self.source_container not in ("Conditions", "Elements"):
            raise ValueError(
                f"Unsupported \"source_container\" \"{self.source_container}\". "
                "Use \"Conditions\" (a skin part) or \"Elements\" (a surface part).")
        self.surface_output_specs = self._ReadFieldSpecs(settings["surface_output_fields"])
        self.volume_output_specs = self._ReadFieldSpecs(settings["volume_output_fields"])

        self.execution_point = settings["execution_point"].GetString()
        if self.execution_point not in _EXECUTION_POINTS:
            raise ValueError(
                f"Unsupported execution point \"{self.execution_point}\". "
                f"Supported: {', '.join(_EXECUTION_POINTS)}.")
        self.output_interval = settings["output_interval"].GetInt()
        if self.output_interval < 1:
            raise ValueError(f"\"output_interval\" must be >= 1 [ output_interval = {self.output_interval} ].")

        self.last_scalars = {}

    @staticmethod
    def _ReadFieldSpecs(fields: Kratos.Parameters):
        return [
            (fields[i]["nim_key"].GetString(),
             fields[i]["variable_name"].GetString(),
             fields[i]["data_location"].GetString())
            for i in range(fields.size())
        ]

    def SetClient(self, client) -> None:
        """Replaces the NimClient - the test seam (a client with a stub
        transport pins the payload without a running NIM)."""
        self._client = client

    def GetClient(self):
        return self._client

    def ExecuteInitializeSolutionStep(self) -> None:
        if self.execution_point == "initialize_solution_step":
            self._ExecuteIfDue()

    def ExecuteFinalizeSolutionStep(self) -> None:
        if self.execution_point == "finalize_solution_step":
            self._ExecuteIfDue()

    def _ExecuteIfDue(self) -> None:
        if self.model_part.ProcessInfo[Kratos.STEP] % self.output_interval == 0:
            self.RunInference()

    def _GatherSurface(self):
        """(points, triangles, surface node row-index) of the triangle skin.

        Pure Kratos + numpy: rows follow the model part's node iteration
        order (ascending ids, the application-wide row contract) and the
        triangles index into those rows. Non-triangle geometries are
        rejected - tessellate a volume part's skin into triangle conditions
        first (mesh_bridge, or SkinDetectionProcess + a simplex skin).
        """
        node_ids = numpy.fromiter(
            (node.Id for node in self.model_part.Nodes), dtype=numpy.int64)
        points = numpy.array(
            [[node.X, node.Y, node.Z] for node in self.model_part.Nodes])
        container = (self.model_part.Conditions if self.source_container == "Conditions"
                     else self.model_part.Elements)
        connectivity = []
        for entity in container:
            geometry = entity.GetGeometry()
            if len(geometry) != 3:
                raise ValueError(
                    f"The NIM takes a TRIANGLE surface; {self.source_container[:-1]} "
                    f"{entity.Id} has {len(geometry)} nodes. Use a triangle skin "
                    "part (source_container \"Conditions\"), or tessellate first "
                    "through the mesh bridge.")
            connectivity.append([node.Id for node in geometry])
        if not connectivity:
            raise ValueError(
                f"\"{self.model_part.Name}\" has no {self.source_container} to send.")
        triangles = numpy.searchsorted(node_ids, numpy.array(connectivity, dtype=numpy.int64))
        surface_rows = numpy.unique(triangles)
        return points, triangles, surface_rows

    def _ScatterByNearest(self, coordinates, rows, specs, response) -> None:
        from scipy.spatial import cKDTree

        nodes = list(self.model_part.Nodes)
        query_points = numpy.array([[nodes[i].X, nodes[i].Y, nodes[i].Z] for i in rows])
        tree = cKDTree(numpy.asarray(coordinates, dtype=numpy.float64))
        _, nearest = tree.query(query_points)

        for nim_key, variable_name, data_location in specs:
            if nim_key not in response:
                raise ValueError(
                    f"The NIM response has no \"{nim_key}\" (it has: "
                    f"{', '.join(sorted(response))}). Check the endpoint - surface "
                    "keys need \"infer\" or \"infer/surface\", volume keys \"infer\" "
                    "or \"infer/volume\".")
            values = numpy.asarray(response[nim_key], dtype=numpy.float64)[nearest]
            variable = Kratos.KratosGlobals.GetVariable(variable_name)
            historical = data_location == "node_historical"
            for row, value in zip(rows, values):
                node = nodes[row]
                # Kratos.Vector silently misreads strided views - keep it contiguous
                payload = (Kratos.Vector(numpy.ascontiguousarray(value))
                           if numpy.ndim(value) else float(value))
                if historical:
                    node.SetSolutionStepValue(variable, payload)
                else:
                    node.SetValue(variable, payload)

    def RunInference(self) -> None:
        points, triangles, surface_rows = self._GatherSurface()
        response = self._client.Infer(
            nim_client.MakeBinaryStlBytes(points, triangles),
            self.form_fields, self.endpoint)

        for specs, coordinates_key, rows in (
                (self.surface_output_specs, "surface_coordinates", surface_rows),
                (self.volume_output_specs, "coordinates",
                 numpy.arange(self.model_part.NumberOfNodes()))):
            if not specs:
                continue
            if coordinates_key not in response:
                raise ValueError(
                    f"The NIM response has no \"{coordinates_key}\" to map from (it "
                    f"has: {', '.join(sorted(response))}). Surface outputs need the "
                    "\"infer\" or \"infer/surface\" endpoint, volume outputs "
                    "\"infer\" or \"infer/volume\".")
            self._ScatterByNearest(response[coordinates_key], rows, specs, response)

        self.last_scalars = {
            key: float(numpy.ravel(value)[0])
            for key, value in response.items()
            if numpy.asarray(value).size == 1
        }
        for key, value in sorted(self.last_scalars.items()):
            Kratos.Logger.PrintInfo("NimInferenceProcess", f"{key} = {value:.6g}")
