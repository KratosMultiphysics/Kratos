import KratosMultiphysics
import KratosMultiphysics.MeshingApplication as Meshing


def Factory(model, settings):
    return MmgModeler(model, settings)


class MmgModeler(KratosMultiphysics.Modeler):
    """Modeler that remeshes a ModelPart using the MMG library.

    The ``library`` parameter selects the MMG variant:
    - ``"MMG2D"``       — 2D planar remeshing
    - ``"MMG3D"``       — 3D volumetric remeshing (default)
    - ``"MMGSurface"``  — 3D surface (triangulated) remeshing

    All remaining parameters are forwarded to MmgProcess.

    Example
    -------
    .. code-block:: json

        {
            "modeler_name"    : "MmgModeler",
            "model_part_name" : "FluidModelPart",
            "library"         : "MMG3D",
            "echo_level"      : 0,
            "force_sizes" : {
                "force_min" : true,
                "minimal_size" : 0.01
            }
        }
    """

    def __init__(self, model, settings):
        super().__init__(model, settings)
        settings.ValidateAndAssignDefaults(self._GetDefaultSettings())
        self._model = model
        self._settings = settings

    def SetupModelPart(self):
        super().SetupModelPart()

        library = self._settings["library"].GetString()

        process_settings = KratosMultiphysics.Parameters(self._settings)
        process_settings.RemoveValue("library")

        model_part_name = process_settings["model_part_name"].GetString()
        process_settings.RemoveValue("model_part_name")
        model_part = self._model[model_part_name]

        if library == "MMG2D":
            Meshing.MmgProcess2D(model_part, process_settings).Execute()
        elif library == "MMGSurface":
            Meshing.MmgProcess3DSurfaces(model_part, process_settings).Execute()
        else:  # default: MMG3D
            Meshing.MmgProcess3D(model_part, process_settings).Execute()

    def _GetDefaultSettings(self):
        return KratosMultiphysics.Parameters("""{
            "model_part_name" : "",
            "library"         : "MMG3D"
        }""")

    def __str__(self):
        return "MmgModeler"
