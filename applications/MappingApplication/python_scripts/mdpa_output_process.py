import KratosMultiphysics as KM
from KratosMultiphysics.kratos_utilities import GenerateVariableListFromInput

def Factory(settings, model):
    if not isinstance(settings, KM.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return MdpaOutputProcess(model, settings["Parameters"])


class MdpaOutputProcess(KM.OutputProcess):
    def __init__(self, model, settings):
        super().__init__()

        self.settings = settings

        default_settings = KM.Parameters('''{
            "model_part_name" : "UNSPECIFIED",
            "mdpa_file_name"  : "UNSPECIFIED",
            "io_flags" : []
        }''')

        settings.ValidateAndAssignDefaults(default_settings)

        self.mdpa_file_name = settings["mdpa_file_name"].GetString()

        self.model_part = model.GetModelPart(settings["model_part_name"].GetString())

    def IsOutputStep(self):
        return True

    def PrintOutput(self):
        KM.ModelPartIO(self.mdpa_file_name, GetIOFlags(self.settings["io_flags"])).WriteModelPart(self.model_part)

def GetIOFlags(settings):
    io_flags_dict = {
        "mesh_only"    : KM.IO.MESH_ONLY,
        "scientific_precision" : KM.IO.SCIENTIFIC_PRECISION
    }
    io_flags = KM.IO.WRITE|KM.IO.SKIP_TIMER
    for flag_name in settings.GetStringArray():
        io_flags |= io_flags_dict[flag_name]

    return io_flags
