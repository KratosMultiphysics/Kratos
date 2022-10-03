import KratosMultiphysics as Kratos
import KratosMultiphysics.kratos_utilities as kratos_utils
from  KratosMultiphysics.deprecation_management import DeprecationManager
import KratosMultiphysics.ShapeOptimizationApplication as KratosSOA

def Factory(settings, model):
    if not isinstance(settings, Kratos.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return UniversalFileFormatOutputProcess(model, settings["Parameters"])


class UniversalFileFormatOutputProcess(Kratos.OutputProcess):
    def __init__(self, model, settings: Kratos.Parameters):
        super().__init__()

        default_parameters = Kratos.Parameters("""{
            "model_part_name"                           : "PLEASE_SPECIFY_MODEL_PART_NAME",
            "output_file_name"                          : "output",
            "write_conditions"                          : false,
            "output_variables_list"                     : [],
            "extrapolate_gauss_point_variables_to_nodes": []
        }""")

        settings.ValidateAndAssignDefaults(default_parameters)

        model_part_name = settings["model_part_name"].GetString()
        self.model_part = model[model_part_name]

        self.interpolation_process = None
        if len(settings["extrapolate_gauss_point_variables_to_nodes"].GetStringArray()) > 0:
            default_interpolation_process_parameters = Kratos.Parameters("""(
                "echo_level"                 : 0,
                "area_average"               : true,
                "average_variable"           : "NODAL_AREA",
                "list_of_variables"          : [],
                "extrapolate_non_historical" : true
            }""")
            default_interpolation_process_parameters["list_of_variables"].SetStringArray(settings["extrapolate_gauss_point_variables_to_nodes"].GetStringArray())
            self.interpolation_process = Kratos.IntegrationValuesExtrapolationToNodesProcess(self.model_part, default_interpolation_process_parameters)

            current_variables_list = settings["output_variables_list"].GetStringArray()
            for var_name in settings["extrapolate_gauss_point_variables_to_nodes"].GetStringArray():
                current_variables_list.append(var_name)

            settings["output_variables_list"].SetStringArray(current_variables_list)

        self.unv_process = KratosSOA.UniversalFileIO(
            self.model_part,
            settings["output_file_name"].GetString(),
            settings["write_conditions"].GetBool(),
            settings["output_variables_list"]
        )

        self.counter = 0

    def PrintOutput(self):
        self.counter += 1
        self.unv_process.InitializeLogging()
        self.unv_process.LogNodalResults(self.counter)

    def IsOutputStep(self):
        return True

