# Importing the Kratos Library
import KratosMultiphysics as KM
import KratosMultiphysics.IgaApplication as IGA


def Factory(settings, model):
    if not isinstance(settings, KM.Parameters):
        raise Exception("Expected input shall be a Parameters object, encapsulating a json string")
    return IgaContactProcessGapSbmALM(model, settings["Parameters"])


class IgaContactProcessGapSbmALM(KM.Process):
    def __init__(self, model, settings):
        super().__init__()

        default_parameters = KM.Parameters(
            """{
                "echo_level" : 0,
                "analysis_model_part_name" : "ModelPart",
                "contact_sub_model_part_name" : "contact",
                "shape_function_derivatives_order" : 1,
                "numbered_considered_neighbours" : 1,
                "alm_parameters" : {
                    "penalty" : 1.0e4,
                    "scale_factor" : 1.0
                },
                "contact_parameters" : {
                    "slave_model_part" : {
                        "sub_model_part_name" : "",
                        "layer_name" : "",
                        "property_id" : 1
                    },
                    "master_model_part" : {
                        "sub_model_part_name" : "",
                        "layer_name" : "",
                        "property_id" : 1
                    }
                }
            }"""
        )

        self.settings = settings
        self.settings.ValidateAndAssignDefaults(default_parameters)

        self.model = model
        self.analysis_model_part = model[self.settings["analysis_model_part_name"].GetString()]
        self.contact_process = IGA.IgaContactProcessGapSbm(model, self.settings)

    def ExecuteInitialize(self):
        self.contact_process.Execute()
        self._AssignALMParameters()

    def ExecuteInitializeSolutionStep(self):
        self.contact_process.ExecuteInitializeSolutionStep()
        self._AssignALMParameters()

    def _AssignALMParameters(self):
        process_info = self.analysis_model_part.ProcessInfo
        process_info[KM.INITIAL_PENALTY] = self.settings["alm_parameters"]["penalty"].GetDouble()
        process_info[KM.SCALE_FACTOR] = self.settings["alm_parameters"]["scale_factor"].GetDouble()
