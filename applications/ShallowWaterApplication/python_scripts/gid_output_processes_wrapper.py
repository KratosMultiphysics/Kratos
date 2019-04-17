import KratosMultiphysics as KM
import KratosMultiphysics.ShallowWaterApplication as SW

def Factory(settings, model):
    if not isinstance(settings, KM.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return GidOutputProcessesWrapper(model, settings["Parameters"])

## This process sets the value of a scalar variable using the AssignScalarVariableProcess.
class GidOutputProcessesWrapper(KM.Process):

    template_settings = KM.Parameters("""
    {
        "Parameters" : {
            "model_part_name"        : "unknown",
            "output_name"            : "file",
            "postprocess_parameters" : {
                "result_file_configuration" : {
                    "gidpost_flags"              : {},
                    "file_label"                 : "time",
                    "output_control_type"        : "step",
                    "output_frequency"           : 1.0,
                    "flush_after_output"         : false,
                    "body_output"                : true,
                    "node_output"                : false,
                    "skin_output"                : false,
                    "plane_output"               : [],
                    "nodal_results"              : [],
                    "nodal_nonhistorical_results": [],
                    "nodal_flags_results"        : [],
                    "elemental_conditional_flags_results": [],
                    "gauss_point_results"        : [],
                    "additional_list_files"      : []
                }
            }
        }
    }
    """)

    def __init__(self, model, settings):

        KM.Process.__init__(self)

        default_settings = KM.Parameters("""
        {
            "file_label"          : "time",
            "output_control_type" : "step",
            "output_frequency"    : 1.0,
            "gidpost_flags"       : {},
            "post_process_list"   : [{
                "model_part_name"       : "unknown_name",
                "output_name"           : "file_name",
                "body_output"           : true,
                "node_output"           : false,
                "nodal_results"         : [],
                "nodal_flags_results"   : [],
                "gauss_point_results"   : [],
                "nodal_nonhistorical_results": [],
                "elemental_conditional_flags_results": []
            }],
            "id_renumbering_process" : {}
        }
        """)
        settings.ValidateAndAssignDefaults(default_settings)

        self.template_settings["Parameters"]["postprocess_parameters"]["result_file_configuration"]["file_label"] = settings["file_label"]
        self.template_settings["Parameters"]["postprocess_parameters"]["result_file_configuration"]["output_control_type"] = settings["output_control_type"]
        self.template_settings["Parameters"]["postprocess_parameters"]["result_file_configuration"]["output_frequency"] = settings["output_frequency"]
        self.template_settings["Parameters"]["postprocess_parameters"]["result_file_configuration"]["gidpost_flags"] = settings["gidpost_flags"]

        self.output_processes = []
        import KratosMultiphysics.gid_output_process as gid_output

        for item in settings["post_process_list"]:
            item.ValidateAndAssignDefaults(default_settings["post_process_list"][0])
            item_settings = self.BuildOutputParameters(item)            
            self.output_processes.append(gid_output.Factory(item_settings, model))

        unique_id_settings = KM.Parameters("""{}""")
        unique_id_settings.AddValue("Parameters", settings["id_renumbering_process"])
        import KratosMultiphysics.ShallowWaterApplication.id_renumbering_process as id_renumbering
        self.unique_id_process = id_renumbering.Factory(unique_id_settings, model)

        self.computing_model_part = model[settings["post_process_list"][0]["model_part_name"].GetString()]
        self.topographic_model_part = model[settings["post_process_list"][1]["model_part_name"].GetString()]

    def ExecuteInitialize(self):
        for process in self.output_processes:
            process.ExecuteInitialize()

    def ExecuteBeforeSolutionLoop(self):
        self.unique_id_process.ExecuteBeforeOutputStep()
        for process in self.output_processes:
            process.ExecuteBeforeSolutionLoop()
        self.unique_id_process.ExecuteAfterOutputStep()

    def ExecuteInitializeSolutionStep(self):
        for process in self.output_processes:
            process.ExecuteInitializeSolutionStep()

    def ExecuteFinalizeSolutionStep(self):
        for process in self.output_processes:
            process.ExecuteFinalizeSolutionStep()

    def ExecuteBeforeOutputStep(self):
        for process in self.output_processes:
            process.ExecuteBeforeOutputStep()

    def IsOutputStep(self):
        is_output_step = False
        for process in self.output_processes:
            if process.IsOutputStep():
                return True
        return False

    def PrintOutput(self):
        self.unique_id_process.ExecuteBeforeOutputStep()
        time = self.computing_model_part.ProcessInfo[KM.TIME]
        step = self.computing_model_part.ProcessInfo[KM.STEP]
        self.topographic_model_part.ProcessInfo[KM.TIME] = time
        self.topographic_model_part.ProcessInfo[KM.STEP] = step
        for process in self.output_processes:
            process.PrintOutput()
        self.unique_id_process.ExecuteAfterOutputStep()

    def ExecuteAfterOutputStep(self):
        for process in self.output_processes:
            process.ExecuteAfterOutputStep()

    def ExecuteFinalize(self):
        for process in self.output_processes:
            process.ExecuteFinalize()

    def Check(self):
        for process in self.output_processes:
            process.Check()

    def BuildOutputParameters(self, item_settings):
        sub_process_settings = self.template_settings.Clone()
        sub_process_settings["Parameters"]["model_part_name"] = item_settings["model_part_name"]
        sub_process_settings["Parameters"]["output_name"] = item_settings["output_name"]
        sub_process_settings["Parameters"]["postprocess_parameters"]["result_file_configuration"]["body_output"] = item_settings["body_output"]
        sub_process_settings["Parameters"]["postprocess_parameters"]["result_file_configuration"]["node_output"] = item_settings["node_output"]
        sub_process_settings["Parameters"]["postprocess_parameters"]["result_file_configuration"]["nodal_results"] = item_settings["nodal_results"]
        sub_process_settings["Parameters"]["postprocess_parameters"]["result_file_configuration"]["nodal_flags_results"] = item_settings["nodal_flags_results"]
        sub_process_settings["Parameters"]["postprocess_parameters"]["result_file_configuration"]["gauss_point_results"] = item_settings["gauss_point_results"]
        sub_process_settings["Parameters"]["postprocess_parameters"]["result_file_configuration"]["nodal_nonhistorical_results"] = item_settings["nodal_nonhistorical_results"]
        sub_process_settings["Parameters"]["postprocess_parameters"]["result_file_configuration"]["elemental_conditional_flags_results"] = item_settings["elemental_conditional_flags_results"]
        return sub_process_settings
