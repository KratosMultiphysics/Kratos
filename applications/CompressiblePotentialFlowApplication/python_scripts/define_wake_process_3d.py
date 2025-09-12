import KratosMultiphysics
import KratosMultiphysics.CompressiblePotentialFlowApplication as CPFApp
import time as time

def Factory(settings, Model):
    if(not isinstance(settings, KratosMultiphysics.Parameters)):
        raise Exception(
            "expected input shall be a Parameters object, encapsulating a json string")

    return DefineWakeProcess3D(Model, settings["Parameters"])

class DefineWakeProcess3D(KratosMultiphysics.Process):
    def __init__(self, Model, settings):
        KratosMultiphysics.Process.__init__(self)

        # Check default settings
        default_settings = KratosMultiphysics.Parameters(r'''{
            "trailing_edge_model_part_name" : "",
            "body_model_part_name"          : "",
            "check_wake_condition_tolerance": 1e-1,
            "wake_process_cpp_parameters":    {
                "wake_normal"                      : [0.0,0.0,1.0],
                "wake_stl_file_name"               : "",
                "wake_dr_translation"              : [0.0,0.0,0.0],
                "tolerance"                        : 1e-9,
                "visualize_wake_vtk"               : false,
                "upper_surface_model_part_name"    : "",
                "lower_surface_model_part_name"    : "",
                "root_points_model_part_name"      : "",
                "tip_points_model_part_name"       : "",
                "blunt_te_surface_model_part_name" : "",
                "shed_wake_from_trailing_edge"     : false,
                "shed_wake_length"                 : 12.5,
                "shed_wake_element_size"           : 0.2,
                "shed_wake_grow_factor"            : 1.05,
                "shed_wake_projection_root_edge"   : 0.0,
                "echo_level"                       : 0
            }
        }''')
        settings.ValidateAndAssignDefaults(default_settings)

        self.model = Model

        self.check_wake_condition_tolerance = settings["check_wake_condition_tolerance"].GetDouble()

        self.wake_process_cpp_parameters = settings["wake_process_cpp_parameters"]

        self.echo_level = self.wake_process_cpp_parameters["echo_level"].GetInt() - 1

        body_model_part_name = settings["body_model_part_name"].GetString()
        if body_model_part_name == "":
            err_msg = "Empty model_part_name in DefineWakeProcess3D\n"
            err_msg += "Please specify the model part that contains the body nodes"
            raise Exception(err_msg)
        self.body_model_part = Model[body_model_part_name]

        self.root_model_part = self.body_model_part.GetRootModelPart()

        trailing_edge_model_part_name = settings["trailing_edge_model_part_name"].GetString()
        if trailing_edge_model_part_name == "":
            err_msg = "Empty trailing_edge_model_part_name in DefineWakeProcess3D\n"
            err_msg += "Please specify the model part that contains the trailing edge nodes"
            raise Exception(err_msg)
        self.trailing_edge_model_part = Model[trailing_edge_model_part_name]

    def ExecuteInitialize(self):
        start_time = time.time()

        CPFApp.Define3DWakeProcess(self.trailing_edge_model_part, self.body_model_part,
                                   self.wake_process_cpp_parameters).ExecuteInitialize()
        
        exe_time = time.time() - start_time
        KratosMultiphysics.Logger.PrintInfo(
            'DefineWakeProcess3D', 'Executing Define3DWakeProcess took ', round(exe_time, 2), ' sec')
        KratosMultiphysics.Logger.PrintInfo(
            'DefineWakeProcess3D', 'Executing Define3DWakeProcess took ', round(exe_time/60, 2), ' min')

    def ExecuteFinalizeSolutionStep(self):
        if not self.root_model_part.HasSubModelPart("wake_elements_model_part"):
            raise Exception("root model part does not have a wake_elements_model_part")
        else:
            self.wake_sub_model_part = self.root_model_part.GetSubModelPart("wake_elements_model_part")

        CPFApp.PotentialFlowUtilities.CheckIfWakeConditionsAreFulfilled3D(self.wake_sub_model_part, self.check_wake_condition_tolerance, self.echo_level)
        CPFApp.PotentialFlowUtilities.ComputePotentialJump3D(self.wake_sub_model_part)