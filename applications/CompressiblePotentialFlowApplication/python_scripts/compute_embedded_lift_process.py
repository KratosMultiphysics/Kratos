import KratosMultiphysics
import KratosMultiphysics.CompressiblePotentialFlowApplication as CPFApp

def Factory(settings, Model):
    if( not isinstance(settings,KratosMultiphysics.Parameters) ):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ComputeEmbeddedLiftProcess(Model, settings["Parameters"])

##all the processes python processes should be derived from "python_process"
class ComputeEmbeddedLiftProcess(KratosMultiphysics.Process):
    def __init__(self, Model, settings ):
        KratosMultiphysics.Process.__init__(self)
        default_parameters = KratosMultiphysics.Parameters("""
            {

                "model_part_name" : "please specify the main model part",
                "reference_area": 1
            }  """)

        settings.ValidateAndAssignDefaults(default_parameters)
        self.main_model_part=Model.GetModelPart(settings["model_part_name"].GetString()).GetRootModelPart()
        self.result_force=KratosMultiphysics.Vector(3)
        self.process=CPFApp.ComputeEmbeddedLiftProcess(self.main_model_part,self.result_force)

    def ExecuteFinalizeSolutionStep(self):
        print("wip_compute_lift_level_set_process")
        self.process.Execute()
        self.lift_coefficient = self.result_force[1]
        self.drag_coefficient = self.result_force[0]

        KratosMultiphysics.Logger.PrintInfo('ComputeEmbeddedLiftProcess',' Cl = ', self.lift_coefficient)
        KratosMultiphysics.Logger.PrintInfo('ComputeEmbeddedLiftProcess',' Cd = ', self.drag_coefficient)

        self.main_model_part.ProcessInfo.SetValue(CPFApp.LIFT_COEFFICIENT, self.lift_coefficient)
        self.main_model_part.ProcessInfo.SetValue(CPFApp.DRAG_COEFFICIENT, self.drag_coefficient)