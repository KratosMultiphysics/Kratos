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
                "reference_area": 1.0
            }  """)

        settings.ValidateAndAssignDefaults(default_parameters)
        self.main_model_part=Model.GetModelPart(settings["model_part_name"].GetString()).GetRootModelPart()
        self.result_force=KratosMultiphysics.Vector(3)
        self.reference_area =  settings["reference_area"].GetDouble()

    def ExecuteFinalizeSolutionStep(self):
        if (self.main_model_part.ProcessInfo.GetValue(KratosMultiphysics.DOMAIN_SIZE)==2):
            CPFApp.ComputeEmbeddedLiftProcess2D(self.main_model_part,self.result_force).Execute()
        else:
            raise(Exception("Dimension of the problem is not 2. Only 2D cases are currently supported."))
        self.lift_coefficient = self.result_force[1]/self.reference_area
        self.drag_coefficient = self.result_force[0]/self.reference_area

        KratosMultiphysics.Logger.PrintInfo('ComputeEmbeddedLiftProcess',' Cl = ', self.lift_coefficient)
        KratosMultiphysics.Logger.PrintInfo('ComputeEmbeddedLiftProcess',' Cd = ', self.drag_coefficient)

        self.main_model_part.ProcessInfo.SetValue(CPFApp.LIFT_COEFFICIENT, self.lift_coefficient)
        self.main_model_part.ProcessInfo.SetValue(CPFApp.DRAG_COEFFICIENT, self.drag_coefficient)