import KratosMultiphysics
import KratosMultiphysics.PfemMeltingApplication as PfemM

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ApplyLaserProcess(Model, settings["Parameters"])

class ApplyLaserProcess(KratosMultiphysics.Process):

    def __init__(self, Model, settings):
         
        KratosMultiphysics.Process.__init__(self)

        # Check the default values
        default_settings = KratosMultiphysics.Parameters( """
        {
            "model_part_name" : "CHOOSE_FLUID_MODELPART_NAME"
        }  """ )


        # Get the fluid model part from the Model container
        self.fluid_model_part = Model[settings["model_part_name"].GetString()]
        
        # Set the Boussinesq force process
        self.ApplyLaserProcess = PfemM.ApplyLaserProcess(self.fluid_model_part, settings)

    def ExecuteInitialize(self):
        self.ApplyLaserProcess.ExecuteInitialize()

    def ExecuteInitializeSolutionStep(self):
        self.ApplyLaserProcess.ExecuteInitializeSolutionStep()
