import KratosMultiphysics
import KratosMultiphysics.FluidDynamicsApplication as KratosFluid
import KratosMultiphysics.ConvectionDiffusionApplication as CD

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ApplyBoussinesqConcentrationFieldProcess(Model, settings["Parameters"])

class ApplyBoussinesqConcentrationFieldProcess(KratosMultiphysics.Process):

    def __init__(self, Model, settings):

        KratosMultiphysics.Process.__init__(self)

        # Check the default values
        default_settings = KratosMultiphysics.Parameters( """
        {
            "model_part_name" : "CHOOSE_FLUID_MODELPART_NAME",
            "base_fluid_density" : 1000.0,
            "particles_density" : 1200.0,
            "gravity" : [0.0,0.0,0.0],
            "modify_pressure" : false,
            "modify_density": false,
            "r0": [0.0, 0.0, 0.0]
        }  """ )

        settings.ValidateAndAssignDefaults(default_settings)

        # Get the fluid model part from the Model container
        self.fluid_model_part = Model[settings["model_part_name"].GetString()]

        # Save the ambient temperature in the fluid model part ProcessInfo
        # self.fluid_model_part = Model[settings["model_part_name"].GetString()]
        # base_density = settings["base_fluid_density"].GetDouble()
        # particles_density = settings["particles_density"].GetDouble()

        # Set the Boussinesq force process
        self.BoussinesqCDProcess = CD.BoussinesqConcentrationFieldProcess(self.fluid_model_part, settings)

    def ExecuteInitialize(self):
        self.BoussinesqCDProcess.ExecuteInitialize()

    def ExecuteInitializeSolutionStep(self):
        self.BoussinesqCDProcess.ExecuteInitializeSolutionStep()
