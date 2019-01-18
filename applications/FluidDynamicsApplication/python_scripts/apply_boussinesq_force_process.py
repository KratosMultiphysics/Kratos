import KratosMultiphysics
import KratosMultiphysics.FluidDynamicsApplication as KratosFluid

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ApplyBoussinesqForceProcess(Model, settings["Parameters"])

class ApplyBoussinesqForceProcess(KratosMultiphysics.Process):

    def __init__(self, Model, settings):

        KratosMultiphysics.Process.__init__(self)

        # Check the default values
        default_settings = KratosMultiphysics.Parameters( """
        {
            "model_part_name" : "CHOOSE_FLUID_MODELPART_NAME",
            "gravity" : [0.0,0.0,0.0],
            "ambient_temperature" : 273.0
        }  """ )

        # Note: if the thermal expansion coefficient is not provided, it is computed as 
        # 1/AMBIENT_TEMPERATURE which is the usual assumption for perfect gases.
        if (settings.Has("thermal_expansion_coefficient")):
            default_settings.AddEmptyValue("thermal_expansion_coefficient").SetDouble(0.0)

        settings.ValidateAndAssignDefaults(default_settings)

        # Get the fluid model part from the Model container
        self.fluid_model_part = Model[settings["model_part_name"].GetString()]

        # Save the ambient temperature in the fluid model part ProcessInfo
        ambient_temperature = settings["ambient_temperature"].GetDouble()
        self.fluid_model_part.ProcessInfo.SetValue(KratosMultiphysics.AMBIENT_TEMPERATURE, ambient_temperature)

        # Set the Boussinesq force process
        self.BoussinesqForceProcess = KratosFluid.BoussinesqForceProcess(self.fluid_model_part, settings)

    def ExecuteInitialize(self):
        self.BoussinesqForceProcess.ExecuteInitialize()

    def ExecuteInitializeSolutionStep(self):
        self.BoussinesqForceProcess.ExecuteInitializeSolutionStep()
