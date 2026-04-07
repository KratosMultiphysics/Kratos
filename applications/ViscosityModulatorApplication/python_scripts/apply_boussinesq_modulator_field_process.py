import KratosMultiphysics
import KratosMultiphysics.FluidDynamicsApplication as KratosFluid
import KratosMultiphysics.ViscosityModulatorApplication as CD

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ApplyBoussinesqModulatorFieldProcess(Model, settings["Parameters"])

class ApplyBoussinesqModulatorFieldProcess(KratosMultiphysics.Process):

    def __init__(self, Model, settings):

        KratosMultiphysics.Process.__init__(self)

        # Check the default values
        default_settings = KratosMultiphysics.Parameters( """
        {
            "model_part_name" : "CHOOSE_FLUID_MODELPART_NAME",
            "base_fluid_density" : 1.0,
            "max_density": 2.0,
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
        self.BoussinesqCDProcess = CD.BoussinesqModulatorFieldProcess(self.fluid_model_part, settings)

    def ExecuteInitialize(self):
        self.BoussinesqCDProcess.ExecuteInitialize()

    def ExecuteInitializeSolutionStep(self):
        # NOTE: When used with CoupledFluidThermalSolverBoussinesq, this method
        # should NOT be called — the solver owns the Boussinesq process and drives
        # it explicitly inside the Picard coupling loop via Execute().
        # If this process is listed under "processes" in ProjectParameters.json
        # while also using CoupledFluidThermalSolverBoussinesq, REMOVE IT from
        # the process list to avoid overwriting the body force with stale φ.
        #
        # For backward compatibility with the original (lagged) coupling, this
        # method is intentionally left active. Simply do not include this process
        # in the JSON when switching to the new solver.
        self.BoussinesqCDProcess.ExecuteInitializeSolutionStep()
