import KratosMultiphysics
import KratosMultiphysics.FluidDynamicsApplication as KratosFluid

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ApplyTaitEquationProcess(Model, settings["Parameters"])

class ApplyTaitEquationProcess(KratosMultiphysics.Process):

    def __init__(self, Model, settings):

        KratosMultiphysics.Process.__init__(self)

        # Check the default values
        default_settings = KratosMultiphysics.Parameters("""
            {
                "help"                     : "This proces sets a given scalar value for the density in all the nodes of FluidModelPart",
                "interval"                 : [0.0,1e+30],
                "model_part_name"          : "",
                "rho_0"                    : 0,                  
                "p_0"                      : 0,                    
                "theta"                    : 0,
                "k_0"                      : 0 
            }
            """
            )

        settings.ValidateAndAssignDefaults(default_settings)

        # getting the ModelPart from the Model
        self.model_part = Model[settings["model_part_name"].GetString()]
        self.settings = settings

    def ExecuteInitializeSolutionStep(self):    
        # Get the parameters
        rho_0 = self.settings["rho_0"].GetDouble()
        p_0 = self.settings["p_0"].GetDouble()
        k_0 = self.settings["k_0"].GetDouble()
        theta = self.settings["theta"].GetDouble()

        for node in self.model_part.Nodes:
            #c = node.GetSolutionStepValue(KratosMultiphysics.SOUND_VELOCITY)
            #print(c)
            p = node.GetSolutionStepValue(KratosMultiphysics.PRESSURE)
            # Calculate the new density magnitude
            mod_rho = rho_0*((p-p_0)/k_0 + 1)**(1/theta)
            # Calculate the new c magnitude
            modified_c = (k_0*theta*(mod_rho/rho_0)**(theta - 1)*(1/rho_0))**(1/2)
            # Setting new solution on the node 
            node.SetSolutionStepValue(KratosMultiphysics.DENSITY,0,mod_rho)
            node.SetSolutionStepValue(KratosMultiphysics.SOUND_VELOCITY,0,modified_c)