import KratosMultiphysics
import KratosMultiphysics.FluidDynamicsApplication as KratosFluid

def Factory(settings, Model):
    if(type(settings) != KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return ExplosionPressureProcess(Model, settings["Parameters"])

class ExplosionPressureProcess(KratosMultiphysics.Process):

    def __init__(self, Model, settings):

        KratosMultiphysics.Process.__init__(self)

        # Check the default values
        default_settings = KratosMultiphysics.Parameters("""
            {
                "help"                     : "This proces sets a given pressure value for all nodes inside circle",
                "interval"                 : [0.0,0.0],
                "model_part_name"	        : "",           
                "p"                        : 0,
                "rho" 		            : 0,
                "radius"                   : 0,
                "p_0"                      : 0,
                "k_0"			    : 0,
                "theta"		    : 0

            }
            """
            )

        settings.ValidateAndAssignDefaults(default_settings)

        # getting the ModelPart from the Model
        model_part = Model[settings["model_part_name"].GetString()]

        p = settings["p"].GetDouble()
        rho = settings["rho"].GetDouble()
        radius = settings["radius"].GetDouble()
        p_0 = settings["p_0"].GetDouble()
        theta = settings["theta"].GetDouble()
        k_0 = settings["k_0"].GetDouble()

        for node in model_part.Nodes:
            x = node.X
            y = node.Y
            if (x**2 + y**2 < radius**2):
                mod_rho = rho*((p-p_0)/k_0 + 1)**(1/theta)
                modified_c = (k_0*theta*(mod_rho/rho)**(theta - 1)*(1/rho))**(1/2)
                node.SetSolutionStepValue(KratosMultiphysics.SOUND_VELOCITY,0,modified_c)
                node.SetSolutionStepValue(KratosMultiphysics.DENSITY,0,mod_rho)
            else:
                c = (k_0*theta**(theta - 1)*(1/rho))**(1/2)
                node.SetSolutionStepValue(KratosMultiphysics.SOUND_VELOCITY,0,c)
                node.SetSolutionStepValue(KratosMultiphysics.DENSITY,0,rho)
                
           
