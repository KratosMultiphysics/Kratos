# Importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.DEMApplication as DEM

def Factory(settings, Model):
    if not isinstance(settings, KratosMultiphysics.Parameters):
        raise Exception("expected input shall be a Parameters object, encapsulating a json string")
    return MultiaxialControlModuleGeneralized2DProcess(Model, settings["Parameters"])

## All the processes python should be derived from "Process"
class MultiaxialControlModuleGeneralized2DProcess(KratosMultiphysics.Process):

    """This process assigns a given value (vector) to the nodes belonging a certain submodelpart

    Only the member variables listed below should be accessed directly.

    Public member variables:
    Model -- the container of the different model parts.
    settings -- Kratos parameters containing solver settings.
    """

    def __init__(self, Model, settings):

        """ The default constructor of the class

        Keyword arguments:
        self -- It signifies an instance of a class.
        Model -- the container of the different model parts.
        settings -- Kratos parameters containing solver settings.
        """

        KratosMultiphysics.Process.__init__(self)

        '''default_settings = KratosMultiphysics.Parameters("""
        {
            "help"                 : "This process applies a given stress to a sample by imposing displacements. The controlled, adaptable displacement is meant to achieve the target stress",
            "dem_model_part_name"      : "please_specify_model_part_name",
            "fem_model_part_name"      : "please_specify_model_part_name",
            "Parameters" : {
                "Parameters"    : {
                    "control_module_delta_time": 2.0e-8,
                    "perturbation_tolerance": 1.0e-4,
                    "perturbation_period": 10,
                    "max_reaction_rate_factor": 10.0,
                    "stiffness_averaging_time_interval": 2.0e-6,
                    "velocity_averaging_time_interval": 2.0e-4,
                    "reaction_averaging_time_interval": 6.0e-8,
                    "output_interval": 0
                },
                "list_of_actuators" : [{
                    "Parameters"    : {
                        "actuator_name": "Z",
                        "initial_velocity" : 0.0,
                        "compression_length" : 1.0,
                        "young_modulus" : 7.0e9
                    },
                    "list_of_dem_boundaries": [{
                        "model_part_name" : "dems",
                        "outer_normal": [0.0,0.0,1.0]
                    }],
                    "list_of_fem_boundaries": [],
                    "target_stress_table": {
                        "input_variable": "TIME",
                        "output_variable": "TARGET_STRESS",
                        "data": [
                            [0.0, 0.0],
                            [2.0E-2, 0.0]
                        ]
                    }
                },{
                    "Parameters"    : {
                        "actuator_name": "Radial",
                        "initial_velocity" : 0.0,
                        "compression_length" : 0.0505,
                        "young_modulus" : 7.0e9
                    },
                    "list_of_dem_boundaries": [],
                    "list_of_fem_boundaries": [{
                        "model_part_name" : "1",
                        "outer_normal": [0.0,0.0,0.0]
                    }],
                    "target_stress_table": {
                        "input_variable": "TIME",
                        "output_variable": "TARGET_STRESS",
                        "data": [
                            [0.0, 0.0],
                            [2.0E-2, -1.0e9]
                        ]
                    }
                }]
            }
        }
        """
        )'''

        #settings.ValidateAndAssignDefaults(default_settings)

        dem_model_part = Model[settings["dem_model_part_name"].GetString()]
        fem_model_part = Model[settings["fem_model_part_name"].GetString()]
        self.cplusplus_version_of_this_process = DEM.MultiaxialControlModuleGeneralized2DUtilities(dem_model_part, fem_model_part, settings)

    def ExecuteInitialize(self):
        self.cplusplus_version_of_this_process.ExecuteInitialize()

    def ExecuteInitializeSolutionStep(self):
        self.cplusplus_version_of_this_process.ExecuteInitializeSolutionStep()

    def ExecuteFinalizeSolutionStep(self):
        self.cplusplus_version_of_this_process.ExecuteFinalizeSolutionStep()

