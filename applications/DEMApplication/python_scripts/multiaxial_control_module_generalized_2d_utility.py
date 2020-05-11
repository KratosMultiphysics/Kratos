import KratosMultiphysics
import KratosMultiphysics.DEMApplication as Dem


class MultiaxialControlModuleGeneralized2DUtility(object):
    def __init__(self, dem_model_part, dem_fem_boundary_model_part):

        self.dem_model_part = dem_model_part
        self.dem_fem_boundary_model_part = dem_fem_boundary_model_part

        project_parameters_file_name = "sp_2d_rigid_fem_parameters.json"

        with open(project_parameters_file_name,'r') as parameter_file:
            project_parameters = KratosMultiphysics.Parameters(parameter_file.read())

        self.parameters = project_parameters["multiaxial_control_module_generalized_2d_utility"]

        # TODO
        # Blind test
        # compression_length = 0.009144
        # face_area = 0.088343
        # alternate_axis_loading = True
        # limit_velocity = -15.0
        # Negative target_stress means compression.

        '''
        self.parameters = KratosMultiphysics.Parameters( """
        {
            "Parameters"    : {
                "control_module_delta_time": 4.0e-7,
                "velocity_factor" : 1.0,
                "stress_tolerance": 1.0e-2,
                "perturbation_period": 10,
                "update_stiffness": true,
                "stiffness_alpha": 1.0,
                "start_time" : 0.0,
                "stress_averaging_time": 1.0e-5
            },
            "list_of_actuators" : [{
                "Parameters"    : {
                    "actuator_name": "Z",
                    "target_stress_table_id": 1,
                    "initial_velocity" : 0.0,
                    "limit_velocity" : -20.0,
                    "compression_length" : 1.0,
                    "young_modulus" : 7.0e9
                },
                "list_of_dem_boundaries": [{
                    "model_part_name" : "Parts_dems",
                    "outer_normal": [0.0,0.0,1.0]
                }],
                "list_of_fem_boundaries": []
            },{
                "Parameters"    : {
                    "actuator_name": "X",
                    "target_stress_table_id": 1,
                    "initial_velocity" : 0.0,
                    "limit_velocity" : -5.0,
                    "compression_length" : 0.1524,
                    "young_modulus" : 21.0e13
                },
                "list_of_dem_boundaries": [],
                "list_of_fem_boundaries": [{
                    "model_part_name" : "1",
                    "outer_normal": [-1.0,0.0,0.0]
                    },{
                    "model_part_name" : "2",
                    "outer_normal": [1.0,0.0,0.0]
                }]
            },{
                "Parameters"    : {
                    "actuator_name": "Y",
                    "target_stress_table_id": 2,
                    "initial_velocity" : 0.0,
                    "limit_velocity" : -5.0,
                    "compression_length" : 0.1524,
                    "young_modulus" : 21.0e13
                },
                "list_of_dem_boundaries": [],
                "list_of_fem_boundaries": [{
                    "model_part_name" : "3",
                    "outer_normal": [0.0,-1.0,0.0]
                    },{
                    "model_part_name" : "4",
                    "outer_normal": [0.0,1.0,0.0]
                }]
            }]
        }  """ )
        '''
        self.cm_utility = Dem.MultiaxialControlModuleGeneralized2DUtilities(self.dem_model_part, 
                                                                            self.dem_fem_boundary_model_part, 
                                                                            self.parameters)

    def ExecuteInitialize(self):
        self.cm_utility.ExecuteInitialize()

    def ExecuteInitializeSolutionStep(self):
        self.cm_utility.ExecuteInitializeSolutionStep()

    def ExecuteFinalizeSolutionStep(self):
        self.cm_utility.ExecuteFinalizeSolutionStep()