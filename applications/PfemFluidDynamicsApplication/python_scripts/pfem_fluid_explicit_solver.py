from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
#import os
#import kratos core and applications
import KratosMultiphysics
import KratosMultiphysics.PfemFluidDynamicsApplication as KratosPfemFluid
#import KratosMultiphysics.SolidMechanicsApplication as KratosSolidMechanics

import pfem_fluid_solver as BaseSolver

# Check that KratosMultiphysics was imported in the main script
KratosMultiphysics.CheckForPreviousImport()

def CreateSolver(main_model_part, custom_settings):
    return PfemFluidExplicitSolver(main_model_part, custom_settings)

class PfemFluidExplicitSolver(BaseSolver.PfemFluidSolver):

    def __init__(self, main_model_part, custom_settings):
        #TODO: shall obtain the computing_model_part from the MODEL once the object is implemented
        self.main_model_part = main_model_part          
        ##settings string in json format
        explicit_solver_settings = KratosMultiphysics.Parameters("""
        {  
            "echo_level": 1,
            "buffer_size": 3,
            "solver_type": "pfem_fluid_explicit_solver",
             "model_import_settings": {
                "input_type": "mdpa",
                "input_filename": "unknown_name",
                "input_file_label": 0
            },
            "dofs"                : [],
            "stabilization_factor": 1.0,
            "reform_dofs_at_each_step": true,
            "line_search": false,
            "compute_reactions": false,
            "compute_contact_forces": false,
            "block_builder": false,
            "clear_storage": false,
            "component_wise": false,
            "move_mesh_flag": true,
            "predictor_corrector": true,
            "time_order": 2,
            "maximum_velocity_iterations": 1,
            "maximum_pressure_iterations": 7,
            "velocity_tolerance": 1e-5,
            "pressure_tolerance": 1e-5,
            "pressure_linear_solver_settings":  {
                "solver_type"                    : "AMGCL",
                "max_iteration"                  : 5000,
                "tolerance"                      : 1e-9,
                "provide_coordinates"            : false,
                "scaling"                        : false,
                "smoother_type"                  : "damped_jacobi",
                "krylov_type"                    : "cg",
                "coarsening_type"                : "aggregation",
                "verbosity"                      : 0
            },
            "velocity_linear_solver_settings": {
                "solver_type"                    : "BICGSTABSolver",
                "max_iteration"                  : 5000,
                "tolerance"                      : 1e-9,
                "preconditioner_type"            : "None",
                "scaling"                        : false
            },
            "solving_strategy_settings":{
               "time_step_prediction_level": 0,
               "max_delta_time": 1.0e-5,
               "fraction_delta_time": 0.9,
               "rayleigh_damping": false,
               "rayleigh_alpha": 0.0,
               "rayleigh_beta" : 0.0
            },
            "bodies_list": [
                {"body_name":"body1",
                "parts_list":["Part1"]
                },
                {"body_name":"body2",
                "parts_list":["Part2","Part3"]
                }
            ],
            "problem_domain_sub_model_part_list": ["fluid_model_part"],
            "processes_sub_model_part_list": [""]
        } 
        """)


        ##overwrite the default settings with user-provided parameters
        self.settings = custom_settings
        self.settings.ValidateAndAssignDefaults(explicit_solver_settings)


         # Construct the base solver.
        super(PfemFluidExplicitSolver,self).__init__(main_model_part,custom_settings)

        print("Construction of Pfem Fluid Explicit Solver finished.")

    def Initialize(self):

        print("::[Pfem Fluid Explicit Solver]:: -START-")

        # Get the computing model part
        self.computing_model_part = self.GetComputingModelPart()

        #mechanical_scheme = KratosSolidMechanics.ExplicitCentralDifferencesScheme(self.settings["max_delta_time"].GetDouble(),
        #                                                                          self.settings["fraction_delta_time"].GetDouble(),
        #                                                                          self.settings["time_step_prediction_level"].GetDouble(),
        #                                                                          self.settings["rayleigh_damping"].GetBool())

        mechanical_scheme = KratosPfemFluid.FirstOrderForwardEulerScheme(1.0e-4,
                                                                         1.0,
                                                                         0,
                                                                         0)
        import linear_solver_factory
        linear_solver = linear_solver_factory.ConstructSolver(self.settings["velocity_linear_solver_settings"])
        #self.fluid_solver = KratosPfemFluid.ExplicitStrategy(self.computing_model_part,
        #                                                     mechanical_scheme,
        #                                                     linear_solver,
        #                                                     self.settings["compute_reactions"].GetBool(),
        #                                                     self.settings["reform_dofs_at_each_step"].GetBool(),
        #                                                     self.settings["move_mesh_flag"].GetBool())

        self.fluid_solver = KratosPfemFluid.ExplicitStrategy(self.computing_model_part,
                                                             mechanical_scheme,
                                                             linear_solver,
                                                             False,
                                                             True,
                                                             True)   
        # Set echo_level
        self.fluid_solver.SetEchoLevel(self.settings["echo_level"].GetInt())

        # Set initialize flag
        if( self.main_model_part.ProcessInfo[KratosMultiphysics.IS_RESTARTED] == True ):
            self.mechanical_solver.SetInitializePerformedFlag(True)

        # Check if everything is assigned correctly
        self.fluid_solver.Check()


        print("::[Pfem Fluid  Explicit Solver]:: -END- ")
