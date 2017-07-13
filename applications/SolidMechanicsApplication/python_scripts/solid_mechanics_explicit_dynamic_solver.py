from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
#import kratos core and applications
import KratosMultiphysics
import KratosMultiphysics.SolidMechanicsApplication as KratosSolid

# Check that KratosMultiphysics was imported in the main script
KratosMultiphysics.CheckForPreviousImport()

# Import the mechanical solver base class
import solid_mechanics_solver

def CreateSolver(main_model_part, custom_settings):
    return ExplicitMechanicalSolver(main_model_part, custom_settings)

class ExplicitMechanicalSolver(solid_mechanics_solver.MechanicalSolver):
    
    ##constructor. the constructor shall only take care of storing the settings 
    ##and the pointer to the main_model part. This is needed since at the point of constructing the 
    ##model part is still not filled and the variables are not yet allocated
    ##
    ##real construction shall be delayed to the function "Initialize" which 
    ##will be called once the model is already filled
    def __init__(self, main_model_part, custom_settings): 
        
        #TODO: shall obtain the computing_model_part from the MODEL once the object is implemented
        self.main_model_part = main_model_part    
        
        #TODO: remove unnecessary fields for the Explicit solver from the defaults
        ##settings string in json format
        default_settings = KratosMultiphysics.Parameters("""
        {
            "solver_type": "solid_mechanics_explicit_dynamic_solver",
            "echo_level": 0,
            "buffer_size": 2,
            "solution_type": "Dynamic",
            "time_integration_method": "Explicit",
            "scheme_type": "CentralDifferences",
            "model_import_settings": {
                "input_type": "mdpa",
                "input_filename": "unknown_name",
                "input_file_label": 0
            },
            "rotation_dofs": false,
            "pressure_dofs": false,
            "stabilization_factor": 1.0,
            "reform_dofs_at_each_step": false,
            "compute_reactions": true,
            "move_mesh_flag": true,
            "clear_storage": false,
            "max_delta_time": 1.0e-5,
            "fraction_delta_time": 0.9,
            "time_step_prediction_level": 0,
            "rayleigh_damping": false,
            "linear_solver_settings":{
                "solver_type": "SuperLUSolver",
                "max_iteration": 500,
                "tolerance": 1e-9,
                "scaling": false,
                "verbosity": 1
            },
            "bodies_list": [],
            "problem_domain_sub_model_part_list": ["solid_model_part"],
            "processes_sub_model_part_list": [""]
        }
        """)

        ##overwrite the default settings with user-provided parameters
        self.settings = custom_settings
        self.settings.ValidateAndAssignDefaults(default_settings)
        
        #construct the linear solver
        import linear_solver_factory
        self.linear_solver = linear_solver_factory.ConstructSolver(self.settings["linear_solver_settings"])
        
        print("Construction of Explicit Dynamics Mechanical Solver finished")
   

    def Initialize(self):

        print("::[Mechanical Solver]:: -START-")
        
        # Get the solid computing model part
        self.computing_model_part = self.GetComputingModelPart()
        
        # Solution scheme creation
        mechanical_scheme = self._GetSolutionScheme(self.settings["max_delta_time"].GetDouble(), 
                                                    self.settings["fraction_delta_time"].GetDouble(),
                                                    self.settings["time_step_prediction_level"].GetInt(),
                                                    self.settings["rayleigh_damping"].GetBool())
        
        # Mechanical solver creation
        self._CreateMechanicalSolver(mechanical_scheme,
                                     self.settings["compute_reactions"].GetBool(),
                                     self.settings["reform_dofs_at_each_step"].GetBool(),
                                     self.settings["move_mesh_flag"].GetBool())

        # Set echo_level
        self.mechanical_solver.SetEchoLevel(self.settings["echo_level"].GetInt())

        # Set initialize flag
        if( self.main_model_part.ProcessInfo[KratosMultiphysics.IS_RESTARTED] == True ):
            self.mechanical_solver.SetInitializePerformedFlag(True)

        # Check if everything is assigned correctly
        self.Check();

        print("::[Mechanical Solver]:: -END- ")
        

    def AddVariables(self):
        
        solid_mechanics_solver.MechanicalSolver.AddVariables(self)

        if self.settings["time_integration_method"].GetString() == "Explicit":
            # Add specific variables for the explicit time integration scheme
            self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NODAL_MASS)
            self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.FORCE_RESIDUAL)
            self.main_model_part.AddNodalSolutionStepVariable(KratosSolid.MIDDLE_VELOCITY)

                    
        print("::[Mechanical Solver]:: Explicit Variables ADDED")


    #### Specific internal functions ####
    def _GetSolutionScheme(self, max_delta_time, fraction_delta_time, time_step_prediction_level, rayleigh_damping):

        if self.settings["scheme_type"].GetString() == "CentralDifferences":

            mechanical_scheme = KratosSolid.ExplicitCentralDifferencesScheme(max_delta_time, 
                                                                             fraction_delta_time, 
                                                                             time_step_prediction_level, 
                                                                             rayleigh_damping)
        
        else:
            raise(self.settings["scheme_type"].GetString()," not implemented yet.")
          
        return mechanical_scheme
        
    def _CreateMechanicalSolver(self, mechanical_scheme, compute_reactions, reform_step_dofs, move_mesh_flag):

        self.mechanical_solver = KratosSolid.ExplicitStrategy(self.computing_model_part, 
                                                              mechanical_scheme, 
                                                              self.linear_solver, 
                                                              compute_reactions, 
                                                              reform_step_dofs, 
                                                              move_mesh_flag)

        self.mechanical_solver.SetRebuildLevel(0) # 1 to recompute the mass matrix in each explicit step 



