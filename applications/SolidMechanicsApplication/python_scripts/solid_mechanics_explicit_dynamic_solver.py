from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
#import kratos core and applications
import KratosMultiphysics
import KratosMultiphysics.SolidMechanicsApplication as SolidMechanicsApplication

# Check that KratosMultiphysics was imported in the main script
KratosMultiphysics.CheckForPreviousImport()

# Import the implicit solver (the explicit one is derived from it)
import solid_mechanics_implicit_dynamic_solver

def CreateSolver(main_model_part, custom_settings):
    return ExplicitMechanicalSolver(main_model_part, custom_settings)

class ExplicitMechanicalSolver(solid_mechanics_implicit_dynamic_solver.ImplicitMechanicalSolver):
    
    ##constructor. the constructor shall only take care of storing the settings 
    ##and the pointer to the main_model part. This is needed since at the point of constructing the 
    ##model part is still not filled and the variables are not yet allocated
    ##
    ##real construction shall be delayed to the function "Initialize" which 
    ##will be called once the model is already filled
    def __init__(self, main_model_part, custom_settings): 
        
        #TODO: shall obtain the compute_model_part from the MODEL once the object is implemented
        self.main_model_part = main_model_part    
        
        #TODO: remove unnecessary fields for the Explicit solver from the defaults
        ##settings string in json format
        default_settings = KratosMultiphysics.Parameters("""
        {
            "solver_type": "solid_mechanics_explicit_dynamic_solver",
            "echo_level": 0,
            "solution_type": "Dynamic",
            "time_integration_method": "Explicit",
            "scheme_type": "CentralDifferences",
            "model_import_settings": {
                "input_type": "mdpa",
                "input_filename": "unknown_name"
            },
            "rotation_dofs": false,
            "pressure_dofs": false,
            "stabilization_factor": 1.0,
            "reform_dofs_at_each_iteration": false,
            "compute_reactions": true,
            "move_mesh_flag": true,
            "max_delta_time": 1.0e-5,
            "fraction_delta_time": 0.9,
            "time_step_prediction_level": 0,
            "rayleigh_damping": false,
            "linear_solver_settings":{
                "solver_type": "Super LU",
                "max_iteration": 500,
                "tolerance": 1e-9,
                "scaling": false,
                "verbosity": 1
            },
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
        
        print("Construction of MechanicalSolver finished")

    def AddVariables(self):
        
        # Add displacements
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        # Add dynamic variables
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.ACCELERATION)
        # Add reactions for the displacements
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION)
        # Add nodal force variables
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.INTERNAL_FORCE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.EXTERNAL_FORCE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.CONTACT_FORCE)
        # Add specific variables for the problem conditions
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.POSITIVE_FACE_PRESSURE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NEGATIVE_FACE_PRESSURE)
        self.main_model_part.AddNodalSolutionStepVariable(SolidMechanicsApplication.POINT_LOAD)
        self.main_model_part.AddNodalSolutionStepVariable(SolidMechanicsApplication.LINE_LOAD)
        self.main_model_part.AddNodalSolutionStepVariable(SolidMechanicsApplication.SURFACE_LOAD)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VOLUME_ACCELERATION)

        if self.settings["time_integration_method"].GetString() == "Explicit":
            # Add specific variables for the explicit time integration scheme
            self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NODAL_MASS)
            self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.FORCE_RESIDUAL)
            self.main_model_part.AddNodalSolutionStepVariable(SolidMechanicsApplication.MIDDLE_VELOCITY)
            
        if self.settings["rotation_dofs"].GetBool():
            # Add specific variables for the problem (rotation dofs)
            self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.ROTATION)
            self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.TORQUE)
            self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.ANGULAR_VELOCITY)
            self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.ANGULAR_ACCELERATION)
            
        if self.settings["pressure_dofs"].GetBool():
            # Add specific variables for the problem (pressure dofs)
            self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.PRESSURE)
            self.main_model_part.AddNodalSolutionStepVariable(SolidMechanicsApplication.PRESSURE_REACTION)
                    
        print("::[Mechanical Solver]:: Variables ADDED")
    
    def Initialize(self):

        print("::[Mechanical Solver]:: -START-")
        
        # Get the solid_computational_model_part 
        self.compute_model_part = self.GetComputeModelPart()
        
        # Solution scheme creation
        mechanical_scheme = self._GetSolutionScheme(self.settings["max_delta_time"].GetDouble(), 
                                                    self.settings["fraction_delta_time"].GetDouble(),
                                                    self.settings["time_step_prediction_level"].GetInt(),
                                                    self.settings["rayleigh_damping"].GetBool())
        
        # Mechanical solver creation
        self._CreateMechanicalSolver(mechanical_scheme,
                                     self.settings["compute_reactions"].GetBool(),
                                     self.settings["reform_dofs_at_each_iteration"].GetBool(),
                                     self.settings["move_mesh_flag"].GetBool())

        # Set the stabilization factor
        self.main_model_part.ProcessInfo[KratosMultiphysics.STABILIZATION_FACTOR] = self.settings["stabilization_factor"].GetDouble()

        # Set echo_level
        self.mechanical_solver.SetEchoLevel(self.settings["echo_level"].GetInt())

        # Check if everything is assigned correctly
        self.Check();

        print("::[Mechanical Solver]:: -END- ")
        
    #### Specific internal functions ####
    def _GetSolutionScheme(self, max_delta_time, fraction_delta_time, time_step_prediction_level, rayleigh_damping):

        if self.settings["scheme_type"].GetString() == "CentralDifferences":

            mechanical_scheme = SolidMechanicsApplication.ExplicitCentralDifferencesScheme(max_delta_time, 
                                                                                           fraction_delta_time, 
                                                                                           time_step_prediction_level, 
                                                                                           rayleigh_damping)
        
        else:
            raise(self.settings["scheme_type"].GetString()," not implemented yet.")
          
        return mechanical_scheme
        
    def _CreateMechanicalSolver(self, mechanical_scheme, compute_reactions, reform_step_dofs, move_mesh_flag):

        self.mechanical_solver = SolidMechanicsApplication.ExplicitStrategy(self.compute_model_part, 
                                                                            mechanical_scheme, 
                                                                            self.linear_solver, 
                                                                            compute_reactions, 
                                                                            reform_step_dofs, 
                                                                            move_mesh_flag)

        self.mechanical_solver.SetRebuildLevel(0) # 1 to recompute the mass matrix in each explicit step 



