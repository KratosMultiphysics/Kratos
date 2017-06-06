from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
#import kratos core and applications
import KratosMultiphysics
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication

# Check that KratosMultiphysics was imported in the main script
KratosMultiphysics.CheckForPreviousImport()

# Import the implicit solver (the explicit one is derived from it)
import structural_mechanics_solver

def CreateSolver(main_model_part, custom_settings):
    return ImplicitMechanicalSolver(main_model_part, custom_settings)

class ImplicitMechanicalSolver(structural_mechanics_solver.MechanicalSolver):
    
    ##constructor. the constructor shall only take care of storing the settings 
    ##and the pointer to the main_model part. This is needed since at the point of constructing the 
    ##model part is still not filled and the variables are not yet allocated
    ##
    ##real construction shall be delayed to the function "Initialize" which 
    ##will be called once the model is already filled
    def __init__(self, main_model_part, custom_settings): 
        
        #TODO: shall obtain the computing_model_part from the MODEL once the object is implemented
        self.main_model_part = main_model_part    
        
        ##settings string in json format
        default_settings = KratosMultiphysics.Parameters("""
        {
            "solver_type": "structural_mechanics_implicit_dynamic_solver",
            "model_import_settings": {
                "input_type": "mdpa",
                "input_filename": "unknown_name"
            },
            "material_import_settings" :{
                "materials_filename": ""
            },
            "echo_level": 0,
            "buffer_size": 2,
            "solution_type": "Dynamic",
            "scheme_type": "Newmark",
            "damp_factor_m" : -0.1,
            "time_integration_method": "Implicit",
            "analysis_type": "Non-Linear",
            "rotation_dofs": false,
            "pressure_dofs": false,
            "stabilization_factor": 1.0,
            "reform_dofs_at_each_step": false,
            "line_search": false,
            "compute_reactions": true,
            "compute_contact_forces": false,
            "block_builder": false,
            "clear_storage": false,
            "move_mesh_flag": true,
            "convergence_criterion": "Residual_criteria",
            "displacement_relative_tolerance": 1.0e-4,
            "displacement_absolute_tolerance": 1.0e-9,
            "residual_relative_tolerance": 1.0e-4,
            "residual_absolute_tolerance": 1.0e-4,
            "max_iteration": 10,
            "split_factor": 10.0,
            "max_number_splits": 3,
            "linear_solver_settings":{
                "solver_type": "Super LU",
                "max_iteration": 500,
                "tolerance": 1e-9,
                "scaling": false,
                "verbosity": 1
            },
            "processes_sub_model_part_list": [""],
            "problem_domain_sub_model_part_list": ["solid_model_part"]
        }
        """)
        
        ##overwrite the default settings with user-provided parameters
        self.settings = custom_settings
        self.settings.ValidateAndAssignDefaults(default_settings)
        
        #construct the linear solver
        import linear_solver_factory
        self.linear_solver = linear_solver_factory.ConstructSolver(self.settings["linear_solver_settings"])
        
        print("Construction of MechanicalSolver finished")

    def Initialize(self):

        print("::[Mechanical Solver]:: -START-")
        
        # Get the solid computing model part
        self.computing_model_part = self.GetComputingModelPart()
        
        # Builder and solver creation
        builder_and_solver = self._GetBuilderAndSolver(self.settings["block_builder"].GetBool())
        
        # Solution scheme creation
        mechanical_scheme = self._GetSolutionScheme(self.settings["scheme_type"].GetString())
        
        # Get the convergence criterion
        mechanical_convergence_criterion = self._GetConvergenceCriterion()
        
        # Mechanical solver creation
        self._CreateMechanicalSolver(mechanical_scheme,
                                     mechanical_convergence_criterion,
                                     builder_and_solver,
                                     self.settings["max_iteration"].GetInt(),
                                     self.settings["compute_reactions"].GetBool(),
                                     self.settings["reform_dofs_at_each_step"].GetBool(),
                                     self.settings["move_mesh_flag"].GetBool(),
                                     self.settings["line_search"].GetBool()
                                     )

        # Set echo_level
        self.mechanical_solver.SetEchoLevel(self.settings["echo_level"].GetInt())

        # Set initialize flag
        if( self.main_model_part.ProcessInfo[KratosMultiphysics.IS_RESTARTED] == True ):
            self.mechanical_solver.SetInitializePerformedFlag(True)
        
        # Check if everything is assigned correctly
        self.Check();

        print("::[Mechanical Solver]:: -END- ")

    def AddVariables(self):
        
        structural_mechanics_solver.MechanicalSolver.AddVariables(self)
            
        if self.settings["rotation_dofs"].GetBool():
            # Add specific variables for the problem (rotation dofs)
            self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.StructuralMechanicsApplication.POINT_TORQUE)
   
        print("::[Mechanical Solver]:: Variables ADDED")

    
    def _GetSolutionScheme(self, scheme_type):

        if(scheme_type == "Newmark") or (scheme_type == "Bossak"):
            if(scheme_type == "Newmark"):
                self.settings.AddEmptyValue("damp_factor_m")  
                self.settings.AddEmptyValue("dynamic_factor")
                self.settings["damp_factor_m"].SetDouble(0.0)
                self.settings["dynamic_factor"].SetDouble(1.0)
                                                                            
            elif(scheme_type == "Bossak"):
                self.settings.AddEmptyValue("damp_factor_m")  
                self.settings.AddEmptyValue("dynamic_factor")
                self.settings["damp_factor_m"].SetDouble(-0.01)
                self.settings["dynamic_factor"].SetDouble(1.0)    
            
            # Creating the implicit solution scheme:  
            #self.main_model_part.ProcessInfo[KratosStructural.RAYLEIGH_ALPHA] = 0.0
            #self.main_model_part.ProcessInfo[KratosStructural.RAYLEIGH_BETA ] = 0.0
            
            mechanical_scheme = KratosMultiphysics.ResidualBasedBossakDisplacementScheme(self.settings["damp_factor_m"].GetDouble())

        elif(scheme_type == "Relaxation"):
          #~ self.main_model_part.GetSubModelPart(self.settings["volume_model_part_name"].GetString()).AddNodalSolutionStepVariable(DISPLACEMENT)  
            
            self.settings.AddEmptyValue("damp_factor_f")  
            self.settings.AddEmptyValue("dynamic_factor_m")
            self.settings["damp_factor_f"].SetDouble(-0.3)
            self.settings["dynamic_factor_m"].SetDouble(10.0) 
            
            mechanical_scheme = KratosMultiphysics.StructuralMechanicsApplication.ResidualBasedRelaxationScheme(self.settings["damp_factor_f"].GetDouble(),
                                                                                            self.settings["dynamic_factor_m"].GetDouble())
                                
        return mechanical_scheme

        
    def _CreateMechanicalSolver(self, mechanical_scheme, mechanical_convergence_criterion, builder_and_solver, max_iters, compute_reactions, reform_step_dofs, move_mesh_flag, line_search):
        if(line_search):
            self.mechanical_solver = KratosMultiphysics.LineSearchStrategy(self.computing_model_part, 
                                                                            mechanical_scheme, 
                                                                            self.linear_solver, 
                                                                            mechanical_convergence_criterion, 
                                                                            builder_and_solver, 
                                                                            max_iters, 
                                                                            compute_reactions, 
                                                                            reform_step_dofs, 
                                                                            move_mesh_flag)


        else:
            self.mechanical_solver = KratosMultiphysics.ResidualBasedNewtonRaphsonStrategy(self.computing_model_part, 
                                                                                            mechanical_scheme, 
                                                                                            self.linear_solver, 
                                                                                            mechanical_convergence_criterion, 
                                                                                            builder_and_solver, 
                                                                                            max_iters, 
                                                                                            compute_reactions, 
                                                                                            reform_step_dofs, 
                                                                                            move_mesh_flag)
