from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
#import kratos core and applications
import KratosMultiphysics
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication

# Check that KratosMultiphysics was imported in the main script
KratosMultiphysics.CheckForPreviousImport()

# Import the implicit solver (the explicit one is derived from it)
import structural_mechanics_solver

def CreateSolver(main_model_part, custom_settings):
    return StaticMechanicalSolver(main_model_part, custom_settings)

class StaticMechanicalSolver(structural_mechanics_solver.MechanicalSolver):
    
    
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
            "solver_type": "structural_mechanics_static_solver",
            "echo_level": 0,
            "buffer_size": 2,
            "solution_type": "Static",
            "analysis_type": "Non-Linear",
            "model_import_settings": {
                "input_type": "mdpa",
                "input_filename": "unknown_name"
            },
            "material_import_settings" :{
                "materials_filename": ""
            },
            "rotation_dofs": false,
            "pressure_dofs": false,
            "stabilization_factor": 1.0,
            "reform_dofs_at_each_step": false,
            "line_search": false,
            "compute_reactions": true,
            "compute_contact_forces": false,
            "block_builder": true,
            "clear_storage": false,
            "move_mesh_flag": true,
            "convergence_criterion": "Residual_criteria",
            "displacement_relative_tolerance": 1.0e-4,
            "displacement_absolute_tolerance": 1.0e-9,
            "residual_relative_tolerance": 1.0e-4,
            "residual_absolute_tolerance": 1.0e-4,
            "max_iteration": 10,
            "linear_solver_settings":{
                "solver_type": "SuperLUSolver",
                "max_iteration": 500,
                "tolerance": 1e-9,
                "scaling": false,
                "verbosity": 1
            },
            "arc_length_settings": {
                "Ide": 5,
                "factor_delta_lmax": 1.00,
                "max_iteration": 20,
                "max_recursive": 50,
                "toler": 1.0E-10,
                "norm": 1.0E-7,
                "MaxLineSearchIterations": 20,
                "tolls": 0.000001,
                "amp": 1.618,
                "etmxa": 5,
                "etmna": 0.1
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
        
        echo_level =  self.settings["echo_level"].GetInt()
        print(self.settings["echo_level"].GetInt())
        
        print("Construction of MechanicalSolver finished")
        
    def Initialize(self):

        print("::[Mechanical Solver]:: -START-")
        
        # Get the solid computing model part
        self.computing_model_part = self.GetComputingModelPart()
                
        # Solution scheme choice
        mechanical_scheme = self._GetSolutionScheme(self.settings["analysis_type"].GetString())
        
        # Get the convergence choice
        mechanical_convergence_criterion = self._GetConvergenceCriterion()
        
        # Builder and solver choice
        builder_and_solver = self._GetBuilderAndSolver(self.settings["block_builder"].GetBool())

        # Mechanical solver choice
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
        self.Check()

        print("::[Mechanical Solver]:: -END- ")
        
    def AddVariables(self):
        
        structural_mechanics_solver.MechanicalSolver.AddVariables(self)

        if self.settings["rotation_dofs"].GetBool():
            # Add specific variables for the problem (rotation dofs)
            self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.StructuralMechanicsApplication.POINT_TORQUE)
        
        if self.settings["analysis_type"].GetString() == "Arc-Length":
            self.main_model_part.ProcessInfo[KratosMultiphysics.StructuralMechanicsApplication.LAMBDA] = 0.00;
   
        print("::[Mechanical Solver]:: Variables ADDED")
    
    def Solve(self):
        KratosMultiphysics.VariableUtils().CheckVariableKeys()
        if self.settings["clear_storage"].GetBool():
            self.Clear()
        
        self.mechanical_solver.Solve()
        
        if self.settings["analysis_type"].GetString() == "Arc-Length":
            echo_level =  self.settings["echo_level"].GetInt()
            lambda_value = self.main_model_part.ProcessInfo[KratosMultiphysics.StructuralMechanicsApplication.LAMBDA]
            if echo_level > 0:
                print("LAMBDA: ", lambda_value)
            self._ChangeCondition(lambda_value)
    
    def _IncreasePointLoad(forcing_nodes_list, Load):
        for node in forcing_nodes_list:
            node.SetSolutionStepValue(KratosMultiphysics.StructuralMechanicsApplication.POINT_LOAD, 0, Load)

    def _IncreaseDisplacement(forcing_nodes_list, disp):
        for node in forcing_nodes_list:
            node.SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT, 0, disp)

    def _ChangeCondition(self, lambda_value):
        echo_level =  self.settings["echo_level"].GetInt()
        force_x = 0.0
        force_y = 0.0
        force_z = 0.0
        for node in self.main_model_part.Nodes:
            new_load = node.GetSolutionStepValue(KratosMultiphysics.StructuralMechanicsApplication.POINT_LOAD, 0) * lambda_value;
            force_x += new_load[0]
            force_y += new_load[1]
            force_z += new_load[2]
            node.SetSolutionStepValue(KratosMultiphysics.StructuralMechanicsApplication.POINT_LOAD, 0, new_load)
        
        if (echo_level > 0):
            print("*********************** ")
            print("The total load applied: ")
            print("POINT_LOAD_X: ", force_x)
            print("POINT_LOAD_Y: ", force_y)
            print("POINT_LOAD_Z: ", force_z)
            print("*********************** ")
    
    def _GetSolutionScheme(self, analysis_type):
        mechanical_scheme = KratosMultiphysics.ResidualBasedIncrementalUpdateStaticScheme()
                                
        return mechanical_scheme
        
    def _CreateMechanicalSolver(self, mechanical_scheme, mechanical_convergence_criterion, builder_and_solver, max_iters, compute_reactions, reform_step_dofs, move_mesh_flag, line_search):
        
        if(line_search):
            self.mechanical_solver = KratosMultiphysics.LineSearchStrategy(
                                                        self.computing_model_part, 
                                                        mechanical_scheme, 
                                                        self.linear_solver, 
                                                        mechanical_convergence_criterion, 
                                                        builder_and_solver, 
                                                        max_iters, 
                                                        compute_reactions, 
                                                        reform_step_dofs, 
                                                        move_mesh_flag)

        else:
            if self.settings["analysis_type"].GetString() == "Linear":
                self.mechanical_solver = KratosMultiphysics.ResidualBasedLinearStrategy(
                                                                        self.computing_model_part, 
                                                                        mechanical_scheme, 
                                                                        self.linear_solver, 
                                                                        builder_and_solver, 
                                                                        compute_reactions, 
                                                                        reform_step_dofs, 
                                                                        False, 
                                                                        move_mesh_flag)
                
            elif self.settings["analysis_type"].GetString() == "Arc-Length":
                Ide = self.settings["arc_length_settings"]["Ide"].GetInt()
                max_iteration = self.settings["arc_length_settings"]["max_iteration"].GetInt()
                max_recursive = self.settings["arc_length_settings"]["max_recursive"].GetInt()
                factor_delta_lmax = self.settings["arc_length_settings"]["factor_delta_lmax"].GetDouble()
                self.mechanical_solver = KratosMultiphysics.StructuralMechanicsApplication.ResidualBasedArcLengthStrategy(
                                                                        self.computing_model_part, 
                                                                        mechanical_scheme, 
                                                                        self.linear_solver, 
                                                                        mechanical_convergence_criterion, 
                                                                        Ide,
                                                                        max_iteration,
                                                                        max_recursive,
                                                                        factor_delta_lmax,
                                                                        compute_reactions, 
                                                                        reform_step_dofs, 
                                                                        move_mesh_flag)

            elif self.settings["analysis_type"].GetString() == "Formfinding":
                print("::[Mechanical Solver]::StructuralMechanicsApplication::Formfinding ")
                self.mechanical_solver = StructuralMechanicsApplication.FormfindingUpdatedReferenceStrategy(
                                                                        self.computing_model_part, 
                                                                        mechanical_scheme, 
                                                                        self.linear_solver, 
                                                                        mechanical_convergence_criterion, 
                                                                        builder_and_solver, 
                                                                        max_iters, 
                                                                        compute_reactions, 
                                                                        reform_step_dofs, 
                                                                        move_mesh_flag)

            else:
                self.mechanical_solver = KratosMultiphysics.ResidualBasedNewtonRaphsonStrategy(
                                                                        self.computing_model_part, 
                                                                        mechanical_scheme, 
                                                                        self.linear_solver, 
                                                                        mechanical_convergence_criterion, 
                                                                        builder_and_solver, 
                                                                        max_iters, 
                                                                        compute_reactions, 
                                                                        reform_step_dofs, 
                                                                        move_mesh_flag)
