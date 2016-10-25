from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
#import kratos core and applications
import KratosMultiphysics
import KratosMultiphysics.SolidMechanicsApplication
import KratosMultiphysics.StructuralMechanicsApplication
import KratosMultiphysics.ContactStructuralMechanicsApplication
import KratosMultiphysics.FSIApplication

# Check that KratosMultiphysics was imported in the main script
KratosMultiphysics.CheckForPreviousImport()

# Import the implicit solver (the explicit one is derived from it)
import structural_mechanics_static_solver

def CreateSolver(main_model_part, custom_settings):
    return StaticMechanicalSolver(main_model_part, custom_settings)

class StaticMechanicalSolver(structural_mechanics_static_solver.StaticMechanicalSolver):
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
            "solver_type": "contact_structural_mechanics_static_solver",
            "echo_level": 0,
            "buffer_size": 2,
            "solution_type": "Static",
            "analysis_type": "Non-Linear",
            "model_import_settings": {
                "input_type": "mdpa",
                "input_filename": "unknown_name"
            },
            "rotation_dofs": false,
            "pressure_dofs": false,
            "stabilization_factor": 1.0,
            "reform_dofs_at_each_step": false,
            "line_search": false,
            "compute_reactions": true,
            "compute_contact_forces": false,
            "compute_mortar_contact": true,
            "block_builder": false,
            "clear_storage": false,
            "component_wise": false,
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
            "accelerate_convergence": false,
            "convergence_accelerator":{
                "solver_type"       : "Relaxation",
                "acceleration_type" : "Aitken",
                "w_0"               :  0.35,
                "max_nl_iterations" :  10,
                "nl_tol"            :  1.0e-5,
                "update_system"     :  false,
                "buffer_size"       :  10
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
        
    def AddVariables(self):
        
        structural_mechanics_static_solver.StaticMechanicalSolver.AddVariables(self)
            
        if  self.settings["compute_mortar_contact"].GetBool():
            # Add normal
            self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NORMAL)
            # Add lagrange multiplier
            self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VECTOR_LAGRANGE_MULTIPLIER)
            # Add weighted gap
            self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.ContactStructuralMechanicsApplication.WEIGHTED_GAP)
            # Add weighted slip
            self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.ContactStructuralMechanicsApplication.WEIGHTED_SLIP)
            # Add weighted friction
            self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.ContactStructuralMechanicsApplication.WEIGHTED_FRICTION)
            # Add normal augmentation factor
            self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.ContactStructuralMechanicsApplication.NORMAL_AUGMENTATION_FACTOR)
            # Add tangent augmentation factor
            self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.ContactStructuralMechanicsApplication.TANGENT_AUGMENTATION_FACTOR)
            # Auxiliar active
            self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.ContactStructuralMechanicsApplication.AUXILIAR_ACTIVE)
            # Auxiliar slip
            self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.ContactStructuralMechanicsApplication.AUXILIAR_SLIP)
            # ACTIVE/INACTIVE
            self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.ContactStructuralMechanicsApplication.IS_ACTIVE_SET)
            # Active check factor
            self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.ContactStructuralMechanicsApplication.ACTIVE_CHECK_FACTOR)
        if self.settings["analysis_type"].GetString() == "Arc-Length":
            self.main_model_part.ProcessInfo[KratosMultiphysics.StructuralMechanicsApplication.LAMBDA] = 0.00;
   
        print("::[Mechanical Solver]:: Variables ADDED")
    
    def AddDofs(self):

        structural_mechanics_static_solver.StaticMechanicalSolver.AddDofs(self)
        
        if  self.settings["compute_mortar_contact"].GetBool():
            for node in self.main_model_part.Nodes:
                node.AddDof(KratosMultiphysics.VECTOR_LAGRANGE_MULTIPLIER_X);
                node.AddDof(KratosMultiphysics.VECTOR_LAGRANGE_MULTIPLIER_Y);
                node.AddDof(KratosMultiphysics.VECTOR_LAGRANGE_MULTIPLIER_Z);

        print("::[Mechanical Solver]:: DOF's ADDED")
    
    def _GetSolutionScheme(self, analysis_type, component_wise, compute_contact_forces):
        if(analysis_type == "Linear"):
            mechanical_scheme = KratosMultiphysics.ResidualBasedIncrementalUpdateStaticScheme()
            
        elif(analysis_type == "Arc-Length"):
            mechanical_scheme = KratosMultiphysics.ResidualBasedIncrementalUpdateStaticScheme()
            
        elif(analysis_type == "Non-Linear" ):
            self.settings.AddEmptyValue("damp_factor_m")  
            self.settings.AddEmptyValue("dynamic_factor")
            self.settings["damp_factor_m"].SetDouble(0.0)
            self.settings["dynamic_factor"].SetDouble(0.0) # Quasi-static scheme
            
            if component_wise:
                mechanical_scheme = KratosMultiphysics.SolidMechanicsApplication.ComponentWiseBossakScheme(
                                                              self.settings["damp_factor_m"].GetDouble(), 
                                                              self.settings["dynamic_factor"].GetDouble())
            else:
                if compute_contact_forces:
                    raise Exception("TODO: change for one that works with contact change")
                    #mechanical_scheme = ResidualBasedContactBossakScheme(self.settings["damp_factor_m"].GetDouble(), 
                                                                         #self.settings["dynamic_factor"].GetDouble())
                elif  self.settings["compute_mortar_contact"].GetBool():
                    mechanical_scheme = KratosMultiphysics.ContactStructuralMechanicsApplication.ResidualBasedIncrementalUpdateStaticContactScheme()
                else:
                    mechanical_scheme = KratosMultiphysics.ResidualBasedIncrementalUpdateStaticScheme()
                                
        return mechanical_scheme
    
    def _GetConvergenceCriterion(self):
        # Creation of an auxiliar Kratos parameters object to store the convergence settings
        conv_params = KratosMultiphysics.Parameters("{}")
        conv_params.AddValue("convergence_criterion",self.settings["convergence_criterion"])
        conv_params.AddValue("rotation_dofs",self.settings["rotation_dofs"])
        conv_params.AddValue("echo_level",self.settings["echo_level"])
        conv_params.AddValue("component_wise",self.settings["component_wise"])
        conv_params.AddValue("displacement_relative_tolerance",self.settings["displacement_relative_tolerance"])
        conv_params.AddValue("displacement_absolute_tolerance",self.settings["displacement_absolute_tolerance"])
        conv_params.AddValue("residual_relative_tolerance",self.settings["residual_relative_tolerance"])
        conv_params.AddValue("residual_absolute_tolerance",self.settings["residual_absolute_tolerance"])
        
        # Construction of the class convergence_criterion
        import convergence_criteria_factory
        convergence_criterion = convergence_criteria_factory.convergence_criterion(conv_params)
        
        if  self.settings["compute_mortar_contact"].GetBool():
            Mortar = KratosMultiphysics.ContactStructuralMechanicsApplication.MortarConvergenceCriteria()
            Mortar.SetEchoLevel(self.settings["echo_level"].GetInt())

            convergence_criterion.mechanical_convergence_criterion = KratosMultiphysics.AndCriteria(Mortar, convergence_criterion.mechanical_convergence_criterion)
        
        return convergence_criterion.mechanical_convergence_criterion
        
    def _CreateMechanicalSolver(self, mechanical_scheme, mechanical_convergence_criterion, builder_and_solver, max_iters, compute_reactions, reform_step_dofs, move_mesh_flag, component_wise, line_search):
        
        if(component_wise):
            self.mechanical_solver = KratosMultiphysics.SolidMechanicsApplication.ComponentWiseNewtonRaphsonStrategy(
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
            if(line_search):
                self.mechanical_solver = KratosMultiphysics.SolidMechanicsApplication.ResidualBasedNewtonRaphsonLineSearchStrategy(
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

                else:
                    if  self.settings["accelerate_convergence"].GetBool():
                        import convergence_accelerator_factory     
                        self.coupling_utility = convergence_accelerator_factory.CreateConvergenceAccelerator(self.settings["convergence_accelerator"])
                        max_nl_it = self.settings["convergence_accelerator"]["max_nl_iterations"].GetInt()
                        nl_tol = self.settings["convergence_accelerator"]["nl_tol"].GetDouble()
                        update_system = self.settings["convergence_accelerator"]["update_system"].GetBool()
                        self.mechanical_solver = KratosMultiphysics.ContactStructuralMechanicsApplication.ResidualBasedNewtonRaphsonContactAcceleratedStrategy(
                                                                                self.computing_model_part, 
                                                                                mechanical_scheme, 
                                                                                self.linear_solver, 
                                                                                mechanical_convergence_criterion, 
                                                                                builder_and_solver, 
                                                                                max_iters, 
                                                                                compute_reactions, 
                                                                                reform_step_dofs, 
                                                                                move_mesh_flag,
                                                                                self.coupling_utility,
                                                                                max_nl_it,
                                                                                nl_tol,
                                                                                update_system
                                                                                )
                    elif  self.settings["compute_mortar_contact"].GetBool():
                        split_factor   = self.settings["split_factor"].GetDouble()
                        max_number_splits = self.settings["max_number_splits"].GetInt()
                        self.mechanical_solver = KratosMultiphysics.ContactStructuralMechanicsApplication.ResidualBasedNewtonRaphsonContactStrategy(
                                                                                self.computing_model_part, 
                                                                                mechanical_scheme, 
                                                                                self.linear_solver, 
                                                                                mechanical_convergence_criterion, 
                                                                                builder_and_solver, 
                                                                                max_iters, 
                                                                                compute_reactions, 
                                                                                reform_step_dofs, 
                                                                                move_mesh_flag,
                                                                                split_factor,
                                                                                max_number_splits
                                                                                )
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
