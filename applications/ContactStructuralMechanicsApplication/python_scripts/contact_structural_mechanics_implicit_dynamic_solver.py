from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
#import kratos core and applications
import KratosMultiphysics
import KratosMultiphysics.SolidMechanicsApplication as SolidMechanicsApplication
import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication
import KratosMultiphysics.ContactStructuralMechanicsApplication as ContactStructuralMechanicsApplication

# Check that KratosMultiphysics was imported in the main script
KratosMultiphysics.CheckForPreviousImport()

# Import the implicit solver (the explicit one is derived from it)
import structural_mechanics_implicit_dynamic_solver

def CreateSolver(main_model_part, custom_settings):
    return ImplicitMechanicalSolver(main_model_part, custom_settings)

class ImplicitMechanicalSolver(structural_mechanics_implicit_dynamic_solver.ImplicitMechanicalSolver):
    
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
            "solver_type": "contact_structural_mechanics_implicit_dynamic_solver",
            "model_import_settings": {
                "input_type": "mdpa",
                "input_filename": "unknown_name"
            },
            "echo_level": 0,
            "buffer_size": 2,
            "solution_type": "Dynamic",
            "rayleigh_alpha": 0.0,
            "rayleigh_beta":  0.0,
            "scheme_type": "Newmark",
            "damp_factor_m" : -0.1,
            "time_integration_method": "Implicit",
            "analysis_type": "Non-Linear",
            "rotation_dofs": false,
            "pressure_dofs": false,
            "damp_factor_m":-0.1,
            "rayleigh_alpha":0.0,
            "rayleigh_beta" :0.0,
            "stabilization_factor": 1.0,
            "reform_dofs_at_each_step": true,
            "line_search": false,
            "implex": false,
            "compute_reactions": true,
            "compute_contact_forces": false,
            "mortar_type": "",
            "block_builder": true,
            "clear_storage": false,
            "component_wise": false,
            "move_mesh_flag": true,
            "contact_tolerance": 0.0e0,
            "fancy_convergence_criterion": true,
            "convergence_criterion": "Residual_criteria",
            "displacement_relative_tolerance": 1.0e-4,
            "displacement_absolute_tolerance": 1.0e-9,
            "residual_relative_tolerance": 1.0e-4,
            "residual_absolute_tolerance": 1.0e-9,
            "max_iteration": 10,
            "ensure_contact": false,
            "split_factor": 10.0,
            "max_number_splits": 3,
            "rescale_factor": false,
            "path_following_penalty": false,
            "linear_solver_settings":{
                "solver_type": "SuperLUSolver",
                "max_iteration": 500,
                "tolerance": 1e-9,
                "scaling": false,
                "verbosity": 1
            },
            "processes_sub_model_part_list": [""],
            "problem_domain_sub_model_part_list": ["solid_model_part"]
        }
        """)
        
        ## Overwrite the default settings with user-provided parameters
        self.settings = custom_settings
        self.settings.ValidateAndAssignDefaults(default_settings)
        
        # Setting reactions true by default
        self.settings["compute_reactions"].SetBool(True)
        
        # Construct the linear solver
        import linear_solver_factory
        self.linear_solver = linear_solver_factory.ConstructSolver(self.settings["linear_solver_settings"])
        
        self.echo_level =  self.settings["echo_level"].GetInt()
        print(self.echo_level)
        
        print("Construction of MechanicalSolver finished")

    def AddVariables(self):
        
        structural_mechanics_implicit_dynamic_solver.ImplicitMechanicalSolver.AddVariables(self)
            
        if  self.settings["mortar_type"].GetString() != "":
            self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NORMAL)                                           # Add normal
            if  self.settings["mortar_type"].GetString() == "ALMContactFrictionless":
                self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NORMAL_CONTACT_STRESS)                        # Add normal contact stress
                self.main_model_part.AddNodalSolutionStepVariable(ContactStructuralMechanicsApplication.WEIGHTED_SCALAR_RESIDUAL)  # Add normal contact gap
                self.main_model_part.AddNodalSolutionStepVariable(ContactStructuralMechanicsApplication.WEIGHTED_GAP)              # Add normal contact gap
                self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NODAL_H)                                      # Add nodal size variable
            elif  self.settings["mortar_type"].GetString() == "ScalarMeshTying":
                self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.SCALAR_LAGRANGE_MULTIPLIER)                   # Add scalar LM
                self.main_model_part.AddNodalSolutionStepVariable(ContactStructuralMechanicsApplication.WEIGHTED_SCALAR_RESIDUAL)  # Add scalar LM residual
            elif  self.settings["mortar_type"].GetString() == "ComponentsMeshTying":
                self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VECTOR_LAGRANGE_MULTIPLIER)                   # Add vector LM
                self.main_model_part.AddNodalSolutionStepVariable(ContactStructuralMechanicsApplication.WEIGHTED_VECTOR_RESIDUAL)  # Add vector LM residual residual
   
        print("::[Mechanical Solver]:: Variables ADDED")
        
    def AddDofs(self):

        structural_mechanics_implicit_dynamic_solver.ImplicitMechanicalSolver.AddDofs(self)
        
        if  self.settings["mortar_type"].GetString() == "ALMContactFrictionless":
            for node in self.main_model_part.Nodes:
                node.AddDof(KratosMultiphysics.NORMAL_CONTACT_STRESS, ContactStructuralMechanicsApplication.WEIGHTED_GAP)
        elif  self.settings["mortar_type"].GetString() == "ScalarMeshTying":
            for node in self.main_model_part.Nodes:
                node.AddDof(KratosMultiphysics.SCALAR_LAGRANGE_MULTIPLIER,ContactStructuralMechanicsApplication.WEIGHTED_SCALAR_RESIDUAL)
        elif  self.settings["mortar_type"].GetString() == "ComponentsMeshTying":
            for node in self.main_model_part.Nodes:
                node.AddDof(KratosMultiphysics.VECTOR_LAGRANGE_MULTIPLIER_X, ContactStructuralMechanicsApplication.WEIGHTED_VECTOR_RESIDUAL_X)
                node.AddDof(KratosMultiphysics.VECTOR_LAGRANGE_MULTIPLIER_Y, ContactStructuralMechanicsApplication.WEIGHTED_VECTOR_RESIDUAL_Y)
                node.AddDof(KratosMultiphysics.VECTOR_LAGRANGE_MULTIPLIER_Z, ContactStructuralMechanicsApplication.WEIGHTED_VECTOR_RESIDUAL_Z)
                    
        print("::[Mechanical Solver]:: DOF's ADDED")
    
    def Initialize(self):
        structural_mechanics_implicit_dynamic_solver.ImplicitMechanicalSolver.Initialize(self)
    
    def _GetConvergenceCriterion(self):
        if "Contact" in self.settings["convergence_criterion"].GetString():
            D_RT = self.settings["displacement_relative_tolerance"].GetDouble()
            D_AT = self.settings["displacement_absolute_tolerance"].GetDouble()
            R_RT = self.settings["residual_relative_tolerance"].GetDouble()
            R_AT = self.settings["residual_absolute_tolerance"].GetDouble()
            contact_tolerance = self.settings["contact_tolerance"].GetDouble()
            fancy_convergence_criterion = self.settings["fancy_convergence_criterion"].GetBool()
            ensure_contact = self.settings["ensure_contact"].GetBool()
            echo_level = self.settings["echo_level"].GetInt()
            
            if (fancy_convergence_criterion == True):
                table = ContactStructuralMechanicsApplication.BprinterUtility()
            
            if(self.settings["convergence_criterion"].GetString() == "Contact_Displacement_criterion"):
                if (fancy_convergence_criterion == True):
                    convergence_criterion = ContactStructuralMechanicsApplication.DisplacementLagrangeMultiplierContactCriteria(D_RT, D_AT, D_RT, D_AT, ensure_contact, table)
                else:
                    convergence_criterion = ContactStructuralMechanicsApplication.DisplacementLagrangeMultiplierContactCriteria(D_RT, D_AT, D_RT, D_AT, ensure_contact)
                convergence_criterion.SetEchoLevel(echo_level)
                
            elif(self.settings["convergence_criterion"].GetString() == "Contact_Residual_criterion"):
                if (fancy_convergence_criterion == True):
                    convergence_criterion = ContactStructuralMechanicsApplication.DisplacementLagrangeMultiplierResidualContactCriteria(R_RT, R_AT, R_RT, R_AT, ensure_contact, table)
                else:
                    convergence_criterion = ContactStructuralMechanicsApplication.DisplacementLagrangeMultiplierResidualContactCriteria(R_RT, R_AT, R_RT, R_AT, ensure_contact)
                convergence_criterion.SetEchoLevel(echo_level)
                    
            elif(self.settings["convergence_criterion"].GetString() == "Contact_And_criterion"):
                if (fancy_convergence_criterion == True):
                    Displacement = ContactStructuralMechanicsApplication.DisplacementLagrangeMultiplierContactCriteria(D_RT, D_AT, D_RT, D_AT, ensure_contact, table)
                    Residual = ContactStructuralMechanicsApplication.DisplacementLagrangeMultiplierResidualContactCriteria(R_RT, R_AT, R_RT, R_AT, ensure_contact, table)
                else:
                    Displacement = ContactStructuralMechanicsApplication.DisplacementLagrangeMultiplierContactCriteria(D_RT, D_AT, D_RT, D_AT, ensure_contact)
                    Residual = ContactStructuralMechanicsApplication.DisplacementLagrangeMultiplierResidualContactCriteria(R_RT, R_AT, R_RT, R_AT, ensure_contact)
                
                Displacement.SetEchoLevel(echo_level)
                Residual.SetEchoLevel(echo_level)
                convergence_criterion = KratosMultiphysics.AndCriteria(Residual, Displacement)
                
            elif(self.settings["convergence_criterion"].GetString() == "Contact_Or_criterion"):
                if (fancy_convergence_criterion == True):
                    Displacement = ContactStructuralMechanicsApplication.DisplacementLagrangeMultiplierContactCriteria(D_RT, D_AT, D_RT, D_AT,ensure_contact, table)
                    Residual = ContactStructuralMechanicsApplication.DisplacementLagrangeMultiplierResidualContactCriteria(R_RT, R_AT, R_RT, R_AT, ensure_contact, table)
                else:
                    Displacement = ContactStructuralMechanicsApplication.DisplacementLagrangeMultiplierContactCriteria(D_RT, D_AT, D_RT, D_AT,ensure_contact)
                    Residual = ContactStructuralMechanicsApplication.DisplacementLagrangeMultiplierResidualContactCriteria(R_RT, R_AT, R_RT, R_AT, ensure_contact)
                
                Displacement.SetEchoLevel(echo_level)
                Residual.SetEchoLevel(echo_level)
                convergence_criterion = KratosMultiphysics.OrCriteria(Residual, Displacement)
            
            # Adding the mortar criteria
            if  (self.settings["mortar_type"].GetString() == "ALMContactFrictionless"):
                if (fancy_convergence_criterion == True):
                    Mortar = ContactStructuralMechanicsApplication.ALMFrictionlessMortarConvergenceCriteria(contact_tolerance, table)
                else:
                    Mortar = ContactStructuralMechanicsApplication.ALMFrictionlessMortarConvergenceCriteria(contact_tolerance)
            elif  (self.settings["mortar_type"].GetString() == "ALMContactFrictional"):
                if (fancy_convergence_criterion == True):
                    Mortar = ContactStructuralMechanicsApplication.ALMFrictionalMortarConvergenceCriteria(contact_tolerance, table)
                else:
                    Mortar = ContactStructuralMechanicsApplication.ALMFrictionalMortarConvergenceCriteria(contact_tolerance)
            elif ("MeshTying" in self.settings["mortar_type"].GetString()):
                if (fancy_convergence_criterion == True):
                    Mortar = ContactStructuralMechanicsApplication.MeshTyingMortarConvergenceCriteria(table)
                else:
                    Mortar = ContactStructuralMechanicsApplication.MeshTyingMortarConvergenceCriteria()
            
            Mortar.SetEchoLevel(self.echo_level)

            if (fancy_convergence_criterion == True):
                convergence_criterion = ContactStructuralMechanicsApplication.MortarAndConvergenceCriteria(convergence_criterion, Mortar, table)
            else:
                convergence_criterion = ContactStructuralMechanicsApplication.MortarAndConvergenceCriteria(convergence_criterion, Mortar)
            convergence_criterion.SetEchoLevel(self.echo_level)
            convergence_criterion.SetActualizeRHSFlag(True)
            
            return convergence_criterion
        
        else: # Standard criteria (same as solid and structural mechanics application)
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
        
            # Adding the mortar criteria
            if  (self.settings["mortar_type"].GetString() == "ALMContactFrictionless"):
                Mortar = ContactStructuralMechanicsApplication.ALMFrictionlessMortarConvergenceCriteria()
                Mortar.SetEchoLevel(self.echo_level)
                convergence_criterion.mechanical_convergence_criterion = KratosMultiphysics.AndCriteria( convergence_criterion.mechanical_convergence_criterion, Mortar)
                (convergence_criterion.mechanical_convergence_criterion).SetActualizeRHSFlag(True)
            elif  (self.settings["mortar_type"].GetString() == "ALMContactFrictional"):
                Mortar = ContactStructuralMechanicsApplication.ALMFrictionalMortarConvergenceCriteria()
                Mortar.SetEchoLevel(self.echo_level)
                convergence_criterion.mechanical_convergence_criterion = KratosMultiphysics.AndCriteria( convergence_criterion.mechanical_convergence_criterion, Mortar)
                (convergence_criterion.mechanical_convergence_criterion).SetActualizeRHSFlag(True)
            elif ("MeshTying" in self.settings["mortar_type"].GetString()):
                Mortar = ContactStructuralMechanicsApplication.MeshTyingMortarConvergenceCriteria()
                Mortar.SetEchoLevel(self.echo_level)
                convergence_criterion.mechanical_convergence_criterion = KratosMultiphysics.AndCriteria( convergence_criterion.mechanical_convergence_criterion, Mortar)
                (convergence_criterion.mechanical_convergence_criterion).SetActualizeRHSFlag(True)
            
            return convergence_criterion.mechanical_convergence_criterion
    
    def _CreateMechanicalSolver(self, mechanical_scheme, mechanical_convergence_criterion, builder_and_solver, max_iters, compute_reactions, reform_step_dofs, move_mesh_flag, component_wise, line_search, implex):
        
        if(component_wise):
            self.mechanical_solver = SolidMechanicsApplication.ComponentWiseNewtonRaphsonStrategy(
                                                                            self.computing_model_part, 
                                                                            mechanical_scheme, 
                                                                            self.linear_solver, 
                                                                            mechanical_convergence_criterion, 
                                                                            builder_and_solver, 
                                                                            max_iters, 
                                                                            compute_reactions, 
                                                                            reform_step_dofs, 
                                                                            move_mesh_flag
                                                                            )
        else:
            if(line_search):
                if(implex):
                    self.mechanical_solver = SolidMechanicsApplication.ResidualBasedNewtonRaphsonLineSearchImplexStrategy(self.computing_model_part, 
                                                                                                            mechanical_scheme, 
                                                                                                            self.linear_solver, 
                                                                                                            mechanical_convergence_criterion, 
                                                                                                            builder_and_solver, 
                                                                                                            max_iters, 
                                                                                                            compute_reactions, 
                                                                                                            reform_step_dofs, 
                                                                                                            move_mesh_flag
                                                                                                            )
                else:
                    newton_parameters = KratosMultiphysics.Parameters("""{}""")
                    newton_parameters.AddValue("rescale_factor",self.settings["rescale_factor"])
                    self.mechanical_solver = ContactStructuralMechanicsApplication.LineSearchContactStrategy(
                                                                                self.computing_model_part, 
                                                                                mechanical_scheme, 
                                                                                self.linear_solver, 
                                                                                mechanical_convergence_criterion, 
                                                                                builder_and_solver, 
                                                                                max_iters, 
                                                                                compute_reactions, 
                                                                                reform_step_dofs, 
                                                                                move_mesh_flag,
                                                                                newton_parameters
                                                                                )

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
                                                                            move_mesh_flag
                                                                            )
                    
                else:
                    if  self.settings["mortar_type"].GetString() != "":
                        newton_parameters = KratosMultiphysics.Parameters("""{}""")
                        newton_parameters.AddValue("split_factor",self.settings["split_factor"])
                        newton_parameters.AddValue("max_number_splits",self.settings["max_number_splits"])
                        newton_parameters.AddValue("rescale_factor",self.settings["rescale_factor"])
                        newton_parameters.AddValue("path_following_penalty",self.settings["path_following_penalty"])
                        self.mechanical_solver = ContactStructuralMechanicsApplication.ResidualBasedNewtonRaphsonContactStrategy(
                                                                                self.computing_model_part, 
                                                                                mechanical_scheme, 
                                                                                self.linear_solver, 
                                                                                mechanical_convergence_criterion, 
                                                                                builder_and_solver, 
                                                                                max_iters, 
                                                                                compute_reactions, 
                                                                                reform_step_dofs, 
                                                                                move_mesh_flag,
                                                                                newton_parameters
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
