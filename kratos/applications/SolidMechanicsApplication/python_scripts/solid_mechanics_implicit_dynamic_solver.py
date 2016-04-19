from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
#import kratos core and applications
import KratosMultiphysics
import KratosMultiphysics.SolidMechanicsApplication as SolidMechanicsApplication

# Check that KratosMultiphysics was imported in the main script
KratosMultiphysics.CheckForPreviousImport()

def CreateSolver(main_model_part, custom_settings):
    return ImplicitMechanicalSolver(main_model_part, custom_settings)

class ImplicitMechanicalSolver:
    
    
    ##constructor. the constructor shall only take care of storing the settings 
    ##and the pointer to the main_model part. This is needed since at the point of constructing the 
    ##model part is still not filled and the variables are not yet allocated
    ##
    ##real construction shall be delayed to the function "Initialize" which 
    ##will be called once the model is already filled
    def __init__(self, main_model_part, custom_settings): 
        
        #TODO: shall obtain the compute_model_part from the MODEL once the object is implemented
        self.main_model_part = main_model_part    
        
        ##settings string in json format
        default_settings = KratosMultiphysics.Parameters("""
        {
            "solver_type": "solid_mechanics_implicit_dynamic_solver",
            "model_import_settings": {
                "input_type": "mdpa",
                "input_filename": "unknown_name"
            },
            "echo_level": 0,
            "time_integration_method": "Implicit",
            "analysis_type": "nonlinear",
            "rotation_dofs": false,
            "pressure_dofs": false,
            "stabilization_factor": 1.0,
            "reform_dofs_at_each_iteration": false,
            "line_search": false,
            "compute_reactions": true,
            "compute_contact_forces": false,
            "block_builder": false,
            "component_wise": false,
            "move_mesh_flag": true,
            "solution_type": "Dynamic",
            "scheme_type": "Newmark",
            "convergence_criterion": "Residual_criteria",
            "displacement_relative_tolerance": 1.0e-4,
            "displacement_absolute_tolerance": 1.0e-9,
            "residual_relative_tolerance": 1.0e-4,
            "residual_absolute_tolerance": 1.0e-4,
            "max_iteration": 10,
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
        
    def GetMinimumBufferSize(self):
        return 2;

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

    def ImportModelPart(self):
        
        print("::[Mechanical Solver]:: Model reading starts.")
        if(self.settings["model_import_settings"]["input_type"].GetString() == "mdpa"):
            
            # Here it would be the place to import restart data if required
            KratosMultiphysics.ModelPartIO(self.settings["model_import_settings"]["input_filename"].GetString()).ReadModelPart(self.main_model_part)
            print("    Import input model part.")
            
            # Auxiliary Kratos parameters object to be called by the CheckAndPepareModelProcess
            aux_params = KratosMultiphysics.Parameters("{}")
            aux_params.AddValue("problem_domain_sub_model_part_list",self.settings["problem_domain_sub_model_part_list"])
            aux_params.AddValue("processes_sub_model_part_list",self.settings["processes_sub_model_part_list"])
            
            # CheckAndPrepareModelProcess creates the solid_computational_model_part
            import check_and_preparemodel_process
            check_and_preparemodel_process.CheckAndPrepareModelProcess(self.main_model_part, aux_params).Execute()
            
            # Constitutive law import
            import constitutive_law_python_utility as constitutive_law_utils
            constitutive_law = constitutive_law_utils.ConstitutiveLawUtility(self.main_model_part, 
                                                                             self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]);
            constitutive_law.Initialize();
            print("    Constitutive law initialized.")
            
        else:
            raise Exception("Other input options are not yet implemented.")
        
        current_buffer_size = self.main_model_part.GetBufferSize()
        if(self.GetMinimumBufferSize() > current_buffer_size):
            self.main_model_part.SetBufferSize( self.GetMinimumBufferSize() )
                
        print ("::[Mechanical Solver]:: Model reading finished.")


    def AddDofs(self):

        for node in self.main_model_part.Nodes:
            # adding dofs
            node.AddDof(KratosMultiphysics.DISPLACEMENT_X, KratosMultiphysics.REACTION_X);
            node.AddDof(KratosMultiphysics.DISPLACEMENT_Y, KratosMultiphysics.REACTION_Y);
            node.AddDof(KratosMultiphysics.DISPLACEMENT_Z, KratosMultiphysics.REACTION_Z);
            
        if self.settings["rotation_dofs"].GetBool():
            for node in model_part.Nodes:
                node.AddDof(KratosMultiphysics.ROTATION_X, KratosMultiphysics.TORQUE_X);
                node.AddDof(KratosMultiphysics.ROTATION_Y, KratosMultiphysics.TORQUE_Y);
                node.AddDof(KratosMultiphysics.ROTATION_Z, KratosMultiphysics.TORQUE_Z);
                
        if self.settings["rotation_dofs"].GetBool():                
            for node in model_part.Nodes:
                node.AddDof(KratosMultiphysics.PRESSURE, SolidMechanicsApplication.PRESSURE_REACTION);

        print("::[Mechanical Solver]:: DOF's ADDED")

    
    def Initialize(self):

        print("::[Mechanical Solver]:: -START-")
        
        # Get the solid_computational_model_part 
        self.compute_model_part = self.GetComputeModelPart()
        
        # Builder and solver creation
        builder_and_solver = self._GetBuilderAndSolver(self.settings["component_wise"].GetBool(), 
                                                       self.settings["block_builder"].GetBool())
        
        # Solution scheme creation
        mechanical_scheme = self._GetSolutionScheme(self.settings["scheme_type"].GetString(), 
                                                    self.settings["component_wise"].GetBool(),
                                                    self.settings["compute_contact_forces"].GetBool())
        
        # Get the convergence criterion
        mechanical_convergence_criterion = self._GetConvergenceCriterion()
        
        # Mechanical solver creation
        self._CreateMechanicalSolver(mechanical_scheme,
                                     mechanical_convergence_criterion,
                                     builder_and_solver,
                                     self.settings["max_iteration"].GetInt(),
                                     self.settings["compute_reactions"].GetBool(),
                                     self.settings["reform_dofs_at_each_iteration"].GetBool(),
                                     self.settings["move_mesh_flag"].GetBool(),
                                     self.settings["component_wise"].GetBool(),
                                     self.settings["line_search"].GetBool())

        # Set the stabilization factor
        self.main_model_part.ProcessInfo[KratosMultiphysics.STABILIZATION_FACTOR] = self.settings["stabilization_factor"].GetDouble()

        # Set echo_level
        self.mechanical_solver.SetEchoLevel(self.settings["echo_level"].GetInt())

        # Check if everything is assigned correctly
        self.Check();

        print("::[Mechanical Solver]:: -END- ")
        
    def GetComputeModelPart(self):
        return self.main_model_part.GetSubModelPart("solid_computational_model_part")
        
    def GetOutputVariables(self):
        pass
        
    def ComputeDeltaTime(self):
        pass
        
    def SaveRestart(self):
        pass #one should write the restart file here
        
    def Solve(self):
        self.mechanical_solver.Solve()

    def SetEchoLevel(self, level):
        self.solver.SetEchoLevel(level)

    def Clear(self):
        self.solver.Clear()
        
    def Check(self):
        self.mechanical_solver.Check()
        
    #### Specific internal functions ####
    
    def _GetBuilderAndSolver(self, component_wise, block_builder):
        # Creating the builder and solver
        if(component_wise):
            builder_and_solver = SolidMechanicsApplication.ComponentWiseBuilderAndSolver(self.linear_solver)
        else:
            if(block_builder):
                # To keep matrix blocks in builder
                builder_and_solver = KratosMultiphysics.ResidualBasedBlockBuilderAndSolver(self.linear_solver)
            else:
                builder_and_solver = KratosMultiphysics.ResidualBasedEliminationBuilderAndSolver(self.linear_solver)
        
        return builder_and_solver
        
    def _GetSolutionScheme(self, scheme_type, component_wise, compute_contact_forces):

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
        if (scheme_type == "Newmark" or scheme_type == "Bossak"):
            #~ self.main_model_part.ProcessInfo[SolidMechanicsApplication.RAYLEIGH_ALPHA] = 0.0
            #~ self.main_model_part.ProcessInfo[SolidMechanicsApplication.RAYLEIGH_BETA ] = 0.0
          
            if(component_wise):
                mechanical_scheme = SolidMechanicsApplication.ComponentWiseBossakScheme(self.settings["damp_factor_m"].GetDouble(), 
                                                                                        self.settings["dynamic_factor"].GetDouble())
            else:
                if(compute_contact_forces):
                    mechanical_scheme = SolidMechanicsApplication.ResidualBasedContactBossakScheme(self.settings["damp_factor_m"].GetDouble(),
                                                                                                   self.settings["dynamic_factor"].GetDouble())
                else:
                    mechanical_scheme = SolidMechanicsApplication.ResidualBasedBossakScheme(self.settings["damp_factor_m"].GetDouble(),
                                                                                            self.settings["dynamic_factor"].GetDouble())

        elif(scheme_type == "Relaxation"):
          #~ self.main_model_part.GetSubModelPart(self.settings["volume_model_part_name"].GetString()).AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)  
            
            import KratosMultiphysics.StructuralMechanicsApplication as StructuralMechanicsApplication
            self.settings.AddEmptyValue("damp_factor_f")  
            self.settings.AddEmptyValue("dynamic_factor_m")
            self.settings["damp_factor_f"].SetDouble(-0.3)
            self.settings["dynamic_factor_m"].SetDouble(10.0) 
            
            mechanical_scheme = StructuralMechanicsApplication.ResidualBasedRelaxationScheme(self.settings["damp_factor_f"].GetDouble(),
                                                                                             self.settings["dynamic_factor_m"].GetDouble())
                                
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
        import convergence_criteria_utility
        convergence_criterion = convergence_criteria_utility.convergence_criterion(conv_params)
            
        return convergence_criterion.mechanical_convergence_criterion
        
    def _CreateMechanicalSolver(self, mechanical_scheme, mechanical_convergence_criterion, builder_and_solver, max_iters, compute_reactions, reform_step_dofs, move_mesh_flag, component_wise, line_search):
        if(component_wise):
            self.mechanical_solver = SolidMechanicsApplication.ComponentWiseNewtonRaphsonStrategy(self.compute_model_part, 
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
                self.mechanical_solver = SolidMechanicsApplication.ResidualBasedNewtonRaphsonLineSearchStrategy(self.compute_model_part, 
                                                                                                                mechanical_scheme, 
                                                                                                                self.linear_solver, 
                                                                                                                mechanical_convergence_criterion, 
                                                                                                                builder_and_solver, 
                                                                                                                max_iters, 
                                                                                                                compute_reactions, 
                                                                                                                reform_step_dofs, 
                                                                                                                move_mesh_flag)

            else:
                self.mechanical_solver = KratosMultiphysics.ResidualBasedNewtonRaphsonStrategy(self.compute_model_part, 
                                                                                               mechanical_scheme, 
                                                                                               self.linear_solver, 
                                                                                               mechanical_convergence_criterion, 
                                                                                               builder_and_solver, 
                                                                                               max_iters, 
                                                                                               compute_reactions, 
                                                                                               reform_step_dofs, 
                                                                                               move_mesh_flag)
