from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
#import kratos core and applications
import KratosMultiphysics
import KratosMultiphysics.SolidMechanicsApplication as SolMechApp

# Check that KratosMultiphysics was imported in the main script
KratosMultiphysics.CheckForPreviousImport()

def CreateSolver(main_model_part, custom_settings):
    return MechanicalSolver(main_model_part, custom_settings)

class MechanicalSolver:
    
    
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
            "solver_type": "solid_mechanics_NEW_SOLVER",
            "model_import_settings": {
                "input_type": "mdpa",
                "input_filename": "unknown_name"
            },
            "echo_level": 0,
            "max_delta_time": 1,
            "fraction_delta_time": 0.9,
            "time_integration_method": "Implicit",
            "explicit_integration_scheme": "CentralDifferences",
            "time_step_prediction_level": "Automatic",
            "rayleigh_damping": false,
            "rotation_dofs": false,
            "pressure_dofs": false,
            "stabilization_factor": 1.0,
            "reform_dofs_at_each_iteration": false,
            "line_search": false,
            "implex": false,
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
            "skin_parts": [""],
            "volume_model_part_name": "volume_model_part"
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
        return 3;

    def AddVariables(self):
        
        #~ self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NODAL_MASS) #MSI, i included the variable becouse i calculate energy
        # add displacements
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        # add dynamic variables
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.ACCELERATION)
        # add reactions for the displacements
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION)
        # add nodal force variables
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.INTERNAL_FORCE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.EXTERNAL_FORCE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.CONTACT_FORCE)
        # add specific variables for the problem conditions
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.POSITIVE_FACE_PRESSURE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NEGATIVE_FACE_PRESSURE)
        self.main_model_part.AddNodalSolutionStepVariable(SolMechApp.POINT_LOAD)
        self.main_model_part.AddNodalSolutionStepVariable(SolMechApp.LINE_LOAD)
        self.main_model_part.AddNodalSolutionStepVariable(SolMechApp.SURFACE_LOAD)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VOLUME_ACCELERATION)
            
        if self.settings["rotation_dofs"].GetBool():
            # add specific variables for the problem (rotation dofs)
            self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.ROTATION)
            self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.TORQUE)
            self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.ANGULAR_VELOCITY)
            self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.ANGULAR_ACCELERATION)
        if self.settings["pressure_dofs"].GetBool():
            # add specific variables for the problem (pressure dofs)
            self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.PRESSURE)
            self.main_model_part.AddNodalSolutionStepVariable(SolMechApp.PRESSURE_REACTION)
                    
        print("::[Mechanical Solver]:: Variables ADDED")

    def ImportModelPart(self):
        
        print("::[Mechanical Solver]:: Model reading starts.")
        if(self.settings["model_import_settings"]["input_type"].GetString() == "mdpa"):
            
            # Here it would be the place to import restart data if required
            KratosMultiphysics.ModelPartIO(self.settings["model_import_settings"]["input_filename"].GetString()).ReadModelPart(self.main_model_part)
            print("    Import input model part.")
            
            # Here we shall check that the input read has the shape we like
            aux_params = KratosMultiphysics.Parameters("{}")
            aux_params.AddValue("volume_model_part_name",self.settings["volume_model_part_name"])
            aux_params.AddValue("skin_parts",self.settings["skin_parts"])
            
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
                node.AddDof(KratosMultiphysics.PRESSURE, SolMechApp.PRESSURE_REACTION);

        print("::[Mechanical Solver]:: DOF's ADDED")

    
    def Initialize(self):

        print("::[Mechanical Solver]:: -START-")
        
        # Builder and solver creation
        builder_and_solver = self._GetBuilderAndSolver(self.settings["component_wise"].GetBool(), 
                                                       self.settings["block_builder"].GetBool())
        
        # Solution scheme creation
        mechanical_scheme = self._GetSolutionSchemeImplicit(self.settings["scheme_type"].GetString(), 
                                                            self.settings["component_wise"].GetBool(),
                                                            self.settings["compute_contact_forces"].GetBool())

        # Convergence criterion creation
        mechanical_convergence_criterion = self._GetConvergenceCriterion(self.settings["convergence_criterion"].GetString(),
                                                                         self.settings["rotation_dofs"].GetBool(),
                                                                         self.settings["echo_level"].GetInt(),
                                                                         self.settings["component_wise"].GetBool())  
        
        # Mechanical solver creation
        self._CreateMechanicalSolver(mechanical_scheme,
                                     mechanical_convergence_criterion,
                                     builder_and_solver,
                                     self.settings["max_iteration"].GetInt(),
                                     self.settings["compute_reactions"].GetBool(),
                                     self.settings["reform_dofs_at_each_iteration"].GetBool(),
                                     self.settings["move_mesh_flag"].GetBool(),
                                     self.settings["component_wise"].GetBool(),
                                     self.settings["line_search"].GetBool(),
                                     self.settings["time_integration_method"].GetString())

        # Set the stabilization factor
        self.main_model_part.ProcessInfo[KratosMultiphysics.STABILIZATION_FACTOR] = self.settings["stabilization_factor"].GetDouble()

        # Set echo_level
        self.mechanical_solver.SetEchoLevel(self.settings["echo_level"].GetInt())

        # Check if everything is assigned correctly
        self.Check();

        print("::[Mechanical Solver]:: -END- ")
        
    def GetComputeModelPart(self):
        return self.main_model_part.GetSubModelPart(self.settings["volume_model_part_name"].GetString())
        
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
        # creating the builder and solver
        if(component_wise):
            builder_and_solver = SolMechApp.ComponentWiseBuilderAndSolver(self.linear_solver)
        else:
            if(block_builder):
                # to keep matrix blocks in builder
                builder_and_solver = KratosMultiphysics.ResidualBasedBlockBuilderAndSolver(self.linear_solver)
            else:
                builder_and_solver = KratosMultiphysics.ResidualBasedEliminationBuilderAndSolver(self.linear_solver)
        
        return builder_and_solver
        
    def _GetSolutionSchemeImplicit(self, scheme_type, component_wise, compute_contact_forces):
        # creating the implicit solution scheme:
        if(scheme_type == "Newmark"):
          damp_factor_m = 0.0;
          dynamic_factor = 1;
          
          self.main_model_part.ProcessInfo[SolMechApp.RAYLEIGH_ALPHA] = 0.0
          self.main_model_part.ProcessInfo[SolMechApp.RAYLEIGH_BETA ] = 0.0
          
          if(component_wise):
              mechanical_scheme = SolMechApp.ComponentWiseBossakScheme(damp_factor_m, 
                                                                       dynamic_factor)
          else:
              if(compute_contact_forces):
                  mechanical_scheme = SolMechApp.ResidualBasedContactBossakScheme(damp_factor_m,
                                                                                  dynamic_factor)
              else:
                  mechanical_scheme = SolMechApp.ResidualBasedBossakScheme(damp_factor_m,
                                                                           dynamic_factor)
        elif(scheme_type == "Bossak"):
          damp_factor_m = -0.01;
          dynamic_factor = 1;
          
          self.main_model_part.ProcessInfo[SolMechApp.RAYLEIGH_ALPHA] = 0.0
          self.main_model_part.ProcessInfo[SolMechApp.RAYLEIGH_BETA ] = 0.0
          
          if(component_wise):
              mechanical_scheme = SolMechApp.ComponentWiseBossakScheme(damp_factor_m, 
                                                                       dynamic_factor)
          else:
              if(compute_contact_forces):
                  mechanical_scheme = SolMechApp.ResidualBasedContactBossakScheme(damp_factor_m,
                                                                                  dynamic_factor)
              else:
                  mechanical_scheme = SolMechApp.ResidualBasedBossakScheme(damp_factor_m,
                                                                           dynamic_factor)
        elif(scheme_type == "Relaxation"):
          #~ self.main_model_part.GetSubModelPart(self.settings["volume_model_part_name"].GetString()).AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)  
            
          import KratosMultiphysics.StructuralMechanicsApplication as StructMechApp
          damp_factor_f = -0.3;
          damp_factor_m = 10.0;
          mechanical_scheme = StructMechApp.ResidualBasedRelaxationScheme(damp_factor_f,
                                                                          damp_factor_m)
                                
        return mechanical_scheme
    
    def _GetConvergenceCriterion(self, convergence_criterion_type, rotation_dofs, echo_level, component_wise):

        D_RT = self.settings["displacement_relative_tolerance"].GetDouble()
        D_AT = self.settings["displacement_absolute_tolerance"].GetDouble()
        R_RT = self.settings["residual_relative_tolerance"].GetDouble()
        R_AT = self.settings["residual_absolute_tolerance"].GetDouble()

        if(rotation_dofs):
            if(convergence_criterion_type == "Displacement_criteria"):
                mechanical_convergence_criterion = SolMechApp.DisplacementCriteria(D_RT, D_AT)
            elif(convergence_criterion_type == "Residual_criteria"):
                mechanical_convergence_criterion = KratosMultiphysics.ResidualCriteria(R_RT, R_AT)
            elif(convergence_criterion_type == "And_criteria"):
                Displacement = SolMechApp.DisplacementCriteria(D_RT, D_AT)
                Residual = KratosMultiphysics.ResidualCriteria(R_RT, R_AT)
                mechanical_convergence_criterion = KratosMultiphysics.AndCriteria(Residual, Displacement)
            elif(convergence_criterion_type == "Or_criteria"):
                Displacement = SolMechApp.DisplacementCriteria(D_RT, D_AT)
                Residual = KratosMultiphysics.ResidualCriteria(R_RT, R_AT)
                mechanical_convergence_criterion = KratosMultiphysics.OrCriteria(Residual, Displacement)
            #~ elif(convergence_criterion_type == "Mixed_criteria"):
                #~ Displacement = KratosMultiphysics.MixedElementCriteria(D_RT, D_AT)
                #~ Residual = KratosMultiphysics.ResidualCriteria(R_RT, R_AT)
                #~ mechanical_convergence_criterion = KratosMultiphysics.AndCriteria(Residual, Displacement)
        else:
            if(echo_level > 1):
                print("::[Mechanical Solver]:: CONVERGENCE CRITERION : ", convergence_criterion_type)

            if(convergence_criterion_type == "Displacement_criteria"):
                mechanical_convergence_criterion = SolMechApp.DisplacementConvergenceCriterion(D_RT, D_AT)
            elif(convergence_criterion_type == "Residual_criteria"):
                if(component_wise):
                    mechanical_convergence_criterion = SolMechApp.ComponentWiseResidualConvergenceCriterion(R_RT, R_AT)
                else:
                    mechanical_convergence_criterion = KratosMultiphysics.ResidualCriteria(R_RT, R_AT)
            elif(convergence_criterion_type == "And_criteria"):
                Displacement = SolMechApp.DisplacementConvergenceCriterion(D_RT, D_AT)
                if(component_wise):
                    Residual = SolMechApp.ComponentWiseResidualConvergenceCriterion(R_RT, R_AT)
                else:
                    Residual = KratosMultiphysics.ResidualCriteria(R_RT, R_AT)
                mechanical_convergence_criterion = KratosMultiphysics.AndCriteria(Residual, Displacement)
            elif(convergence_criterion_type == "Or_criteria"):
                Displacement = SolMechApp.DisplacementConvergenceCriterion(D_RT, D_AT)
                if(self.component_wise):
                    Residual = SolMechApp.ComponentWiseResidualConvergenceCriterion(R_RT, R_AT)
                else:
                    Residual = KratosMultiphysics.ResidualCriteria(R_RT, R_AT)
                mechanical_convergence_criterion = KratosMultiphysics.OrCriteria(Residual, Displacement)
            #~ elif(convergence_criterion_type == "Mixed_criteria"):
                #~ Displacement = KratosMultiphysics.MixedElementConvergeCriteria(D_RT, D_AT)
                #~ Residual = KratosMultiphysics.ResidualCriteria(R_RT, R_AT)
                #~ mechanical_convergence_criterion = KratosMultiphysics.AndCriteria(Residual, Displacement)
            
        return mechanical_convergence_criterion
        
    def _CreateMechanicalSolver(self, mechanical_scheme, mechanical_convergence_criterion, builder_and_solver, max_iters, compute_reactions, reform_step_dofs, move_mesh_flag, component_wise, line_search, time_integration_method):
        if (time_integration_method == "Implicit"):
            if(component_wise):
                self.mechanical_solver = SolMechApp.ComponentWiseNewtonRaphsonStrategy(self.main_model_part, 
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
                    self.mechanical_solver = SolMechApp.ResidualBasedNewtonRaphsonLineSearchStrategy(self.main_model_part, 
                                                                                                    mechanical_scheme, 
                                                                                                    self.linear_solver, 
                                                                                                    mechanical_convergence_criterion, 
                                                                                                    builder_and_solver, 
                                                                                                    max_iters, 
                                                                                                    compute_reactions, 
                                                                                                    reform_step_dofs, 
                                                                                                    move_mesh_flag)
                else:
                    self.mechanical_solver = KratosMultiphysics.ResidualBasedNewtonRaphsonStrategy(self.main_model_part, 
                                                                                                   mechanical_scheme, 
                                                                                                   self.linear_solver, 
                                                                                                   mechanical_convergence_criterion, 
                                                                                                   builder_and_solver, 
                                                                                                   max_iters, 
                                                                                                   compute_reactions, 
                                                                                                   reform_step_dofs, 
                                                                                                   move_mesh_flag)
