from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
#import kratos core and applications
import KratosMultiphysics
import KratosMultiphysics.SolidMechanicsApplication as SolidMechanicsApplication
#~ import KratosMultiphysics.ExternalSolversApplication as ExtSolvApp

#~ from KratosMultiphysics import *
#~ from KratosMultiphysics.SolidMechanicsApplication import *
#~ from KratosMultiphysics.ExternalSolversApplication import *

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
            "solution_type": "Static",
            "scheme_type": "Linear",
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
        self.main_model_part.AddNodalSolutionStepVariable(SolidMechanicsApplication.IMPOSED_DISPLACEMENT)       # Required variable which was not included
        self.main_model_part.AddNodalSolutionStepVariable(SolidMechanicsApplication.IMPOSED_ROTATION)           # Required variable which was not included
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.POSITIVE_FACE_PRESSURE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NEGATIVE_FACE_PRESSURE)
        self.main_model_part.AddNodalSolutionStepVariable(SolidMechanicsApplication.POINT_LOAD)
        self.main_model_part.AddNodalSolutionStepVariable(SolidMechanicsApplication.LINE_LOAD)
        self.main_model_part.AddNodalSolutionStepVariable(SolidMechanicsApplication.SURFACE_LOAD)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VOLUME_ACCELERATION)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DENSITY)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.YOUNG_MODULUS)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.POISSON_RATIO)
            
        if self.settings["rotation_dofs"].GetBool():
            # add specific variables for the problem (rotation dofs)
            self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.ROTATION)
            self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.TORQUE)
            self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.ANGULAR_VELOCITY)
            self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.ANGULAR_ACCELERATION)
        if self.settings["pressure_dofs"].GetBool():
            # add specific variables for the problem (pressure dofs)
            self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.PRESSURE)
            self.main_model_part.AddNodalSolutionStepVariable(SolidMechanicsApplication.PRESSURE_REACTION)
        if self.settings["time_integration_method"].GetString() == "Explicit":
            self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NODAL_MASS)
            self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.FORCE_RESIDUAL)
            self.main_model_part.AddNodalSolutionStepVariable(SolidMechanicsApplication.MIDDLE_VELOCITY)
                    
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
            
            ######### TO BE REVIEWED AS SOON AS THE REPLACE SETTINGS ARE STATED ###########
            #~ #here we replace the dummy elements we read with proper elements
            #~ self.settings.AddEmptyValue("element_replace_settings")
            #~ if(self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] == 3):
                #~ self.settings["element_replace_settings"] = KratosMultiphysics.Parameters("""
                    #~ {
                    #~ "element_name":"VMS3D4N",
                    #~ "condition_name": "MonolithicWallCondition3D"
                    #~ }
                    #~ """)
            #~ elif(self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] == 2):
                #~ self.settings["element_replace_settings"] = KratosMultiphysics.Parameters("""
                    #~ {
                    #~ "element_name":"VMS2D3N",
                    #~ "condition_name": "MonolithicWallCondition2D"
                    #~ }
                    #~ """)
            #~ else:
                #~ raise Exception("domain size is not 2 or 3")
            #~ 
            #~ KratosMultiphysics.ReplaceElementsAndConditionsProcess(self.main_model_part, self.settings["element_replace_settings"]).Execute()
            #~ 
            #~ import check_and_preparemodel_process
            #~ check_and_preparemodel_process.CheckAndPrepareModelProcess(self.main_model_part, aux_params).Execute()
            
            ###### TODO: This manner does not allow to set different materials (unique material model)...
            # Set density, Young modulus and Poisson ratio to the nodes
            for el in self.main_model_part.Elements:
                density = el.Properties.GetValue(KratosMultiphysics.DENSITY)
                young_modulus = el.Properties.GetValue(KratosMultiphysics.YOUNG_MODULUS)
                poisson_ratio = el.Properties.GetValue(KratosMultiphysics.POISSON_RATIO)
                break
            
            KratosMultiphysics.VariableUtils().SetScalarVar(KratosMultiphysics.DENSITY, density, self.main_model_part.Nodes)
            KratosMultiphysics.VariableUtils().SetScalarVar(KratosMultiphysics.YOUNG_MODULUS, young_modulus, self.main_model_part.Nodes)
            KratosMultiphysics.VariableUtils().SetScalarVar(KratosMultiphysics.POISSON_RATIO, poisson_ratio, self.main_model_part.Nodes)
            print("    Density, Young modulus and Poisson ratio set.")
            
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
                
        # Default settings
        echo_level = self.settings["echo_level"].GetInt()
        domain_size = self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]
        buffer_size = 3 #default buffer_size

        pressure_dofs = self.settings["pressure_dofs"].GetBool()
        rotation_dofs = self.settings["rotation_dofs"].GetBool()
        stabilization_factor = self.settings["stabilization_factor"].GetDouble()

        # Definition of the solvers
        solution_type = self.settings["solution_type"].GetString()
        scheme_type = self.settings["scheme_type"].GetString()
        time_integration_method = self.settings["time_integration_method"].GetString()
        explicit_integration_method = self.settings["explicit_integration_scheme"].GetString()
        
        arlequin = 0 #?¿¿?¿¿??¿?¿?¿¿?

        # Definition of the convergence criteria
        rel_disp_tol = self.settings["displacement_relative_tolerance"].GetDouble()
        abs_disp_tol = self.settings["displacement_absolute_tolerance"].GetDouble()
        rel_res_tol = self.settings["residual_relative_tolerance"].GetDouble()
        abs_res_tol = self.settings["residual_absolute_tolerance"].GetDouble()
        max_iters = self.settings["max_iteration"].GetInt()
        convergence_criterion_type = self.settings["convergence_criterion"].GetString()

        # Definition of the default builder_and_solver:
        block_builder = self.settings["block_builder"].GetBool()
        builder_and_solver = SolidMechanicsApplication.ResidualBasedBuilderAndSolver(self.linear_solver)

        # Definition of the component wise calculation "computation is slower"
        #(it affects to strategy, builder_and_solver, scheme and convergence_criterion)
        component_wise = self.settings["component_wise"].GetBool()

        # Definition of computing flags
        compute_reactions = self.settings["compute_reactions"].GetBool()
        compute_contact_forces = self.settings["compute_contact_forces"].GetBool()
        line_search = self.settings["line_search"].GetBool()
        implex = self.settings["implex"].GetBool()
        reform_step_dofs = self.settings["reform_dofs_at_each_iteration"].GetBool()

        # Definition of the (x_n+1 = x_n+dx) in each iteration -> strategy base class includes de MoveMesh method
        move_mesh_flag = self.settings["move_mesh_flag"].GetBool()
        
        max_delta_time = self.settings["max_delta_time"].GetDouble()
        fraction_delta_time = self.settings["fraction_delta_time"].GetDouble()
        time_step_prediction_level = self.settings["time_step_prediction_level"].GetString()
        #~ self.time_step_prediction_level = 0 ## IN THE ORIGINAL SOLVER IT WAS AN INTEGER...
        rayleigh_damping = self.settings["rayleigh_damping"].GetBool()

        print("::[Mechanical Solver]:: -START-")
        
        ##### ?¿?¿??¿¿?¿ ######
        if(time_integration_method == "Explicit"):
            if(explicit_integration_scheme == "CentralDifferences"):
                buffer_size = 2 
                #could be 1 becouse it uses a extra variable - MIDDLE_VELOCITY for previous time step velocity and avoids extra buffer_size. However
                #the prediction/application of B.C. need last step values 
            if(arlequin == 1): #### ARLEQUIN IS NEVER = 1, IS PREVIOUSLY SET TO 0.
                buffer_size = 2
        ##### ¿?¿?¿?¿?¿¿ ######

        # Builder and solver creation
        builder_and_solver = self.GetBuilderAndSolver(component_wise, block_builder)
        
        # Solution scheme creation
        if(time_integration_method == "Explicit"):
            mechanical_scheme = self.GetSolutionSchemeExplicit(explicit_integration_scheme, 
                                                               max_delta_time, 
                                                               fraction_delta_time, 
                                                               time_step_prediction_level, 
                                                               rayleigh_damping)
        elif(time_integration_method == "Implicit"):
            mechanical_scheme = self.GetSolutionSchemeImplicit(scheme_type, 
                                                               component_wise,
                                                               compute_contact_forces)

        # Convergence criterion creation
        mechanical_convergence_criterion = self.GetConvergenceCriterion(convergence_criterion_type,
                                                                        rotation_dofs,
                                                                        echo_level,
                                                                        component_wise)  
        # Mechanical solver creation
        self.CreateMechanicalSolver(mechanical_scheme,
                                    mechanical_convergence_criterion,
                                    builder_and_solver,
                                    max_iters,
                                    compute_reactions,
                                    reform_step_dofs,
                                    move_mesh_flag,
                                    component_wise,
                                    line_search,
                                    time_integration_method)

        # Set the stabilization factor
        self.main_model_part.ProcessInfo[KratosMultiphysics.STABILIZATION_FACTOR] = stabilization_factor

        # Set echo_level
        self.mechanical_solver.SetEchoLevel(echo_level)

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
    
    def GetBuilderAndSolver(self, component_wise, block_builder):
        # creating the builder and solver
        if(component_wise):
            builder_and_solver = SolidMechanicsApplication.ComponentWiseBuilderAndSolver(self.linear_solver)
        else:
            if(block_builder):
                # to keep matrix blocks in builder
                builder_and_solver = KratosMultiphysics.ResidualBasedBlockBuilderAndSolver(self.linear_solver)
            else:
                builder_and_solver = SolidMechanicsApplication.ResidualBasedBuilderAndSolver(self.linear_solver)
        
        return builder_and_solver
        
    def GetSolutionSchemeExplicit(self, explicit_integration_scheme, max_delta_time, fraction_delta_time, time_step_prediction_level, rayleigh_damping):
        # creating the explicit solution scheme:
        if(explicit_integration_scheme == "CentralDifferences"):
            mechanical_scheme = SolidMechanicsApplication.ExplicitCentralDifferencesScheme(max_delta_time, 
                                                                            fraction_delta_time, 
                                                                            time_step_prediction_level,
                                                                            rayleigh_damping)
        else:
            raise(explicit_integration_scheme," not implemented yet.")
                                                                                
        return mechanical_scheme
        
    def GetSolutionSchemeImplicit(self, scheme_type, component_wise, compute_contact_forces):
        # creating the implicit solution scheme:
        if(scheme_type == "Newmark"):
          damp_factor_m = 0.0;
          dynamic_factor = 1;
          
          self.main_model_part.ProcessInfo[SolidMechanicsApplication.RAYLEIGH_ALPHA] = 0.0
          self.main_model_part.ProcessInfo[SolidMechanicsApplication.RAYLEIGH_BETA ] = 0.0
          
          if(component_wise):
              mechanical_scheme = SolidMechanicsApplication.ComponentWiseBossakScheme(damp_factor_m, 
                                                                       dynamic_factor)
          else:
              if(compute_contact_forces):
                  mechanical_scheme = SolidMechanicsApplication.ResidualBasedContactBossakScheme(damp_factor_m,
                                                                                  dynamic_factor)
              else:
                  mechanical_scheme = SolidMechanicsApplication.ResidualBasedBossakScheme(damp_factor_m,
                                                                           dynamic_factor)
        elif(scheme_type == "Bossak"):
          damp_factor_m = -0.01;
          dynamic_factor = 1;
          
          self.main_model_part.ProcessInfo[SolidMechanicsApplication.RAYLEIGH_ALPHA] = 0.0
          self.main_model_part.ProcessInfo[SolidMechanicsApplication.RAYLEIGH_BETA ] = 0.0
          
          if(component_wise):
              mechanical_scheme = SolidMechanicsApplication.ComponentWiseBossakScheme(damp_factor_m, 
                                                                       dynamic_factor)
          else:
              if(compute_contact_forces):
                  mechanical_scheme = SolidMechanicsApplication.ResidualBasedContactBossakScheme(damp_factor_m,
                                                                                  dynamic_factor)
              else:
                  mechanical_scheme = SolidMechanicsApplication.ResidualBasedBossakScheme(damp_factor_m,
                                                                           dynamic_factor)
        elif(scheme_type == "Relaxation"):
          #~ self.main_model_part.GetSubModelPart(self.settings["volume_model_part_name"].GetString()).AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)  
            
          import KratosMultiphysics.StructuralMechanicsApplication as StructMechApp
          damp_factor_f = -0.3;
          damp_factor_m = 10.0;
          mechanical_scheme = StructMechApp.ResidualBasedRelaxationScheme(damp_factor_f,
                                                                          damp_factor_m)
                                
        return mechanical_scheme
    
    def GetConvergenceCriterion(self, convergence_criterion_type, rotation_dofs, echo_level, component_wise):

        D_RT = self.settings["displacement_relative_tolerance"].GetDouble()
        D_AT = self.settings["displacement_absolute_tolerance"].GetDouble()
        R_RT = self.settings["residual_relative_tolerance"].GetDouble()
        R_AT = self.settings["residual_absolute_tolerance"].GetDouble()

        if(rotation_dofs):
            if(convergence_criterion_type == "Displacement_criteria"):
                mechanical_convergence_criterion = SolidMechanicsApplication.DisplacementCriteria(D_RT, D_AT)
            elif(convergence_criterion_type == "Residual_criteria"):
                mechanical_convergence_criterion = KratosMultiphysics.ResidualCriteria(R_RT, R_AT)
            elif(convergence_criterion_type == "And_criteria"):
                Displacement = SolidMechanicsApplication.DisplacementCriteria(D_RT, D_AT)
                Residual = KratosMultiphysics.ResidualCriteria(R_RT, R_AT)
                mechanical_convergence_criterion = KratosMultiphysics.AndCriteria(Residual, Displacement)
            elif(convergence_criterion_type == "Or_criteria"):
                Displacement = SolidMechanicsApplication.DisplacementCriteria(D_RT, D_AT)
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
                mechanical_convergence_criterion = SolidMechanicsApplication.DisplacementConvergenceCriterion(D_RT, D_AT)
            elif(convergence_criterion_type == "Residual_criteria"):
                if(component_wise):
                    mechanical_convergence_criterion = SolidMechanicsApplication.ComponentWiseResidualConvergenceCriterion(R_RT, R_AT)
                else:
                    mechanical_convergence_criterion = KratosMultiphysics.ResidualCriteria(R_RT, R_AT)
            elif(convergence_criterion_type == "And_criteria"):
                Displacement = SolidMechanicsApplication.DisplacementConvergenceCriterion(D_RT, D_AT)
                if(component_wise):
                    Residual = SolidMechanicsApplication.ComponentWiseResidualConvergenceCriterion(R_RT, R_AT)
                else:
                    Residual = KratosMultiphysics.ResidualCriteria(R_RT, R_AT)
                mechanical_convergence_criterion = KratosMultiphysics.AndCriteria(Residual, Displacement)
            elif(convergence_criterion_type == "Or_criteria"):
                Displacement = SolidMechanicsApplication.DisplacementConvergenceCriterion(D_RT, D_AT)
                if(self.component_wise):
                    Residual = SolidMechanicsApplication.ComponentWiseResidualConvergenceCriterion(R_RT, R_AT)
                else:
                    Residual = KratosMultiphysics.ResidualCriteria(R_RT, R_AT)
                mechanical_convergence_criterion = KratosMultiphysics.OrCriteria(Residual, Displacement)
            #~ elif(convergence_criterion_type == "Mixed_criteria"):
                #~ Displacement = KratosMultiphysics.MixedElementConvergeCriteria(D_RT, D_AT)
                #~ Residual = KratosMultiphysics.ResidualCriteria(R_RT, R_AT)
                #~ mechanical_convergence_criterion = KratosMultiphysics.AndCriteria(Residual, Displacement)
            
        return mechanical_convergence_criterion
        
    def CreateMechanicalSolver(self, mechanical_scheme, mechanical_convergence_criterion, builder_and_solver, max_iters, compute_reactions, reform_step_dofs, move_mesh_flag, component_wise, line_search, time_integration_method):
        if (time_integration_method == "Implicit"):
            if(component_wise):
                self.mechanical_solver = SolidMechanicsApplication.ComponentWiseNewtonRaphsonStrategy(self.main_model_part, 
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
                    self.mechanical_solver = SolidMechanicsApplication.ResidualBasedNewtonRaphsonLineSearchStrategy(self.main_model_part, 
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
        elif(time_integration_method == "Explicit"):
            self.mechanical_solver = SolidMechanicsApplication.ExplicitStrategy(self.main_model_part, 
                                                                mechanical_scheme, 
                                                                self.linear_solver, 
                                                                compute_reactions, 
                                                                reform_step_dofs, 
                                                                move_mesh_flag)
            self.mechanical_solver.SetRebuildLevel(0) # 1 to recompute the mass matrix in each explicit step 


