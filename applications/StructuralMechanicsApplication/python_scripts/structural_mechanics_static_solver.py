from __future__ import print_function, absolute_import, division  # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
#import kratos core and applications
import KratosMultiphysics
import KratosMultiphysics.SolidMechanicsApplication
import KratosMultiphysics.StructuralMechanicsApplication

# Check that KratosMultiphysics was imported in the main script
KratosMultiphysics.CheckForPreviousImport()

# Import the implicit solver (the explicit one is derived from it)
import solid_mechanics_static_solver

def IncreasePointLoad(forcing_nodes_list, Load):
    for node in forcing_nodes_list:
        node.SetSolutionStepValue(KratosMultiphysics.SolidMechanicsApplication.POINT_LOAD_X, 0, Load[0])
        node.SetSolutionStepValue(KratosMultiphysics.SolidMechanicsApplication.POINT_LOAD_Y, 0, Load[1])
        node.SetSolutionStepValue(KratosMultiphysics.SolidMechanicsApplication.POINT_LOAD_Z, 0, Load[2])

def IncreaseDisplacement(forcing_nodes_list, disp):
    for node in forcing_nodes_list:
        node.SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_X, 0, disp[0])
        node.SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Y, 0, disp[1])
        node.SetSolutionStepValue(KratosMultiphysics.DISPLACEMENT_Z, 0, disp[2])

def ChangeCondition(model_part, echo_level):
    force_x = 0.0
    force_y = 0.0
    force_z = 0.0
    for node in model_part.Nodes:
        new_load = node.GetSolutionStepValue(KratosMultiphysics.SolidMechanicsApplication.POINT_LOAD, 0) * model_part.ProcessInfo[KratosMultiphysics.StructuralMechanicsApplication.LAMBDA];
        force_x += new_load[0]
        force_y += new_load[1]
        force_z += new_load[2]
        node.SetSolutionStepValue(KratosMultiphysics.SolidMechanicsApplication.POINT_LOAD, 0, new_load)
    
    if (echo_level > 0):
        print("*********************** ")
        print("The total load applied: ")
        print("POINT_LOAD_X: ", force_x)
        print("POINT_LOAD_Y: ", force_y)
        print("POINT_LOAD_Z: ", force_z)
        print("*********************** ")

def CreateSolver(main_model_part, custom_settings):
    return StaticStructuralSolver(main_model_part, custom_settings)

class StaticStructuralSolver(solid_mechanics_static_solver.StaticMechanicalSolver):
    
    
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
            "solver_type": "structural_mechanics_static_solver",
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
            "reform_dofs_at_each_iteration": false,
            "line_search": false,
            "compute_reactions": true,
            "compute_contact_forces": false,
            "compute_mortar_contact": false,
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
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.SolidMechanicsApplication.POINT_LOAD)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.SolidMechanicsApplication.LINE_LOAD)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.SolidMechanicsApplication.SURFACE_LOAD)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VOLUME_ACCELERATION)
            
        if self.settings["rotation_dofs"].GetBool():
            # Add specific variables for the problem (rotation dofs)
            self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.ROTATION)
            self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.TORQUE)
            self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.StructuralMechanicsApplication.POINT_TORQUE)
            self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.ANGULAR_VELOCITY)
            self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.ANGULAR_ACCELERATION)
        if self.settings["pressure_dofs"].GetBool():
            # Add specific variables for the problem (pressure dofs)
            self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.PRESSURE)
            self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.SolidMechanicsApplication.PRESSURE_REACTION)
            
        if  self.settings["compute_mortar_contact"].GetBool():
            # Add normal
            self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NORMAL)
            # Add lagrange multiplier
            self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VECTOR_LAGRANGE_MULTIPLIER)
            # Add weighted gap
            self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.StructuralMechanicsApplication.WEIGHTED_GAP)
            
        if self.settings["analysis_type"].GetString() == "Arc-Length":
            self.main_model_part.ProcessInfo[KratosMultiphysics.StructuralMechanicsApplication.LAMBDA] = 0.00;
   
        print("::[Mechanical Solver]:: Variables ADDED")
        
    def AddDofs(self):

        for node in self.main_model_part.Nodes:
            # adding dofs
            node.AddDof(KratosMultiphysics.DISPLACEMENT_X, KratosMultiphysics.REACTION_X);
            node.AddDof(KratosMultiphysics.DISPLACEMENT_Y, KratosMultiphysics.REACTION_Y);
            node.AddDof(KratosMultiphysics.DISPLACEMENT_Z, KratosMultiphysics.REACTION_Z);
            
        if self.settings["rotation_dofs"].GetBool():
            for node in self.main_model_part.Nodes:
                node.AddDof(KratosMultiphysics.ROTATION_X, KratosMultiphysics.TORQUE_X);
                node.AddDof(KratosMultiphysics.ROTATION_Y, KratosMultiphysics.TORQUE_Y);
                node.AddDof(KratosMultiphysics.ROTATION_Z, KratosMultiphysics.TORQUE_Z);
                
        if self.settings["pressure_dofs"].GetBool():                
            for node in self.main_model_part.Nodes:
                node.AddDof(KratosMultiphysics.PRESSURE, KratosSolid.PRESSURE_REACTION);
        
        if  self.settings["compute_mortar_contact"].GetBool():
            for node in self.main_model_part.Nodes:
                node.AddDof(KratosMultiphysics.VECTOR_LAGRANGE_MULTIPLIER_X);
                node.AddDof(KratosMultiphysics.VECTOR_LAGRANGE_MULTIPLIER_Y);
                node.AddDof(KratosMultiphysics.VECTOR_LAGRANGE_MULTIPLIER_Z);

        print("::[Mechanical Solver]:: DOF's ADDED")

    def GetComputeModelPart(self):
        return self.main_model_part
    
    def Solve(self):
        
        if self.settings["clear_storage"].GetBool():
            self.Clear()
            
        self.mechanical_solver.Solve()
        
        if self.settings["analysis_type"].GetString() == "Arc-Length":
            echo_level =  self.settings["echo_level"].GetInt()
            if echo_level > 0:
                print("LAMBDA: ", self.main_model_part.ProcessInfo[KratosMultiphysics.StructuralMechanicsApplication.LAMBDA])
            ChangeCondition(self.main_model_part, echo_level)
    
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
                    mechanical_scheme = KratosMultiphysics.StructuralMechanicsApplication.ResidualBasedIncrementalUpdateStaticContactScheme()
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
            Mortar = KratosMultiphysics.StructuralMechanicsApplication.MortarConvergenceCriteria()
            Mortar.SetEchoLevel(self.settings["echo_level"].GetInt())

            convergence_criterion.mechanical_convergence_criterion = KratosMultiphysics.AndCriteria(Mortar, convergence_criterion.mechanical_convergence_criterion)
        
        return convergence_criterion.mechanical_convergence_criterion
        
    def _CreateMechanicalSolver(self, mechanical_scheme, mechanical_convergence_criterion, builder_and_solver, max_iters, compute_reactions, reform_step_dofs, move_mesh_flag, component_wise, line_search):
        
        if(component_wise):
            self.mechanical_solver = KratosMultiphysics.SolidMechanicsApplication.ComponentWiseNewtonRaphsonStrategy(
                                                                            self.compute_model_part, 
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
                                                                            self.compute_model_part, 
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
                                                                            self.compute_model_part, 
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
                                                                            self.compute_model_part, 
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
                    if  self.settings["compute_mortar_contact"].GetBool():
                        split_factor   = self.settings["split_factor"].GetDouble()
                        max_number_splits = self.settings["max_number_splits"].GetInt()
                        self.mechanical_solver = KratosMultiphysics.StructuralMechanicsApplication.ResidualBasedNewtonRaphsonContactStrategy(
                                                                                self.compute_model_part, 
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
                                                                                self.compute_model_part, 
                                                                                mechanical_scheme, 
                                                                                self.linear_solver, 
                                                                                mechanical_convergence_criterion, 
                                                                                builder_and_solver, 
                                                                                max_iters, 
                                                                                compute_reactions, 
                                                                                reform_step_dofs, 
                                                                                move_mesh_flag)
