from __future__ import print_function, absolute_import, division # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
#import kratos core and applications
import KratosMultiphysics
import KratosMultiphysics.SolidMechanicsApplication
import KratosMultiphysics.PoromechanicsApplication as KratosPoro
import KratosMultiphysics.FluidDynamicsApplication as KratosCFD

# Check that KratosMultiphysics was imported in the main script
KratosMultiphysics.CheckForPreviousImport()


def CreateSolver(main_model_part, custom_settings):

    return FractureUPwSolver(main_model_part, custom_settings)


class FractureUPwSolver(object):

    def __init__(self, main_model_part, custom_settings):

        self.main_model_part = main_model_part

        ##settings string in json format
        default_settings = KratosMultiphysics.Parameters("""
        {
            "solver_type": "poromechanics_fracture_U_Pw_solver",
            "model_part_name": "PorousDomain",
            "domain_size": 2,
            "start_time": 0.0,
            "time_step": 0.1,
            "model_import_settings":{
                "input_type": "mdpa",
                "input_filename": "unknown_name"
            },
            "buffer_size": 2,
            "echo_level": 0,
            "reform_dofs_at_each_step": false,
            "clear_storage": false,
            "compute_reactions": false,
            "move_mesh_flag": false,
            "nodal_smoothing": false,
            "periodic_interface_conditions": false,
            "solution_type": "quasi_static",
            "scheme_type": "Newmark",
            "newmark_beta": 0.25,
            "newmark_gamma": 0.5,
            "newmark_theta": 0.5,
            "rayleigh_m": 0.0,
            "rayleigh_k": 0.0,
            "strategy_type": "newton_raphson",
            "convergence_criterion": "Displacement_criterion",
            "displacement_relative_tolerance": 1.0e-4,
            "displacement_absolute_tolerance": 1.0e-9,
            "residual_relative_tolerance": 1.0e-4,
            "residual_absolute_tolerance": 1.0e-9,
            "max_iteration": 15,
            "desired_iterations": 4,
            "max_radius_factor": 20.0,
            "min_radius_factor": 0.5,
            "block_builder": true,
            "nonlocal_damage": false,
            "characteristic_length": 0.05,
            "search_neighbours_step": false,
            "linear_solver_settings":{
                "solver_type": "AMGCL",
                "tolerance": 1.0e-6,
                "max_iteration": 100,
                "scaling": false,
                "verbosity": 0,
                "preconditioner_type": "ILU0Preconditioner",
                "smoother_type": "ilu0",
                "krylov_type": "gmres",
                "coarsening_type": "aggregation"
            },
            "problem_domain_sub_model_part_list": [""],
            "processes_sub_model_part_list": [""],
            "body_domain_sub_model_part_list": [""],
            "loads_sub_model_part_list": [],
            "loads_variable_list": []
        }
        """)

        # Overwrite the default settings with user-provided parameters
        self.settings = custom_settings
        self.settings.ValidateAndAssignDefaults(default_settings)

        # Construct the linear solver
        import linear_solver_factory
        self.linear_solver = linear_solver_factory.ConstructSolver(self.settings["linear_solver_settings"])

        print("Construction of FractureUPwSolver finished")

    def AddVariables(self):

        ## Solid Variables
        # Add displacements
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        # Add reactions for the displacements
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION)
        # Add dynamic variables
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.ACCELERATION)
        # Add variables for the solid conditions
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.FORCE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.FACE_LOAD)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NORMAL_CONTACT_STRESS)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.TANGENTIAL_CONTACT_STRESS)
        ## Fluid Variables
        # Add water pressure
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.WATER_PRESSURE)
        # Add reactions for the water pressure
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION_WATER_PRESSURE)
        # Add dynamic variables
        self.main_model_part.AddNodalSolutionStepVariable(KratosPoro.DT_WATER_PRESSURE)
        # Add variables for the water conditions
        self.main_model_part.AddNodalSolutionStepVariable(KratosPoro.NORMAL_FLUID_FLUX)
        ## Other variables
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VOLUME_ACCELERATION)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.PERIODIC_PAIR_INDEX)
        if(self.settings["nodal_smoothing"].GetBool() == True):
            self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NODAL_AREA)
            self.main_model_part.AddNodalSolutionStepVariable(KratosPoro.NODAL_CAUCHY_STRESS_TENSOR)
            self.main_model_part.AddNodalSolutionStepVariable(KratosPoro.NODAL_DAMAGE_VARIABLE)
            self.main_model_part.AddNodalSolutionStepVariable(KratosPoro.NODAL_JOINT_AREA)
            self.main_model_part.AddNodalSolutionStepVariable(KratosPoro.NODAL_JOINT_WIDTH)
            self.main_model_part.AddNodalSolutionStepVariable(KratosPoro.NODAL_JOINT_DAMAGE)


        print("Variables correctly added")

    def GetMinimumBufferSize(self):
        return 2

    def AddDofs(self):

        for node in self.main_model_part.Nodes:
            ## Solid dofs
            node.AddDof(KratosMultiphysics.DISPLACEMENT_X,KratosMultiphysics.REACTION_X)
            node.AddDof(KratosMultiphysics.DISPLACEMENT_Y,KratosMultiphysics.REACTION_Y)
            node.AddDof(KratosMultiphysics.DISPLACEMENT_Z,KratosMultiphysics.REACTION_Z)
            ## Fluid dofs
            node.AddDof(KratosMultiphysics.WATER_PRESSURE,KratosMultiphysics.REACTION_WATER_PRESSURE)

        if(self.settings["solution_type"].GetString() == "Dynamic"):
            for node in self.main_model_part.Nodes:
                # adding VELOCITY as dofs
                node.AddDof(KratosMultiphysics.VELOCITY_X)
                node.AddDof(KratosMultiphysics.VELOCITY_Y)
                node.AddDof(KratosMultiphysics.VELOCITY_Z)
                # adding ACCELERATION as dofs
                node.AddDof(KratosMultiphysics.ACCELERATION_X)
                node.AddDof(KratosMultiphysics.ACCELERATION_Y)
                node.AddDof(KratosMultiphysics.ACCELERATION_Z)

        print("DOFs correctly added")

    def ImportModelPart(self):

        if(self.settings["model_import_settings"]["input_type"].GetString() == "mdpa"):

            # Read ModelPart
            KratosMultiphysics.ModelPartIO(self.settings["model_import_settings"]["input_filename"].GetString()).ReadModelPart(self.main_model_part)

            # Create computing_model_part, set constitutive law and buffer size
            self._ExecuteAfterReading()

        else:
            raise Exception("Other input options are not yet implemented.")

        print ("Model reading finished")

    def Initialize(self):

        self.computing_model_part = self.GetComputingModelPart()

        # Set ProcessInfo variables
        if(self.settings["nodal_smoothing"].GetBool() == True):
            self.main_model_part.ProcessInfo.SetValue(KratosPoro.NODAL_SMOOTHING, True)
        else:
            self.main_model_part.ProcessInfo.SetValue(KratosPoro.NODAL_SMOOTHING, False)

        # Builder and solver creation
        builder_and_solver = self._ConstructBuilderAndSolver(self.settings["block_builder"].GetBool())

        # Solution scheme creation
        scheme = self._ConstructScheme(self.settings["scheme_type"].GetString(),
                                         self.settings["solution_type"].GetString())

        # Get the convergence criterion
        convergence_criterion = self._ConstructConvergenceCriterion(self.settings["convergence_criterion"].GetString())

        # Solver creation
        self.Solver = self._ConstructSolver(builder_and_solver,
                                            scheme,
                                            convergence_criterion,
                                            self.settings["strategy_type"].GetString())

        # Set echo_level
        self.Solver.SetEchoLevel(self.settings["echo_level"].GetInt())

        # Check if everything is assigned correctly
        self.Solver.Check()

        print ("Initialization FractureUPwSolver finished")

    def GetComputingModelPart(self):
        return self.main_model_part.GetSubModelPart(self.computing_model_part_name)

    def GetOutputVariables(self):
        pass

    def ComputeDeltaTime(self):
        pass

    def SaveRestart(self):
        pass #one should write the restart file here

    def Solve(self):
        if self.settings["clear_storage"].GetBool():
            self.Clear()

        self.Solver.Solve()

    # solve :: sequencial calls

    def InitializeStrategy(self):
        if self.settings["clear_storage"].GetBool():
            self.Clear()

        self.Solver.Initialize()

    def InitializeSolutionStep(self):
        self.Solver.InitializeSolutionStep()

    def Predict(self):
        self.Solver.Predict()

    def SolveSolutionStep(self):
        self.Solver.SolveSolutionStep()

    def FinalizeSolutionStep(self):
        self.Solver.FinalizeSolutionStep()

    # solve :: sequencial calls

    def SetEchoLevel(self, level):

        self.Solver.SetEchoLevel(level)

    def Clear(self):

        self.Solver.Clear()

    def Check(self):

        self.Solver.Check()

    #### Specific internal functions ####

    def _ExecuteAfterReading(self):

        self.computing_model_part_name = "porous_computational_model_part"

        # Create list of sub sub model parts (it is a copy of the standard lists with a different name)
        import json

        self.body_domain_sub_sub_model_part_list = []
        for i in range(self.settings["body_domain_sub_model_part_list"].size()):
            self.body_domain_sub_sub_model_part_list.append("sub_"+self.settings["body_domain_sub_model_part_list"][i].GetString())
        self.body_domain_sub_sub_model_part_list = KratosMultiphysics.Parameters(json.dumps(self.body_domain_sub_sub_model_part_list))

        self.loads_sub_sub_model_part_list = []
        for i in range(self.settings["loads_sub_model_part_list"].size()):
            self.loads_sub_sub_model_part_list.append("sub_"+self.settings["loads_sub_model_part_list"][i].GetString())
        self.loads_sub_sub_model_part_list = KratosMultiphysics.Parameters(json.dumps(self.loads_sub_sub_model_part_list))

        # Auxiliary parameters object for the CheckAndPepareModelProcess
        params = KratosMultiphysics.Parameters("{}")
        params.AddEmptyValue("computing_model_part_name").SetString(self.computing_model_part_name)
        params.AddValue("problem_domain_sub_model_part_list",self.settings["problem_domain_sub_model_part_list"])
        params.AddValue("processes_sub_model_part_list",self.settings["processes_sub_model_part_list"])
        params.AddValue("body_domain_sub_model_part_list",self.settings["body_domain_sub_model_part_list"])
        params.AddValue("body_domain_sub_sub_model_part_list",self.body_domain_sub_sub_model_part_list)
        params.AddValue("loads_sub_model_part_list",self.settings["loads_sub_model_part_list"])
        params.AddValue("loads_sub_sub_model_part_list",self.loads_sub_sub_model_part_list)
        # CheckAndPrepareModelProcess creates the porous_computational_model_part
        import check_and_prepare_model_process_poro
        check_and_prepare_model_process_poro.CheckAndPrepareModelProcess(self.main_model_part, params).Execute()

        # Constitutive law import
        import poromechanics_constitutivelaw_utility
        poromechanics_constitutivelaw_utility.SetConstitutiveLaw(self.main_model_part)

        self.main_model_part.SetBufferSize( self.settings["buffer_size"].GetInt() )
        minimum_buffer_size = self.GetMinimumBufferSize()
        if(minimum_buffer_size > self.main_model_part.GetBufferSize()):
            self.main_model_part.SetBufferSize( minimum_buffer_size )

    def _ConstructBuilderAndSolver(self, block_builder):

        # Creating the builder and solver
        if(self.settings["periodic_interface_conditions"].GetBool() == True):
            builder_and_solver = KratosCFD.ResidualBasedBlockBuilderAndSolverPeriodic(self.linear_solver,KratosMultiphysics.PERIODIC_PAIR_INDEX)
        else:
            if(block_builder):
                builder_and_solver = KratosMultiphysics.ResidualBasedBlockBuilderAndSolver(self.linear_solver)
            else:
                builder_and_solver = KratosMultiphysics.ResidualBasedEliminationBuilderAndSolver(self.linear_solver)

        return builder_and_solver

    def _ConstructScheme(self, scheme_type, solution_type):

        if(scheme_type == "Newmark"):
            beta = self.settings["newmark_beta"].GetDouble()
            gamma = self.settings["newmark_gamma"].GetDouble()
            theta = self.settings["newmark_theta"].GetDouble()
            rayleigh_m = self.settings["rayleigh_m"].GetDouble()
            rayleigh_k = self.settings["rayleigh_k"].GetDouble()
            if(solution_type == "quasi_static"):
                if(rayleigh_m<1.0e-20 and rayleigh_k<1.0e-20):
                    scheme = KratosPoro.NewmarkQuasistaticUPwScheme(beta,gamma,theta)
                else:
                    scheme = KratosPoro.NewmarkQuasistaticDampedUPwScheme(beta,gamma,theta,rayleigh_m,rayleigh_k)
            else:
                scheme = KratosPoro.NewmarkDynamicUPwScheme(beta,gamma,theta,rayleigh_m,rayleigh_k)
        else:
            raise Exception("Apart from Newmark, other scheme_type are not available.")

        return scheme

    def _ConstructConvergenceCriterion(self, convergence_criterion):

        D_RT = self.settings["displacement_relative_tolerance"].GetDouble()
        D_AT = self.settings["displacement_absolute_tolerance"].GetDouble()
        R_RT = self.settings["residual_relative_tolerance"].GetDouble()
        R_AT = self.settings["residual_absolute_tolerance"].GetDouble()
        echo_level = self.settings["echo_level"].GetInt()

        if(convergence_criterion == "Displacement_criterion"):
            convergence_criterion = KratosMultiphysics.DisplacementCriteria(D_RT, D_AT)
            convergence_criterion.SetEchoLevel(echo_level)
        elif(convergence_criterion == "Residual_criterion"):
            convergence_criterion = KratosMultiphysics.ResidualCriteria(R_RT, R_AT)
            convergence_criterion.SetEchoLevel(echo_level)
        elif(convergence_criterion == "And_criterion"):
            Displacement = KratosMultiphysics.DisplacementCriteria(D_RT, D_AT)
            Displacement.SetEchoLevel(echo_level)
            Residual = KratosMultiphysics.ResidualCriteria(R_RT, R_AT)
            Residual.SetEchoLevel(echo_level)
            convergence_criterion = KratosMultiphysics.AndCriteria(Residual, Displacement)
        elif(convergence_criterion == "Or_criterion"):
            Displacement = KratosMultiphysics.DisplacementCriteria(D_RT, D_AT)
            Displacement.SetEchoLevel(echo_level)
            Residual = KratosMultiphysics.ResidualCriteria(R_RT, R_AT)
            Residual.SetEchoLevel(echo_level)
            convergence_criterion = KratosMultiphysics.OrCriteria(Residual, Displacement)

        return convergence_criterion

    def _ConstructSolver(self, builder_and_solver, scheme, convergence_criterion, strategy_type):

        nonlocal_damage = self.settings["nonlocal_damage"].GetBool()
        max_iters = self.settings["max_iteration"].GetInt()
        compute_reactions = self.settings["compute_reactions"].GetBool()
        reform_step_dofs = self.settings["reform_dofs_at_each_step"].GetBool()
        move_mesh_flag = self.settings["move_mesh_flag"].GetBool()

        if strategy_type == "newton_raphson":
            self.strategy_params = KratosMultiphysics.Parameters("{}")
            self.strategy_params.AddValue("loads_sub_model_part_list",self.loads_sub_sub_model_part_list)
            self.strategy_params.AddValue("loads_variable_list",self.settings["loads_variable_list"])
            if nonlocal_damage:
                self.strategy_params.AddValue("body_domain_sub_model_part_list",self.body_domain_sub_sub_model_part_list)
                self.strategy_params.AddValue("characteristic_length",self.settings["characteristic_length"])
                self.strategy_params.AddValue("search_neighbours_step",self.settings["search_neighbours_step"])
                solver = KratosPoro.PoromechanicsNewtonRaphsonNonlocalStrategy(self.computing_model_part,
                                                                               scheme,
                                                                               self.linear_solver,
                                                                               convergence_criterion,
                                                                               builder_and_solver,
                                                                               self.strategy_params,
                                                                               max_iters,
                                                                               compute_reactions,
                                                                               reform_step_dofs,
                                                                               move_mesh_flag)
            else:
                solver = KratosPoro.PoromechanicsNewtonRaphsonStrategy(self.computing_model_part,
                                                                       scheme,
                                                                       self.linear_solver,
                                                                       convergence_criterion,
                                                                       builder_and_solver,
                                                                       self.strategy_params,
                                                                       max_iters,
                                                                       compute_reactions,
                                                                       reform_step_dofs,
                                                                       move_mesh_flag)
        else:
            # Arc-Length strategy
            self.strategy_params = KratosMultiphysics.Parameters("{}")
            self.strategy_params.AddValue("desired_iterations",self.settings["desired_iterations"])
            self.strategy_params.AddValue("max_radius_factor",self.settings["max_radius_factor"])
            self.strategy_params.AddValue("min_radius_factor",self.settings["min_radius_factor"])
            self.strategy_params.AddValue("loads_sub_model_part_list",self.loads_sub_sub_model_part_list)
            self.strategy_params.AddValue("loads_variable_list",self.settings["loads_variable_list"])
            if nonlocal_damage:
                self.strategy_params.AddValue("body_domain_sub_model_part_list",self.body_domain_sub_sub_model_part_list)
                self.strategy_params.AddValue("characteristic_length",self.settings["characteristic_length"])
                self.strategy_params.AddValue("search_neighbours_step",self.settings["search_neighbours_step"])
                solver = KratosPoro.PoromechanicsRammArcLengthNonlocalStrategy(self.computing_model_part,
                                                                               scheme,
                                                                               self.linear_solver,
                                                                               convergence_criterion,
                                                                               builder_and_solver,
                                                                               self.strategy_params,
                                                                               max_iters,
                                                                               compute_reactions,
                                                                               reform_step_dofs,
                                                                               move_mesh_flag)
            else:
                solver = KratosPoro.PoromechanicsRammArcLengthStrategy(self.computing_model_part,
                                                                       scheme,
                                                                       self.linear_solver,
                                                                       convergence_criterion,
                                                                       builder_and_solver,
                                                                       self.strategy_params,
                                                                       max_iters,
                                                                       compute_reactions,
                                                                       reform_step_dofs,
                                                                       move_mesh_flag)

        return solver

    def _CheckConvergence(self):

        IsConverged = self.Solver.IsConverged()

        return IsConverged

    def _UpdateLoads(self):

        self.Solver.UpdateLoads()
