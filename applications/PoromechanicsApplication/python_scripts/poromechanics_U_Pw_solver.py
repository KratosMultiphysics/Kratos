# Importing the Kratos Library
import KratosMultiphysics
from KratosMultiphysics.python_solver import PythonSolver

# Import applications
import KratosMultiphysics.FluidDynamicsApplication as KratosCFD
import KratosMultiphysics.PoromechanicsApplication as KratosPoro
import KratosMultiphysics.StructuralMechanicsApplication as KratosStructural

def CreateSolver(model, custom_settings):
    return UPwSolver(model, custom_settings)

class UPwSolver(PythonSolver):
    '''Solver for the solution of displacement-pore pressure coupled problems.'''

    def __init__(self, model, custom_settings):

        self._validate_settings_in_baseclass=True # To be removed eventually

        super(UPwSolver,self).__init__(model, custom_settings)

        self.min_buffer_size = 2

        # Either retrieve the model part from the model or create a new one
        model_part_name = self.settings["model_part_name"].GetString()

        if model_part_name == "":
            raise Exception('Please specify a model_part name!')

        if self.model.HasModelPart(model_part_name):
            self.main_model_part = self.model.GetModelPart(model_part_name)
        else:
            self.main_model_part = self.model.CreateModelPart(model_part_name,self.min_buffer_size)

        self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE,
                                                  self.settings["domain_size"].GetInt())

        KratosMultiphysics.Logger.PrintInfo("UPwSolver", "Construction of UPwSolver finished.")

    @classmethod
    def GetDefaultParameters(cls):
        this_defaults = KratosMultiphysics.Parameters("""{
            "solver_type": "poromechanics_U_Pw_solver",
            "model_part_name": "PorousModelPart",
            "domain_size": 2,
            "start_time": 0.0,
            "time_step": 0.1,
            "model_import_settings":{
                "input_type": "mdpa",
                "input_filename": "unknown_name"
            },
            "material_import_settings" :{
                "materials_filename": ""
            },
            "buffer_size": 2,
            "echo_level": 0,
            "reform_dofs_at_each_step": false,
            "clear_storage": false,
            "compute_reactions": false,
            "move_mesh_flag": false,
            "nodal_smoothing": false,
            "gp_to_nodal_variable_list": [],
            "gp_to_nodal_variable_extrapolate_non_historical": false,
            "periodic_interface_conditions": false,
            "solution_type": "implicit_quasi_static",
            "scheme_type": "Newmark",
            "newmark_beta": 0.25,
            "newmark_gamma": 0.5,
            "newmark_theta": 0.5,
            "calculate_alpha_beta"       : false,
            "omega_1"                    : 1.0,
            "omega_n"                    : 10.0,
            "xi_1"                       : 1.0,
            "xi_n"                       : 0.05,
            "rayleigh_alpha": 0.0,
            "rayleigh_beta": 0.0,
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
                "solver_type": "amgcl",
                "tolerance": 1.0e-6,
                "max_iteration": 100,
                "scaling": false,
                "verbosity": 0,
                "preconditioner_type": "ilu0",
                "smoother_type": "ilu0",
                "krylov_type": "gmres",
                "coarsening_type": "aggregation"
            },
            "problem_domain_sub_model_part_list": [""],
            "processes_sub_model_part_list": [""],
            "body_domain_sub_model_part_list": [""],
            "loads_sub_model_part_list": [],
            "loads_variable_list": []
        }""")

        this_defaults.AddMissingParameters(super(UPwSolver, cls).GetDefaultParameters())
        return this_defaults

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

        self.main_model_part.AddNodalSolutionStepVariable(KratosPoro.NODAL_DAMAGE_VARIABLE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosPoro.NODAL_JOINT_AREA)
        self.main_model_part.AddNodalSolutionStepVariable(KratosPoro.NODAL_JOINT_WIDTH)
        self.main_model_part.AddNodalSolutionStepVariable(KratosPoro.NODAL_JOINT_DAMAGE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NODAL_AREA)
        self.main_model_part.AddNodalSolutionStepVariable(KratosPoro.NODAL_EFFECTIVE_STRESS_TENSOR)
        self.main_model_part.AddNodalSolutionStepVariable(KratosPoro.NODAL_WATER_PRESSURE_GRADIENT)
        self.main_model_part.AddNodalSolutionStepVariable(KratosPoro.INITIAL_STRESS_TENSOR)

        # Add variables from gp to nodal variable list
        if self.settings["gp_to_nodal_variable_extrapolate_non_historical"].GetBool()==False:
            for i in range(self.settings["gp_to_nodal_variable_list"].size()):
                variable_name = self.settings["gp_to_nodal_variable_list"][i].GetString()
                variable = KratosMultiphysics.KratosGlobals.GetVariable(variable_name)
                self.main_model_part.AddNodalSolutionStepVariable(variable)

        KratosMultiphysics.Logger.PrintInfo("UPwSolver", "Variables added correctly.")

    def ImportModelPart(self):
        # we can use the default implementation in the base class
        self._ImportModelPart(self.main_model_part,self.settings["model_import_settings"])

    def PrepareModelPart(self):

        # Set ProcessInfo variables
        self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.TIME,
                                                  self.settings["start_time"].GetDouble())
        self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.DELTA_TIME,
                                                  self.settings["time_step"].GetDouble())
        self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.STEP, 0)
        self.main_model_part.ProcessInfo.SetValue(KratosPoro.TIME_UNIT_CONVERTER, 1.0)
        if(self.settings["nodal_smoothing"].GetBool() == True):
            self.main_model_part.ProcessInfo.SetValue(KratosPoro.NODAL_SMOOTHING, True)
        else:
            self.main_model_part.ProcessInfo.SetValue(KratosPoro.NODAL_SMOOTHING, False)

        if not self.main_model_part.ProcessInfo[KratosMultiphysics.IS_RESTARTED]:
            ## Executes the check and prepare model process (Create computing_model_part and set constitutive law)
            self._ExecuteCheckAndPrepare()
            ## Set buffer size
            self._SetBufferSize()

        KratosMultiphysics.Logger.PrintInfo("UPwSolver", "Model reading finished.")

    def AddDofs(self):
        ## Solid dofs
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_X, KratosMultiphysics.REACTION_X,self.main_model_part)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_Y, KratosMultiphysics.REACTION_Y,self.main_model_part)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_Z, KratosMultiphysics.REACTION_Z,self.main_model_part)
        ## Fluid dofs
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.WATER_PRESSURE, KratosMultiphysics.REACTION_WATER_PRESSURE,self.main_model_part)

        if(self.settings["solution_type"].GetString() == "implicit_dynamic"):
            for node in self.main_model_part.Nodes:
                # adding VELOCITY as dofs
                node.AddDof(KratosMultiphysics.VELOCITY_X)
                node.AddDof(KratosMultiphysics.VELOCITY_Y)
                node.AddDof(KratosMultiphysics.VELOCITY_Z)
                # adding ACCELERATION as dofs
                node.AddDof(KratosMultiphysics.ACCELERATION_X)
                node.AddDof(KratosMultiphysics.ACCELERATION_Y)
                node.AddDof(KratosMultiphysics.ACCELERATION_Z)

        KratosMultiphysics.Logger.PrintInfo("UPwSolver", "DOFs added correctly.")

    def GetMinimumBufferSize(self):
        return self.min_buffer_size

    def Initialize(self):
        self.computing_model_part = self.GetComputingModelPart()

        # Fill the previous steps of the buffer with the initial conditions
        self._FillBuffer()

        # Construct the linear solver
        self.linear_solver = self._ConstructLinearSolver()

        # Builder and solver creation
        builder_and_solver = self._ConstructBuilderAndSolver(self.settings["block_builder"].GetBool())

        # Solution scheme creation
        self.scheme = self._ConstructScheme(self.settings["scheme_type"].GetString(),
                                         self.settings["solution_type"].GetString())

        # Get the convergence criterion
        self.convergence_criterion = self._ConstructConvergenceCriterion(self.settings["convergence_criterion"].GetString())

        # Solver creation
        self.solver = self._ConstructSolver(builder_and_solver,
                                            self.settings["strategy_type"].GetString())

        # Set echo_level
        self.SetEchoLevel(self.settings["echo_level"].GetInt())

        # Initialize Strategy
        if self.settings["clear_storage"].GetBool():
            self.Clear()

        self.solver.Initialize()

        # Check if everything is assigned correctly
        self.Check()

        # Check and construct gp_to_nodal_variable process
        self._CheckAndConstructGPtoNodalVariableExtrapolationProcess()

        KratosMultiphysics.Logger.PrintInfo("UPwSolver", "Solver initialization finished.")

    def GetComputingModelPart(self):
        return self.main_model_part.GetSubModelPart(self.computing_model_part_name)

    def ComputeDeltaTime(self):
        return self.main_model_part.ProcessInfo[KratosMultiphysics.DELTA_TIME]

    def Clear(self):
        self.solver.Clear()

    def Check(self):
        self.solver.Check()

    def SetEchoLevel(self, level):
        self.solver.SetEchoLevel(level)

    def AdvanceInTime(self, current_time):
        dt = self.ComputeDeltaTime()
        new_time = current_time + dt

        self.main_model_part.CloneTimeStep(new_time)
        self.main_model_part.ProcessInfo[KratosMultiphysics.STEP] += 1

        return new_time

    def InitializeSolutionStep(self):
        self.solver.InitializeSolutionStep()

    def Predict(self):
        self.solver.Predict()

    def SolveSolutionStep(self):
        is_converged = self.solver.SolveSolutionStep()
        return is_converged

    def FinalizeSolutionStep(self):
        self.solver.FinalizeSolutionStep()

        if self.gp_to_nodal_variable:
            self.integration_values_extrapolation_to_nodes_process.ExecuteBeforeSolutionLoop()
            self.integration_values_extrapolation_to_nodes_process.ExecuteFinalizeSolutionStep()

    def Finalize(self):
        if self.gp_to_nodal_variable:
            self.integration_values_extrapolation_to_nodes_process.ExecuteFinalize()

    def Solve(self):
        message = "".join([
            "Calling UPwSolver.Solve() method, which is deprecated\n",
            "Please call the individual methods instead:\n",
            "solver.InitializeSolutionStep()\n",
            "solver.Predict()\n",
            "solver.SolveSolutionStep()\n",
            "solver.FinalizeSolutionStep()\n"]
        )
        KratosMultiphysics.Logger.PrintWarning("UPwSolver",message)

        if self.settings["clear_storage"].GetBool():
            self.Clear()

        self.InitializeSolutionStep()
        self.Predict()
        self.SolveSolutionStep()
        self.FinalizeSolutionStep()

    #### Specific internal functions ####

    def _ExecuteCheckAndPrepare(self):

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
        from KratosMultiphysics.PoromechanicsApplication import check_and_prepare_model_process_poro
        check_and_prepare_model_process_poro.CheckAndPrepareModelProcess(self.main_model_part, params).Execute()

        # Constitutive law import
        materials_imported = self.import_constitutive_laws()
        if materials_imported:
            KratosMultiphysics.Logger.PrintInfo("UPwSolver", "Constitutive law was successfully imported via json.")
        else:
            KratosMultiphysics.Logger.PrintInfo("UPwSolver", "Constitutive law was not successfully imported.")

    def import_constitutive_laws(self):
        materials_filename = self.settings["material_import_settings"]["materials_filename"].GetString()
        if (materials_filename != ""):
            # Add constitutive laws and material properties from json file to model parts.
            material_settings = KratosMultiphysics.Parameters("""{"Parameters": {"materials_filename": ""}} """)
            material_settings["Parameters"]["materials_filename"].SetString(materials_filename)
            KratosMultiphysics.ReadMaterialsUtility(material_settings, self.model)
            materials_imported = True
        else:
            materials_imported = False
        return materials_imported

    def _SetBufferSize(self):
        required_buffer_size = self.settings["buffer_size"].GetInt()
        if required_buffer_size < self.GetMinimumBufferSize():
            required_buffer_size = self.GetMinimumBufferSize()
        current_buffer_size = self.main_model_part.GetBufferSize()
        buffer_size = max(current_buffer_size, required_buffer_size)
        self.main_model_part.SetBufferSize(buffer_size)

    def _FillBuffer(self):
        buffer_size = self.main_model_part.GetBufferSize()
        time = self.main_model_part.ProcessInfo[KratosMultiphysics.TIME]
        delta_time = self.main_model_part.ProcessInfo[KratosMultiphysics.DELTA_TIME]
        step = self.main_model_part.ProcessInfo[KratosMultiphysics.STEP]

        step = step - (buffer_size-1)*1
        self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.STEP, step)
        time = time - (buffer_size-1)*delta_time
        self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.TIME, time)
        for i in range(buffer_size-1):
            step = step + 1
            self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.STEP, step)
            time = time + delta_time
            self.main_model_part.CloneTimeStep(time)

    def _ConstructLinearSolver(self):
        import KratosMultiphysics.python_linear_solver_factory as linear_solver_factory
        return linear_solver_factory.ConstructSolver(self.settings["linear_solver_settings"])

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

        self.main_model_part.ProcessInfo.SetValue(KratosPoro.VELOCITY_COEFFICIENT, 1.0)
        self.main_model_part.ProcessInfo.SetValue(KratosPoro.DT_PRESSURE_COEFFICIENT, 1.0)

        if(scheme_type == "Newmark"):
            beta = self.settings["newmark_beta"].GetDouble()
            gamma = self.settings["newmark_gamma"].GetDouble()
            theta = self.settings["newmark_theta"].GetDouble()
            rayleigh_alpha = self.settings["rayleigh_alpha"].GetDouble()
            rayleigh_beta = self.settings["rayleigh_beta"].GetDouble()
            if self.settings["calculate_alpha_beta"].GetBool():
                omega_1 = self.settings["omega_1"].GetDouble()
                omega_n = self.settings["omega_n"].GetDouble()
                xi_1 = self.settings["xi_1"].GetDouble()
                xi_n = self.settings["xi_n"].GetDouble()
                rayleigh_beta = 2.0*(xi_n*omega_n-xi_1*omega_1)/(omega_n*omega_n-omega_1*omega_1)
                rayleigh_alpha = 2.0*xi_1*omega_1-rayleigh_beta*omega_1*omega_1
                KratosMultiphysics.Logger.PrintInfo("::[UPwSolver]:: Scheme Information")
                KratosMultiphysics.Logger.PrintInfo("::[UPwSolver]:: omega_1: ",omega_1)
                KratosMultiphysics.Logger.PrintInfo("::[UPwSolver]:: omega_n: ",omega_n)
                KratosMultiphysics.Logger.PrintInfo("::[UPwSolver]:: xi_1: ",xi_1)
                KratosMultiphysics.Logger.PrintInfo("::[UPwSolver]:: xi_n: ",xi_n)
                KratosMultiphysics.Logger.PrintInfo("::[UPwSolver]:: Alpha and Beta output")
                KratosMultiphysics.Logger.PrintInfo("::[UPwSolver]:: rayleigh_alpha: ",rayleigh_alpha)
                KratosMultiphysics.Logger.PrintInfo("::[UPwSolver]:: rayleigh_beta: ",rayleigh_beta)
            
            self.main_model_part.ProcessInfo.SetValue(KratosStructural.RAYLEIGH_ALPHA,rayleigh_alpha)
            self.main_model_part.ProcessInfo.SetValue(KratosStructural.RAYLEIGH_BETA,rayleigh_beta)
            if(solution_type == "implicit_quasi_static"):
                if(rayleigh_alpha<1.0e-20 and rayleigh_beta<1.0e-20):
                    scheme = KratosPoro.NewmarkQuasistaticUPwScheme(beta,gamma,theta)
                else:
                    scheme = KratosPoro.NewmarkQuasistaticDampedUPwScheme(beta,gamma,theta)
            else:
                scheme = KratosPoro.NewmarkDynamicUPwScheme(beta,gamma,theta)
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

    def _ConstructSolver(self, builder_and_solver, strategy_type):

        self.main_model_part.ProcessInfo.SetValue(KratosPoro.IS_CONVERGED, True)
        self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.NL_ITERATION_NUMBER, 1)

        nonlocal_damage = self.settings["nonlocal_damage"].GetBool()
        max_iters = self.settings["max_iteration"].GetInt()
        compute_reactions = self.settings["compute_reactions"].GetBool()
        reform_step_dofs = self.settings["reform_dofs_at_each_step"].GetBool()
        move_mesh_flag = self.settings["move_mesh_flag"].GetBool()

        self.strategy_params = KratosMultiphysics.Parameters("{}")
        self.strategy_params.AddValue("loads_sub_model_part_list",self.loads_sub_sub_model_part_list)
        self.strategy_params.AddValue("loads_variable_list",self.settings["loads_variable_list"])

        if strategy_type == "newton_raphson":
            if nonlocal_damage:
                self.strategy_params.AddValue("body_domain_sub_model_part_list",self.body_domain_sub_sub_model_part_list)
                self.strategy_params.AddValue("characteristic_length",self.settings["characteristic_length"])
                self.strategy_params.AddValue("search_neighbours_step",self.settings["search_neighbours_step"])
                solving_strategy = KratosPoro.PoromechanicsNewtonRaphsonNonlocalStrategy(self.computing_model_part,
                                                                               self.scheme,
                                                                               self.convergence_criterion,
                                                                               builder_and_solver,
                                                                               self.strategy_params,
                                                                               max_iters,
                                                                               compute_reactions,
                                                                               reform_step_dofs,
                                                                               move_mesh_flag)
            else:
                solving_strategy = KratosPoro.PoromechanicsNewtonRaphsonStrategy(self.computing_model_part,
                                                                       self.scheme,
                                                                       self.convergence_criterion,
                                                                       builder_and_solver,
                                                                       self.strategy_params,
                                                                       max_iters,
                                                                       compute_reactions,
                                                                       reform_step_dofs,
                                                                       move_mesh_flag)
        else:
            # Arc-Length strategy
            self.main_model_part.ProcessInfo.SetValue(KratosPoro.ARC_LENGTH_LAMBDA,1.0)
            self.main_model_part.ProcessInfo.SetValue(KratosPoro.ARC_LENGTH_RADIUS_FACTOR,1.0)

            self.strategy_params.AddValue("desired_iterations",self.settings["desired_iterations"])
            self.strategy_params.AddValue("max_radius_factor",self.settings["max_radius_factor"])
            self.strategy_params.AddValue("min_radius_factor",self.settings["min_radius_factor"])
            if nonlocal_damage:
                self.strategy_params.AddValue("body_domain_sub_model_part_list",self.body_domain_sub_sub_model_part_list)
                self.strategy_params.AddValue("characteristic_length",self.settings["characteristic_length"])
                self.strategy_params.AddValue("search_neighbours_step",self.settings["search_neighbours_step"])
                solving_strategy = KratosPoro.PoromechanicsRammArcLengthNonlocalStrategy(self.computing_model_part,
                                                                               self.scheme,
                                                                               self.convergence_criterion,
                                                                               builder_and_solver,
                                                                               self.strategy_params,
                                                                               max_iters,
                                                                               compute_reactions,
                                                                               reform_step_dofs,
                                                                               move_mesh_flag)
            else:
                solving_strategy = KratosPoro.PoromechanicsRammArcLengthStrategy(self.computing_model_part,
                                                                       self.scheme,
                                                                       self.convergence_criterion,
                                                                       builder_and_solver,
                                                                       self.strategy_params,
                                                                       max_iters,
                                                                       compute_reactions,
                                                                       reform_step_dofs,
                                                                       move_mesh_flag)

        return solving_strategy

    def _CheckConvergence(self):

        IsConverged = self.solver.IsConverged()

        return IsConverged

    def _UpdateLoads(self):

        self.solver.UpdateLoads()

    def _CheckAndConstructGPtoNodalVariableExtrapolationProcess(self):
        if self.settings["gp_to_nodal_variable_list"].size() > 0:
            self.gp_to_nodal_variable = True
            # Create extrapolation process
            extrapolation_parameters = KratosMultiphysics.Parameters("""
            {
                "model_part_name"            : "",
                "echo_level"                 : 0,
                "average_variable"           : "NODAL_AREA",
                "area_average"               : true,
                "list_of_variables"          : [],
                "extrapolate_non_historical" : false
            }
            """)
            for i in range(self.settings["gp_to_nodal_variable_list"].size()):
                var_name = self.settings["gp_to_nodal_variable_list"][i].GetString()
                extrapolation_parameters["list_of_variables"].Append(var_name)
            extrapolation_parameters["model_part_name"].SetString(self.settings["model_part_name"].GetString())
            extrapolation_parameters["extrapolate_non_historical"].SetBool(self.settings["gp_to_nodal_variable_extrapolate_non_historical"].GetBool())
            self.integration_values_extrapolation_to_nodes_process = KratosMultiphysics.IntegrationValuesExtrapolationToNodesProcess(self.main_model_part, extrapolation_parameters)
        else:
            self.gp_to_nodal_variable = False
            self.integration_values_extrapolation_to_nodes_process = KratosMultiphysics.Process()