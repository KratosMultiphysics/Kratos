# Importing the Kratos Library
import KratosMultiphysics
from KratosMultiphysics.python_solver import PythonSolver

# Import applications
import KratosMultiphysics.GeoMechanicsApplication as KratosGeo
import KratosMultiphysics.StructuralMechanicsApplication as KratosStructure

# Import base class file
import KratosMultiphysics.GeoMechanicsApplication.geomechanics_solver as GeoSolver

def CreateSolver(model, custom_settings):
    return UPwSolver(model, custom_settings)

class UPwSolver(GeoSolver.GeoMechanicalSolver):
    '''Solver for the solution of displacement-pore pressure coupled problems.'''

    def __init__(self, model, custom_settings):
        super().__init__(model, custom_settings)

        KratosMultiphysics.Logger.PrintInfo("GeoMechanics_U_Pw_Solver", "Construction of Solver finished.")

    @classmethod
    def GetDefaultParameters(cls):
        this_defaults = KratosMultiphysics.Parameters("""{
            "solver_type": "geomechanics_U_Pw_solver",
            "model_part_name": "PorousDomain",
            "domain_size": 2,
            "start_time": 0.0,
            "model_import_settings":{
                "input_type": "mdpa",
                "input_filename": "unknown_name"
            },
            "material_import_settings" :{
                "materials_filename": ""
            },
            "time_stepping": {
                "end_time" : 1.0,
                "time_step": 0.1
            },
            "buffer_size": 2,
            "echo_level": 0,
            "reform_dofs_at_each_step": false,
            "clear_storage": false,
            "compute_reactions": false,
            "move_mesh_flag": false,
            "nodal_smoothing": false,
            "reset_displacements":  false,
            "solution_type": "quasi_static",
            "scheme_type": "Newmark",
            "newmark_beta": 0.25,
            "newmark_gamma": 0.5,
            "newmark_theta": 0.5,
            "rayleigh_m": 0.0,
            "rayleigh_k": 0.0,
            "strategy_type": "newton_raphson",
            "convergence_criterion": "Displacement_criterion",
            "water_pressure_relative_tolerance": 1.0e-4,
            "water_pressure_absolute_tolerance": 1.0e-9,
            "displacement_relative_tolerance": 1.0e-4,
            "displacement_absolute_tolerance": 1.0e-9,
            "residual_relative_tolerance": 1.0e-4,
            "residual_absolute_tolerance": 1.0e-9,
            "desired_iterations"         : 4,
            "max_radius_factor"          : 20.0,
            "min_radius_factor"          : 0.5,
            "max_iterations"             : 15,
            "min_iterations"             : 6,
            "number_cycles"              : 5,
            "increase_factor"            : 2.0,
            "reduction_factor"           : 0.5,
            "realised_factor"            : 1.0,
            "calculate_reactions"        : true,
            "max_line_search_iterations" : 5,
            "first_alpha_value"          : 0.5,
            "second_alpha_value"         : 1.0,
            "min_alpha"                  : 0.1,
            "max_alpha"                  : 2.0,
            "line_search_tolerance"      : 0.5,
            "rotation_dofs"              : false,
            "block_builder"              : true,
            "search_neighbours_step"     : false,
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
        }""")

        this_defaults.AddMissingParameters(super().GetDefaultParameters())
        return this_defaults

    def PrepareModelPart(self):

        # Set ProcessInfo variables
        self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.TIME,
                                                  self.settings["start_time"].GetDouble())
        self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.DELTA_TIME,
                                                  self.settings["time_stepping"]["time_step"].GetDouble())
        self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.STEP, 0)

        self.main_model_part.ProcessInfo.SetValue(KratosGeo.TIME_UNIT_CONVERTER, 1.0)

        self.main_model_part.ProcessInfo.SetValue(KratosGeo.NODAL_SMOOTHING,
                                                  self.settings["nodal_smoothing"].GetBool())

        if not self.main_model_part.ProcessInfo[KratosMultiphysics.IS_RESTARTED]:
            ## Executes the check and prepare model process (Create computing_model_part and set constitutive law)
            self._ExecuteCheckAndPrepare()
            ## Set buffer size
            self._SetBufferSize()

        if not self.model.HasModelPart(self.settings["model_part_name"].GetString()):
            self.model.AddModelPart(self.main_model_part)

        KratosMultiphysics.Logger.PrintInfo("GeoMechanics_U_Pw_Solver", "Model reading finished.")

    def AddDofs(self):
        ## displacement dofs
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_X, KratosMultiphysics.REACTION_X,self.main_model_part)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_Y, KratosMultiphysics.REACTION_Y,self.main_model_part)
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.DISPLACEMENT_Z, KratosMultiphysics.REACTION_Z,self.main_model_part)

        for node in self.main_model_part.Nodes:
            # adding TOTAL_DISPLACEMENT as dofs
            node.AddDof(KratosGeo.TOTAL_DISPLACEMENT_X)
            node.AddDof(KratosGeo.TOTAL_DISPLACEMENT_Y)
            node.AddDof(KratosGeo.TOTAL_DISPLACEMENT_Z)

        ## Fluid dofs
        KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.WATER_PRESSURE, KratosMultiphysics.REACTION_WATER_PRESSURE,self.main_model_part)

        if self.settings["rotation_dofs"].GetBool():
            KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.ROTATION_X, KratosMultiphysics.REACTION_MOMENT_X,self.main_model_part)
            KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.ROTATION_Y, KratosMultiphysics.REACTION_MOMENT_Y,self.main_model_part)
            KratosMultiphysics.VariableUtils().AddDof(KratosMultiphysics.ROTATION_Z, KratosMultiphysics.REACTION_MOMENT_Z,self.main_model_part)

        if (self.settings["solution_type"].GetString() == "Dynamic"):
            KratosMultiphysics.Logger.PrintInfo("GeoMechanics_U_Pw_Solver", "Dynamic analysis.")
            for node in self.main_model_part.Nodes:
                # adding VELOCITY as dofs
                node.AddDof(KratosMultiphysics.VELOCITY_X)
                node.AddDof(KratosMultiphysics.VELOCITY_Y)
                node.AddDof(KratosMultiphysics.VELOCITY_Z)
                # adding ACCELERATION as dofs
                node.AddDof(KratosMultiphysics.ACCELERATION_X)
                node.AddDof(KratosMultiphysics.ACCELERATION_Y)
                node.AddDof(KratosMultiphysics.ACCELERATION_Z)
            # if self.settings["rotation_dofs"].GetBool():
            #     for node in self.main_model_part.Nodes:
            #         # adding ANGULAR_VELOCITY as dofs
            #         node.AddDof(KratosMultiphysics.ANGULAR_VELOCITY_X)
            #         node.AddDof(KratosMultiphysics.ANGULAR_VELOCITY_Y)
            #         node.AddDof(KratosMultiphysics.ANGULAR_VELOCITY_Z)
            #         # adding ANGULAR_ACCELERATION as dofs
            #         node.AddDof(KratosMultiphysics.ANGULAR_ACCELERATION_X)
            #         node.AddDof(KratosMultiphysics.ANGULAR_ACCELERATION_Y)
            #         node.AddDof(KratosMultiphysics.ANGULAR_ACCELERATION_Z)

        KratosMultiphysics.Logger.PrintInfo("GeoMechanics_U_Pw_Solver", "DOFs added correctly.")


    def Initialize(self):
        KratosMultiphysics.Logger.PrintInfo("::[GeoMechanics_U_Pw_Solver]:: ", "Initialisation ...")

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

        KratosMultiphysics.Logger.PrintInfo("GeoMechanics_U_Pw_Solver", "solver.Initialize is set successfully")

        # Check if everything is assigned correctly
        self.Check()

        KratosMultiphysics.Logger.PrintInfo("GeoMechanics_U_Pw_Solver", "Solver initialization finished.")

    def InitializeSolutionStep(self):
        self.solver.InitializeSolutionStep()

    def Predict(self):
        self.solver.Predict()

    def SolveSolutionStep(self):
        is_converged = self.solver.SolveSolutionStep()
        return is_converged

    def FinalizeSolutionStep(self):
        self.solver.FinalizeSolutionStep()

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
        from KratosMultiphysics.GeoMechanicsApplication import check_and_prepare_model_process_geo
        check_and_prepare_model_process_geo.CheckAndPrepareModelProcess(self.main_model_part, params).Execute()

        # NOTE: We do this here in case the model is empty, so the properties can be assigned
        if not self.model.HasModelPart(self.main_model_part.Name):
            self.model.AddModelPart(self.main_model_part)

        # Import constitutive laws.
        materials_imported = self.import_constitutive_laws()
        if materials_imported:
            KratosMultiphysics.Logger.PrintInfo("::[GeoMechanicalSolver]:: ", "Constitutive law was successfully imported.")
        else:
            # KratosMultiphysics.Logger.PrintInfo("::[GeoMechanicalSolver]:: ", "Constitutive law was not imported.")
            raise Exception("::[GeoMechanicalSolver]:: ", "Constitutive law was not imported.")

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
        if (block_builder):
            builder_and_solver = KratosMultiphysics.ResidualBasedBlockBuilderAndSolver(self.linear_solver)
        else:
            builder_and_solver = KratosMultiphysics.ResidualBasedEliminationBuilderAndSolver(self.linear_solver)

        return builder_and_solver

    def _ConstructScheme(self, scheme_type, solution_type):

        self.main_model_part.ProcessInfo.SetValue(KratosGeo.VELOCITY_COEFFICIENT, 1.0)
        self.main_model_part.ProcessInfo.SetValue(KratosGeo.DT_PRESSURE_COEFFICIENT, 1.0)

        if (scheme_type.lower() == "newmark"):
            beta = self.settings["newmark_beta"].GetDouble()
            gamma = self.settings["newmark_gamma"].GetDouble()
            theta = self.settings["newmark_theta"].GetDouble()
            rayleigh_m = self.settings["rayleigh_m"].GetDouble()
            rayleigh_k = self.settings["rayleigh_k"].GetDouble()
            self.main_model_part.ProcessInfo.SetValue(KratosStructure.RAYLEIGH_ALPHA, rayleigh_m)
            self.main_model_part.ProcessInfo.SetValue(KratosStructure.RAYLEIGH_BETA, rayleigh_k)
            KratosMultiphysics.Logger.PrintInfo("GeoMechanics_U_Pw_Solver, solution_type", solution_type)
            if (solution_type.lower() == "quasi-static" or solution_type.lower() == "quasi_static"):
                if (rayleigh_m < 1.0e-20 and rayleigh_k < 1.0e-20):
                    KratosMultiphysics.Logger.PrintInfo("GeoMechanics_U_Pw_Solver, scheme", "Quasi-UnDamped.")
                    scheme = KratosGeo.NewmarkQuasistaticUPwScheme(beta,gamma,theta)
                else:
                    KratosMultiphysics.Logger.PrintInfo("GeoMechanics_U_Pw_Solver, scheme", "Quasi-Damped.")
                    scheme = KratosGeo.NewmarkQuasistaticDampedUPwScheme(beta,gamma,theta)
            elif (solution_type.lower() == "dynamic"):
                KratosMultiphysics.Logger.PrintInfo("GeoMechanics_U_Pw_Solver, scheme", "Dynamic.")
                scheme = KratosGeo.NewmarkDynamicUPwScheme(beta,gamma,theta)
            elif (solution_type.lower() == "k0-procedure" or solution_type.lower() == "k0_procedure"):
                if (rayleigh_m < 1.0e-20 and rayleigh_k < 1.0e-20):
                    KratosMultiphysics.Logger.PrintInfo("GeoMechanics_U_Pw_Solver, scheme", "Quasi-UnDamped.")
                    scheme = KratosGeo.NewmarkQuasistaticUPwScheme(beta,gamma,theta)
                else:
                    KratosMultiphysics.Logger.PrintInfo("GeoMechanics_U_Pw_Solver, scheme", "Quasi-Damped.")
                    scheme = KratosGeo.NewmarkQuasistaticDampedUPwScheme(beta,gamma,theta)
            else:
              raise Exception("Undefined solution type", solution_type)
        elif (scheme_type.lower() == "backward_euler"or solution_type.lower() == "backward-euler"):
            if (solution_type.lower() == "quasi-static" or solution_type.lower() == "quasi_static"):
                KratosMultiphysics.Logger.PrintInfo("GeoMechanics_U_Pw_Solver, scheme", "Backward Euler.")
                scheme = KratosGeo.BackwardEulerQuasistaticUPwScheme()
            elif (solution_type.lower() == "k0-procedure" or solution_type.lower() == "k0_procedure"):
                KratosMultiphysics.Logger.PrintInfo("GeoMechanics_U_Pw_Solver, scheme", "Backward Euler.")
                scheme = KratosGeo.BackwardEulerQuasistaticUPwScheme()
            else:
              raise Exception("Undefined/uncompatible solution type with Backward Euler", solution_type)
        else:
            raise Exception("Apart from Newmark, other scheme_type are not available.")

        return scheme

    def _ConstructConvergenceCriterion(self, convergence_criterion):

        D_RT = self.settings["displacement_relative_tolerance"].GetDouble()
        D_AT = self.settings["displacement_absolute_tolerance"].GetDouble()
        R_RT = self.settings["residual_relative_tolerance"].GetDouble()
        R_AT = self.settings["residual_absolute_tolerance"].GetDouble()
        echo_level = self.settings["echo_level"].GetInt()

        if(convergence_criterion.lower() == "displacement_criterion"):
            convergence_criterion = KratosMultiphysics.DisplacementCriteria(D_RT, D_AT)
            convergence_criterion.SetEchoLevel(echo_level)
        elif(convergence_criterion.lower() == "residual_criterion"):
            convergence_criterion = KratosMultiphysics.ResidualCriteria(R_RT, R_AT)
            convergence_criterion.SetEchoLevel(echo_level)
        elif(convergence_criterion.lower() == "and_criterion"):
            Displacement = KratosMultiphysics.DisplacementCriteria(D_RT, D_AT)
            Displacement.SetEchoLevel(echo_level)
            Residual = KratosMultiphysics.ResidualCriteria(R_RT, R_AT)
            Residual.SetEchoLevel(echo_level)
            convergence_criterion = KratosMultiphysics.AndCriteria(Residual, Displacement)
        elif(convergence_criterion.lower() == "or_criterion"):
            Displacement = KratosMultiphysics.DisplacementCriteria(D_RT, D_AT)
            Displacement.SetEchoLevel(echo_level)
            Residual = KratosMultiphysics.ResidualCriteria(R_RT, R_AT)
            Residual.SetEchoLevel(echo_level)
            convergence_criterion = KratosMultiphysics.OrCriteria(Residual, Displacement)
        elif(convergence_criterion.lower() == "water_pressure_criterion"):
            convergence_criterion = KratosMultiphysics.MixedGenericCriteria([(KratosMultiphysics.WATER_PRESSURE, D_RT, D_AT)])
            convergence_criterion.SetEchoLevel(echo_level)
        elif(convergence_criterion.lower() == "displacement_and_water_pressure_criterion"):
            convergence_criterion = KratosMultiphysics.MixedGenericCriteria([(KratosMultiphysics.DisplacementCriteria(D_RT, D_AT)),(KratosMultiphysics.WATER_PRESSURE, D_RT, D_AT)])
            convergence_criterion.SetEchoLevel(echo_level)
        else:
            err_msg =  "The requested convergence criterion \"" + convergence_criterion + "\" is not available!\n"
            err_msg += "Available options are: \"displacement_criterion\", \"residual_criterion\", \"and_criterion\", \"or_criterion\", \"water_pressure_criterion\", \"displacement_and_water_pressure_criterion\""
            raise Exception(err_msg)

        return convergence_criterion

    def _ConstructSolver(self, builder_and_solver, strategy_type):

        self.main_model_part.ProcessInfo.SetValue(KratosGeo.IS_CONVERGED, True)
        self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.NL_ITERATION_NUMBER, 1)

        max_iters = self.settings["max_iterations"].GetInt()
        compute_reactions = self.settings["compute_reactions"].GetBool()
        reform_step_dofs = self.settings["reform_dofs_at_each_step"].GetBool()
        move_mesh_flag = self.settings["move_mesh_flag"].GetBool()

        if strategy_type.lower() == "newton_raphson":
            self.strategy_params = KratosMultiphysics.Parameters("{}")
            self.strategy_params.AddValue("loads_sub_model_part_list",self.loads_sub_sub_model_part_list)
            self.strategy_params.AddValue("loads_variable_list",self.settings["loads_variable_list"])
            solving_strategy = KratosGeo.GeoMechanicsNewtonRaphsonStrategy(self.computing_model_part,
                                                                           self.scheme,
                                                                           self.linear_solver,
                                                                           self.convergence_criterion,
                                                                           builder_and_solver,
                                                                           self.strategy_params,
                                                                           max_iters,
                                                                           compute_reactions,
                                                                           reform_step_dofs,
                                                                           move_mesh_flag)
        elif strategy_type.lower() == "line_search":
            self.strategy_params = KratosMultiphysics.Parameters("{}")
            self.strategy_params.AddValue("max_iteration",             self.settings["max_iterations"])
            self.strategy_params.AddValue("compute_reactions",          self.settings["compute_reactions"])
            self.strategy_params.AddValue("max_line_search_iterations", self.settings["max_line_search_iterations"])
            self.strategy_params.AddValue("first_alpha_value",          self.settings["first_alpha_value"])
            self.strategy_params.AddValue("second_alpha_value",         self.settings["second_alpha_value"])
            self.strategy_params.AddValue("min_alpha",                  self.settings["min_alpha"])
            self.strategy_params.AddValue("max_alpha",                  self.settings["max_alpha"])
            self.strategy_params.AddValue("line_search_tolerance",      self.settings["line_search_tolerance"])
            solving_strategy = KratosMultiphysics.LineSearchStrategy(self.computing_model_part,
                                                                     self.scheme,
                                                                     self.linear_solver,
                                                                     self.convergence_criterion,
                                                                     self.strategy_params)

        elif strategy_type.lower() == "arc_length":
            # Arc-Length strategy
            self.main_model_part.ProcessInfo.SetValue(KratosGeo.ARC_LENGTH_LAMBDA,        1.0)
            self.main_model_part.ProcessInfo.SetValue(KratosGeo.ARC_LENGTH_RADIUS_FACTOR, 1.0)
            self.strategy_params = KratosMultiphysics.Parameters("{}")
            self.strategy_params.AddValue("desired_iterations",self.settings["desired_iterations"])
            self.strategy_params.AddValue("max_radius_factor",self.settings["max_radius_factor"])
            self.strategy_params.AddValue("min_radius_factor",self.settings["min_radius_factor"])
            self.strategy_params.AddValue("loads_sub_model_part_list",self.loads_sub_sub_model_part_list)
            self.strategy_params.AddValue("loads_variable_list",self.settings["loads_variable_list"])
            solving_strategy = KratosGeo.GeoMechanicsRammArcLengthStrategy(self.computing_model_part,
                                                                           self.scheme,
                                                                           self.linear_solver,
                                                                           self.convergence_criterion,
                                                                           builder_and_solver,
                                                                           self.strategy_params,
                                                                           max_iters,
                                                                           compute_reactions,
                                                                           reform_step_dofs,
                                                                           move_mesh_flag)
        else:
            raise Exception("Undefined strategy type", strategy_type)

        return solving_strategy

    def _CheckConvergence(self):

        IsConverged = self.solver.IsConverged()

        return IsConverged

    def _UpdateLoads(self):

        self.solver.UpdateLoads()
