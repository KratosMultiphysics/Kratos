from __future__ import print_function, absolute_import, division # makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics
from KratosMultiphysics.python_solver import PythonSolver

# Import applications
import KratosMultiphysics.FluidDynamicsApplication
import KratosMultiphysics.ConvectionDiffusionApplication as KratosConvDiff
import KratosMultiphysics.FluidTransportApplication as KratosFluidTransport

def CreateSolver(main_model_part, custom_settings):

    return FluidTransportSolver(main_model_part, custom_settings)


class FluidTransportSolver(PythonSolver):

    def __init__(self, model, custom_settings):

        self._validate_settings_in_baseclass=True # To be removed eventually

        super(FluidTransportSolver,self).__init__(model, custom_settings)

        # There is only a single rank in OpenMP, we always print
        self._is_printing_rank = True

        self.min_buffer_size = 2

        # Either retrieve the model part from the model or create a new one
        model_part_name = self.settings["model_part_name"].GetString()

        if model_part_name == "":
            raise Exception('Please specify a model_part name!')

        if self.model.HasModelPart(model_part_name):
            self.main_model_part = self.model.GetModelPart(model_part_name)
        else:
            self.main_model_part = self.model.CreateModelPart(model_part_name)

        self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE,
                                                  self.settings["domain_size"].GetInt())

        KratosMultiphysics.Logger.PrintInfo("FluidTransportSolver", "Construction of FluidTransportSolver finished.")

    @classmethod
    def GetDefaultSettings(cls):
        this_defaults = KratosMultiphysics.Parameters("""{
            "solver_type": "fluid_transport_solver",
            "model_part_name": "FluidTransportDomain",
            "domain_size": 2,
            "start_time": 0.0,
            "time_step": 0.1,
            "model_import_settings":{
                "input_type": "mdpa",
                "input_filename": "unknown_name",
                "input_file_label": 0
            },
            "buffer_size":                        2,
            "echo_level":                         0,
            "clear_storage":                      false,
            "compute_reactions":                  false,
            "move_mesh_flag":                     false,
            "reform_dofs_at_each_step":           false,
            "block_builder":                      true,
            "solution_type":                      "Steady",
            "scheme_type":                        "Implicit",
            "newmark_theta":                      0.5,
            "strategy_type":                      "Linear",
            "convergence_criterion":              "And_criterion",
            "displacement_relative_tolerance":    1.0E-4,
            "displacement_absolute_tolerance":    1.0E-9,
            "residual_relative_tolerance":        1.0E-4,
            "residual_absolute_tolerance":        1.0E-9,
            "max_iteration":                      15,
            "linear_solver_settings":             {
                "solver_type":   "ExternalSolversApplication.super_lu",
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
            "pfem2_convection_settings"    : {
                "use_pfem2_convection"         : false
	        }
        }""")

        this_defaults.AddMissingParameters(super(FluidTransportSolver, cls).GetDefaultSettings())
        return this_defaults

    def AddVariables(self):

        ## ConvectionDiffusionSettings
        thermal_settings = KratosMultiphysics.ConvectionDiffusionSettings()
        thermal_settings.SetReactionVariable(KratosMultiphysics.ABSORPTION_COEFFICIENT)
        thermal_settings.SetDiffusionVariable(KratosMultiphysics.CONDUCTIVITY)

        if(self.settings["solution_type"].GetString() == "Steady"):
            thermal_settings.SetUnknownVariable(KratosMultiphysics.TEMPERATURE)
        elif(self.settings["scheme_type"].GetString() == "Implicit"):
            thermal_settings.SetUnknownVariable(KratosFluidTransport.PHI_THETA)
        else:
            thermal_settings.SetUnknownVariable(KratosMultiphysics.TEMPERATURE)

        thermal_settings.SetSpecificHeatVariable(KratosMultiphysics.SPECIFIC_HEAT)
        thermal_settings.SetDensityVariable(KratosMultiphysics.DENSITY)
        thermal_settings.SetVolumeSourceVariable(KratosMultiphysics.HEAT_FLUX)
        thermal_settings.SetSurfaceSourceVariable(KratosMultiphysics.FACE_HEAT_FLUX)
        thermal_settings.SetMeshVelocityVariable(KratosMultiphysics.MESH_VELOCITY)
        thermal_settings.SetVelocityVariable(KratosMultiphysics.VELOCITY)

        if self.settings["pfem2_convection_settings"]["use_pfem2_convection"].GetBool():
            thermal_settings.SetProjectionVariable(KratosConvDiff.PROJECTED_SCALAR1) # Required by PFEM2 convection
            self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.YP)

        ## ConvectionDiffusionSettings Variable
        self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.CONVECTION_DIFFUSION_SETTINGS, thermal_settings)

        ## Convection Variables
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.MESH_VELOCITY)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NORMAL)
        self.main_model_part.AddNodalSolutionStepVariable(KratosConvDiff.PROJECTED_SCALAR1) # Required by PFEM2 convection
        self.main_model_part.AddNodalSolutionStepVariable(KratosConvDiff.DELTA_SCALAR1) # Required by PFEM2 convection
        self.main_model_part.AddNodalSolutionStepVariable(KratosConvDiff.MEAN_SIZE) # Required by PFEM2 convection

        # Add thermal variables
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.TEMPERATURE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION_FLUX)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.HEAT_FLUX)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.FACE_HEAT_FLUX)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NODAL_AREA)
        self.main_model_part.AddNodalSolutionStepVariable(KratosFluidTransport.PHI_THETA) # Phi variable refering to the n+theta step
        self.main_model_part.AddNodalSolutionStepVariable(KratosFluidTransport.NODAL_PHI_GRADIENT)
        self.main_model_part.AddNodalSolutionStepVariable(KratosFluidTransport.NODAL_ANALYTIC_SOLUTION)

        print("Variables correctly added")

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

        if not self.main_model_part.ProcessInfo[KratosMultiphysics.IS_RESTARTED]:
            ## Executes the check and prepare model process (Create computing_model_part and set constitutive law)
            self._ExecuteCheckAndPrepare()
            ## Set buffer size
            self._SetBufferSize()

        if self._is_printing_rank:
            KratosMultiphysics.Logger.PrintInfo("FluidTransportSolver", "Model reading finished.")

    def AddDofs(self):

        for node in self.main_model_part.Nodes:
            ## Fluid dofs
            #node.AddDof(KratosMultiphysics.VELOCITY_X,KratosMultiphysics.REACTION_X)
            #node.AddDof(KratosMultiphysics.VELOCITY_Y,KratosMultiphysics.REACTION_Y)
            #node.AddDof(KratosMultiphysics.VELOCITY_Z,KratosMultiphysics.REACTION_Z)

            ## Thermal dofs

            if(self.settings["solution_type"].GetString() == "Steady"):
                node.AddDof(KratosMultiphysics.TEMPERATURE, KratosMultiphysics.REACTION_FLUX)
            elif(self.settings["scheme_type"].GetString() == "Implicit"):
                node.AddDof(KratosFluidTransport.PHI_THETA, KratosMultiphysics.REACTION_FLUX)
            else:
                node.AddDof(KratosMultiphysics.TEMPERATURE, KratosMultiphysics.REACTION_FLUX)

        if self._is_printing_rank:
            KratosMultiphysics.Logger.PrintInfo("FluidTransportSolver", "DOFs added correctly.")

    def GetMinimumBufferSize(self):
        return self.min_buffer_size

    def Initialize(self):

        # Set ProcessInfo variables
        self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.NL_ITERATION_NUMBER, 1)

        # Get the computing model parts
        self.computing_model_part = self.GetComputingModelPart()

        # Fill the previous steps of the buffer with the initial conditions
        self._FillBuffer()

        # Construct the linear solver
        self.linear_solver = self._ConstructLinearSolver()

        # Builder and solver creation
        builder_and_solver = self._ConstructBuilderAndSolver(self.settings["block_builder"].GetBool())

        # Solution scheme creation
        self.scheme = self._ConstructScheme(self.settings["solution_type"].GetString())

        # Get the convergence criterion
        self.convergence_criterion = self._ConstructConvergenceCriterion(self.settings["convergence_criterion"].GetString())

        # Solver creation
        self.Solver = self._ConstructSolver(builder_and_solver,
                                            self.scheme,
                                            self.settings["strategy_type"].GetString())

        # Set echo_level
        self.Solver.SetEchoLevel(self.settings["echo_level"].GetInt())

        # Initialize Strategy
        if self.settings["clear_storage"].GetBool():
            self.Clear()

        self.domain_size = self.main_model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE]

        # Calculate Nodal Area
        self.nodal_area_process = KratosMultiphysics.CalculateNodalAreaProcess(self.main_model_part, self.domain_size)
        self.nodal_area_process.Execute()

        # KratosMultiphysics.BodyNormalCalculationUtils().CalculateBodyNormals(self.main_model_part, self.domain_size)

        KratosMultiphysics.BodyNormalCalculationUtils().CalculateBodyNormals(self.main_model_part, self.domain_size)

        # Check if everything is assigned correctly
        self.Solver.Check()
        self._GetParticlesStage().Check()

        self._GetParticlesStage().ExecuteBeforeSolutionLoop()

        KratosMultiphysics.Logger.PrintInfo("FluidTransportSolver", "Solver initialization finished.")

    def GetComputingModelPart(self):
        return self.main_model_part.GetSubModelPart(self.computing_model_part_name)

    def ComputeDeltaTime(self):
        return self.main_model_part.ProcessInfo[KratosMultiphysics.DELTA_TIME]

    def Clear(self):
        self.Solver.Clear()

    def Check(self):
        self.Solver.Check()

    def AdvanceInTime(self, current_time):
        dt = self.ComputeDeltaTime()
        new_time = current_time + dt

        self.main_model_part.CloneTimeStep(new_time)
        self.main_model_part.ProcessInfo[KratosMultiphysics.STEP] += 1

        return new_time

    def InitializeSolutionStep(self):
        self._GetParticlesStage().ExecuteInitializeSolutionStep()
        self.Solver.InitializeSolutionStep()

    def Predict(self):
        self.Solver.Predict()

    def SolveSolutionStep(self):
        is_converged = self.Solver.SolveSolutionStep()
        return is_converged

    def FinalizeSolutionStep(self):
        self.Solver.FinalizeSolutionStep()
        self._GetParticlesStage().ExecuteFinalizeSolutionStep()

    def Solve(self):
        message = "".join([
            "Calling FluidTransportSolver.Solve() method, which is deprecated\n",
            "Please call the individual methods instead:\n",
            "solver.InitializeSolutionStep()\n",
            "solver.Predict()\n",
            "solver.SolveSolutionStep()\n",
            "solver.FinalizeSolutionStep()\n"]
        )
        KratosMultiphysics.Logger.PrintWarning("FluidTransportSolver",message)

        if self.settings["clear_storage"].GetBool():
            self.Clear()

        self.InitializeSolutionStep()
        self.Predict()
        self.SolveSolutionStep()
        self.FinalizeSolutionStep()

    #### Specific internal functions ####

    def _ExecuteCheckAndPrepare(self):

        self.computing_model_part_name = "fluid_transport_computing_domain"

        # Auxiliary Kratos parameters object to be called by the CheckAndPepareModelProcess
        aux_params = KratosMultiphysics.Parameters("{}")
        aux_params.AddEmptyValue("computing_model_part_name").SetString(self.computing_model_part_name)
        aux_params.AddValue("problem_domain_sub_model_part_list",self.settings["problem_domain_sub_model_part_list"])
        aux_params.AddValue("processes_sub_model_part_list",self.settings["processes_sub_model_part_list"])

        # CheckAndPrepareModelProcess creates the solid_computational_model_part
        from KratosMultiphysics.FluidTransportApplication import check_and_prepare_model_process_fluid_transport
        check_and_prepare_model_process_fluid_transport.CheckAndPrepareModelProcess(self.main_model_part, aux_params).Execute()

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
        linear_solver = linear_solver_factory.ConstructSolver(self.settings["linear_solver_settings"])
        return linear_solver

    def _ConstructBuilderAndSolver(self, block_builder):

        # Creating the builder and solver
        if(block_builder):
            builder_and_solver = KratosMultiphysics.ResidualBasedBlockBuilderAndSolver(self.linear_solver)
        else:
            builder_and_solver = KratosMultiphysics.ResidualBasedEliminationBuilderAndSolver(self.linear_solver)

        return builder_and_solver

    def _ConstructScheme(self, solution_type):

        # Creating the builder and solver
        if(solution_type == "Steady"):
            scheme = KratosMultiphysics.ResidualBasedIncrementalUpdateStaticScheme()
        elif(self.settings["scheme_type"].GetString() == "Implicit"):
            theta = self.settings["newmark_theta"].GetDouble()
            scheme = KratosFluidTransport.GeneralizedNewmarkGN11Scheme(theta)
        else:
            theta = 1.0
            scheme = KratosFluidTransport.ExplicitForwardEulerScheme(theta)

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

    def _ConstructSolver(self, builder_and_solver, scheme, strategy_type):

        compute_reactions = self.settings["compute_reactions"].GetBool()
        reform_step_dofs = self.settings["reform_dofs_at_each_step"].GetBool()
        move_mesh_flag = self.settings["move_mesh_flag"].GetBool()

        if strategy_type == "Newton-Raphson":

            solver = KratosMultiphysics.ResidualBasedNewtonRaphsonStrategy(self.computing_model_part,
                                                                            self.scheme,
                                                                            self.linear_solver,
                                                                            self.convergence_criterion,
                                                                            builder_and_solver,
                                                                            self.settings["max_iteration"].GetInt(),
                                                                            compute_reactions,
                                                                            reform_step_dofs,
                                                                            move_mesh_flag)
        else:
            compute_norm_dx_flag = False

            solver = KratosMultiphysics.ResidualBasedLinearStrategy(self.computing_model_part,
                                                                        self.scheme,
                                                                        self.linear_solver,
                                                                        builder_and_solver,
                                                                        compute_reactions,
                                                                        reform_step_dofs,
                                                                        compute_norm_dx_flag,
                                                                        move_mesh_flag)

        return solver

    def _GetParticlesStage(self):
        if not hasattr(self, '_particles_stage'):
            self._particles_stage = self._CreateParticlesStage()
        return self._particles_stage

    def _CreateParticlesStage(self):
        if self.settings["pfem2_convection_settings"]["use_pfem2_convection"].GetBool():
            convection_settings = KratosMultiphysics.Parameters("""{"Parameters" : {}}""")
            convection_settings["Parameters"] = self.settings["pfem2_convection_settings"].Clone()
            convection_settings["Parameters"].RemoveValue("use_pfem2_convection")
            convection_settings["Parameters"].AddValue("model_part_name", self.settings["model_part_name"])
            convection_settings["Parameters"].AddValue("crank_nicolson_theta", self.settings["newmark_theta"])
            import KratosMultiphysics.FluidTransportApplication.pfem2_fluid_transport_process as module
            return module.Factory(convection_settings, self.model)
        else:
            return KratosMultiphysics.Process()
