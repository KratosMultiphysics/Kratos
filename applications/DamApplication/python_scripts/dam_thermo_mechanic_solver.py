from __future__ import print_function, absolute_import, division # makes KratosMultiphysics backward compatible with python 2.6 and 2.7
#import kratos core and applications
import KratosMultiphysics
import KratosMultiphysics.ConvectionDiffusionApplication as KratosConvDiff
import KratosMultiphysics.SolidMechanicsApplication as KratosSolid
import KratosMultiphysics.PoromechanicsApplication as KratosPoro
import KratosMultiphysics.DamApplication as KratosDam
import json

# Check that KratosMultiphysics was imported in the main script
KratosMultiphysics.CheckForPreviousImport()


def CreateSolver(main_model_part, custom_settings):

    return DamThermoMechanicSolver(main_model_part, custom_settings)


class DamThermoMechanicSolver(object):

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
            "solver_type": "dam_thermo_mechanic_solver",
            "model_import_settings":{
                "input_type": "mdpa",
                "input_filename": "unknown_name",
                "input_file_label": 0
            },
            "echo_level": 0,
            "buffer_size": 2,
            "processes_sub_model_part_list": [""],
            "thermal_solver_settings":{
                "echo_level": 0,
                "reform_dofs_at_each_step": false,
                "clear_storage": false,
                "compute_reactions": false,
                "move_mesh_flag": true,
                "compute_norm_dx_flag": false,
                "theta_scheme": 1.0,
                "block_builder": true,
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
                "thermal_loads_sub_model_part_list": []
            },
            "mechanical_solver_settings":{
                "echo_level": 0,
                "reform_dofs_at_each_step": false,
                "clear_storage": false,
                "compute_reactions": false,
                "move_mesh_flag": true,
                "solution_type": "Quasi-Static",
                "scheme_type": "Newmark",
                "rayleigh_m": 0.0,
                "rayleigh_k": 0.0,
                "strategy_type": "Newton-Raphson",
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
                "body_domain_sub_model_part_list": [],
                "mechanical_loads_sub_model_part_list": [],
                "loads_sub_model_part_list": [],
                "loads_variable_list": []
            }
        }
        """)

        # Overwrite the default settings with user-provided parameters
        self.settings = custom_settings
        self.settings.ValidateAndAssignDefaults(default_settings)

        # Construct the linear solver
        import linear_solver_factory
        self.thermal_linear_solver = linear_solver_factory.ConstructSolver(self.settings["thermal_solver_settings"]["linear_solver_settings"])
        self.mechanical_linear_solver = linear_solver_factory.ConstructSolver(self.settings["mechanical_solver_settings"]["linear_solver_settings"])

        print("Construction of DamThermoMechanicSolver finished")

    def AddVariables(self):

        ## Mechanical Variables
        # Add displacements
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        # Add reactions for the displacements
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.REACTION)
        # Add dynamic variables
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.ACCELERATION)
        # Add variables for the solid conditions
        self.main_model_part.AddNodalSolutionStepVariable(KratosSolid.POINT_LOAD)
        self.main_model_part.AddNodalSolutionStepVariable(KratosSolid.FORCE_LOAD)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.POSITIVE_FACE_PRESSURE)
        # Add volume acceleration
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VOLUME_ACCELERATION)
        # Add variables for post-processing
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.NODAL_AREA)
        self.main_model_part.AddNodalSolutionStepVariable(KratosDam.NODAL_CAUCHY_STRESS_TENSOR)
        self.main_model_part.AddNodalSolutionStepVariable(KratosDam.Vi_POSITIVE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosDam.Viii_POSITIVE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosPoro.NODAL_JOINT_WIDTH)
        self.main_model_part.AddNodalSolutionStepVariable(KratosPoro.NODAL_JOINT_AREA)
        self.main_model_part.AddNodalSolutionStepVariable(KratosDam.NODAL_YOUNG_MODULUS)

        ## Thermal variables
        thermal_settings = KratosMultiphysics.ConvectionDiffusionSettings()
        thermal_settings.SetDiffusionVariable(KratosMultiphysics.CONDUCTIVITY)
        thermal_settings.SetUnknownVariable(KratosMultiphysics.TEMPERATURE)
        thermal_settings.SetSpecificHeatVariable(KratosMultiphysics.SPECIFIC_HEAT)
        thermal_settings.SetDensityVariable(KratosMultiphysics.DENSITY)
        thermal_settings.SetVolumeSourceVariable(KratosMultiphysics.HEAT_FLUX)
        thermal_settings.SetSurfaceSourceVariable(KratosMultiphysics.FACE_HEAT_FLUX)
        self.main_model_part.ProcessInfo.SetValue(KratosMultiphysics.CONVECTION_DIFFUSION_SETTINGS, thermal_settings)

        # Add thermal variables
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.CONDUCTIVITY)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.TEMPERATURE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.SPECIFIC_HEAT)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DENSITY)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.HEAT_FLUX)
        self.main_model_part.AddNodalSolutionStepVariable(KratosMultiphysics.FACE_HEAT_FLUX)
        self.main_model_part.AddNodalSolutionStepVariable(KratosDam.DAM_SURFACE_NODE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosDam.NODAL_REFERENCE_TEMPERATURE)

        self.main_model_part.AddNodalSolutionStepVariable(KratosDam.INITIAL_NODAL_CAUCHY_STRESS_TENSOR)

        # This Variable is used for computing heat source according Azenha Formulation
        self.main_model_part.AddNodalSolutionStepVariable(KratosDam.ALPHA_HEAT_SOURCE)
        self.main_model_part.AddNodalSolutionStepVariable(KratosDam.TIME_ACTIVATION)


        print("Variables correctly added")

    def GetMinimumBufferSize(self):
        return 2

    def AddDofs(self):

        for node in self.main_model_part.Nodes:
            ## Solid dofs
            node.AddDof(KratosMultiphysics.DISPLACEMENT_X,KratosMultiphysics.REACTION_X)
            node.AddDof(KratosMultiphysics.DISPLACEMENT_Y,KratosMultiphysics.REACTION_Y)
            node.AddDof(KratosMultiphysics.DISPLACEMENT_Z,KratosMultiphysics.REACTION_Z)
            ## Thermal dofs
            node.AddDof(KratosMultiphysics.TEMPERATURE)

        if(self.settings["mechanical_solver_settings"]["solution_type"].GetString() == "Dynamic"):
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

        # Set ProcessInfo variables
        self.main_model_part.ProcessInfo.SetValue(KratosConvDiff.THETA, self.settings["thermal_solver_settings"]["theta_scheme"].GetDouble())

        # Get the computing model parts
        self.thermal_computing_model_part = self.GetComputingThermalModelPart()
        self.mechanical_computing_model_part = self.GetComputingModelPart()

        # Builder and solver creation
        thermal_builder_and_solver = self._ConstructBuilderAndSolver(self.settings["thermal_solver_settings"]["block_builder"].GetBool(),self.thermal_linear_solver)
        mechanical_builder_and_solver = self._ConstructBuilderAndSolver(self.settings["mechanical_solver_settings"]["block_builder"].GetBool(),self.mechanical_linear_solver)

        # Solution scheme creation
        thermal_scheme = KratosMultiphysics.ResidualBasedIncrementalUpdateStaticScheme()
        mechanical_scheme = self._ConstructScheme(self.settings["mechanical_solver_settings"]["scheme_type"].GetString(),
                                         self.settings["mechanical_solver_settings"]["solution_type"].GetString())

        # Get the convergence criterion
        convergence_criterion = self._ConstructConvergenceCriterion(self.settings["mechanical_solver_settings"]["convergence_criterion"].GetString())

        # Solver creation
        self.Thermal_Solver = KratosMultiphysics.ResidualBasedLinearStrategy(self.thermal_computing_model_part,
                                                                             thermal_scheme,
                                                                             self.thermal_linear_solver,
                                                                             thermal_builder_and_solver,
                                                                             self.settings["thermal_solver_settings"]["compute_reactions"].GetBool(),
                                                                             self.settings["thermal_solver_settings"]["reform_dofs_at_each_step"].GetBool(),
                                                                             self.settings["thermal_solver_settings"]["compute_norm_dx_flag"].GetBool(),
                                                                             self.settings["thermal_solver_settings"]["move_mesh_flag"].GetBool())
        self.Mechanical_Solver = self._ConstructSolver(mechanical_builder_and_solver,
                                                       mechanical_scheme,
                                                       convergence_criterion,
                                                       self.settings["mechanical_solver_settings"]["strategy_type"].GetString())

        # Set echo_level
        self.Thermal_Solver.SetEchoLevel(self.settings["thermal_solver_settings"]["echo_level"].GetInt())
        self.Mechanical_Solver.SetEchoLevel(self.settings["mechanical_solver_settings"]["echo_level"].GetInt())

        # Check if everything is assigned correctly
        self.Thermal_Solver.Check()
        self.Mechanical_Solver.Check()

        print ("Initialization DamThermoMechanicSolver finished")

    def GetComputingModelPart(self):
        return self.main_model_part.GetSubModelPart(self.mechanical_model_part_name)

    def GetComputingThermalModelPart(self):
        return self.main_model_part.GetSubModelPart(self.thermal_model_part_name)

    def GetOutputVariables(self):
        pass

    def ComputeDeltaTime(self):
        pass

    def SaveRestart(self):
        pass #one should write the restart file here

    def Solve(self):
        if self.settings["thermal_solver_settings"]["clear_storage"].GetBool():
            self.Thermal_Solver.Clear()
        if self.settings["mechanical_solver_settings"]["clear_storage"].GetBool():
            self.Mechanical_Solver.Clear()

        self.Thermal_Solver.Solve()
        self.Mechanical_Solver.Solve()

    # solve :: sequencial calls

    def InitializeStrategy(self):
        if self.settings["thermal_solver_settings"]["clear_storage"].GetBool():
            self.Thermal_Solver.Clear()
        if self.settings["mechanical_solver_settings"]["clear_storage"].GetBool():
            self.Mechanical_Solver.Clear()

        self.Thermal_Solver.Initialize()
        self.Mechanical_Solver.Initialize()

    def InitializeSolutionStep(self):
        self.Thermal_Solver.InitializeSolutionStep()
        self.Mechanical_Solver.InitializeSolutionStep()

    def Predict(self):
        self.Thermal_Solver.Predict()
        self.Mechanical_Solver.Predict()

    def SolveSolutionStep(self):
        self.Thermal_Solver.SolveSolutionStep()
        self.Mechanical_Solver.SolveSolutionStep()

    def FinalizeSolutionStep(self):
        self.Thermal_Solver.FinalizeSolutionStep()
        self.Mechanical_Solver.FinalizeSolutionStep()

    # solve :: sequencial calls

    def SetEchoLevel(self, level):
        self.Thermal_Solver.SetEchoLevel(level)
        self.Mechanical_Solver.SetEchoLevel(level)

    def Clear(self):
        self.Thermal_Solver.Clear()
        self.Mechanical_Solver.Clear()

    def Check(self):
        self.Thermal_Solver.Check()
        self.Mechanical_Solver.Check()

    #### Specific internal functions ####

    def _ExecuteAfterReading(self):

        self.thermal_model_part_name = "thermal_computing_domain"
        self.mechanical_model_part_name = "mechanical_computing_domain"

        # Create list of sub sub model parts (it is a copy of the standard lists with a different name)
        thermal_loads_sub_sub_model_part_list = []
        for i in range(self.settings["thermal_solver_settings"]["problem_domain_sub_model_part_list"].size()):
            thermal_loads_sub_sub_model_part_list.append("sub_"+self.settings["thermal_solver_settings"]["problem_domain_sub_model_part_list"][i].GetString())
        thermal_loads_sub_sub_model_part_list = KratosMultiphysics.Parameters(json.dumps(thermal_loads_sub_sub_model_part_list))

        self.body_domain_sub_sub_model_part_list = []
        for i in range(self.settings["mechanical_solver_settings"]["body_domain_sub_model_part_list"].size()):
            self.body_domain_sub_sub_model_part_list.append("sub_"+self.settings["mechanical_solver_settings"]["body_domain_sub_model_part_list"][i].GetString())
        self.body_domain_sub_sub_model_part_list = KratosMultiphysics.Parameters(json.dumps(self.body_domain_sub_sub_model_part_list))

        self.loads_sub_sub_model_part_list = []
        for i in range(self.settings["mechanical_solver_settings"]["loads_sub_model_part_list"].size()):
            self.loads_sub_sub_model_part_list.append("sub_"+self.settings["mechanical_solver_settings"]["loads_sub_model_part_list"][i].GetString())
        self.loads_sub_sub_model_part_list = KratosMultiphysics.Parameters(json.dumps(self.loads_sub_sub_model_part_list))

        # Auxiliary Kratos parameters object to be called by the CheckAndPepareModelProcessDamThermal
        aux_params = KratosMultiphysics.Parameters("{}")
        aux_params.AddEmptyValue("thermal_model_part_name").SetString(self.thermal_model_part_name)
        aux_params.AddValue("thermal_domain_sub_model_part_list",self.settings["thermal_solver_settings"]["problem_domain_sub_model_part_list"])
        aux_params.AddValue("thermal_loads_sub_model_part_list",self.settings["thermal_solver_settings"]["thermal_loads_sub_model_part_list"])
        aux_params.AddValue("thermal_domain_sub_sub_model_part_list",thermal_loads_sub_sub_model_part_list)

        aux_params.AddEmptyValue("mechanical_model_part_name").SetString(self.mechanical_model_part_name)
        aux_params.AddValue("mechanical_domain_sub_model_part_list",self.settings["mechanical_solver_settings"]["problem_domain_sub_model_part_list"])
        aux_params.AddValue("mechanical_loads_sub_model_part_list",self.settings["mechanical_solver_settings"]["mechanical_loads_sub_model_part_list"])
        aux_params.AddValue("body_domain_sub_model_part_list",self.settings["mechanical_solver_settings"]["body_domain_sub_model_part_list"])
        aux_params.AddValue("body_domain_sub_sub_model_part_list",self.body_domain_sub_sub_model_part_list)
        aux_params.AddValue("loads_sub_model_part_list",self.settings["mechanical_solver_settings"]["loads_sub_model_part_list"])
        aux_params.AddValue("loads_sub_sub_model_part_list",self.loads_sub_sub_model_part_list)

        # CheckAndPrepareModelProcessDamMechanical creates the solid_computational_model_part
        import check_and_prepare_model_process_dam_thermal
        check_and_prepare_model_process_dam_thermal.CheckAndPrepareModelProcessDamThermal(self.main_model_part, aux_params).Execute()

        # Constitutive law import
        import dam_constitutive_law_utility
        dam_constitutive_law_utility.SetConstitutiveLaw(self.main_model_part)

        self.main_model_part.SetBufferSize( self.settings["buffer_size"].GetInt() )
        minimum_buffer_size = self.GetMinimumBufferSize()
        if(minimum_buffer_size > self.main_model_part.GetBufferSize()):
            self.main_model_part.SetBufferSize( minimum_buffer_size )

    def _ConstructBuilderAndSolver(self, block_builder, linear_solver):

        # Creating the builder and solver
        if(block_builder):
            builder_and_solver = KratosMultiphysics.ResidualBasedBlockBuilderAndSolver(linear_solver)
        else:
            builder_and_solver = KratosMultiphysics.ResidualBasedEliminationBuilderAndSolver(linear_solver)

        return builder_and_solver

    def _ConstructScheme(self, scheme_type, solution_type):

        rayleigh_m = self.settings["mechanical_solver_settings"]["rayleigh_m"].GetDouble()
        rayleigh_k = self.settings["mechanical_solver_settings"]["rayleigh_k"].GetDouble()

        if(solution_type == "Quasi-Static"):
            if(rayleigh_m<1.0e-20 and rayleigh_k<1.0e-20):
                scheme =  KratosDam.IncrementalUpdateStaticSmoothingScheme()
            else:
                scheme =  KratosDam.IncrementalUpdateStaticDampedSmoothingScheme(rayleigh_m,rayleigh_k)
        else:
            if(scheme_type == "Newmark"):
                damp_factor_m = 0.0
            else:
                damp_factor_m = -0.01
            scheme = KratosDam.BossakDisplacementSmoothingScheme(damp_factor_m,rayleigh_m,rayleigh_k)

        return scheme

    def _ConstructConvergenceCriterion(self, convergence_criterion):

        D_RT = self.settings["mechanical_solver_settings"]["displacement_relative_tolerance"].GetDouble()
        D_AT = self.settings["mechanical_solver_settings"]["displacement_absolute_tolerance"].GetDouble()
        R_RT = self.settings["mechanical_solver_settings"]["residual_relative_tolerance"].GetDouble()
        R_AT = self.settings["mechanical_solver_settings"]["residual_absolute_tolerance"].GetDouble()
        echo_level = self.settings["mechanical_solver_settings"]["echo_level"].GetInt()

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

        nonlocal_damage = self.settings["mechanical_solver_settings"]["nonlocal_damage"].GetBool()
        max_iters = self.settings["mechanical_solver_settings"]["max_iteration"].GetInt()
        compute_reactions = self.settings["mechanical_solver_settings"]["compute_reactions"].GetBool()
        reform_step_dofs = self.settings["mechanical_solver_settings"]["reform_dofs_at_each_step"].GetBool()
        move_mesh_flag = self.settings["mechanical_solver_settings"]["move_mesh_flag"].GetBool()

        if strategy_type == "Newton-Raphson":
            if nonlocal_damage:
                self.strategy_params = KratosMultiphysics.Parameters("{}")
                self.strategy_params.AddValue("loads_sub_model_part_list",self.loads_sub_sub_model_part_list)
                self.strategy_params.AddValue("loads_variable_list",self.settings["mechanical_solver_settings"]["loads_variable_list"])
                self.strategy_params.AddValue("body_domain_sub_model_part_list",self.body_domain_sub_sub_model_part_list)
                self.strategy_params.AddValue("characteristic_length",self.settings["mechanical_solver_settings"]["characteristic_length"])
                self.strategy_params.AddValue("search_neighbours_step",self.settings["mechanical_solver_settings"]["search_neighbours_step"])
                solver = KratosPoro.PoromechanicsNewtonRaphsonNonlocalStrategy(self.mechanical_computing_model_part,
                                                                               scheme,
                                                                               self.mechanical_linear_solver,
                                                                               convergence_criterion,
                                                                               builder_and_solver,
                                                                               self.strategy_params,
                                                                               max_iters,
                                                                               compute_reactions,
                                                                               reform_step_dofs,
                                                                               move_mesh_flag)
            else:
                self.main_model_part.ProcessInfo.SetValue(KratosPoro.IS_CONVERGED, True)
                solver = KratosMultiphysics.ResidualBasedNewtonRaphsonStrategy(self.mechanical_computing_model_part,
                                                                                    scheme,
                                                                                    self.mechanical_linear_solver,
                                                                                    convergence_criterion,
                                                                                    builder_and_solver,
                                                                                    max_iters,
                                                                                    compute_reactions,
                                                                                    reform_step_dofs,
                                                                                    move_mesh_flag)
        else:
            # Arc-Length strategy
            self.strategy_params = KratosMultiphysics.Parameters("{}")
            self.strategy_params.AddValue("desired_iterations",self.settings["mechanical_solver_settings"]["desired_iterations"])
            self.strategy_params.AddValue("max_radius_factor",self.settings["mechanical_solver_settings"]["max_radius_factor"])
            self.strategy_params.AddValue("min_radius_factor",self.settings["mechanical_solver_settings"]["min_radius_factor"])
            self.strategy_params.AddValue("loads_sub_model_part_list",self.loads_sub_sub_model_part_list)
            self.strategy_params.AddValue("loads_variable_list",self.settings["mechanical_solver_settings"]["loads_variable_list"])
            if nonlocal_damage:
                self.strategy_params.AddValue("body_domain_sub_model_part_list",self.body_domain_sub_sub_model_part_list)
                self.strategy_params.AddValue("characteristic_length",self.settings["mechanical_solver_settings"]["characteristic_length"])
                self.strategy_params.AddValue("search_neighbours_step",self.settings["mechanical_solver_settings"]["search_neighbours_step"])
                solver = KratosPoro.PoromechanicsRammArcLengthNonlocalStrategy(self.mechanical_computing_model_part,
                                                                               scheme,
                                                                               self.mechanical_linear_solver,
                                                                               convergence_criterion,
                                                                               builder_and_solver,
                                                                               self.strategy_params,
                                                                               max_iters,
                                                                               compute_reactions,
                                                                               reform_step_dofs,
                                                                               move_mesh_flag)
            else:
                solver = KratosPoro.PoromechanicsRammArcLengthStrategy(self.mechanical_computing_model_part,
                                                                       scheme,
                                                                       self.mechanical_linear_solver,
                                                                       convergence_criterion,
                                                                       builder_and_solver,
                                                                       self.strategy_params,
                                                                       max_iters,
                                                                       compute_reactions,
                                                                       reform_step_dofs,
                                                                       move_mesh_flag)

        return solver

    def _CheckConvergence(self):

        IsConverged = self.Mechanical_Solver.IsConverged()

        return IsConverged

    def _UpdateLoads(self):

        self.Mechanical_Solver.UpdateLoads()
