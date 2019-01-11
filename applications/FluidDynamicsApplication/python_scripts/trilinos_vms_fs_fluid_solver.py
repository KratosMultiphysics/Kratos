from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# -*- coding: utf-8 -*-
from KratosMultiphysics import *
from KratosMultiphysics.mpi import *
from KratosMultiphysics.FluidDynamicsApplication import *
from KratosMultiphysics.MetisApplication import *
from KratosMultiphysics.TrilinosApplication import *

def AddVariables(model_part, config=None):
    model_part.AddNodalSolutionStepVariable(VELOCITY)
    model_part.AddNodalSolutionStepVariable(FRACT_VEL)
    model_part.AddNodalSolutionStepVariable(MESH_VELOCITY)
    model_part.AddNodalSolutionStepVariable(PRESSURE)
    model_part.AddNodalSolutionStepVariable(PRESSURE_OLD_IT)
    model_part.AddNodalSolutionStepVariable(PRESS_PROJ)
    model_part.AddNodalSolutionStepVariable(CONV_PROJ)
    model_part.AddNodalSolutionStepVariable(DIVPROJ)
    model_part.AddNodalSolutionStepVariable(NODAL_AREA)
    model_part.AddNodalSolutionStepVariable(BODY_FORCE)
    model_part.AddNodalSolutionStepVariable(DENSITY)
    model_part.AddNodalSolutionStepVariable(VISCOSITY)
    model_part.AddNodalSolutionStepVariable(FLAG_VARIABLE)
    model_part.AddNodalSolutionStepVariable(IS_STRUCTURE)
    model_part.AddNodalSolutionStepVariable(PARTITION_INDEX)
    model_part.AddNodalSolutionStepVariable(REACTION)
    model_part.AddNodalSolutionStepVariable(NORMAL)
    model_part.AddNodalSolutionStepVariable(Y_WALL)
    model_part.AddNodalSolutionStepVariable(PATCH_INDEX)

    if config is not None:
        if hasattr(config, "TurbulenceModel"):
            if config.TurbulenceModel == "Spalart-Allmaras":
                model_part.AddNodalSolutionStepVariable(TURBULENT_VISCOSITY)
                model_part.AddNodalSolutionStepVariable(MOLECULAR_VISCOSITY)
                model_part.AddNodalSolutionStepVariable(TEMP_CONV_PROJ)
                model_part.AddNodalSolutionStepVariable(DISTANCE)
    mpi.world.barrier()
    if mpi.rank == 0:
        print("variables for the trilinos fractional step solver added correctly")


def AddDofs(model_part, config=None):

    for node in model_part.Nodes:
        # adding dofs
        node.AddDof(PRESSURE)
        node.AddDof(VELOCITY_X, REACTION_X)
        node.AddDof(VELOCITY_Y, REACTION_Y)
        node.AddDof(VELOCITY_Z, REACTION_Z)

    if config is not None:
        if hasattr(config, "TurbulenceModel"):
            if config.TurbulenceModel == "Spalart-Allmaras":
                for node in model_part.Nodes:
                    node.AddDof(TURBULENT_VISCOSITY)

    mpi.world.barrier()
    if mpi.rank == 0:
        print("dofs for the trilinos fractional step solver added correctly")


class IncompressibleFluidSolver:

    def __init__(self, model_part, domain_size):

        # neighbour search
        number_of_avg_elems = 10
        number_of_avg_nodes = 10
        self.neighbour_search = FindNodalNeighboursProcess(
            model_part, number_of_avg_elems, number_of_avg_nodes)

        self.model_part = model_part
        self.domain_size = domain_size

        # assignation of parameters to be used
        self.vel_toll = float(0.001)
        self.press_toll = float(0.001)
        self.max_vel_its = 6
        self.max_press_its = 3
        self.time_order = 2
        self.compute_reactions = False
        self.ReformDofAtEachIteration = False
        self.CalculateNormDxFlag = True
        self.predictor_corrector = False
        self.use_dt_in_stabilization = False

        self.echo_level = 0

        self.step = 1
        self.projections_are_initialized = False

        self.dynamic_tau = 0.001
        self.activate_tau2 = False

        self.Comm = CreateCommunicator()

#
        velocity_nit_max = 100
        velocity_linear_tol = 1e-6

        import MonolithicMultiLevelSolver
        self.velocity_linear_solver = MonolithicMultiLevelSolver.LinearSolver(
            velocity_linear_tol, velocity_nit_max)

        pressure_nit_max = 100
        pressure_linear_tol = 1e-6

        import PressureMultiLevelSolver
        self.pressure_linear_solver = PressureMultiLevelSolver.MultilevelLinearSolver(
            pressure_linear_tol, pressure_nit_max)

        self.use_slip_conditions = False

        self.use_spalart_allmaras = False
        self.use_des = False
        self.Cdes = 1.0
        self.wall_nodes = list()
        self.spalart_allmaras_linear_solver = None

        self.divergence_clearance_steps = 0

    def Initialize(self):
        (self.neighbour_search).Execute()

        self.model_part.ProcessInfo.SetValue(DYNAMIC_TAU, self.dynamic_tau)

        # check if slip conditions are defined
        slip_cond_count = 0
        if self.use_slip_conditions == False:
            for cond in self.model_part.Conditions:
                if cond.GetValue(IS_STRUCTURE) != 0.0:
                    slip_cond_count += 1
                    break
        all_results = mpi.allgather(mpi.world, slip_cond_count)
        for count in all_results:
            slip_cond_count += count
        if slip_cond_count > 0:
            self.use_slip_conditions = True

        # if we use slip conditions, calculate normals on the boundary
        if self.use_slip_conditions:
            self.normal_util = NormalCalculationUtils()
            self.normal_util.CalculateOnSimplex(
                self.model_part,
                self.domain_size,
                IS_STRUCTURE)

        MoveMeshFlag = False

        self.solver_settings = TrilinosFractionalStepSettingsPeriodic(
            self.Comm,
            self.model_part,
            self.domain_size,
            self.time_order,
            self.use_slip_conditions,
            MoveMeshFlag,
            self.ReformDofAtEachIteration,
            PATCH_INDEX)

        self.solver_settings.SetStrategy(TrilinosStrategyLabel.Velocity,
                                         self.velocity_linear_solver,
                                         self.vel_toll,
                                         self.max_vel_its)

        self.solver_settings.SetStrategy(TrilinosStrategyLabel.Pressure,
                                         self.pressure_linear_solver,
                                         self.press_toll,
                                         self.max_press_its)

        if self.use_spalart_allmaras:
            for node in self.wall_nodes:
                node.SetValue(IS_VISITED, 1.0)
                node.SetSolutionStepValue(DISTANCE, 0, 0.0)

            if(self.domain_size == 2):
                self.redistance_utils = ParallelDistanceCalculator2D()
            else:
                self.redistance_utils = ParallelDistanceCalculator3D()

            max_levels = 100
            max_distance = 1000
            self.redistance_utils.CalculateDistancesLagrangianSurface(
                self.model_part,
                DISTANCE,
                NODAL_AREA,
                max_levels,
                max_distance)

            non_linear_tol = 0.001
            max_it = 10
            reform_dofset = self.ReformDofAtEachIteration
            time_order = self.time_order

            if self.spalart_allmaras_linear_solver is None:
                turb_aztec_parameters = ParameterList()
                turb_aztec_parameters.set("AZ_solver", "AZ_gmres")
                turb_aztec_parameters.set("AZ_kspace", 100)
                turb_aztec_parameters.set("AZ_output", "AZ_none")

                turb_preconditioner_type = "ILU"
                turb_preconditioner_parameters = ParameterList()
                turb_overlap_level = 0
                turb_nit_max = 1000
                turb_linear_tol = 1e-9

                self.spalart_allmaras_linear_solver = AztecSolver(
                    turb_aztec_parameters,
                    turb_preconditioner_type,
                    turb_preconditioner_parameters,
                    turb_linear_tol,
                    turb_nit_max,
                    turb_overlap_level)
                self.spalart_allmaras_linear_solver.SetScalingType(
                    AztecScalingType.LeftScaling)

            turbulence_model = TrilinosSpalartAllmarasTurbulenceModel(
                self.Comm,
                self.model_part,
                self.spalart_allmaras_linear_solver,
                self.domain_size,
                non_linear_tol,
                max_it,
                reform_dofset,
                time_order)

            if self.use_des:
                self.turbulence_model.ActivateDES(self.Cdes)

            self.solver_settings.SetTurbulenceModel(
                TrilinosTurbulenceModelLabel.SpalartAllmaras,
                self.spalart_allmaras_linear_solver,
                non_linear_tol,
                max_it)

        self.solver_settings.SetEchoLevel(self.echo_level)

        self.solver = TrilinosFSStrategy(self.model_part,
                                         self.solver_settings,
                                         self.predictor_corrector,
                                         PATCH_INDEX)

# generating the slip conditions
# self.create_slip_conditions.Execute()
# (self.solver).SetSlipProcess(self.create_slip_conditions);
# self.slip_conditions_initialized = True
# (self.solver).SetEchoLevel(self.echo_level)
# (self.solver).SetEchoLevel(self.echo_level)
        print("finished initialization of the fluid strategy")

    def Solve(self):
        if self.divergence_clearance_steps > 0:
            if mpi.rank == 0:
                print("Calculating divergence-free initial condition")
            ## initialize with a Stokes solution step
            #stokes_aztec_parameters = ParameterList()
            #stokes_aztec_parameters.set("AZ_solver", "AZ_gmres")
            #stokes_aztec_parameters.set("AZ_kspace", 100)
            #stokes_aztec_parameters.set("AZ_output", "AZ_none")

            #stokes_preconditioner_type = "ILU"
            #stokes_preconditioner_parameters = ParameterList()
            #stokes_overlap_level = 0
            #stokes_nit_max = 1000
            #stokes_linear_tol = 1e-9

            #stokes_linear_solver = AztecSolver(stokes_aztec_parameters,
                                               #stokes_preconditioner_type,
                                               #stokes_preconditioner_parameters,
                                               #stokes_linear_tol,
                                               #stokes_nit_max,
                                               #stokes_overlap_level)
            #stokes_linear_solver.SetScalingType(AztecScalingType.LeftScaling)
            #stokes_process = TrilinosStokesInitializationProcess(
                #self.Comm,
                #self.model_part,
                #stokes_linear_solver,
                #self.domain_size,
                #PATCH_INDEX)
            ## copy periodic conditions to Stokes problem
            #stokes_process.SetConditions(self.model_part.Conditions)
            ## execute Stokes process
            #stokes_process.Execute()
            #stokes_process = None

            #for node in self.model_part.Nodes:
                #node.SetSolutionStepValue(PRESSURE, 0, 0.0)



            dt = self.model_part.ProcessInfo.GetValue(DELTA_TIME)
            self.model_part.ProcessInfo.SetValue(DELTA_TIME, 100.0 * dt)

            while self.divergence_clearance_steps > 0:
                if self.divergence_clearance_steps > 1:
                    for node in self.model_part.Nodes:
                        node.SetSolutionStepValue(PRESSURE, 0, 0.0)
                        vel = node.GetSolutionStepValue(VELOCITY, 0)
                        node.SetSolutionStepValue(VELOCITY, 1, vel)
                        node.SetSolutionStepValue(VELOCITY, 2, vel)

                self.divergence_clearance_steps -= 1

                (self.solver).Solve()

            self.model_part.ProcessInfo.SetValue(DELTA_TIME, dt)

            self.divergence_clearance_steps = 0
            if mpi.rank == 0:
                print("Finished divergence clearance")

        if(self.ReformDofAtEachIteration):
            self.slip_conditions_initialized = False
            (self.neighbour_search).Execute()
# if(self.slip_conditions_initialized == False):
# self.create_slip_conditions.Execute()
# (self.solver).SetSlipProcess(self.create_slip_conditions);
# self.slip_conditions_initialized = True

        (self.solver).Solve()

        if(self.compute_reactions):
            self.solver.CalculateReactions()

    def Clear(self):
        (self.solver).Clear()
        self.slip_conditions_initialized = True

    def WriteRestartFile(self, FileName):
        backupfile = open(FileName + ".py", 'w')
        backupfile.write("from KratosMultiphysics import *\n")
        backupfile.write("def Restart(NODES):\n")

        import restart_utilities
        restart_utilities.PrintRestart_ScalarVariable_PyFormat(
            VELOCITY_X,
            "VELOCITY_X",
            self.model_part.Nodes,
            backupfile)
        restart_utilities.PrintRestart_ScalarVariable_PyFormat(
            VELOCITY_Y,
            "VELOCITY_Y",
            self.model_part.Nodes,
            backupfile)
        restart_utilities.PrintRestart_ScalarVariable_PyFormat(
            VELOCITY_Z,
            "VELOCITY_Z",
            self.model_part.Nodes,
            backupfile)
        restart_utilities.PrintRestart_ScalarVariable_PyFormat(
            PRESSURE, "PRESSURE", self.model_part.Nodes, backupfile)
        restart_utilities.PrintRestart_ScalarVariable_PyFormat(
            DENSITY, "DENSITY", self.model_part.Nodes, backupfile)
        restart_utilities.PrintRestart_ScalarVariable_PyFormat(
            VISCOSITY, "VISCOSITY", self.model_part.Nodes, backupfile)

        restart_utilities.PrintRestartFixity_PyFormat(
            VELOCITY_X,
            "VELOCITY_X",
            self.model_part.Nodes,
            backupfile)
        restart_utilities.PrintRestartFixity_PyFormat(
            VELOCITY_Y,
            "VELOCITY_Y",
            self.model_part.Nodes,
            backupfile)
        restart_utilities.PrintRestartFixity_PyFormat(
            VELOCITY_Z,
            "VELOCITY_Z",
            self.model_part.Nodes,
            backupfile)
        restart_utilities.PrintRestartFixity_PyFormat(
            PRESSURE,
            "PRESSURE",
            self.model_part.Nodes,
            backupfile)

        backupfile.close()

    #
    def activate_smagorinsky(self, C):
        for elem in self.model_part.Elements:
            elem.SetValue(C_SMAGORINSKY, C)

#
#


def CreateSolver(model_part, config):
    fluid_solver = IncompressibleFluidSolver(model_part, config.domain_size)

    # default settings
    fluid_solver.vel_toll = config.vel_toll
    if(hasattr(config, "vel_toll")):
        fluid_solver.vel_toll = config.vel_toll
    if(hasattr(config, "press_toll")):
        fluid_solver.press_toll = config.press_toll
    if(hasattr(config, "max_vel_its")):
        fluid_solver.max_vel_its = config.max_vel_its
    if(hasattr(config, "max_press_its")):
        fluid_solver.max_press_its = config.max_press_its
    if(hasattr(config, "time_order")):
        fluid_solver.time_order = config.time_order
    if(hasattr(config, "compute_reactions")):
        fluid_solver.compute_reactions = config.compute_reactions
    if(hasattr(config, "ReformDofAtEachIteration")):
        fluid_solver.ReformDofAtEachIteration = config.ReformDofAtEachIteration
    if(hasattr(config, "predictor_corrector")):
        fluid_solver.predictor_corrector = config.predictor_corrector
    if(hasattr(config, "echo_level")):
        fluid_solver.echo_level = config.echo_level
    if(hasattr(config, "dynamic_tau")):
        fluid_solver.dynamic_tau = config.dynamic_tau

    # linear solver settings
    import deprecated_trilinos_linear_solver_factory
    if(hasattr(config, "pressure_linear_solver_config")):
        fluid_solver.pressure_linear_solver = deprecated_trilinos_linear_solver_factory.ConstructSolver(
            config.pressure_linear_solver_config)
    if(hasattr(config, "velocity_linear_solver_config")):
        fluid_solver.velocity_linear_solver = deprecated_trilinos_linear_solver_factory.ConstructSolver(
            config.velocity_linear_solver_config)
    if(hasattr(config, "divergence_cleareance_step")):
        fluid_solver.divergence_clearance_steps = config.divergence_cleareance_step

    # RANS or DES settings
    if hasattr(config, "TurbulenceModel"):
        if config.TurbulenceModel == "Spalart-Allmaras":
            fluid_solver.use_spalart_allmaras = True
        elif config.TurbulenceModel == "Smagorinsky-Lilly":
            if hasattr(config, "SmagorinskyConstant"):
                fluid_solver.activate_smagorinsky(config.SmagorinskyConstant)
            else:
                msg = """Fluid solver error: Smagorinsky model requested, but
                         the value for the Smagorinsky constant is
                         undefined."""
                raise Exception(msg)

    return fluid_solver
