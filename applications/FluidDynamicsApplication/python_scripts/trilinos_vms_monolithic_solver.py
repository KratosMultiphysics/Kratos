from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.FluidDynamicsApplication import *
from KratosMultiphysics.mpi import *
from KratosMultiphysics.MetisApplication import *
from KratosMultiphysics.TrilinosApplication import *
# Check that KratosMultiphysics was imported in the main script
CheckForPreviousImport()


def AddVariables(model_part, config=None):
    model_part.AddNodalSolutionStepVariable(VELOCITY)
    model_part.AddNodalSolutionStepVariable(ACCELERATION)
    model_part.AddNodalSolutionStepVariable(MESH_VELOCITY)
    model_part.AddNodalSolutionStepVariable(PRESSURE)
    model_part.AddNodalSolutionStepVariable(IS_STRUCTURE)
    model_part.AddNodalSolutionStepVariable(DISPLACEMENT)
    model_part.AddNodalSolutionStepVariable(VISCOSITY)
    model_part.AddNodalSolutionStepVariable(DENSITY)
    model_part.AddNodalSolutionStepVariable(BODY_FORCE)
    model_part.AddNodalSolutionStepVariable(NODAL_AREA)
    model_part.AddNodalSolutionStepVariable(NODAL_H)
    model_part.AddNodalSolutionStepVariable(ADVPROJ)
    model_part.AddNodalSolutionStepVariable(DIVPROJ)
    model_part.AddNodalSolutionStepVariable(REACTION)
    model_part.AddNodalSolutionStepVariable(REACTION_WATER_PRESSURE)
    model_part.AddNodalSolutionStepVariable(EXTERNAL_PRESSURE)
    model_part.AddNodalSolutionStepVariable(FLAG_VARIABLE)
    model_part.AddNodalSolutionStepVariable(NORMAL)
    model_part.AddNodalSolutionStepVariable(Y_WALL)
    model_part.AddNodalSolutionStepVariable(PARTITION_INDEX)
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
        print("variables for the trilinos monolithic fluid solver added correctly")


def AddDofs(model_part, config=None):
    for node in model_part.Nodes:
        # adding dofs
        node.AddDof(VELOCITY_X, REACTION_X)
        node.AddDof(VELOCITY_Y, REACTION_Y)
        node.AddDof(VELOCITY_Z, REACTION_Z)
        node.AddDof(PRESSURE, REACTION_WATER_PRESSURE)

    if config is not None:
        if hasattr(config, "TurbulenceModel"):
            if config.TurbulenceModel == "Spalart-Allmaras":
                for node in model_part.Nodes:
                    node.AddDof(TURBULENT_VISCOSITY)

    mpi.world.barrier()
    if mpi.rank == 0:
        print("dofs for the trilinos monolithic fluid solver added correctly")


class MonolithicSolver:
    #

    def __init__(self, model_part, domain_size):

        self.model_part = model_part
        self.domain_size = domain_size

        self.alpha = -0.3
        self.move_mesh_strategy = 0

        self.Comm = CreateCommunicator()

        # for Spalart-Allmaras
        self.use_spalart_allmaras = False
        self.wall_nodes = None

        # definition of the convergence criteria
        self.vel_criteria = 1e-3
        self.press_criteria = 1e-3
        self.vel_abs_criteria = 1e-9
        self.press_abs_criteria = 1e-9

        self.max_iter = 20

        # default settings
        self.echo_level = 1
        self.CalculateReactionFlag = False
        self.ReformDofSetAtEachStep = False
        self.CalculateNormDxFlag = True
        self.MoveMeshFlag = False

        self.linear_solver_tol = 1e-6
        self.linear_solver_max_it = 100

        import MonolithicMultiLevelSolver
        self.linear_solver = MonolithicMultiLevelSolver.LinearSolver(
            self.linear_solver_tol, self.linear_solver_max_it)

        if(domain_size == 2):
            estimate_neighbours = 10
            self.guess_row_size = estimate_neighbours * (self.domain_size + 1)
        else:
            estimate_neighbours = 25
            self.guess_row_size = estimate_neighbours * (self.domain_size + 1)
        self.use_spalart_allmaras = False
        self.use_des = False
        self.Cdes = 1.0
        self.wall_nodes = list()
        self.spalart_allmaras_linear_solver = None
        self.turbulence_model = None

        self.divergence_clearance_steps = 0

    #
    def Initialize(self):

        # CHAPUZA to set the non historical value of IS_STRUCTURE correctly...
        # to be improved
        for condition in self.model_part.Conditions:
            if condition.GetValue(IS_STRUCTURE) == 1.0:
                for node in condition.GetNodes():
                    node.SetSolutionStepValue(IS_STRUCTURE, 0, 1.0)
        self.model_part.GetCommunicator().AssembleCurrentData(IS_STRUCTURE)
        mpi.world.barrier()
        for node in self.model_part.Nodes:
            if node.GetSolutionStepValue(IS_STRUCTURE, 0) != 0.0:
                node.SetValue(IS_STRUCTURE, 1.0)
                node.SetSolutionStepValue(IS_STRUCTURE, 0, 1.0)

        # compute normals "correctly"
        self.normal_calculator = NormalCalculationUtils()
        self.normal_calculator.CalculateOnSimplex(
            self.model_part,
            self.domain_size,
            IS_STRUCTURE)

        # If Spalart-Allmaras: Initialize Spalart-Allmaras solver
        if self.use_spalart_allmaras:
            self.activate_spalart_allmaras()

        if self.turbulence_model is not None:
            self.time_scheme = TrilinosPredictorCorrectorVelocityBossakSchemeTurbulent(
                self.alpha, self.move_mesh_strategy, self.domain_size, self.turbulence_model)
        else:  # No turbulence model
            self.time_scheme = TrilinosPredictorCorrectorVelocityBossakSchemeTurbulent(
                self.alpha, self.move_mesh_strategy, self.domain_size, PATCH_INDEX)

        self.time_scheme.Check(self.model_part)

        self.conv_criteria = TrilinosUPCriteria(
            self.vel_criteria,
            self.vel_abs_criteria,
            self.press_criteria,
            self.press_abs_criteria)

        # creating the solution strategy
        #import trilinos_strategy_python_periodic
        #self.solver = trilinos_strategy_python_periodic.SolvingStrategyPeriodic(
        #    self.domain_size,
        #    self.model_part,
        #    self.time_scheme,
        #    self.linear_solver,
        #    self.conv_criteria,
        #    self.CalculateReactionFlag,
        #    self.ReformDofSetAtEachStep,
        #    self.MoveMeshFlag,
        #    self.Comm,
        #    self.guess_row_size,
        #    PATCH_INDEX)
        import trilinos_strategy_python
        self.solver = trilinos_strategy_python.SolvingStrategyPython(
            "standard",
            self.model_part,
            self.time_scheme,
            self.linear_solver,
            self.conv_criteria,
            self.CalculateReactionFlag,
            self.ReformDofSetAtEachStep,
            self.MoveMeshFlag,
            self.Comm,
            self.guess_row_size)
        self.solver.max_iter = self.max_iter

        (self.solver).SetEchoLevel(self.echo_level)

    #
    def Solve(self):
        if self.divergence_clearance_steps > 0:
            if mpi.rank == 0:
                print("Calculating divergence-free initial condition")
            # initialize with a Stokes solution step
            stokes_aztec_parameters = ParameterList()
            stokes_aztec_parameters.set("AZ_solver", "AZ_gmres")
            stokes_aztec_parameters.set("AZ_kspace", 100)
            stokes_aztec_parameters.set("AZ_output", "AZ_none")

            stokes_preconditioner_type = "ILU"
            stokes_preconditioner_parameters = ParameterList()
            stokes_overlap_level = 0
            stokes_nit_max = 1000
            stokes_linear_tol = 1e-9

            stokes_linear_solver = AztecSolver(stokes_aztec_parameters,
                                               stokes_preconditioner_type,
                                               stokes_preconditioner_parameters,
                                               stokes_linear_tol,
                                               stokes_nit_max,
                                               stokes_overlap_level)
            stokes_linear_solver.SetScalingType(AztecScalingType.LeftScaling)
            stokes_process = TrilinosStokesInitializationProcess(
                self.Comm,
                self.model_part,
                stokes_linear_solver,
                self.domain_size,
                PATCH_INDEX)
            # copy periodic conditions to Stokes problem
            stokes_process.SetConditions(self.model_part.Conditions)
            # execute Stokes process
            stokes_process.Execute()
            stokes_process = None

            for node in self.model_part.Nodes:
                node.SetSolutionStepValue(PRESSURE, 0, 0.0)
                node.SetSolutionStepValue(ACCELERATION_X, 0, 0.0)
                node.SetSolutionStepValue(ACCELERATION_Y, 0, 0.0)
                node.SetSolutionStepValue(ACCELERATION_Z, 0, 0.0)
#                vel = node.GetSolutionStepValue(VELOCITY)
#                for i in range(0,2):
#                    node.SetSolutionStepValue(VELOCITY,i,vel)

            self.divergence_clearance_steps = 0
            if mpi.rank == 0:
                print("Finished divergence clearance")

        (self.solver).Solve()

    #
    def SetEchoLevel(self, level):
        (self.solver).SetEchoLevel(level)

    #
    def Clear(self):
        (self.solver).Clear()

    #
    def activate_smagorinsky(self, C):
        for elem in self.model_part.Elements:
            elem.SetValue(C_SMAGORINSKY, C)

    def activate_spalart_allmaras(self):

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
        reform_dofset = self.ReformDofSetAtEachStep
        time_order = 2

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

        self.turbulence_model = TrilinosSpalartAllmarasTurbulenceModel(
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


#
#
def CreateSolver(model_part, config):
    fluid_solver = MonolithicSolver(model_part, config.domain_size)

    if(hasattr(config, "alpha")):
        fluid_solver.alpha = config.alpha

    # definition of the convergence criteria
    if(hasattr(config, "velocity_relative_tolerance")):
        fluid_solver.rel_vel_tol = config.velocity_relative_tolerance
    if(hasattr(config, "velocity_absolute_tolerance")):
        fluid_solver.abs_vel_tol = config.velocity_absolute_tolerance
    if(hasattr(config, "pressure_relative_tolerance")):
        fluid_solver.rel_pres_tol = config.pressure_relative_tolerance
    if(hasattr(config, "pressure_absolute_tolerance")):
        fluid_solver.abs_pres_tol = config.pressure_absolute_tolerance
    if(hasattr(config, "dynamic_tau")):
        fluid_solver.dynamic_tau = config.dynamic_tau
    if(hasattr(config, "oss_switch")):
        fluid_solver.oss_switch = config.oss_switch
    if(hasattr(config, "max_iteration")):
        fluid_solver.max_iter = config.max_iteration
    if(hasattr(config, "echo_level")):
        fluid_solver.echo_level = config.echo_level
    if(hasattr(config, "compute_reactions")):
        fluid_solver.compute_reactions = config.compute_reactions
    if(hasattr(config, "ReformDofSetAtEachStep")):
        fluid_solver.ReformDofSetAtEachStep = config.ReformDofSetAtEachStep
    if(hasattr(config, "divergence_cleareance_step")):
        fluid_solver.divergence_clearance_steps = config.divergence_cleareance_step

    import deprecated_trilinos_linear_solver_factory
    if(hasattr(config, "linear_solver_config")):
        fluid_solver.linear_solver = deprecated_trilinos_linear_solver_factory.ConstructSolver(
            config.linear_solver_config)

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
