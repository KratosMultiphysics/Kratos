from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# -*- coding: utf-8 -*-
from KratosMultiphysics import *
from KratosMultiphysics.mpi import *
from KratosMultiphysics.MetisApplication import *
from KratosMultiphysics.IncompressibleFluidApplication import *
from KratosMultiphysics.TrilinosApplication import *
# Check that KratosMultiphysics was imported in the main script
CheckForPreviousImport()


def AddVariables(model_part):
    model_part.AddNodalSolutionStepVariable(VELOCITY)
    model_part.AddNodalSolutionStepVariable(FRACT_VEL)
    model_part.AddNodalSolutionStepVariable(MESH_VELOCITY)
    model_part.AddNodalSolutionStepVariable(PRESSURE)
    model_part.AddNodalSolutionStepVariable(PRESSURE_OLD_IT)
    model_part.AddNodalSolutionStepVariable(PRESS_PROJ)
    model_part.AddNodalSolutionStepVariable(CONV_PROJ)
    model_part.AddNodalSolutionStepVariable(NODAL_MASS)
    model_part.AddNodalSolutionStepVariable(BODY_FORCE)
    model_part.AddNodalSolutionStepVariable(DENSITY)
    model_part.AddNodalSolutionStepVariable(VISCOSITY)
    model_part.AddNodalSolutionStepVariable(EXTERNAL_PRESSURE)
    model_part.AddNodalSolutionStepVariable(FLAG_VARIABLE)

    model_part.AddNodalSolutionStepVariable(DISPLACEMENT)
    model_part.AddNodalSolutionStepVariable(IS_STRUCTURE)
    model_part.AddNodalSolutionStepVariable(IS_INTERFACE)
    model_part.AddNodalSolutionStepVariable(ARRHENIUS)

    model_part.AddNodalSolutionStepVariable(PARTITION_INDEX)

    mpi.world.barrier()
    print("variables for the trilinos incompressible fluid solver added correctly")


def AddDofs(model_part):

    for node in model_part.Nodes:
        # adding dofs
        node.AddDof(PRESSURE)
        node.AddDof(FRACT_VEL_X)
        node.AddDof(FRACT_VEL_Y)
        node.AddDof(FRACT_VEL_Z)
        node.AddDof(VELOCITY_X)
        node.AddDof(VELOCITY_Y)
        node.AddDof(VELOCITY_Z)
    mpi.world.barrier()

    print("dofs for the trilinos incompressible fluid solver added correctly")


def ReadRestartFile(FileName, nodes):
    NODES = nodes
    aaa = open(FileName)
    for line in aaa:
        exec(line)

# import start.pyinc

# aaa = __import__(FileName)
# aaa.Restart(nodes)


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
        self.CalculateReactions = False
        self.ReformDofAtEachIteration = False
        self.CalculateNormDxFlag = True
        self.laplacian_form = 2
        # 1 = laplacian, 2 = Discrete Laplacian
        self.predictor_corrector = False
        self.use_dt_in_stabilization = False

        self.echo_level = 0

        self.step = 1
        self.projections_are_initialized = False

        self.dynamic_tau = 0.001
        self.activate_tau2 = False

#
        # defining the linear solver
        velocity_aztec_parameters = ParameterList()
        velocity_aztec_parameters.set("AZ_solver", "AZ_gmres")
        velocity_aztec_parameters.set("AZ_kspace", 100)
        velocity_aztec_parameters.set("AZ_output", "AZ_none")
        # velocity_aztec_parameters.set("AZ_output",32);

# preconditioner_type = "Amesos"
# preconditioner_parameters = ParameterList()
# preconditioner_parameters.set("amesos: solver type", "Amesos_Klu");
# preconditioner_type = "Ifpack_PointRelaxation"
# preconditioner_parameters = ParameterList()
# preconditioner_parameters.set("relaxation: type", "Jacobi");

# preconditioner_type = "ILU"
# preconditioner_parameters = ParameterList()
# overlap_level = 1
# nit_max = 1000
# tol = 1e-6

        velocity_preconditioner_type = "ILU"
        velocity_preconditioner_parameters = ParameterList()
        velocity_overlap_level = 0
        velocity_nit_max = 1000
        velocity_linear_tol = 1e-6

        pressure_aztec_parameters = ParameterList()
        pressure_aztec_parameters.set("AZ_solver", "AZ_bicgstab")
        # pressure_preconditioner_type = "IC"
        pressure_preconditioner_type = "ILU"
        # pressure_preconditioner_type = "AZ_none"
        pressure_aztec_parameters.set("AZ_output", 32)
        # pressure_aztec_parameters.set("AZ_output","AZ_none");
        pressure_preconditioner_parameters = ParameterList()
        pressure_overlap_level = 0
        pressure_nit_max = 1000
        pressure_linear_tol = 1e-3

        self.velocity_linear_solver = AztecSolver(
            velocity_aztec_parameters,
            velocity_preconditioner_type,
            velocity_preconditioner_parameters,
            velocity_linear_tol,
            velocity_nit_max,
            velocity_overlap_level)
        # self.velocity_linear_solver.SetScalingType(AztecScalingType.NoScaling)
        # self.pressure_linear_solver.SetScalingType(AztecScalingType.NoScaling)
        self.velocity_linear_solver.SetScalingType(
            AztecScalingType.LeftScaling)
        # self.velocity_linear_solver.SetScalingType(AztecScalingType.SymmetricScaling)
        # self.pressure_linear_solver.SetScalingType(AztecScalingType.SymmetricScaling)

        self.pressure_linear_solver = AztecSolver(
            pressure_aztec_parameters,
            pressure_preconditioner_type,
            pressure_preconditioner_parameters,
            pressure_linear_tol,
            pressure_nit_max,
            pressure_overlap_level)
        self.pressure_linear_solver.SetScalingType(
            AztecScalingType.LeftScaling)

        # handling slip condition
        self.slip_conditions_initialized = False
        self.create_slip_conditions = GenerateSlipConditionProcess(
            self.model_part, domain_size)

#

# vel_solver_parameters = ParameterList()
# self.velocity_linear_solver =  AmesosSolver("Superludist",vel_solver_parameters);
#
# press_solver_parameters = ParameterList()
# self.pressure_linear_solver =  AmesosSolver("Superludist",press_solver_parameters);

    def Initialize(self):
        (self.neighbour_search).Execute()

        self.model_part.ProcessInfo.SetValue(DYNAMIC_TAU, self.dynamic_tau)
        self.model_part.ProcessInfo.SetValue(
            ACTIVATE_TAU2, self.activate_tau2)

#        self.solver = ResidualBasedFluidStrategyCoupled(self.model_part,self.velocity_linear_solver,self.pressure_linear_solver,self.CalculateReactions,self.ReformDofAtEachIteration,self.CalculateNormDxFlag,self.vel_toll,self.press_toll,self.max_vel_its,self.max_press_its, self.time_order,self.domain_size, self.laplacian_form, self.predictor_corrector)
#        print "in python: okkio using Coupled Strategy"
#        self.solver = ResidualBasedFluidStrategy(self.model_part,self.velocity_linear_solver,self.pressure_linear_solver,self.CalculateReactions,self.ReformDofAtEachIteration,self.CalculateNormDxFlag,self.vel_toll,self.press_toll,self.max_vel_its,self.max_press_its, self.time_order,self.domain_size, self.laplacian_form, self.predictor_corrector)

        self.Comm = CreateCommunicator()
        print(self.Comm)

        solver_configuration = TrilinosFractionalStepConfiguration(
            self.Comm,
            self.model_part,
            self.velocity_linear_solver,
            self.pressure_linear_solver,
            self.domain_size,
            self.laplacian_form,
            self.use_dt_in_stabilization)

        print("solver configuration created correctly")
        self.solver = TrilinosFractionalStepStrategy(
            self.model_part,
            solver_configuration,
            self.ReformDofAtEachIteration,
            self.vel_toll,
            self.press_toll,
            self.max_vel_its,
            self.max_press_its,
            self.time_order,
            self.domain_size,
            self.predictor_corrector)

        # generating the slip conditions
        self.create_slip_conditions.Execute()
        #(self.solver).SetSlipProcess(self.create_slip_conditions);
        self.slip_conditions_initialized = True

        (self.solver).SetEchoLevel(self.echo_level)

# (self.solver).SetEchoLevel(self.echo_level)
        print("finished initialization of the fluid strategy")

    def Solve(self):
        if(self.ReformDofAtEachIteration):
            self.slip_conditions_initialized = False

        if(self.slip_conditions_initialized == False):
            self.create_slip_conditions.Execute()
            #(self.solver).SetSlipProcess(self.create_slip_conditions);
            self.slip_conditions_initialized = True

        (self.solver).Solve()

        if(self.CalculateReactions):
            self.solver.ComputeReactions(REACTION)
        #(self.solver).ApplyFractionalVelocityFixity()
        #(self.solver).InitializeFractionalStep(self.step, self.time_order);
        #(self.solver).InitializeProjections(self.step,self.projections_are_initialized);
        # self.projections_are_initialized = True;

        #(self.solver).AssignInitialStepValues();

        #(self.solver).SolveStep1(self.vel_toll,self.max_vel_its)
        #(self.solver).SolveStep2();
        #(self.solver).ActOnLonelyNodes();
        #(self.solver).SolveStep3();
        #(self.solver).SolveStep4();

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

    def ActivateSpalartAllmaras(self, wall_nodes, DES, CDES=1.0):
        for node in wall_nodes:
            node.SetValue(IS_VISITED, 1.0)

        distance_calculator = BodyDistanceCalculationUtils()
        distance_calculator.CalculateDistances2D(
            self.model_part.Elements, DISTANCE, 100.0)

        non_linear_tol = 0.001
        max_it = 10
        reform_dofset = self.ReformDofAtEachIteration
        time_order = self.time_order

        turb_aztec_parameters = ParameterList()
        turb_aztec_parameters.set("AZ_solver", "AZ_gmres")
        turb_aztec_parameters.set("AZ_kspace", 100)
        turb_aztec_parameters.set("AZ_output", "AZ_none")

        turb_preconditioner_type = "ILU"
        turb_preconditioner_parameters = ParameterList()
        turb_overlap_level = 0
        turb_nit_max = 1000
        turb_linear_tol = 1e-9

        turb_linear_solver = AztecSolver(
            turb_aztec_parameters,
            turb_preconditioner_type,
            turb_preconditioner_parameters,
            turb_linear_tol,
            turb_nit_max,
            turb_overlap_level)
        turb_linear_solver.SetScalingType(AztecScalingType.LeftScaling)

        turbulence_model = TrilinosSpalartAllmarasTurbulenceModel(
            self.Comm,
            self.model_part,
            turb_linear_solver,
            self.domain_size,
            non_linear_tol,
            max_it,
            reform_dofset,
            time_order)
        turbulence_model.AdaptForFractionalStep()
        if(DES):
            turbulence_model.ActivateDES(CDES)

        self.solver.AddInitializeIterationProcess(turbulence_model)
