from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.IncompressibleFluidApplication import *
from KratosMultiphysics.ThermoMechanicalApplication import *
from KratosMultiphysics.MeshingApplication import *
from KratosMultiphysics.FluidDynamicsApplication import *
from KratosMultiphysics.ExternalSolversApplication import *
from KratosMultiphysics.mpi import *
from KratosMultiphysics.TrilinosApplication import *
from KratosMultiphysics.MetisApplication import *

import trilinos_just_levelset_solver

# settings for the convection solver
distance_settings = ConvectionDiffusionSettings()
distance_settings.SetUnknownVariable(DISTANCE)
distance_settings.SetConvectionVariable(VELOCITY)
distance_settings.SetMeshVelocityVariable(MESH_VELOCITY)
# distance_settings.SetVolumeSourceVariable(HEAT_FLUX)
# distance_settings.SetDiffusionVariable(ARRHENIUSAUX)
# distance_settings.SetDensityVariable(AUX_INDEX) # For level set solver
# rho and C are assigned to 1


def AddVariables(model_part):
    model_part.AddNodalSolutionStepVariable(VELOCITY)
    model_part.AddNodalSolutionStepVariable(ACCELERATION)
    model_part.AddNodalSolutionStepVariable(MESH_VELOCITY)
    model_part.AddNodalSolutionStepVariable(PRESSURE)
    model_part.AddNodalSolutionStepVariable(IS_FLUID)
    model_part.AddNodalSolutionStepVariable(IS_STRUCTURE)
    model_part.AddNodalSolutionStepVariable(IS_FREE_SURFACE)
    model_part.AddNodalSolutionStepVariable(IS_INTERFACE)
    model_part.AddNodalSolutionStepVariable(IS_BOUNDARY)
    model_part.AddNodalSolutionStepVariable(DISPLACEMENT)
    model_part.AddNodalSolutionStepVariable(VISCOSITY)
    model_part.AddNodalSolutionStepVariable(DENSITY)
    model_part.AddNodalSolutionStepVariable(BODY_FORCE)
    model_part.AddNodalSolutionStepVariable(NODAL_AREA)
    model_part.AddNodalSolutionStepVariable(NODAL_H)
    model_part.AddNodalSolutionStepVariable(THAWONE)
    model_part.AddNodalSolutionStepVariable(THAWTWO)
    model_part.AddNodalSolutionStepVariable(FLAG_VARIABLE)
    model_part.AddNodalSolutionStepVariable(NORMAL)
    model_part.AddNodalSolutionStepVariable(DISTANCE)
    model_part.AddNodalSolutionStepVariable(IS_SLIP)
    model_part.AddNodalSolutionStepVariable(PRESSURES)
    model_part.AddNodalSolutionStepVariable(MATERIAL)
    model_part.AddNodalSolutionStepVariable(LAST_AIR)
    model_part.AddNodalSolutionStepVariable(NODAL_MASS)
    model_part.AddNodalSolutionStepVariable(NODAL_PAUX)
    # model_part.AddNodalSolutionStepVariable(WET_VOLUME)
    # variables needed for the distance solver
    trilinos_just_levelset_solver.AddVariables(model_part, distance_settings)

    print("variables for the MONOLITHIC_SOLVER_EULERIAN added correctly")


def AddDofs(model_part):
    for node in model_part.Nodes:
        # adding dofs
        node.AddDof(VELOCITY_X)
        node.AddDof(VELOCITY_Y)
        node.AddDof(VELOCITY_Z)
        node.AddDof(PRESSURE)

    trilinos_just_levelset_solver.AddDofs(model_part, distance_settings)
    print("dofs for the monolithic solver added correctly")


class MonolithicSolver:
    #

    def __init__(self, model_part, domain_size):

        self.model_part = model_part
        self.domain_size = domain_size

        self.Comm = CreateCommunicator()
        self.linear_solver = TrilinosLinearSolver()

        self.alpha = -0.3
        self.move_mesh_strategy = 0
        # self.time_scheme = ResidualBasedPredictorCorrectorVelocityBossakSchemeDPGEnriched( self.alpha,self.move_mesh_strategy,self.domain_size )
        self.time_scheme = TrilinosResidualBasedPredictorCorrectorVelocityBossakSchemeDPGEnriched(
            self.alpha, self.move_mesh_strategy, self.domain_size)
        self.time_scheme.Check(self.model_part)

        if(domain_size == 2):
            estimate_neighbours = 10
            self.guess_row_size = estimate_neighbours * (self.domain_size + 1)
            # self.buildertype="ML2Dpress"
        else:
            estimate_neighbours = 25
            self.guess_row_size = estimate_neighbours * (self.domain_size + 1)
            # self.buildertype="ML3Dpress"
        self.buildertype = "standard"
        # aztec_parameters = ParameterList()
        # aztec_parameters.set("AZ_solver","AZ_gmres");
        # aztec_parameters.set("AZ_kspace",200);
        # aztec_parameters.set("AZ_output","AZ_none");
        # aztec_parameters.set("AZ_output",10);
        # preconditioner_type = "ILU"
        # preconditioner_parameters = ParameterList()
        # preconditioner_parameters.set ("fact: drop tolerance", 1e-9);
        # preconditioner_parameters.set ("fact: level-of-fill", 1);
        # overlap_level = 0
        # nit_max = 1000
        # linear_tol = 1e-5
        # self.linear_solver =  AztecSolver(aztec_parameters,preconditioner_type,preconditioner_parameters,linear_tol,nit_max,overlap_level);

        aztec_parameters = ParameterList()
        aztec_parameters.set("AZ_solver", "AZ_gmres")
        aztec_parameters.set("AZ_output", 1)
        aztec_parameters.set("AZ_kspace", 50)

        # settings of the ML solver
        MLList = ParameterList()
        default_settings = EpetraDefaultSetter()
        default_settings.SetDefaults(MLList, "NSSA")
        MLList.set("ML output", 1)
        MLList.set("coarse: max size", 10000)
        MLList.set("max levels", 3)
        MLList.set("aggregation: type", "Uncoupled")
        MLList.set("coarse: type", "Amesos-Superludist")
        tolerance = 1e-4
        max_iterations = 50
        self.linear_solver = MultiLevelSolver(
            aztec_parameters,
            MLList,
            tolerance,
            max_iterations)

        # definition of the convergence criteria
        self.rel_vel_tol = 1e-5
        self.abs_vel_tol = 1e-7
        self.rel_pres_tol = 1e-5
        self.abs_pres_tol = 1e-7

        self.dynamic_tau_levelset = 0.01
        self.dynamic_tau_fluid = 1.0
        self.oss_switch = 0

        # non newtonian setting
        self.regularization_coef = 1000

        self.max_iter = 30

        # default settings
        self.echo_level = 0
        self.CalculateReactionFlag = True
        self.ReformDofSetAtEachStep = True
        self.CalculateNormDxFlag = True
        self.MoveMeshFlag = False
        self.volume_correction = True

# print "Construction monolithic solver finished"
        mpi.world.barrier()
        # Creat Lavel_set solver
        # construct the model part
        if(domain_size == 2):
            raise "error, still not implemented in 2D"
            conv_elem = "SUPGConv2D"
            conv_cond = "Condition2D"
        else:
            conv_elem = "SUPGConv3D"
            conv_cond = "Condition3D"
        self.level_set_model_part = ModelPart("level_set_model_part")
        self.conv_generator = ConnectivityPreserveModeler()
        (self.conv_generator).GenerateModelPart(self.model_part,
                                                self.level_set_model_part, conv_elem, conv_cond)
        (ParallelFillCommunicator(self.level_set_model_part)).Execute()
        if(mpi.rank == 0):
            print("finished initializing the parallel fill communicator for convection_model_part")
        # constructing the convection solver for the distance
        self.trilinos_just_levelset_solver = trilinos_just_levelset_solver.Solver(
            self.level_set_model_part, domain_size, distance_settings)
        self.trilinos_just_levelset_solver.max_iterations = 8

        #
        # properties of the two fluids
        self.rho1 = 2400.0  # applied on the negative part of the domain 1000.0
        self.conductivity1 = 1.0

        self.rho2 = 1.0  # applied to the positive part of the domain#1.0
        self.conductivity2 = 1.0

        self.mu = 3.0e-3
        self.divergence_clearance_performed = False
        #

        # Distance utilities

         #
        if(self.domain_size == 2):
            self.redistance_utils = ParallelDistanceCalculator2D()
        else:
            self.redistance_utils = ParallelDistanceCalculator3D()

        self.redistance_frequency = 1
        self.max_edge_size = self.redistance_utils.FindMaximumEdgeSize(
            self.level_set_model_part)
        self.max_distance = self.max_edge_size * 5.0
        self.max_levels = 25  # self.max_distance/self.min_edge_size

        self.max_ns_iterations = 8
        self.internal_step_counter = 1

        # Slip condition
        self.use_slip_conditions = False

        # volume correction
        self.volume_correction_switch = True
    #

    def ApplyFluidProperties(self):
        # apply density
        mu1 = 1.0 * self.mu / self.rho1
        # mu1 = self.mu
        # mu2 = 0.01*self.mu/self.rho2
        mu2 = mu1

        BiphasicFillingUtilities().ApplyFluidProperties(
            self.model_part, mu1, self.rho1, mu2, self.rho2)
# for node in self.model_part.Nodes:
# dist = node.GetSolutionStepValue(DISTANCE)
# if(dist < 0):
# node.SetSolutionStepValue(DENSITY,0,self.rho1)
# node.SetSolutionStepValue(VISCOSITY,0,mu1)
# else:
# node.SetSolutionStepValue(DENSITY,0,self.rho2)
# node.SetSolutionStepValue(VISCOSITY,0,mu2)
    #

    def Initialize(self):
        # creating the solution strategy
        # self.conv_criteria = VelPrCriteria(self.rel_vel_tol,self.abs_vel_tol,\
                                           # self.rel_pres_tol,self.abs_pres_tol)
        self.conv_criteria = TrilinosUPCriteria(
            self.rel_vel_tol,
            self.abs_vel_tol,
            self.rel_pres_tol,
            self.abs_pres_tol,
            self.Comm)

        # builder_and_solver = ResidualBasedBlockBuilderAndSolver(self.linear_solver)
        # self.solver = ResidualBasedNewtonRaphsonStrategy(self.model_part,self.time_scheme,self.linear_solver,self.conv_criteria,builder_and_solver,self.max_iter,self.CalculateReactionFlag, self.ReformDofSetAtEachStep,self.MoveMeshFlag)
        import trilinos_strategy_python
        self.solver = trilinos_strategy_python.SolvingStrategyPython(
            self.buildertype,
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
        print(">>>>>>>>>>>>>>>", self.oss_switch)
        self.model_part.ProcessInfo.SetValue(
            DYNAMIC_TAU, self.dynamic_tau_fluid)
        self.model_part.ProcessInfo.SetValue(OSS_SWITCH, self.oss_switch)

        # LEvel_set solver initialization
        self.trilinos_just_levelset_solver.dynamic_tau = self.dynamic_tau_levelset
# self.redistance_utils.CalculateDistances(self.model_part,DISTANCE,NODAL_AREA,self.max_levels,self.max_distance)
        self.redistance_utils.CalculateInterfacePreservingDistances(
            self.model_part,
            DISTANCE,
            NODAL_AREA,
            self.max_levels,
            self.max_distance)
        self.trilinos_just_levelset_solver.Initialize()

        self.ApplyFluidProperties()

        self.next_redistance = self.redistance_frequency
        #
        # FOR SLIP
        #
        # Manullay assign!
        for cond in self.model_part.Conditions:
            cond.SetValue(IS_STRUCTURE, 1.0)
        # if we use slip conditions, calculate normals on the boundary
        if (self.use_slip_conditions):
            (FindConditionsNeighboursProcess(
                self.model_part, 3, 20)).ClearNeighbours()
            (FindConditionsNeighboursProcess(self.model_part, 3, 20)).Execute()
            self.normal_util = NormalCalculationUtils()
            self.normal_util.CalculateOnSimplex(
                self.model_part,
                self.domain_size,
                IS_STRUCTURE,
                0.0,
                35.0)  # ,0.0,35.0

        # saving inlet nodes
        self.inlet_nodes = []
        for cond in self.model_part.Conditions:
            if(cond.GetValue(IS_INLET) > 0):
                for node in cond.GetNodes():
                    self.inlet_nodes.append(node)

    #
    def DoRedistance(self):
        mpi.world.barrier()
        # self.redistance_utils.CalculateDistances(self.model_part,DISTANCE,NODAL_AREA,self.max_levels,self.max_distance)
        self.redistance_utils.CalculateInterfacePreservingDistances(
            self.model_part,
            DISTANCE,
            NODAL_AREA,
            self.max_levels,
            self.max_distance)
        mpi.world.barrier()

     #
    def ConvectDistance(self):
        mpi.world.barrier()
        self.level_set_model_part.ProcessInfo = self.model_part.ProcessInfo
        (self.level_set_model_part.ProcessInfo).SetValue(
            CONVECTION_DIFFUSION_SETTINGS, distance_settings)
        (self.level_set_model_part.ProcessInfo).SetValue(
            DYNAMIC_TAU, self.dynamic_tau_levelset)  # self.dynamic_tau
        (self.trilinos_just_levelset_solver).Solve()
        BiphasicFillingUtilities(
        ).DistanceFarRegionCorrection(
            self.model_part,
            self.max_distance)
        mpi.world.barrier()

     #
      #
    def Solve(self, step):
        # at the beginning of the calculations do a div clearance step
        if(self.divergence_clearance_performed == False):
            for node in self.model_part.Nodes:
                node.SetSolutionStepValue(
                    DISTANCE,
                    1,
                    node.GetSolutionStepValue(DISTANCE))
            self.divergence_clearance_performed = True

        # convect distance function
        self.ConvectDistance()
        # recompute distance function as needed
        if(self.internal_step_counter >= self.next_redistance):
            self.DoRedistance()
            # BiphasicFillingUtilities().DistanceFarRegionCorrection(self.model_part,  self.max_distance)
            self.next_redistance = self.internal_step_counter + \
                self.redistance_frequency

        if(self.volume_correction_switch and step > 10):
            net_volume = self.model_part.ProcessInfo[NET_INPUT_MATERIAL]
            BiphasicFillingUtilities().VolumeCorrection(
                self.model_part, net_volume)
        self.ApplyFluidProperties()
        # Recompute normals if necessary
# if(self.ReformDofSetAtEachStep == True):
# if self.use_slip_conditions == True:
# self.normal_util.CalculateOnSimplex(self.model_part,self.domain_size,IS_STRUCTURE,0.0,35.0)#,0.0,35.0
        mpi.world.barrier()
        (self.solver).Solve()
        self.internal_step_counter += 1
        mpi.world.barrier()
    #

    def SetEchoLevel(self, level):
        (self.solver).SetEchoLevel(level)

    #
