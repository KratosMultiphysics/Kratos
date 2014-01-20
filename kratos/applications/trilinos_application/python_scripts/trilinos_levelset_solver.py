from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# -*- coding: utf-8 -*-
# importing the Kratos Library
from KratosMultiphysics import *
# from KratosMultiphysics.ConvectionDiffusionApplication import *
from KratosMultiphysics.ThermoMechanicalApplication import *
from KratosMultiphysics.TrilinosApplication import *
from KratosMultiphysics.MeshingApplication import *
from KratosMultiphysics.IncompressibleFluidApplication import EstimateDt3D
from KratosMultiphysics.mpi import *
# Check that KratosMultiphysics was imported in the main script
CheckForPreviousImport()

import trilinos_pureconvection_solver
import trilinos_thermal_solver
import trilinos_monolithic_solver_eulerian

# settings for the convection solver
distance_settings = ConvectionDiffusionSettings()
distance_settings.SetUnknownVariable(DISTANCE)
distance_settings.SetConvectionVariable(VELOCITY)
distance_settings.SetMeshVelocityVariable(MESH_VELOCITY)

temperature_settings = ConvectionDiffusionSettings()
temperature_settings.SetDensityVariable(DENSITY)
temperature_settings.SetDiffusionVariable(CONDUCTIVITY)
temperature_settings.SetUnknownVariable(TEMPERATURE)
temperature_settings.SetVolumeSourceVariable(HEAT_FLUX)
temperature_settings.SetSurfaceSourceVariable(FACE_HEAT_FLUX)
temperature_settings.SetMeshVelocityVariable(MESH_VELOCITY)
temperature_settings.SetConvectionVariable(VELOCITY)


def AddVariables(model_part):
    print(distance_settings)
    # variables needed for the fluid solver
    trilinos_monolithic_solver_eulerian.AddVariables(model_part)
    model_part.AddNodalSolutionStepVariable(DISTANCE)
    model_part.AddNodalSolutionStepVariable(NODAL_AREA)
    model_part.AddNodalSolutionStepVariable(IS_STRUCTURE)

    # variables needed for the distance solver
    trilinos_pureconvection_solver.AddVariables(model_part, distance_settings)
    trilinos_thermal_solver.AddVariables(model_part, temperature_settings)

    model_part.AddNodalSolutionStepVariable(PARTITION_INDEX)


def AddDofs(model_part):
    trilinos_monolithic_solver_eulerian.AddDofs(model_part)
    trilinos_pureconvection_solver.AddDofs(model_part, distance_settings)
    trilinos_thermal_solver.AddDofs(model_part, temperature_settings)

    print("variables for the convection diffusion solver added correctly")


class TrilinosLevelSetSolver:

    def __init__(self, model_part, domain_size):
        self.model_part = model_part
        self.domain_size = domain_size

        mpi.world.barrier()

        # construct the model part for the convection solver
        if(self.domain_size == 2):
            raise "error, still not implemented in 2D"
            conv_elem = "SUPGConv2D"
            conv_cond = "Condition2D"
        else:
            conv_elem = "SUPGConv3D"
            conv_cond = "Condition3D"
        self.convection_model_part = ModelPart("convection_model_part")
        self.conv_generator = ConnectivityPreserveModeler()
        (self.conv_generator).GenerateModelPart(self.model_part, self.convection_model_part, conv_elem, conv_cond)
        (ParallelFillCommunicator(self.convection_model_part)).Execute()
        if(mpi.rank == 0):
            print("finished initializing the parallel fill communicator for convection_model_part")

        # construct the model part for the t solver
        if(self.domain_size == 2):
            conv_elem = "SUPGConvDiff2D"
            conv_cond = "ThermalFace2D"
        else:
            conv_elem = "SUPGConvDiff3D"
            conv_cond = "ThermalFace3D"
        self.thermal_model_part = ModelPart("thermal_model_part")
        self.conv_generator = ConnectivityPreserveModeler()
        (self.conv_generator).GenerateModelPart(self.model_part, self.thermal_model_part, conv_elem, conv_cond)
        (ParallelFillCommunicator(self.thermal_model_part)).Execute()
        if(mpi.rank == 0):
            print("finished initializing the parallel fill communicator for thermal_model_part")

        # constructing the convection solver for the distance
        self.convection_solver = trilinos_pureconvection_solver.Solver(self.convection_model_part, self.domain_size, distance_settings)
        self.convection_solver.max_iterations = 8

        # constructing the convection solver for the distance
        self.thermal_solver = trilinos_thermal_solver.Solver(self.thermal_model_part, self.domain_size, temperature_settings)

        # constructing the fluid solver
        self.fluid_solver = trilinos_monolithic_solver_eulerian.MonolithicSolver(self.model_part, self.domain_size)
        self.vel_criteria = 1e-3
        self.press_criteria = 1e-5
        self.vel_abs_criteria = 1e-9
        self.press_abs_criteria = 1e-9
        self.fluid_solver.ReformDofSetAtEachStep = False
        self.fluid_solver.alpha = -0.3

        #
        # properties of the two fluids
        self.rho1 = 1000.0  # applied on the negative part of the domain
        self.conductivity1 = 1.0
        self.specific_heat1 = 1.0

        self.rho2 = 1.0  # applied to the positive part of the domain
        self.conductivity2 = 1.0
        self.specific_heat2 = 1.0

        # common properties for the two fluids
        self.mu = 1.0e-3
        self.ambient_temperature = 293.15  # note that temperatures should be given in Kelvins
        self.mould_temperature = self.ambient_temperature
        self.inlet_temperature = 800.0
        self.convection_coefficient = 0.0
        self.emissivity = 0.0
        #

        if(self.domain_size == 2):
            self.redistance_utils = ParallelDistanceCalculator2D()
        else:
            self.redistance_utils = ParallelDistanceCalculator3D()

        self.max_levels = 50
        self.redistance_frequency = 1
        self.max_edge_size = self.redistance_utils.FindMaximumEdgeSize(self.convection_model_part)
        self.max_distance = self.max_edge_size * 3.0

        # assigning the fluid properties
#        conductivity = 0.0;
#        density = 1.0;
#        specific_heat = 1.0;
#        for node in model_part.Nodes:
#            node.SetSolutionStepValue(temperature_settings.GetDiffusionVariable(),0,conductivity);
#            node.SetSolutionStepValue(temperature_settings.GetDensityVariable(),0,density);
#            node.SetSolutionStepValue(SPECIFIC_HEAT,0,specific_heat);

        self.max_ns_iterations = 8
        self.dynamic_tau = 1.00

        self.divergence_clearance_performed = True  # setting to true it will not perform it

        mpi.world.barrier()

        # create utility to estimate Dt
        self.dt_estimator = EstimateDt3D(self.model_part)

        # utility to set to zero variables
        self.variable_utils = VariableUtils()

        self.internal_step_counter = 0

    def EchoSettings(self):
        mpi.world.barrier()
        if(mpi.rank == 0):
            print(" ")
            print("*******************************************")
            print("settings currently used : ")
            print("*******************************************")
            print("rho negative domain (rho1)                        = ", self.rho1)
            print("rho positive domain (rho2)                        = ", self.rho2)
            print("mu                                                = ", self.mu)
            print(" ")
            print("FLUID SOLVER SETTINGS : ")
            print("vel_criteria                                      =", self.vel_criteria)
            print("vel_abs_criteria                                  =", self.vel_abs_criteria)
            print("press_criteria                                    =", self.press_criteria)
            print("press_abs_criteria                                =", self.press_abs_criteria)
            print("max_ns_iterations                                 =", self.max_ns_iterations)
            print("redistance_frequency                              =", self.redistance_frequency)
            print(" ")
            print("THERMAL SOLVER SETTINGS : ")
            print("conductivity negative domain (conductivity1)      = ", self.conductivity1)
            print("conductivity positive domain (conductivity2)      = ", self.conductivity2)
            print("specific_heat negative domain (specific_heat1)    = ", self.specific_heat1)
            print("specific_heat positive domain (specific_heat2)    = ", self.specific_heat2)
            print("inlet temperature                                 = ", self.inlet_temperature)
            print("ambient_temperature                               = ", self.ambient_temperature)
            print("mould_temperature                                 = ", self.mould_temperature)
            print("convection_coefficient                            = ", self.convection_coefficient)
            print("emissivity                                        = ", self.emissivity)
            print(" ")
            print("TABLES USED : ")
            print("...pooyan please fill this ...!!!!!!!!!!!!!! ")
            print(" ")
            print("*******************************************")
            print(" ")

        mpi.world.barrier()

    def ApplyFluidProperties(self):
        # apply density
        mu1 = self.mu / self.rho1
        mu2 = self.mu / self.rho2
        for node in self.model_part.Nodes:
            dist = node.GetSolutionStepValue(DISTANCE)
            if(dist < 0):
                node.SetSolutionStepValue(DENSITY, 0, self.rho1)
                node.SetSolutionStepValue(VISCOSITY, 0, mu1)
            else:
                node.SetSolutionStepValue(DENSITY, 0, self.rho2)
                node.SetSolutionStepValue(VISCOSITY, 0, mu2)

    def Initialize(self):
        self.EchoSettings()

        self.convection_solver.dynamic_tau = self.dynamic_tau
        self.thermal_solver.dynamic_tau = self.dynamic_tau
        self.fluid_solver.dynamic_tau = self.dynamic_tau

        self.fluid_solver.vel_criteria = self.vel_criteria
        self.fluid_solver.press_criteria = self.press_criteria
        self.fluid_solver.vel_abs_criteria = self.vel_abs_criteria
        self.fluid_solver.press_abs_criteria = self.press_abs_criteria
        self.fluid_solver.max_iter = self.max_ns_iterations

        self.redistance_utils.CalculateDistances(self.model_part, DISTANCE, NODAL_AREA, self.max_levels, self.max_distance)
        self.convection_solver.Initialize()
        self.thermal_solver.Initialize()
        self.fluid_solver.Initialize()

        # CHAPUZA to set the non historical value of IS_STRUCTURE correctly... to be improved
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
        self.normal_calculator.CalculateOnSimplex(self.model_part, self.domain_size, IS_STRUCTURE)

        #(BodyNormalCalculationUtils()).CalculateBodyNormals(self.model_part,self.domain_size)

        # for node in self.model_part.Nodes:
           # normal = node.GetSolutionStepValue(NORMAL,0)
           # node.SetValue(NORMAL,normal)

        self.thermal_solver.SetEchoLevel(0)
        self.convection_solver.SetEchoLevel(1)

        self.next_redistance = self.redistance_frequency

        # set the temperature to the mould temperature
        for node in self.model_part.Nodes:
            node.SetSolutionStepValue(TEMPERATURE, 0, self.mould_temperature)
            node.SetSolutionStepValue(TEMPERATURE, 1, self.mould_temperature)

        # build a list of inlet nodes
        self.inlet_nodes = []
        for node in self.model_part.Nodes:
            # set to inlet temperature all of the nodes with fixed velocity
            if(node.IsFixed(VELOCITY_X) and node.IsFixed(VELOCITY_Y) and node.IsFixed(VELOCITY_Z)):
                self.inlet_nodes.append(node)
                node.Fix(TEMPERATURE)
                node.SetSolutionStepValue(TEMPERATURE, 0, self.inlet_temperature)
                node.SetSolutionStepValue(TEMPERATURE, 1, self.inlet_temperature)

            # also set to inlet temperature all of the nodes which fall within the negative areo of the domain
            if(node.GetSolutionStepValue(DISTANCE) < 0.0):
                node.SetSolutionStepValue(TEMPERATURE, 0, self.inlet_temperature)
                node.SetSolutionStepValue(TEMPERATURE, 1, self.inlet_temperature)

        self.ApplyFluidProperties()

#        self.model_part.ProcessInfo.SetValue(DYNAMIC_TAU, self.dynamic_tau);

        # set the thermal properties to the appropriate values
        for prop in self.thermal_model_part.Properties:
            prop.SetValue(AMBIENT_TEMPERATURE, self.ambient_temperature)
            prop.SetValue(EMISSIVITY, self.emissivity)
            prop.SetValue(CONVECTION_COEFFICIENT, self.convection_coefficient)

#        self.model_part.ProcessInfo.SetValue(DYNAMIC_TAU, self.dynamic_tau);

    def EstimateDt(self, CFL, dt_max):
        return self.dt_estimator.EstimateDt(CFL, dt_max)

    def DoDivergenceClearance(self):
        # do a redistance
        self.redistance_utils.CalculateDistances(self.model_part, DISTANCE, NODAL_AREA, self.max_levels, self.max_distance)

        for node in self.model_part.Nodes:
            node.SetSolutionStepValue(TEMPERATURE, 0, self.rho1)

        (self.fluid_solver).Solve()
        mpi.world.barrier()

        zero = Vector(3)
        zero[0] = 0.0
        zero[1] = 0.0
        zero[2] = 0.0

        for node in self.model_part.Nodes:
            # old_pressure = node.GetSolutionStepValue(PRESSURE,1)
            # node.SetSolutionStepValue(PRESSURE,0,old_pressure)

            node.SetSolutionStepValue(ACCELERATION, 0, zero)

            vel = node.GetSolutionStepValue(VELOCITY)
            node.SetSolutionStepValue(VELOCITY, 1, vel)
        mpi.world.barrier()

    def DoRedistance(self):
        # redistance if required
        if(mpi.rank == 0):
            print("beginning recalculation of distances")
        self.redistance_utils.CalculateDistances(self.model_part, DISTANCE, NODAL_AREA, self.max_levels, self.max_distance)
        mpi.world.barrier()
        if(mpi.rank == 0):
            print("finished recalculation of distances")

    def ConvectDistance(self):
        mpi.world.barrier()

        if(mpi.rank == 0):
            print("beginning convection step for the distance function")

        self.convection_model_part.ProcessInfo = self.model_part.ProcessInfo
        (self.convection_model_part.ProcessInfo).SetValue(CONVECTION_DIFFUSION_SETTINGS, distance_settings)
        (self.convection_model_part.ProcessInfo).SetValue(DYNAMIC_TAU, self.dynamic_tau)

        (self.convection_solver).Solve()
        mpi.world.barrier()
        if(mpi.rank == 0):
            print("finished convection step for the distance function")

    def ComputeThermalSolution(self):
        mpi.world.barrier()

        if(mpi.rank == 0):
            print("beginning thermal solution")

        self.thermal_model_part.ProcessInfo = self.model_part.ProcessInfo
        (self.thermal_model_part.ProcessInfo).SetValue(CONVECTION_DIFFUSION_SETTINGS, temperature_settings)
        (self.thermal_model_part.ProcessInfo).SetValue(DYNAMIC_TAU, 0.0)  # self.dynamic_tau)

        for node in self.thermal_model_part.Nodes:
            dist = node.GetSolutionStepValue(DISTANCE)
            if(dist < 0):
                node.SetSolutionStepValue(CONDUCTIVITY, 0, self.conductivity1)
                node.SetSolutionStepValue(SPECIFIC_HEAT, 0, self.specific_heat1)
            else:
                node.SetSolutionStepValue(CONDUCTIVITY, 0, self.conductivity2)
                node.SetSolutionStepValue(SPECIFIC_HEAT, 0, self.specific_heat2)

        (self.thermal_solver).Solve()
        mpi.world.barrier()

        # get rid of overshoots and undershoots
        for node in self.thermal_model_part.Nodes:
            temperature = node.GetSolutionStepValue(TEMPERATURE)
            if(temperature > self.inlet_temperature):
                node.SetSolutionStepValue(TEMPERATURE, 0, self.inlet_temperature)
            elif(temperature < self.ambient_temperature):
                node.SetSolutionStepValue(TEMPERATURE, 0, self.ambient_temperature)

        if(mpi.rank == 0):
            print("finished thermal solution")

    def ComputeFluidSolution(self):
        mpi.world.barrier()
        # snap distance to grid
        # eps = 1e-3*self.max_edge_size;
        # for node in self.model_part.Nodes:
            # dist = node.GetSolutionStepValue(DISTANCE)
            # node.SetValue(DISTANCE,dist)
            # if(abs(dist) < eps):
                # node.SetSolutionStepValue(DISTANCE,0,0.0)

        # apply density
#        mu1 = self.mu/self.rho1
#        mu2 = self.mu/self.rho2
#        for node in self.model_part.Nodes:
#            dist = node.GetSolutionStepValue(DISTANCE)
#            if(dist < 0):
#                node.SetSolutionStepValue(DENSITY,0,self.rho1)
#                node.SetSolutionStepValue(VISCOSITY,0,mu1)
#            else:
#                node.SetSolutionStepValue(DENSITY,0,self.rho2)
#                node.SetSolutionStepValue(VISCOSITY,0,mu2)

        self.ApplyFluidProperties()

        # apply temperature dependent properties ...

        # solve fluid
        (self.convection_model_part.ProcessInfo).SetValue(DYNAMIC_TAU, self.dynamic_tau)
        (self.fluid_solver).Solve()
        mpi.world.barrier()
        if(mpi.rank == 0):
            print("Solution Step finished")

        # for node in self.model_part.Nodes:
            # dist = node.GetValue(DISTANCE)
            # if(abs(dist) < eps):
                # node.SetSolutionStepValue(DISTANCE,0,dist)

    def Solve(self):
        mpi.world.barrier()

        # at the beginning of the calculations do a div clearance step
        if(self.divergence_clearance_performed == False):
            self.DoDivergenceClearance()
            self.divergence_clearance_performed = True

        # convect distance function
        self.ConvectDistance()

        # recompute distance function as needed
        if(self.internal_step_counter >= self.next_redistance):
            self.DoRedistance()
            self.next_redistance = self.internal_step_counter + self.redistance_frequency

        # compute temperature distribution
        self.ComputeThermalSolution()

        # solve the two-fluid problem with the properties computed so far
        self.ComputeFluidSolution()

        self.internal_step_counter += 1

    def GetSolverIterations(self):
        return (self.fluid_solver).solver.iterations_last_solution
