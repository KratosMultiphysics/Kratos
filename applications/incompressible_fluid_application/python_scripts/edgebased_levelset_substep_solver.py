from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.IncompressibleFluidApplication import *
from KratosMultiphysics.MeshingApplication import *

def AddVariables(model_part):
    model_part.AddNodalSolutionStepVariable(VELOCITY)
    model_part.AddNodalSolutionStepVariable(PRESSURE)
    model_part.AddNodalSolutionStepVariable(NORMAL)
    model_part.AddNodalSolutionStepVariable(AUX_INDEX)
    model_part.AddNodalSolutionStepVariable(DISTANCE)
    model_part.AddNodalSolutionStepVariable(PRESS_PROJ)
    model_part.AddNodalSolutionStepVariable(POROSITY)
    model_part.AddNodalSolutionStepVariable(VISCOSITY)

    model_part.AddNodalSolutionStepVariable(DIAMETER)
    model_part.AddNodalSolutionStepVariable(LIN_DARCY_COEF)
    model_part.AddNodalSolutionStepVariable(NONLIN_DARCY_COEF)

    model_part.AddNodalSolutionStepVariable(NODAL_AREA)

    print("variables for the edgebased incompressible fluid substep solver added correctly")


def AddDofs(model_part):
    for node in model_part.Nodes:
        node.AddDof(PRESSURE)
        node.AddDof(VELOCITY_X)
        node.AddDof(VELOCITY_Y)
        node.AddDof(VELOCITY_Z)


class EdgeBasedLevelSetSolver:

    def __init__(self, model_part, domain_size,
                 body_force, viscosity, density):

        print("entered in EdgeBasedLevelSetSubstepSolver python constructor")
        # data of the problem
        self.model_part = model_part
        self.domain_size = domain_size
        self.body_force = body_force
        self.density = density
        self.viscosity = viscosity
        self.stabdt_pressure_factor = 1.0
        self.stabdt_convection_factor = 0.01
        self.use_mass_correction = True
        self.redistance_frequency = 5
        self.step = 0
        self.extrapolation_layers = 5
        self.tau2_factor = 0.0
        self.edge_detection_angle = 45.0
        self.assume_constant_pressure = True
        self.timer = Timer()

        self.use_parallel_distance_calculation = False

        # neighbour search
        number_of_avg_elems = 10
        number_of_avg_nodes = 10
        self.neighbour_search = FindNodalNeighboursProcess(
            model_part, number_of_avg_elems, number_of_avg_nodes)
        (self.neighbour_search).Execute()

        # erase isolated notes
        eliminate_isolated = EliminateIsolatedNodesProcess(model_part)
        eliminate_isolated.Execute()

        # definition of the solvers
#        pDiagPrecond = DiagonalPreconditioner()
#        self.pressure_linear_solver =  CGSolver(1e-3, 5000,pDiagPrecond)
 #       self.pressure_linear_solver =  CGSolver(1e-3, 5000)

#        pDiagPrecond = DiagonalPreconditioner()
#        self.pressure_linear_solver =  BICGSTABSolver(1e-3, 5000,pDiagPrecond)
        self.pressure_linear_solver = BICGSTABSolver(1e-3, 5000)

        # initializing the press proj to -body_force
        press_proj_init = Vector(3)
        press_proj_init[0] = body_force[0] * density
        press_proj_init[1] = body_force[1] * density
        press_proj_init[2] = body_force[2] * density
        for node in self.model_part.Nodes:
            eps = node.GetSolutionStepValue(POROSITY)
            node.SetSolutionStepValue(PRESS_PROJ, 0, press_proj_init * eps)
        print("entered in EdgeBasedLevelSetSubstepSolver initialize")

        self.keep_inlet_nodes = True

    def Initialize(self):
        print("entered in EdgeBasedLevelSetSubstepSolver python constructor")
        # build the edge data structure
        if(self.domain_size == 2):
            self.matrix_container = MatrixContainerC2C2D()
        else:
            self.matrix_container = MatrixContainerC2C3D()
        self.matrix_container.ConstructCSRVector(self.model_part)
        self.matrix_container.BuildCSRData(self.model_part)
        # for 3D problems we need to evaluate the condition's neighbours
        if(self.domain_size == 3):
            self.condition_neighbours_finder = FindConditionsNeighboursProcess(
                self.model_part, self.domain_size, 10)
            self.condition_neighbours_finder.Execute()
        # constructing the solver
        if(self.domain_size == 2):
            if(self.use_parallel_distance_calculation == False):
                self.distance_utils = SignedDistanceCalculationUtils2D()
            else:
                self.distance_utils = ParallelDistanceCalculator2D()

            self.fluid_solver = EdgeBasedLevelSetSubstep2D(
                self.matrix_container,
                self.model_part,
                self.viscosity,
                self.density,
                self.body_force,
                self.use_mass_correction,
                self.edge_detection_angle,
                self.stabdt_pressure_factor,
                self.stabdt_convection_factor,
                self.edge_detection_angle,
                self.assume_constant_pressure)
        else:
            if(self.use_parallel_distance_calculation == False):
                self.distance_utils = SignedDistanceCalculationUtils3D()
            else:
                self.distance_utils = ParallelDistanceCalculator3D()

            self.fluid_solver = EdgeBasedLevelSetSubstep3D(
                self.matrix_container,
                self.model_part,
                self.viscosity,
                self.density,
                self.body_force,
                self.use_mass_correction,
                self.edge_detection_angle,
                self.stabdt_pressure_factor,
                self.stabdt_convection_factor,
                self.edge_detection_angle,
                self.assume_constant_pressure)
#
        self.max_edge_size = self.distance_utils.FindMaximumEdgeSize(
            self.model_part)
        self.distance_size = self.max_edge_size * 3.0
        print("###################### max distance = ", self.distance_size)

        self.fluid_solver.SetShockCapturingCoefficient(0.1)  # 0.7)

#        self.reorder = True
#        self.distance_tools = BodyDistanceCalculationUtils()
        nneg = 0
        npos = 0
        for node in self.model_part.Nodes:
            if(node.GetSolutionStepValue(DISTANCE) < 0.0):
                nneg = nneg + 1
            else:
                npos = npos + 1

#        print "nneg=",nneg;
#        print "npos=",npos

        self.fluid_solver.Initialize()
#        print "in 128"
        nneg = 0
        npos = 0
        for node in self.model_part.Nodes:
            if(node.GetSolutionStepValue(DISTANCE) < 0.0):
                nneg = nneg + 1
            else:
                npos = npos + 1

        # saving inlet nodes
        self.inlet_nodes = []
        for node in self.model_part.Nodes:
            if(node.IsFixed(DISTANCE)):
                self.inlet_nodes.append(node)

#        print "nneg=",nneg;
#        print "npos=",npos
        self.Redistance()

#        for node in self.model_part.Nodes:
#            dist = node.GetSolutionStepValue(DISTANCE)
#            node.SetSolutionStepValue(DISTANCE,1,dist)
#        self.Redistance()
        self.measured_volume = 0.0
        self.expected_volume = self.fluid_solver.ComputeWetVolume()
        self.tot_volume = self.fluid_solver.ComputeTotalVolume()

        self.vol_variation = 0.0
        print("initial wet volume = ", self.expected_volume)

#        print "**********************************************"
        print("finished EdgeBasedLevelSetSubstepSolver initialize")

    #
    #
    def WettenNodes(self):
        for node in self.inlet_nodes:
            if(node.GetSolutionStepValue(DISTANCE) > 0):
                node.SetSolutionStepValue(DISTANCE, 0, -1.0e-3)

    #
    #
    def Redistance(self):
        if(self.keep_inlet_nodes):
            self.WettenNodes()

        if(self.use_parallel_distance_calculation == False):
            self.distance_utils.CalculateDistances(
                self.model_part,
                DISTANCE,
                self.distance_size)
        else:
            # print "max distance", self.distance_size
            # print "max extrapolation layers",self.extrapolation_layers
            self.distance_utils.CalculateDistances(
                self.model_part,
                DISTANCE,
                NODAL_AREA,
                self.extrapolation_layers,
                self.distance_size)

    #
    #
    def FluidOnlySolve(self):
        if (self.extrapolation_layers < 3):
            print("insufficient number of extrapolation layers. Minimum is 3")
            raise ValueError

        print("entered in EdgeBasedLevelSetSolver fluid only solve")
        (self.fluid_solver).ExtrapolateValues(self.extrapolation_layers)

        (self.fluid_solver).SolveStep1()
        (self.fluid_solver).SolveStep2(self.pressure_linear_solver)
        (self.fluid_solver).SolveStep3()

        (self.fluid_solver).ExtrapolateValues(self.extrapolation_layers)
        print("finished EdgeBasedLevelSetSolver fluid only solve")

    #
    #
    def Solve(self):
        self.AuxSolve(True)

    def AuxSolve(self, allow_redistancing):
        if (self.extrapolation_layers < 3):
            print("insufficient number of extrapolation layers. Minimum is 3")
            raise ValueError

        (self.fluid_solver).GatherValues()

        self.timer.Start("Update Fixed Velocity Values")
        (self.fluid_solver).UpdateFixedVelocityValues()
        self.timer.Stop("Update Fixed Velocity Values")

        self.timer.Start("Extrapolate Values")
        (self.fluid_solver).ExtrapolateValues(self.extrapolation_layers)
        self.timer.Stop("Extrapolate Values")

        # convect levelset function
       # self.convection_solver.Solve();
        self.timer.Start("Convect Distance")
        (self.fluid_solver).ConvectDistance()

        if(self.keep_inlet_nodes):
            self.WettenNodes()

        self.timer.Stop("Convect Distance")

        # correct volume
        if(self.use_mass_correction):
            self.timer.Start("MassCorrection")
            self.fluid_solver.ComputeWetVolume()
            self.measured_volume = self.fluid_solver.ComputeWetVolume()
            self.vol_variation = self.fluid_solver.ComputeVolumeVariation()
            self.expected_volume = self.expected_volume + self.vol_variation
            if(self.measured_volume >= 0.9999 * self.tot_volume):
                return 1
            # print "measured volume = ", measured_volume
            # print "vol_variation   = ",vol_variation
            # print "expected volume = ", self.expected_volume
            max_volume_error = 0.99
            # if(self.measured_volume / self.expected_volume < max_volume_error):
            # print "artificial mass correction"
            aaa = self.fluid_solver.ContinuousVolumeCorrection(
                self.expected_volume, self.measured_volume)
            self.timer.Stop("MassCorrection")

        if(self.step == self.redistance_frequency):
            self.timer.Start("Redistance")
            self.Redistance()
            self.timer.Stop("Redistance")
            self.step = 0
            print("redistance was executed")
        self.step += 1

        # solve fluid
        self.timer.Start("Solve Step 1")
        (self.fluid_solver).SolveStep1()
        self.timer.Stop("Solve Step 1")
        self.timer.Start("Solve Step 2")
        status = (self.fluid_solver).SolveStep2(self.pressure_linear_solver)
        self.timer.Stop("Solve Step 2")
        if(status == 0):  # everything went fine
            self.timer.Start("Solve Step 3")
            (self.fluid_solver).SolveStep3()
            self.timer.Stop("Solve Step 3")
        # something went wrong ... restart step and do redistance
        elif(allow_redistancing):
            print("*********************************** TIMESTEP REDUCTION EXECUTED ********************************")
            self.fluid_solver.ReduceTimeStep(
                self.model_part,
                self.model_part.ProcessInfo[TIME])
            self.step = self.redistance_frequency
            self.AuxSolve(False)

# if(self.step == self.redistance_frequency):
# self.Redistance()
# self.step = 0
# print "redistance was executed"
# self.step += 1
    #
#
#    def Solve(self):
#        (self.fluid_solver).ExtrapolateValues(self.extrapolation_layers)
#
# convect levelset function
#        (self.fluid_solver).ConvectDistance()
#
#        convection_success = (self.fluid_solver).CheckDistanceConvection()
#        if(convection_success == False):
# time step reduction is needed
# print "############### distance convection failed!! ###############"
#            return False
#            errrrrrrrrr
#
# solve fluid
#        (self.fluid_solver).SolveStep1();
#        (self.fluid_solver).SolveStep2(self.pressure_linear_solver);
#        (self.fluid_solver).SolveStep3();
#
#        if(self.step == self.redistance_frequency):
#            self.Redistance()
#            self.step = 0
#            print "redistance was executed"
#        self.step += 1
#
#        return True
    #
    #
    def EstimateTimeStep(self, safety_factor, max_Dt):
        dt = (self.fluid_solver).ComputeTimeStep(safety_factor, max_Dt)

        if(dt > max_Dt):
            dt = max_Dt

        # print dt

        return dt

    #
    #
    def CalculateInitialPressureDistribution(self):
        # prova!
        dt_aux = 1e-6
        self.model_part.ProcessInfo.SetValue(DELTA_TIME, dt_aux)
        (self.fluid_solver).SolveStep1()
        aaa = Vector(3)
        aaa[0] = self.body_force[0] * dt_aux
        aaa[1] = self.body_force[1] * dt_aux
        aaa[2] = self.body_force[2] * dt_aux
        for node in self.model_part.Nodes:
            if(node.IsFixed(VELOCITY_X) == False):
                node.SetSolutionStepValue(VELOCITY, 0, aaa)

#        for node in self.model_part.Nodes:
#            print node.GetSolutionStepValue(VELOCITY)

        (self.fluid_solver).SolveStep2(self.pressure_linear_solver)

        zero = Vector(3)
        zero[0] = 0.0
        zero[1] = 0.0
        zero[2] = 0.0
        for node in self.model_part.Nodes:
            if(node.IsFixed(VELOCITY_X) == False):
                node.SetSolutionStepValue(VELOCITY, 0, zero)
        self.model_part.ProcessInfo.SetValue(DELTA_TIME, 0.0)
