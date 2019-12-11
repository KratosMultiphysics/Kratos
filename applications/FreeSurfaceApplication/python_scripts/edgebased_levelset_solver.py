from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7

# Importing the Kratos Library
import KratosMultiphysics as KM

# Import applications
import KratosMultiphysics.FreeSurfaceApplication as KratosFreeSurf



def AddVariables(model_part):
    model_part.AddNodalSolutionStepVariable(KM.VELOCITY)
    model_part.AddNodalSolutionStepVariable(KM.PRESSURE)
    model_part.AddNodalSolutionStepVariable(KM.NORMAL)
    model_part.AddNodalSolutionStepVariable(KM.AUX_INDEX)
    model_part.AddNodalSolutionStepVariable(KM.DISTANCE)
    model_part.AddNodalSolutionStepVariable(KM.PRESS_PROJ)
    model_part.AddNodalSolutionStepVariable(KM.POROSITY)
    model_part.AddNodalSolutionStepVariable(KM.VISCOSITY)
    model_part.AddNodalSolutionStepVariable(KM.DIAMETER)
    model_part.AddNodalSolutionStepVariable(KM.LIN_DARCY_COEF)
    model_part.AddNodalSolutionStepVariable(KM.NONLIN_DARCY_COEF)
    model_part.AddNodalSolutionStepVariable(KM.NODAL_AREA)
    model_part.AddNodalSolutionStepVariable(KM.STRUCTURE_VELOCITY)

    print("variables for the edgebased incompressible fluid solver added correctly")


def AddDofs(model_part):
    for node in model_part.Nodes:
        node.AddDof(KM.PRESSURE)
        node.AddDof(KM.VELOCITY_X)
        node.AddDof(KM.VELOCITY_Y)
        node.AddDof(KM.VELOCITY_Z)


class EdgeBasedLevelSetSolver:

    def __init__(self, model_part, domain_size,
                 body_force, viscosity, density):

        print("entered in EdgeBasedLevelSetSolver python constructor")
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
        self.timer = KM.Timer()

        self.use_parallel_distance_calculation = False
        # 0 = None; 1 = Ergun; 2 = Custom A y B;
        self.compute_porous_resistance_law = 0

        # neighbour search
        number_of_avg_elems = 10
        number_of_avg_nodes = 10
        self.neighbour_search = KM.FindNodalNeighboursProcess(
            model_part, number_of_avg_elems, number_of_avg_nodes)
        (self.neighbour_search).Execute()

        # erase isolated notes
        eliminate_isolated = KM.EliminateIsolatedNodesProcess(model_part)
        eliminate_isolated.Execute()

        # definition of the solvers
#        pDiagPrecond = DiagonalPreconditioner()
#        self.pressure_linear_solver =  CGSolver(1e-3, 5000,pDiagPrecond)
 #       self.pressure_linear_solver =  CGSolver(1e-3, 5000)

#        pDiagPrecond = DiagonalPreconditioner()
#        self.pressure_linear_solver =  BICGSTABSolver(1e-3, 5000,pDiagPrecond)
        self.pressure_linear_solver = KM.BICGSTABSolver(1e-3, 5000)

        # initializing the press proj to -body_force
        press_proj_init = KM.Vector(3)
        press_proj_init[0] = body_force[0] * density
        press_proj_init[1] = body_force[1] * density
        press_proj_init[2] = body_force[2] * density
        for node in self.model_part.Nodes:
            eps = node.GetSolutionStepValue(KM.POROSITY)
            node.SetSolutionStepValue(KM.PRESS_PROJ, 0, press_proj_init * eps)
        print("entered in EdgeBasedLevelSetSolver initialize")

    def Initialize(self):
        print("entered in EdgeBasedLevelSetSolver python constructor")
        # build the edge data structure
        if(self.domain_size == 2):
            self.matrix_container = KratosFreeSurf.MatrixContainer2D()
        else:
            self.matrix_container = KratosFreeSurf.MatrixContainer3D()
        self.matrix_container.ConstructCSRVector(self.model_part)
        self.matrix_container.BuildCSRData(self.model_part)
        # for 3D problems we need to evaluate the condition's neighbours
        if(self.domain_size == 3):
            self.condition_neighbours_finder = KM.FindConditionsNeighboursProcess(
                self.model_part, self.domain_size, 10)
            self.condition_neighbours_finder.Execute()
        # constructing the solver
        if(self.domain_size == 2):
            if(self.use_parallel_distance_calculation == False):
                self.distance_utils = KM.SignedDistanceCalculationUtils2D()
            else:
                self.distance_utils = KM.ParallelDistanceCalculator2D()

            self.fluid_solver = KratosFreeSurf.EdgeBasedLevelSet2D(
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
                self.distance_utils = KM.SignedDistanceCalculationUtils3D()
            else:
                self.distance_utils = KM.ParallelDistanceCalculator3D()

            self.fluid_solver = KratosFreeSurf.EdgeBasedLevelSet3D(
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

        self.fluid_solver.SetShockCapturingCoefficient(0.0)

#        self.reorder = True
#        self.distance_tools = BodyDistanceCalculationUtils()
        nneg = 0
        npos = 0
        for node in self.model_part.Nodes:
            if(node.GetSolutionStepValue(KM.DISTANCE) < 0.0):
                nneg = nneg + 1
            else:
                npos = npos + 1

        print("nneg=", nneg)
        print("npos=", npos)

        self.fluid_solver.Initialize()
        nneg = 0
        npos = 0
        for node in self.model_part.Nodes:
            if(node.GetSolutionStepValue(KM.DISTANCE) < 0.0):
                nneg = nneg + 1
            else:
                npos = npos + 1

        print("nneg=", nneg)
        print("npos=", npos)
        self.Redistance()

#        for node in self.model_part.Nodes:
#            dist = node.GetSolutionStepValue(KM.DISTANCE)
#            node.SetSolutionStepValue(KM.DISTANCE,1,dist)
#        self.Redistance()

        print("**********************************************")
        print("finished EdgeBasedLevelSetSolver initialize")

    #
    #
    def Redistance(self):
        if(self.use_parallel_distance_calculation == False):
            self.distance_utils.CalculateDistances(
                self.model_part,
                KM.DISTANCE,
                self.distance_size)
        else:
            print("max distance", self.distance_size)
            print("max extrapolation layers", self.extrapolation_layers)
            self.distance_utils.CalculateDistances(
                self.model_part,
                KM.DISTANCE,
                KM.NODAL_AREA,
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
        if (self.extrapolation_layers < 3):
            print("insufficient number of extrapolation layers. Minimum is 3")
            raise ValueError

        self.timer.Start("Calculate Porous Resistance Law")
        (self.fluid_solver).CalculatePorousResistanceLaw(
            self.compute_porous_resistance_law)
        self.timer.Stop("Calculate Porous Resistance Law")

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
        self.timer.Stop("Convect Distance")

        if(self.step == self.redistance_frequency):
            self.timer.Start("Redistance")
            self.Redistance()
            self.timer.Stop("Redistance")
            self.step = 0
            # print "redistance was executed"
        self.step += 1

        # solve fluid
        self.timer.Start("Solve Step 1")
        (self.fluid_solver).SolveStep1()
        self.timer.Stop("Solve Step 1")
        self.timer.Start("Solve Step 2")
        (self.fluid_solver).SolveStep2(self.pressure_linear_solver)
        self.timer.Stop("Solve Step 2")
        self.timer.Start("Solve Step 3")
        (self.fluid_solver).SolveStep3()
        self.timer.Stop("Solve Step 3")

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
    def EstimateBoundedTimeStep(self, safety_factor, max_Dt):
        dt = (self.fluid_solver).ComputeBoundedTimeStep(safety_factor, max_Dt)

        if(dt > max_Dt):
            dt = max_Dt

        # print dt

        return dt

    #
    #
    def CalculateInitialPressureDistribution(self):
        # prova!
        dt_aux = 1e-6
        self.model_part.ProcessInfo.SetValue(KM.DELTA_TIME, dt_aux)
        (self.fluid_solver).SolveStep1()
        aaa = KM.Vector(3)
        aaa[0] = self.body_force[0] * dt_aux
        aaa[1] = self.body_force[1] * dt_aux
        aaa[2] = self.body_force[2] * dt_aux
        for node in self.model_part.Nodes:
            if(node.IsFixed(KM.VELOCITY_X) == False):
                node.SetSolutionStepValue(KM.VELOCITY, 0, aaa)

#        for node in self.model_part.Nodes:
#            print node.GetSolutionStepValue(KM.VELOCITY)

        (self.fluid_solver).SolveStep2(self.pressure_linear_solver)

        zero = Vector(3)
        zero[0] = 0.0
        zero[1] = 0.0
        zero[2] = 0.0
        for node in self.model_part.Nodes:
            if(node.IsFixed(KM.VELOCITY_X) == False):
                node.SetSolutionStepValue(KM.VELOCITY, 0, zero)
        self.model_part.ProcessInfo.SetValue(KM.DELTA_TIME, 0.0)
