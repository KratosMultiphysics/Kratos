from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
# importing the Kratos Library
from KratosMultiphysics import *
from KratosMultiphysics.OpenCLApplication import *
from KratosMultiphysics.IncompressibleFluidApplication import *
CheckForPreviousImport()

import time as timer


def AddVariables(model_part):
    model_part.AddNodalSolutionStepVariable(VELOCITY)
    model_part.AddNodalSolutionStepVariable(PRESSURE)
    model_part.AddNodalSolutionStepVariable(NORMAL)
    model_part.AddNodalSolutionStepVariable(AUX_INDEX)
    model_part.AddNodalSolutionStepVariable(PRESS_PROJ)

    print("variables for the edgebased incompressible fluid solver added correctly")


def AddDofs(model_part):
    for node in model_part.Nodes:
        node.AddDof(PRESSURE)
        node.AddDof(VELOCITY_X)
        node.AddDof(VELOCITY_Y)
        node.AddDof(VELOCITY_Z);


class OpenClSolver:

    def __init__(self, model_part, domain_size, body_force, viscosity, density, cl_source_path):

        print("entered in EdgeBasedEulerianSolver python constructor")
        # data of the problem
        self.model_part = model_part
        self.domain_size = domain_size
        self.body_force = body_force
        self.density = density
        self.viscosity = viscosity

        self.use_mass_correction = True

        self.stabdt_pressure_factor = 1.0;
        self.stabdt_convection_factor = 0.01;

        self.redistance_frequency = 5;
        self.step = 0;

        self.extrapolation_layers = 5
        self.tau2_factor = 1.0
        self.edge_detection_angle = 45.0
        self.assume_constant_pressure = False

        # neighbour search
        number_of_avg_elems = 10
        number_of_avg_nodes = 10
        self.neighbour_search = FindNodalNeighboursProcess(model_part, number_of_avg_elems, number_of_avg_nodes)
        (self.neighbour_search).Execute()

        self.cl_source_path = cl_source_path
        self.single_device_flag = True

        self.tot_solve_time = 0.0
        self.cpu_to_gpu_time = 0.0
        self.step1_time = 0.0
        self.step2_time = 0.0
        self.step3_time = 0.0
        self.gpu_to_cpu_time = 0.0
        self.total_solves = 0

    def Initialize(self):

        print(NORMAL)
        print(IS_STRUCTURE)

        if(self.domain_size == 2):
            raise "error! only 3D implemented"

        print("entered in intialize")
        self.device_group = OpenCLDeviceGroup(cl_device_type.CL_DEVICE_TYPE_GPU, self.single_device_flag)
        self.device_group.AddCLSearchPath(self.cl_source_path)

        self.matrix_container = OpenCLMatrixContainer3D(self.device_group)
        self.matrix_container.ConstructCSRVector(self.model_part)
        self.matrix_container.BuildCSRData(self.model_part)
        print("contructed matrix_containers")

        self.condition_neighbours_finder = FindConditionsNeighboursProcess(self.model_part, self.domain_size, 10)
        self.condition_neighbours_finder.Execute()

        # constructing the solver
        self.fluid_solver = OpenCLFluidSolver3D(self.matrix_container, self.model_part, self.viscosity, self.density, self.body_force, self.use_mass_correction, self.edge_detection_angle, self.stabdt_pressure_factor, self.stabdt_convection_factor, self.edge_detection_angle, self.assume_constant_pressure)

        self.fluid_solver.Initialize()

        print("**********************************************")
        print("finished OpenCLNSFluidSolver Initialize")

    #
    #
    def Solve(self):
        t0 = timer.time()
        (self.fluid_solver).LoadDataToGPU();
# (self.fluid_solver).UpdateFixedVelocityValues()

        t1 = timer.time()
        (self.fluid_solver).SolveStep1();

        t2 = timer.time()
        (self.fluid_solver).SolveStep2();

        t3 = timer.time()
        (self.fluid_solver).SolveStep3();

        t4 = timer.time()
        (self.fluid_solver).WriteDataToCPU();

        t5 = timer.time()
        tot = t5 - t0
        print("TOTAL STEP	  time --->", t5 - t0)
        print("CPU->GPU transfer time --->", t1 - t0, "tot % -->", (t1 - t0) / tot)
        print("Step1		  time --->", t2 - t1, "tot % -->", (t2 - t1) / tot)
        print("Step2		  time --->", t3 - t2, "tot % -->", (t3 - t2) / tot)
        print("Step3		  time --->", t4 - t3, "tot % -->", (t4 - t3) / tot)
        print("GPU->CPU transfer time --->", t5 - t4, "tot % -->", (t5 - t4) / tot)
        self.tot_solve_time += tot
        self.cpu_to_gpu_time += t1 - t0
        self.step1_time += t2 - t1
        self.step2_time += t3 - t2
        self.step3_time += t4 - t3
        self.gpu_to_cpu_time += t5 - t4
        self.total_solves += 1

    #
    #
    def EstimateTimeStep(self, safety_factor, max_Dt):
        dt = (self.fluid_solver).ComputeTimeStep(safety_factor, max_Dt);

        if(dt > max_Dt):
            dt = max_Dt

        print(dt)

        return dt

    #
    #
    def PrintTimings(self):
        print("FINAL TIMINGS IN SOLVE (sum of all step timings)")
        print("TOTAL STEP	  time --->", self.tot_solve_time)
        print("CPU->GPU transfer time --->", self.cpu_to_gpu_time, "tot % -->", self.cpu_to_gpu_time / self.tot_solve_time)
        print("Step1		  time --->", self.step1_time, "tot % -->", self.step1_time / self.tot_solve_time)
        print("Step2		  time --->", self.step2_time, "tot % -->", self.step2_time / self.tot_solve_time)
        print("Step3		  time --->", self.step3_time, "tot % -->", self.step3_time / self.tot_solve_time)
        print("GPU->CPU transfer time --->", self.gpu_to_cpu_time, "tot % -->", self.gpu_to_cpu_time / self.tot_solve_time)
        print("total solves performed --->", self.total_solves)
