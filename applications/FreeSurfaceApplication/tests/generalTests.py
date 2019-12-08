import KratosMultiphysics
import KratosMultiphysics.FreeSurfaceApplication as KratosFreeSurface

import KratosMultiphysics.KratosUnittest as KratosUnittest


def GetFilePath(fileName):
    return os.path.dirname(__file__) + "/" + fileName


class KratosFreeSurfaceGeneralTests(KratosUnittest.TestCase):

    def testReservoir2D(self):
        self.domain_size = 2
        self.safety_factor = 0.1
        self.body_force = KratosMultiphysics.Vector(3)
        self.body_force[0] = 0.0
        self.body_force[1] = -9.81
        self.body_force[2] = 0.0
        self.work_folder = "Reservoir2D"
        self.problem_name = "Reservoir2D"
        self.reference_file = "reference_reservoir_2D"
        self.ExecuteReservoirTests()

    def testReservoir3D(self):
        self.domain_size = 3
        self.safety_factor = 0.5
        self.body_force = KratosMultiphysics.Vector(3)
        self.body_force[0] = 0.0
        self.body_force[1] = 0.0
        self.body_force[2] = -9.81
        self.work_folder = "Reservoir3D"
        self.problem_name = "Reservoir3D"
        self.reference_file = "reference_reservoir_3D"
        self.ExecuteReservoirTests()

    def ExecuteReservoirTests(self):
        with KratosUnittest.WorkFolderScope(self.work_folder, __file__):
            self.setUp()
            self.setUpProblem()
            self.runTest()
            self.tearDown()
            self.checkResults()

    def setUp(self):
        self.check_tolerance = 1e-6
        self.print_reference_values = False

    def setUpProblem(self):
        pass

    def runTest(self):
        # Set number of OMP threads
        parallel=KratosMultiphysics.OpenMPUtils()
        parallel.SetNumThreads(1)

        # setting the domain size for the problem to be solved
        domain_size = self.domain_size

        # defining a model part for the fluid and one for the structure
        FreeSurfaceModel = KratosMultiphysics.Model()
        self.main_model_part = FreeSurfaceModel.CreateModelPart("FluidPart")

        # importing the solvers needed
        from KratosMultiphysics.FreeSurfaceApplication.edgebased_levelset_solver import EdgeBasedLevelSetSolver
        EdgeBasedLevelSetSolver.AddVariables(self.main_model_part)

        # introducing input file name
        input_file_name = self.problem_name

        model_part_io_fluid = KratosMultiphysics.ModelPartIO(input_file_name)
        model_part_io_fluid.ReadModelPart(self.main_model_part)

        # setting up the buffer size: SHOULD BE DONE AFTER READING!!!
        self.main_model_part.SetBufferSize(2)

        # adding dofs
        EdgeBasedLevelSetSolver.AddDofs(self.main_model_part)

        # we assume here that all of the internal nodes are marked with a negative distance
        # set the distance of all of the internal nodes to a small value
        small_value = 0.0001
        n_active = 0
        for node in self.main_model_part.Nodes:
            dist = node.GetSolutionStepValue(KratosMultiphysics.DISTANCE)
            if(dist < 0.0):
                n_active = n_active + 1
                node.SetSolutionStepValue(KratosMultiphysics.DISTANCE, 0, -small_value)
            else:
                node.SetSolutionStepValue(KratosMultiphysics.DISTANCE, 0, small_value)

        if(n_active == 0):
            raise "ERROR. At least one node has to be initialized with a distance lesser than 0"

        # make sure that the porosity is not zero on any node (set by default to fluid only)
        for node in self.main_model_part.Nodes:
            if(node.GetSolutionStepValue(KratosMultiphysics.POROSITY) == 0.0):
                node.SetSolutionStepValue(KratosMultiphysics.POROSITY, 0, 1.0)
            if(node.GetSolutionStepValue(KratosMultiphysics.DIAMETER) == 0.0):
                node.SetSolutionStepValue(KratosMultiphysics.DIAMETER, 0, 1.0)

        # constructing the solver
        body_force = KratosMultiphysics.Vector(3)
        body_force = self.body_force
        if(body_force[0] == 0.0 and body_force[1] == 0.0 and body_force[2] == 0.0):
            raise "ERROR. Body Force cannot be a ZERO VECTOR"

        viscosity = 1.0e-06
        density = 1000.0
        fluid_solver = EdgeBasedLevelSetSolver(self.main_model_part, domain_size, body_force, viscosity, density)
        fluid_solver.redistance_frequency = 5.0
        fluid_solver.extrapolation_layers = int(5)
        fluid_solver.stabdt_pressure_factor = 1.0
        fluid_solver.stabdt_convection_factor = 1.0e-02
        fluid_solver.use_mass_correction = True
        fluid_solver.tau2_factor = 1.0
        fluid_solver.edge_detection_angle = 45.0
        fluid_solver.assume_constant_pressure = False
        fluid_solver.compute_porous_resistance_law = 0

        fluid_solver.Initialize()

        # settings to be changed
        max_Dt = 1.0e-02
        final_time = 5.0
        safety_factor = self.safety_factor

        number_of_inital_steps = 10.0
        initial_time_step = 1.0e-05

        original_max_dt = max_Dt

        max_safety_factor = safety_factor

        time = 0.0
        step = 0

        while(time < final_time):

            if(step < number_of_inital_steps):
                max_Dt = initial_time_step
            else:
                max_Dt = original_max_dt
                # progressively increment the safety factor
                # in the steps that follow a reduction of it
                safety_factor = safety_factor * 1.2
                if(safety_factor > max_safety_factor):
                    safety_factor = max_safety_factor

            Dt = fluid_solver.EstimateTimeStep(safety_factor, max_Dt)

            time = time + Dt
            self.main_model_part.CloneTimeStep(time)

            if(step >= 3):
                fluid_solver.Solve()

                check_dt = fluid_solver.EstimateTimeStep(0.95, max_Dt)

                if(check_dt < Dt):

                    # we found a velocity too large! we need to reduce the time step
                    fluid_solver.ReduceTimeStep(self.main_model_part, time)  # this is to set the database to the value at the beginning of the step

                    safety_factor *= 3.00000e-01
                    reduced_dt = fluid_solver.EstimateTimeStep(safety_factor, max_Dt)

                    time = time - Dt + reduced_dt

                    fluid_solver.ReduceTimeStep(self.main_model_part, time)  # this is to set the database to the value at the beginning of the step

                    fluid_solver.Solve()

            step = step + 1

    def tearDown(self):
        pass

    def checkResults(self):
        with KratosUnittest.WorkFolderScope(self.work_folder, __file__):
                if self.print_reference_values:
                    with open(self.reference_file+'.csv','w') as ref_file:
                        ref_file.write("#ID, DISTANCE\n")
                        for node in self.main_model_part.Nodes:
                            dist = node.GetSolutionStepValue(KratosMultiphysics.DISTANCE,0)
                            ref_file.write("{0}, {1}\n".format(node.Id, dist))
                else:
                    with open(self.reference_file+'.csv','r') as reference_file:
                        reference_file.readline() # skip header
                        line = reference_file.readline()

                        for node in self.main_model_part.Nodes:
                            values = [ float(i) for i in line.rstrip('\n ').split(',') ]
                            reference_dist = values[1]

                            distance = node.GetSolutionStepValue(KratosMultiphysics.DISTANCE)
                            self.assertAlmostEqual(reference_dist, distance, delta = self.check_tolerance)

                            line = reference_file.readline()
                        if line != '': # If we did not reach the end of the reference file
                            self.fail("The number of nodes in the mdpa is smaller than the number of nodes in the output file")
