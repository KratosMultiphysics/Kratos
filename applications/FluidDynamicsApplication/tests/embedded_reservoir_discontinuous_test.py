import KratosMultiphysics
import KratosMultiphysics.FluidDynamicsApplication
try:
    import KratosMultiphysics.ExternalSolversApplication
    have_external_solvers = True
except ImportError as e:
    have_external_solvers = False


import KratosMultiphysics.KratosUnittest as UnitTest

import os

class WorkFolderScope:
    def __init__(self, work_folder):
        self.currentPath = os.getcwd()
        self.scope = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)),work_folder))

    def __enter__(self):
        os.chdir(self.scope)

    def __exit__(self, type, value, traceback):
        os.chdir(self.currentPath)

@UnitTest.skipUnless(have_external_solvers,"Missing required application: ExternalSolversApplication")
class EmbeddedReservoirDiscontinuousTest(UnitTest.TestCase):
    def testEmbeddedReservoirDiscontinuous3D(self):
        self.distance = 0.99
        self.slip_level_set = True
        self.work_folder = "EmbeddedReservoirTest"
        self.reference_file = "reference_slip_reservoir_3D"
        self.settings = "EmbeddedReservoirDiscontinuous3DTestParameters.json"
        self.ExecuteEmbeddedReservoirTest()

    def ExecuteEmbeddedReservoirTest(self):
        with WorkFolderScope(self.work_folder):
            self.setUp()
            self.setUpProblem()
            self.runTest()
            self.tearDown()
            self.checkResults()

    def setUp(self):
        self.check_tolerance = 1.0e-8
        self.print_output = False
        self.print_reference_values = False

    def tearDown(self):
        with WorkFolderScope(self.work_folder):
            try:
                os.remove(self.ProjectParameters["solver_settings"]["model_import_settings"]["input_filename"].GetString()+'.time')
            except FileNotFoundError as e:
                pass

    def setUpProblem(self):
        with WorkFolderScope(self.work_folder):
            with open(self.settings, 'r') as parameter_file:
                self.ProjectParameters = KratosMultiphysics.Parameters(parameter_file.read())

            self.Model = KratosMultiphysics.Model()
            import fluid_dynamics_analysis
            self.simulation = fluid_dynamics_analysis.FluidDynamicsAnalysis(self.Model, self.ProjectParameters)

    def setUpDistanceField(self):
        # Get the model part containing the domain
        fluid_model_part = self.simulation._GetSolver().main_model_part

        # Set continuous distance field
        for node in fluid_model_part.Nodes:
            d = self.distance - node.Z
            node.SetSolutionStepValue(KratosMultiphysics.DISTANCE, d)

        # Set discontinuous distance field
        for element in fluid_model_part.Elements:
            i_node = 0
            elem_dist = KratosMultiphysics.Vector(4)
            for node in element.GetNodes():
                elem_dist[i_node] = node.GetSolutionStepValue(KratosMultiphysics.DISTANCE)
                i_node += 1
            element.SetValue(KratosMultiphysics.ELEMENTAL_DISTANCES, elem_dist)

    def setUpInitialCondition(self):
        # Set exact initial solution
        v_zero = KratosMultiphysics.Vector(3,0.0)
        for node in self.simulation._GetSolver().main_model_part.Nodes:
            node.SetSolutionStepValue(KratosMultiphysics.VELOCITY, 0, v_zero)
            node.SetSolutionStepValue(KratosMultiphysics.VELOCITY, 1, v_zero)
            node.SetSolutionStepValue(KratosMultiphysics.VELOCITY, 2, v_zero)
            if node.GetSolutionStepValue(KratosMultiphysics.DISTANCE) < 0.0:
                node.SetSolutionStepValue(KratosMultiphysics.PRESSURE, 0, 1.0e+06)
                node.SetSolutionStepValue(KratosMultiphysics.PRESSURE, 1, 1.0e+06)
                node.SetSolutionStepValue(KratosMultiphysics.PRESSURE, 2, 1.0e+06)
            else:
                node.SetSolutionStepValue(KratosMultiphysics.PRESSURE, 0, 0.0)
                node.SetSolutionStepValue(KratosMultiphysics.PRESSURE, 1, 0.0)
                node.SetSolutionStepValue(KratosMultiphysics.PRESSURE, 2, 0.0)

    def runTest(self):
        with WorkFolderScope(self.work_folder):
            # Set up the test
            self.simulation.Initialize()
            self.setUpDistanceField()
            self.setUpInitialCondition()

            # Run the test
            self.simulation.RunSolutionLoop()

            # Finalize the test
            self.simulation.Finalize()

    def checkResults(self):
        with WorkFolderScope(self.work_folder):
            fluid_model_part = self.simulation._GetSolver().main_model_part
            if self.print_reference_values:
                with open(self.reference_file+'.csv','w') as ref_file:
                    ref_file.write("#ID, PRESSURE, VELOCITY_X, VELOCITY_Y, VELOCITY_Z\n")
                    for node in fluid_model_part.Nodes:
                        pres = node.GetSolutionStepValue(KratosMultiphysics.PRESSURE)
                        v_x = node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_X)
                        v_y = node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_Y)
                        v_z = node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_Z)
                        ref_file.write("{0}, {1}, {2}, {3}, {4}\n".format(node.Id, pres, v_x, v_y, v_z))
            else:
                with open(self.reference_file+'.csv','r') as reference_file:
                    reference_file.readline() # skip header
                    line = reference_file.readline()

                    for node in fluid_model_part.Nodes:
                        values = [ float(i) for i in line.rstrip('\n ').split(',') ]
                        reference_pres = values[1]
                        reference_v_x = values[2]
                        reference_v_y = values[3]
                        reference_v_z = values[4]

                        pres = node.GetSolutionStepValue(KratosMultiphysics.PRESSURE)
                        v_x = node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_X)
                        v_y = node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_Y)
                        v_z = node.GetSolutionStepValue(KratosMultiphysics.VELOCITY_Z)
                        self.assertAlmostEqual(reference_v_x, v_x, delta = self.check_tolerance)
                        self.assertAlmostEqual(reference_v_y, v_y, delta = self.check_tolerance)
                        self.assertAlmostEqual(reference_v_z, v_z, delta = self.check_tolerance)
                        self.assertAlmostEqual(reference_pres, pres, delta = self.check_tolerance)

                        line = reference_file.readline()
                    if line != '': # If we did not reach the end of the reference file
                        self.fail("The number of nodes in the mdpa is smaller than the number of nodes in the output file")

if __name__ == '__main__':
    test = EmbeddedReservoirDiscontinuousTest()
    test.setUp()
    test.distance = 0.99
    test.slip_level_set = True
    test.print_output = False
    test.print_reference_values = False
    test.work_folder = "EmbeddedReservoirDiscontinuousTest"
    test.reference_file = "reference_embedded_reservoir_discontinuous_3D"
    test.settings = "EmbeddedReservoirDiscontinuous3DTestParameters.json"
    test.setUpProblem()
    test.runTest()
    test.tearDown()
    test.checkResults()

