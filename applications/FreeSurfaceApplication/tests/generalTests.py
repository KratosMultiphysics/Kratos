import KratosMultiphysics
import KratosMultiphysics.FreeSurfaceApplication as KratosFreeSurface
from KratosMultiphysics.FreeSurfaceApplication.free_surface_analysis import FreeSurfaceAnalysis

import KratosMultiphysics.KratosUnittest as KratosUnittest
import pathlib


def GetFilePath(fileName):
    return os.path.dirname(__file__) + "/" + fileName


class KratosFreeSurfaceGeneralTests(KratosUnittest.TestCase):

    def ExecuteReservoirTests(self):
        with KratosUnittest.WorkFolderScope(self.work_folder, __file__):
            self.setUp()
            self.setUpProblem()
            self.runTest()
            self.tearDown()
            self.checkResults(self.reference_file)

    def test_Reservoir2D(self):
        self.__RunFreeSurfaceAnalysis("Reservoir2D", "reference_reservoir_2D")

    def test_Reservoir3D(self):
        self.__RunFreeSurfaceAnalysis("Reservoir3D", "reference_reservoir_3D")

    def __RunFreeSurfaceAnalysis(self, test_name: str, reference_name: str):
        with KratosUnittest.WorkFolderScope(test_name, __file__):
            parameters_path = pathlib.Path("ProjectParameters.json")
            model_part_path = pathlib.Path(test_name + ".mdpa")
            with open(parameters_path, 'r') as file:
                parameters = KratosMultiphysics.Parameters(file.read())

            model = KratosMultiphysics.Model()

            analysis = FreeSurfaceAnalysis(model, parameters)
            analysis.Run()

            self.main_model_part = analysis._GetSolver().GetComputingModelPart()
            self.checkResults(reference_name)

    def setUp(self):
        self.check_tolerance = 1e-5
        self.print_reference_values = False

    def checkResults(self, reference_file_name: str):
        if self.print_reference_values:
            with open(reference_file_name + '.csv','w') as ref_file:
                ref_file.write("#ID, DISTANCE\n")
                for node in self.main_model_part.Nodes:
                    dist = node.GetSolutionStepValue(KratosMultiphysics.DISTANCE,0)
                    ref_file.write("{0}, {1}\n".format(node.Id, dist))
        else:
            with open(reference_file_name + '.csv','r') as reference_file:
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


if __name__ == "__main__":
    KratosUnittest.main()