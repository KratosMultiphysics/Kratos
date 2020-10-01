import os

import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as kratos_utils
from KratosMultiphysics.testing.utilities import ReadModelPart

from KratosMultiphysics import IsDistributedRun

if (IsDistributedRun()):
    from KratosMultiphysics.mpi import MPIAssembleUtilities as assemble_utilities
else:
    from KratosMultiphysics import AssembleUtilities as assemble_utilities


def GetFilePath(fileName):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), fileName)

def RemoveFiles(mdpa_name):
    kratos_utils.DeleteFileIfExisting(mdpa_name + ".time")

class TestAssembleUtilities(KratosUnittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.current_model = KratosMultiphysics.Model()
        cls.model_part = cls.current_model.CreateModelPart("Main")
        cls.model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] = 3
        cls.model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DENSITY)
        cls.model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VELOCITY)
        cls.mdpa_name = GetFilePath("auxiliar_files_for_python_unittest/mdpa_files/three_dim_symmetrical_cube")
        ReadModelPart(cls.mdpa_name, cls.model_part)

    @classmethod
    def tearDownClass(cls):
        RemoveFiles(cls.mdpa_name)

    def test_AssembleCurrentDataWithValuesMap(self):
        ## initialize nodal data
        for node in self.model_part.Nodes:
            node_id = node.Id
            node.SetSolutionStepValue(KratosMultiphysics.DENSITY, 0, node_id)
            node.SetSolutionStepValue(KratosMultiphysics.VELOCITY, 0, KratosMultiphysics.Array3([node_id, node_id * 2, node_id * 3]))

        # prepare the assemble maps
        assemble_double_map, assemble_array3_map = self.__generate_maps()

        # assemble based on given nodal maps
        assemble_utilities().AssembleCurrentDataWithValuesMap(self.model_part, KratosMultiphysics.DENSITY, assemble_double_map)
        assemble_utilities().AssembleCurrentDataWithValuesMap(self.model_part, KratosMultiphysics.VELOCITY, assemble_array3_map)

        coefficient = self.model_part.GetCommunicator().TotalProcesses()

        # check for values
        for node in self.model_part.Nodes:
            node_id = node.Id

            if (node_id % 10 == 0):
                self.assertAlmostEqual(
                    node.GetSolutionStepValue(KratosMultiphysics.DENSITY),
                    node_id + coefficient * assemble_double_map[node_id], 12)
                self.assertVectorAlmostEqual(
                    node.GetSolutionStepValue(KratosMultiphysics.VELOCITY),
                    KratosMultiphysics.Array3([node_id, node_id * 2, node_id * 3]) + coefficient * assemble_array3_map[node_id], 12)
            else:
                self.assertAlmostEqual(
                    node.GetSolutionStepValue(KratosMultiphysics.DENSITY),
                    node_id, 12)
                self.assertVectorAlmostEqual(
                    node.GetSolutionStepValue(KratosMultiphysics.VELOCITY),
                    KratosMultiphysics.Array3([node_id, node_id * 2, node_id * 3]), 12)

    def test_AssembleNonHistoricalDataWithValuesMap(self):
        ## initialize nodal data
        for node in self.model_part.Nodes:
            node_id = node.Id
            node.SetValue(KratosMultiphysics.PRESSURE, node_id * 2)
            node.SetValue(KratosMultiphysics.DISPLACEMENT, KratosMultiphysics.Array3([node_id * 2, node_id * 3, node_id * 4]))

        # prepare the assemble maps
        assemble_double_map, assemble_array3_map = self.__generate_maps()

        # assemble based on given nodal maps
        assemble_utilities().AssembleNonHistoricalDataWithValuesMap(self.model_part, KratosMultiphysics.PRESSURE, assemble_double_map)
        assemble_utilities().AssembleNonHistoricalDataWithValuesMap(self.model_part, KratosMultiphysics.DISPLACEMENT, assemble_array3_map)

        coefficient = self.model_part.GetCommunicator().TotalProcesses()

        # check for values
        for node in self.model_part.Nodes:
            node_id = node.Id

            if (node_id % 10 == 0):
                self.assertAlmostEqual(
                    node.GetValue(KratosMultiphysics.PRESSURE),
                    node_id * 2 + coefficient * assemble_double_map[node_id], 12)
                self.assertVectorAlmostEqual(
                    node.GetValue(KratosMultiphysics.DISPLACEMENT),
                    KratosMultiphysics.Array3([node_id * 2, node_id * 3, node_id * 4]) + coefficient * assemble_array3_map[node_id], 12)
            else:
                self.assertAlmostEqual(
                    node.GetValue(KratosMultiphysics.PRESSURE),
                    node_id * 2, 12)
                self.assertVectorAlmostEqual(
                    node.GetValue(KratosMultiphysics.DISPLACEMENT),
                    KratosMultiphysics.Array3([node_id * 2, node_id * 3, node_id * 4]), 12)

    def __generate_maps(self):
        number_of_nodes = self.model_part.GetCommunicator().GlobalNumberOfNodes()
        assemble_double_map = {}
        assemble_array3_map = {}
        for node_id in range(1, number_of_nodes):
            if (node_id % 10 == 0):
                assemble_double_map[node_id] = node_id * 2
                assemble_array3_map[node_id] = KratosMultiphysics.Array3([node_id * 2, node_id * 3, 0.0])
        return assemble_double_map, assemble_array3_map

if __name__ == '__main__':
    KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)
    KratosUnittest.main()