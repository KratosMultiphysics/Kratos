import os

import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as kratos_utils
from KratosMultiphysics.testing.utilities import ReadModelPart
from KratosMultiphysics import AssembleUtilities


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
        cls.mdpa_name = GetFilePath("auxiliar_files_for_python_unittest/mdpa_files/coarse_sphere_with_conditions")
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
        assemble_double_map, assemble_array3_map = TestAssembleUtilities.__generate_maps(self.model_part.GetCommunicator().GlobalNumberOfNodes())

        # assemble based on given nodal maps
        AssembleUtilities.AssembleCurrentDataWithValuesMap(self.model_part, KratosMultiphysics.DENSITY, assemble_double_map)
        AssembleUtilities.AssembleCurrentDataWithValuesMap(self.model_part, KratosMultiphysics.VELOCITY, assemble_array3_map)

        coefficient = self.model_part.GetCommunicator().TotalProcesses()

        # check for values
        for node in self.model_part.Nodes:
            node_id = node.Id

            if (node_id % 3 == 0):
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

    def test_AssembleNonHistoricalNodalDataWithValuesMap(self):
        entities = self.model_part.Nodes
        ## initialize entity data
        TestAssembleUtilities.__fill_entities(entities)

        # prepare the assemble maps
        assemble_double_map, assemble_array3_map = TestAssembleUtilities.__generate_maps(self.model_part.GetCommunicator().GlobalNumberOfNodes())

        # assemble based on given maps
        AssembleUtilities.AssembleNonHistoricalNodalDataWithValuesMap(self.model_part, KratosMultiphysics.PRESSURE, assemble_double_map)
        AssembleUtilities.AssembleNonHistoricalNodalDataWithValuesMap(self.model_part, KratosMultiphysics.DISPLACEMENT, assemble_array3_map)

        # check for values
        self.__check_entity_data(entities, assemble_double_map, assemble_array3_map)

    def test_AssembleElementDataWithValuesMap(self):
        entities = self.model_part.Elements
        ## initialize entity data
        TestAssembleUtilities.__fill_entities(entities)

        # prepare the assemble maps
        assemble_double_map, assemble_array3_map = TestAssembleUtilities.__generate_maps(self.model_part.GetCommunicator().GlobalNumberOfElements())

        # assemble based on given maps
        AssembleUtilities.AssembleElementDataWithValuesMap(self.model_part, KratosMultiphysics.PRESSURE, assemble_double_map)
        AssembleUtilities.AssembleElementDataWithValuesMap(self.model_part, KratosMultiphysics.DISPLACEMENT, assemble_array3_map)

        # check for values
        self.__check_entity_data(entities, assemble_double_map, assemble_array3_map)

    def test_AssembleConditionDataWithValuesMap(self):
        entities = self.model_part.Conditions
        ## initialize entity data
        TestAssembleUtilities.__fill_entities(entities)

        # prepare the assemble maps
        assemble_double_map, assemble_array3_map = TestAssembleUtilities.__generate_maps(self.model_part.GetCommunicator().GlobalNumberOfConditions())

        # assemble based on given maps
        AssembleUtilities.AssembleConditionDataWithValuesMap(self.model_part, KratosMultiphysics.PRESSURE, assemble_double_map)
        AssembleUtilities.AssembleConditionDataWithValuesMap(self.model_part, KratosMultiphysics.DISPLACEMENT, assemble_array3_map)

        # check for values
        self.__check_entity_data(entities, assemble_double_map, assemble_array3_map)

    def __check_entity_data(self, entities, assemble_double_map, assemble_array3_map):
        coefficient = self.model_part.GetCommunicator().TotalProcesses()

        for entity in entities:
            entity_id = entity.Id

            if (entity_id % 3 == 0):
                self.assertAlmostEqual(
                    entity.GetValue(KratosMultiphysics.PRESSURE),
                    entity_id * 2 + coefficient * assemble_double_map[entity_id], 12)
                self.assertVectorAlmostEqual(
                    entity.GetValue(KratosMultiphysics.DISPLACEMENT),
                    KratosMultiphysics.Array3([entity_id * 2, entity_id * 3, entity_id * 4]) + coefficient * assemble_array3_map[entity_id], 12)
            else:
                self.assertAlmostEqual(
                    entity.GetValue(KratosMultiphysics.PRESSURE),
                    entity_id * 2, 12)
                self.assertVectorAlmostEqual(
                    entity.GetValue(KratosMultiphysics.DISPLACEMENT),
                    KratosMultiphysics.Array3([entity_id * 2, entity_id * 3, entity_id * 4]), 12)

    @staticmethod
    def __generate_maps(number_of_entities):
        assemble_double_map = {}
        assemble_array3_map = {}
        for entity_id in range(1, number_of_entities + 1):
            if (entity_id % 3 == 0):
                assemble_double_map[entity_id] = entity_id * 2
                assemble_array3_map[entity_id] = KratosMultiphysics.Array3([entity_id * 2, entity_id * 3, 0.0])
        return assemble_double_map, assemble_array3_map

    @staticmethod
    def __fill_entities(entities):
        for entity in entities:
            entity_id = entity.Id
            entity.SetValue(KratosMultiphysics.PRESSURE, entity_id * 2)
            entity.SetValue(KratosMultiphysics.DISPLACEMENT, KratosMultiphysics.Array3([entity_id * 2, entity_id * 3, entity_id * 4]))


if __name__ == '__main__':
    KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)
    KratosUnittest.main()