from __future__ import print_function, absolute_import, division  # makes these scripts backward compatible with python 2.6 and 2.7

import KratosMultiphysics as KM
import KratosMultiphysics.KratosUnittest as KratosUnittest

from KratosMultiphysics.CoSimulationApplication.processes import create_point_based_entites_process


if KM.IsDistributedRun():
    import KratosMultiphysics.mpi as KratosMPI

class TestCreatePointBasedEntitiesProcess(KratosUnittest.TestCase):

    def setUp(self):
        self.model = KM.Model()
        self.root_model_part = self.model.CreateModelPart("root_mp")
        self.root_model_part.AddNodalSolutionStepVariable(KM.PRESSURE)
        self.root_model_part.AddNodalSolutionStepVariable(KM.PARTITION_INDEX)
        self.dimension = 3

        self.root_model_part.CreateNewProperties(0)

        data_comm = KM.DataCommunicator.GetDefault()
        self.my_pid = data_comm.Rank()
        self.num_nodes = self.my_pid % 5 + 4 # num_nodes in range (4 ... 8)
        if self.my_pid == 4:
            self.num_nodes = 0 # in order to emulate one partition not having local nodes

        smp_nodes_1 = self.root_model_part.CreateSubModelPart("smp_nodes_1")
        smp_nodes_2 = self.root_model_part.CreateSubModelPart("smp_nodes_2")
        smp_nodes_3 = self.root_model_part.CreateSubModelPart("smp_nodes_3")

        scan_sum_num_nodes = data_comm.ScanSum(self.num_nodes)

        if self.num_nodes > 0:
            for i_node, id_node in enumerate(range(scan_sum_num_nodes-self.num_nodes, scan_sum_num_nodes)):
                # this creates the same coords in different ranks, which does not matter for this test
                if i_node < self.num_nodes-2: # last two nodes go to other sub-model-part
                    smp_nodes_1.CreateNewNode(id_node+1, 0.1*i_node, 0.0, 0.0)
                else:
                    smp_nodes_2.CreateNewNode(id_node+1, 1.3, 0.1*i_node, 0.0)

        for node in self.root_model_part.Nodes:
            node.SetSolutionStepValue(KM.PARTITION_INDEX, self.my_pid)

        # adding the first node of each smp to another smp to emulate an "overlapping" interface
        for node in smp_nodes_1.Nodes:
            smp_nodes_3.AddNodes([node.Id])
            break
        for node in smp_nodes_2.Nodes:
            smp_nodes_3.AddNodes([node.Id])
            break

        if KM.IsDistributedRun():
            KratosMPI.ParallelFillCommunicator(self.root_model_part).Execute()

    def test_create_entities_from_one_model_part(self):
        settings = KM.Parameters("""{
            "Parameters" : {
                "root_model_part_name"       : "root_mp",
                "new_sub_model_part_name"    : "smp_with_conditions",
                "sub_model_part_names"       : ["smp_nodes_1"],
                "entity_name"                : "PointCondition2D1N",
                "entity_type"                : "condition",
                "properties_id"              : 0
            }
        }""")

        self.process = self.__CreateProcess(settings)
        self.assertTrue(self.root_model_part.HasSubModelPart("smp_with_conditions"))
        self.assertEqual(self.root_model_part.GetSubModelPart("smp_nodes_1").NumberOfNodes(), self.root_model_part.NumberOfConditions())

        self.__CheckCreatedEntitiesIdAreCorrectlyNumbered(self.root_model_part.Conditions)

    def test_create_entities_from_multiple_model_parts(self):
        settings = KM.Parameters("""{
            "Parameters" : {
                "root_model_part_name"       : "root_mp",
                "new_sub_model_part_name"    : "smp_with_conditions",
                "sub_model_part_names"       : ["smp_nodes_1", "smp_nodes_2"],
                "entity_name"                : "PointCondition2D1N",
                "entity_type"                : "condition",
                "properties_id"              : 0
            }
        }""")

        self.process = self.__CreateProcess(settings)
        self.assertTrue(self.root_model_part.HasSubModelPart("smp_with_conditions"))
        self.assertEqual(self.root_model_part.NumberOfNodes(), self.root_model_part.NumberOfConditions())

        self.__CheckCreatedEntitiesIdAreCorrectlyNumbered(self.root_model_part.Conditions)

    def test_create_entities_from_overlapping_model_parts(self):
        settings = KM.Parameters("""{
            "Parameters" : {
                "root_model_part_name"       : "root_mp",
                "new_sub_model_part_name"    : "smp_with_conditions",
                "sub_model_part_names"       : ["smp_nodes_1", "smp_nodes_2", "smp_nodes_3"],
                "entity_name"                : "PointCondition2D1N",
                "entity_type"                : "condition",
                "properties_id"              : 0
            }
        }""")

        self.process = self.__CreateProcess(settings)
        self.assertTrue(self.root_model_part.HasSubModelPart("smp_with_conditions"))
        self.assertEqual(self.root_model_part.NumberOfNodes(), self.root_model_part.NumberOfConditions())

        self.__CheckCreatedEntitiesIdAreCorrectlyNumbered(self.root_model_part.Conditions)

    # def test_create_entities_with_preexisting_entites(self):
    #     pass

    def __CreateProcess(self, settings):
        return create_point_based_entites_process.Factory(settings, self.model)

    def __CheckCreatedEntitiesIdAreCorrectlyNumbered(self, entities):
        num_local_entites = len(entities)
        local_ids = [entity.Id for entity in entities]

        for entity in entities:
            local_first_id = entity.Id
            break

        scan_sum_num_entities = KM.DataCommunicator.GetDefault().ScanSum(num_local_entites)

        for exp_id, entity in zip(range(scan_sum_num_entities-num_local_entites, scan_sum_num_entities), entities):
            self.assertEqual(exp_id+1, entity.Id)

if __name__ == '__main__':
    KratosUnittest.main()
