import os

import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as kratos_utils
from KratosMultiphysics.testing.utilities import ReadModelPart

def GetFilePath(fileName):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), fileName)

class TestMetisSubModelPartList(KratosUnittest.TestCase):
    def setUp(self):
        self.comm = KratosMultiphysics.Testing.GetDefaultDataCommunicator()
        self.size = self.comm.Size()
        self.rank = self.comm.Rank()

    def tearDown(self):
        with KratosUnittest.WorkFolderScope(self.work_folder, __file__):
            self.comm.Barrier()
            kratos_utils.DeleteDirectoryIfExisting("cube_partitioned")
            self.comm.Barrier()

    def test_ReadTwoSubModelParts(self):
        """Checks that all processor have entities from the given list
           of sub model parts.
        """
        self.work_folder = ""
        self.file_name = "cube"
        current_model = KratosMultiphysics.Model()
        model_part = current_model.CreateModelPart("Main")
        model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] = 3
        settings = KratosMultiphysics.Parameters("""{
            "model_import_settings" : {
                "input_type": "mdpa",
                "input_filename"                             : \"""" + GetFilePath(self.file_name) + """\",
                "partition_in_memory"                        : false,
                "sub_model_part_list" : ["submodelpart_solid", "submodelpart_liquid"]
            },
            "echo_level" : 0
        }""")
        results = [["Main.submodelpart_liquid" , [133, 381, 228]],
                   ["Main.submodelpart_solid" ,  [280, 810, 552]]]
        ReadModelPart(self.file_name, model_part, settings)
        for i_result in results:
            submodel_part_name = i_result[0]
            submodel_part = current_model[submodel_part_name]
            local_number_nodes = submodel_part.GetCommunicator().LocalMesh().NumberOfNodes()
            local_number_elements = submodel_part.GetCommunicator().LocalMesh().NumberOfElements()
            local_number_conditions = submodel_part.GetCommunicator().LocalMesh().NumberOfConditions()
            if self.size<=10: #if too many partitions are used, some may end up with no nodes/elems
                self.assertTrue(local_number_nodes > 0)
                self.assertTrue(local_number_elements > 0)
                self.assertTrue(local_number_conditions > 0)
            total_nodes = submodel_part.GetCommunicator().GlobalNumberOfNodes()
            total_elements = submodel_part.GetCommunicator().GlobalNumberOfElements()
            total_conditions = submodel_part.GetCommunicator().GlobalNumberOfConditions()
            self.assertEqual(total_nodes, i_result[1][0])
            self.assertEqual(total_elements, i_result[1][1])
            self.assertEqual(total_conditions, i_result[1][2])

        total_main_nodes = model_part.GetCommunicator().GlobalNumberOfNodes()
        total_main_elements = model_part.GetCommunicator().GlobalNumberOfElements()
        total_main_conditions = model_part.GetCommunicator().GlobalNumberOfConditions()

        self.assertEqual(total_main_nodes, 413 )
        self.assertEqual(total_main_elements, 1191 )
        self.assertEqual(total_main_conditions, 780 )


    def test_ReadSubSubModelParts(self):
        """Checks that all processor have entities from the given list
           of sub-sub model parts.
        """
        self.work_folder = ""
        self.file_name = "cube"
        current_model = KratosMultiphysics.Model()
        model_part = current_model.CreateModelPart("Main")
        model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] = 3
        settings = KratosMultiphysics.Parameters("""{
            "model_import_settings" : {
                "input_type": "mdpa",
                "input_filename"                             : \"""" + GetFilePath(self.file_name) + """\",
                "partition_in_memory"                        : false,
                "sub_model_part_list" : ["ingate", "mainPart", "submodelpart_solid"]
            },
            "echo_level" : 0
        }""")
        results = [["Main.submodelpart_liquid.ingate" , [81, 188, 110]],
                   ["Main.submodelpart_liquid.mainPart" , [85, 193, 118]],
                   ["Main.submodelpart_solid" , [280,810,552]]]
        ReadModelPart(self.file_name, model_part, settings)
        for i_result in results:
            submodel_part_name = i_result[0]
            submodel_part = current_model[submodel_part_name]
            local_number_nodes = submodel_part.GetCommunicator().LocalMesh().NumberOfNodes()
            local_number_elements = submodel_part.GetCommunicator().LocalMesh().NumberOfElements()
            local_number_conditions = submodel_part.GetCommunicator().LocalMesh().NumberOfConditions()
            if self.size<=10: #if too many partitions are used, some may end up with no nodes/elems
                self.assertTrue(local_number_nodes > 0)
                self.assertTrue(local_number_elements > 0)
                self.assertTrue(local_number_conditions > 0)
            total_nodes = submodel_part.GetCommunicator().GlobalNumberOfNodes()
            total_elements = submodel_part.GetCommunicator().GlobalNumberOfElements()
            total_conditions = submodel_part.GetCommunicator().GlobalNumberOfConditions()
            self.assertEqual(total_nodes, i_result[1][0])
            self.assertEqual(total_elements, i_result[1][1])
            self.assertEqual(total_conditions, i_result[1][2])

        total_main_nodes = model_part.GetCommunicator().GlobalNumberOfNodes()
        total_main_elements = model_part.GetCommunicator().GlobalNumberOfElements()
        total_main_conditions = model_part.GetCommunicator().GlobalNumberOfConditions()
        self.assertEqual(total_main_nodes, 413 )
        self.assertEqual(total_main_elements, 1191 )
        self.assertEqual(total_main_conditions, 780 )

    def test_ReadWithoutSubModelParts(self):
        """Checks that all processor have entities of main model part
           if sub_model_parts_list is empty.
        """
        self.work_folder = ""
        self.file_name = "cube"
        current_model = KratosMultiphysics.Model()
        model_part = current_model.CreateModelPart("Main")
        model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] = 3
        settings = KratosMultiphysics.Parameters("""{
            "model_import_settings" : {
                "input_type": "mdpa",
                "input_filename"                             : \"""" + GetFilePath(self.file_name) + """\",
                "partition_in_memory"                        : false
            },
            "echo_level" : 0
        }""")
        ReadModelPart(self.file_name, model_part, settings)

        local_main_number_nodes = model_part.GetCommunicator().LocalMesh().NumberOfNodes()
        local_main_number_elements = model_part.GetCommunicator().LocalMesh().NumberOfElements()
        local_main_number_conditions = model_part.GetCommunicator().LocalMesh().NumberOfConditions()

        self.assertTrue(local_main_number_nodes > 0)
        self.assertTrue(local_main_number_elements > 0)
        self.assertTrue(local_main_number_conditions > 0)

        total_main_nodes = model_part.GetCommunicator().GlobalNumberOfNodes()
        total_main_elements = model_part.GetCommunicator().GlobalNumberOfElements()
        total_main_conditions = model_part.GetCommunicator().GlobalNumberOfConditions()
        self.assertEqual(total_main_nodes, 413 )
        self.assertEqual(total_main_elements, 1191 )
        self.assertEqual(total_main_conditions, 780 )

    def test_NodesAreNotBeingReordered(self):
        """Checks that all processor have entities of main model part
           if sub_model_parts_list is empty.
        """
        self.work_folder = ""
        self.file_name = "cube"
        current_model = KratosMultiphysics.Model()
        model_part = current_model.CreateModelPart("Main")
        model_part.ProcessInfo[KratosMultiphysics.DOMAIN_SIZE] = 3
        settings = KratosMultiphysics.Parameters("""{
            "model_import_settings" : {
                "input_type": "mdpa",
                "input_filename"                             : \"""" + GetFilePath(self.file_name) + """\",
                "partition_in_memory"                        : false,
                "sub_model_part_list" : ["ingate", "mainPart", "submodelpart_solid"]
            },
            "echo_level" : 0
        }""")
        ReadModelPart(self.file_name, model_part, settings)
        results = {1:[0.0,	0.0508782,	0.0514686],
            20:[0.0,	0.0176281,	0.0138362],
            50:[-0.0248201,	0.025,	0.1],
            80:[-0.01241,	0.0375,	0.0],
            100:[0.0244074,	0.025,	0.0],
            120:[0.0244074,	0.0,	0.05],
            140:[0.0244074,	0.0725731,	0.0483441],
            170:[-0.0373201,	0.0619352,	0.0700256],
            200:[-0.0373201,	0.1125,	0.034375],
            220:[-0.000206333,	0.1125,	-0.0125],
            240:[0.0369074,	0.1125,	0.1125],
            260:[0.0369074,	-0.0125,	0.01875],
            280:[0.0369074,	0.0947413,	0.0296595],
            300:[-0.0248201,	0.025,	0.1],
            330:[-0.000103589,	0.0983335,	0.1125],
            360:[-0.0128035,	0.0,	0.00872692],
            393:[-0.0159504,	0.0707544,	-0.0125],
            413:[-0.01241,	0.1,	0.0375]}

        for node in model_part.Nodes:
            if node in results:
                self.assertAlmostEqual(results.get(node.Id)[0], node.X, 7)
                self.assertAlmostEqual(results.get(node.Id)[1], node.Y, 7)
                self.assertAlmostEqual(results.get(node.Id)[2], node.Z, 7)

if __name__ == '__main__':
    KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)
    KratosUnittest.main()
