import os

import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as kratos_utils

if KratosMultiphysics.IsDistributedRun():
    from KratosMultiphysics.mpi import distributed_import_model_part_utility

def GetFilePath(fileName):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), fileName)

class TestMetisSubModelPartList(KratosUnittest.TestCase):
    def setUp(self):
        self.comm = KratosMultiphysics.DataCommunicator.GetDefault()
        self.size = self.comm.Size()
        self.rank = self.comm.Rank()

    def tearDown(self):
        with KratosUnittest.WorkFolderScope(self.work_folder, __file__):
            self.comm.Barrier()
            if self.rank == 0:
                kratos_utils.DeleteFileIfExisting(self.file_name + ".time")
            kratos_utils.DeleteFileIfExisting(self.file_name + "_" + str(self.rank) + ".mdpa")
            self.comm.Barrier()

    def ReadModelPart(self, model_part, settings):
        if KratosMultiphysics.IsDistributedRun():
            model_part.AddNodalSolutionStepVariable(KratosMultiphysics.PARTITION_INDEX)
            model_part_import_util = distributed_import_model_part_utility.DistributedImportModelPartUtility(
                model_part, settings)
            model_part_import_util.ImportModelPart()
            model_part_import_util.CreateCommunicators()
        else:
            import_flags = KratosMultiphysics.ModelPartIO.READ | KratosMultiphysics.ModelPartIO.SKIP_TIMER
            mdpa_file_name = settings["model_import_settings"]["input_filename"].GetString()
            KratosMultiphysics.ModelPartIO(mdpa_file_name, import_flags).ReadModelPart(model_part)

    def test_ReadTwoSubModelParts(self):
        self.work_folder = ""
        self.file_name = "cube"
        current_model = KratosMultiphysics.Model()
        model_part = current_model.CreateModelPart("Main")
        settings = KratosMultiphysics.Parameters("""{
            "model_import_settings" : {
                "input_type": "mdpa",
                "input_filename"                             : \"""" + GetFilePath(self.file_name) + """\",
                "partition_in_memory"                        : false,
                "partition_sub_model_parts_list"             : true,
                "submodelpart_list" : ["submodelpart_solid", "submodelpart_liquid"]
            },
            "echo_level" : 0
        }""")
        results = {"Main.submodelpart_liquid" : [289, 924, 464],
                   "Main.submodelpart_solid" :  [563, 1652, 1118]}
        self.ReadModelPart(model_part, settings)
        for submodel_part_name in results:
            submodel_part = current_model[submodel_part_name]
            local_number_nodes = submodel_part.GetCommunicator().LocalMesh().NumberOfNodes()
            local_number_elements = submodel_part.GetCommunicator().LocalMesh().NumberOfElements()
            local_number_conditions = submodel_part.GetCommunicator().LocalMesh().NumberOfConditions()
            self.assertTrue(local_number_nodes > 0)
            self.assertTrue(local_number_elements > 0)
            self.assertTrue(local_number_conditions > 0)
            total_nodes = submodel_part.GetCommunicator().GetDataCommunicator().SumAll(local_number_nodes)
            total_elements =submodel_part.GetCommunicator().GetDataCommunicator().SumAll(local_number_elements)
            total_conditions = submodel_part.GetCommunicator().GetDataCommunicator().SumAll(local_number_conditions)
            self.assertEqual(total_nodes, results.get(submodel_part_name)[0])
            self.assertEqual(total_elements, results.get(submodel_part_name)[1])
            self.assertEqual(total_conditions, results.get(submodel_part_name)[2])

        local_main_number_nodes = model_part.GetCommunicator().LocalMesh().NumberOfNodes()
        local_main_number_elements = model_part.GetCommunicator().LocalMesh().NumberOfElements()
        local_main_number_conditions = model_part.GetCommunicator().LocalMesh().NumberOfConditions()

        total_main_nodes = model_part.GetCommunicator().GetDataCommunicator().SumAll(local_main_number_nodes)
        total_main_elements =model_part.GetCommunicator().GetDataCommunicator().SumAll(local_main_number_elements)
        total_main_conditions = model_part.GetCommunicator().GetDataCommunicator().SumAll(local_main_number_conditions)
        self.assertEqual(total_main_nodes, 852 )
        self.assertEqual(total_main_elements, 2576 )
        self.assertEqual(total_main_conditions, 1582 )


    def test_ReadSubSubModelParts(self):
        self.work_folder = ""
        self.file_name = "cube"
        current_model = KratosMultiphysics.Model()
        model_part = current_model.CreateModelPart("Main")
        settings = KratosMultiphysics.Parameters("""{
            "model_import_settings" : {
                "input_type": "mdpa",
                "input_filename"                             : \"""" + GetFilePath(self.file_name) + """\",
                "partition_in_memory"                        : false,
                "partition_sub_model_parts_list"             : true,
                "submodelpart_list" : ["ingate", "mainPart", "submodelpart_solid"]
            },
            "echo_level" : 0
        }""")
        results = {"Main.submodelpart_liquid.ingate" : [183, 537, 268],
                   "Main.submodelpart_liquid.mainPart" : [137, 387, 196],
                   "Main.submodelpart_solid" : [563,1652,1118]}
        self.ReadModelPart(model_part, settings)
        for submodel_part_name in results:
            submodel_part = current_model[submodel_part_name]
            local_number_nodes = submodel_part.GetCommunicator().LocalMesh().NumberOfNodes()
            local_number_elements = submodel_part.GetCommunicator().LocalMesh().NumberOfElements()
            local_number_conditions = submodel_part.GetCommunicator().LocalMesh().NumberOfConditions()
            self.assertTrue(local_number_nodes > 0)
            self.assertTrue(local_number_elements > 0)
            self.assertTrue(local_number_conditions > 0)
            total_nodes = submodel_part.GetCommunicator().GetDataCommunicator().SumAll(local_number_nodes)
            total_elements =submodel_part.GetCommunicator().GetDataCommunicator().SumAll(local_number_elements)
            total_conditions = submodel_part.GetCommunicator().GetDataCommunicator().SumAll(local_number_conditions)
            self.assertEqual(total_nodes, results.get(submodel_part_name)[0])
            self.assertEqual(total_elements, results.get(submodel_part_name)[1])
            self.assertEqual(total_conditions, results.get(submodel_part_name)[2])

        local_main_number_nodes = model_part.GetCommunicator().LocalMesh().NumberOfNodes()
        local_main_number_elements = model_part.GetCommunicator().LocalMesh().NumberOfElements()
        local_main_number_conditions = model_part.GetCommunicator().LocalMesh().NumberOfConditions()

        total_main_nodes = model_part.GetCommunicator().GetDataCommunicator().SumAll(local_main_number_nodes)
        total_main_elements =model_part.GetCommunicator().GetDataCommunicator().SumAll(local_main_number_elements)
        total_main_conditions = model_part.GetCommunicator().GetDataCommunicator().SumAll(local_main_number_conditions)
        self.assertEqual(total_main_nodes, 852 )
        self.assertEqual(total_main_elements, 2576 )
        self.assertEqual(total_main_conditions, 1582 )




if __name__ == '__main__':
    KratosUnittest.main()