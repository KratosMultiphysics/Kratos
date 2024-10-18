import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as KratosUtilities
from KratosMultiphysics.testing.utilities import ReadModelPart

# STL imports
import pathlib


def GetFilePath(file_name):
    return pathlib.Path(__file__).absolute().parent / file_name

class TestCombineModelPartModeler(KratosUnittest.TestCase):

    def setUp(self):
        self.comm = KratosMultiphysics.Testing.GetDefaultDataCommunicator()
        self.file_name = "combine_model_part_modeler"
        self.work_folder = "test_files/combine_model_part_modeler"
        self.size = self.comm.Size()
        self.rank = self.comm.Rank()

    def tearDown(self):
        with KratosUnittest.WorkFolderScope(self.work_folder, __file__):
            self.comm.Barrier()
            if self.rank == 0:
                KratosUtilities.DeleteFileIfExisting(self.file_name + ".time")
                KratosUtilities.DeleteDirectoryIfExisting(self.file_name + "_partitioned")
            self.comm.Barrier()

    def test_CheckModelPartsHaveSameRoot(self):
        self.model = KratosMultiphysics.Model()
        self.model.CreateModelPart("main_model_part.submodel_1")
        self.model.CreateModelPart("not_main_model_part.submodel_2")
        settings = KratosMultiphysics.Parameters('''{
            "combined_model_part_name" : "thermal_model_part",
            "model_part_list" : [
                {
                    "origin_model_part": "main_model_part.submodel_1",
                    "destination_model_part": "thermal_model_part.submodel_1"
                },
                {
                    "origin_model_part": "not_main_model_part.submodel_2",
                    "destination_model_part": "thermal_model_part.submodel_2"
                }
            ]
        }''')
        modeler = KratosMultiphysics.CombineModelPartModeler(self.model, settings)
        expected_error_msg = "The origin model part \"not_main_model_part.submodel_2\" does not have"
        expected_error_msg += " the same root as the rest of origin model parts: \"main_model_part\"."
        with self.assertRaisesRegex(RuntimeError, expected_error_msg):
            modeler.SetupModelPart()

    def test_CombineModelParts(self):
        model_part = self._CreateModelPart()
        material_settings = KratosMultiphysics.Parameters('''{
            "Parameters" : {
                "materials_filename" : "test_files/combine_model_part_modeler/material.json"
            }
        }''')
        full_path = GetFilePath(material_settings["Parameters"]["materials_filename"].GetString())
        material_settings["Parameters"]["materials_filename"].SetString(str(full_path))
        KratosMultiphysics.ReadMaterialsUtility(material_settings,model_part.GetModel())
        settings = KratosMultiphysics.Parameters('''{
            "combined_model_part_name" : "thermal_model_part",
            "condition_name" : "SurfaceCondition3D3N",
            "model_part_list" : [
                {
                    "origin_model_part": "main_model_part.submodelpart_solid.chiller",
                    "destination_model_part": "thermal_model_part.submodelpart_solid.chiller"
                },
                {
                    "origin_model_part": "main_model_part.submodelpart_solid.contact_solid2liquid.contact_chiller_62-2liquid",
                    "destination_model_part": "thermal_model_part.submodelpart_solid.contact_solid2liquid.contact_chiller_62-2liquid"
                },
                {
                    "origin_model_part": "main_model_part.submodelpart_solid.contact_solid2air.contact_chiller_62-2air",
                    "destination_model_part": "thermal_model_part.submodelpart_solid.contact_solid2air.contact_chiller_62-2air"
                },
                {
                    "origin_model_part": "main_model_part.submodelpart_liquid.contact_liquid2air",
                    "destination_model_part": "thermal_model_part.submodelpart_liquid.contact_liquid2air"
                },
                {
                    "origin_model_part": "main_model_part.submodelpart_liquid.contact_liquid2solid.contact_liquid2chiller_62",
                    "destination_model_part": "thermal_model_part.submodelpart_liquid.contact_liquid2solid.contact_liquid2chiller_62"
                },
                {
                    "origin_model_part": "main_model_part.submodelpart_liquid.mainPart.mainPart_57",
                    "destination_model_part": "thermal_model_part.submodelpart_liquid.mainPart.mainPart_57"
                },
                {
                    "origin_model_part": "main_model_part.submodelpart_liquid.contact_liquid2solid.contact_liquid2mold_158",
                    "destination_model_part": "thermal_model_part.submodelpart_liquid.contact_liquid2air.contact_liquid2mold_158"
                }
            ]
        }''')
        modeler = KratosMultiphysics.CombineModelPartModeler(self.model, settings)
        modeler.SetupModelPart()

        self.assertEqual(self.model["thermal_model_part"].GetCommunicator().GlobalNumberOfNodes(), 79)
        self.assertEqual(self.model["thermal_model_part"].GetCommunicator().GlobalNumberOfElements(), 160)
        self.assertEqual(self.model["thermal_model_part"].GetCommunicator().GlobalNumberOfConditions(), 108)
        self.assertEqual(len(self.model["thermal_model_part"].Properties), 3)
        self.assertEqual(self.model["thermal_model_part"].GetBufferSize(), 3)

        #Solid
        self.assertEqual(self.model["thermal_model_part.submodelpart_solid.chiller.chiller_62"].GetCommunicator().GlobalNumberOfNodes(), 27)
        self.assertEqual(self.model["thermal_model_part.submodelpart_solid.chiller.chiller_62"].GetCommunicator().GlobalNumberOfElements(), 46)
        self.assertEqual(self.model["thermal_model_part.submodelpart_solid.chiller.chiller_62"].GetCommunicator().GlobalNumberOfConditions(), 10)
        self.assertEqual(self.model["thermal_model_part.submodelpart_solid.chiller"].GetCommunicator().GlobalNumberOfNodes(), 27)
        self.assertEqual(self.model["thermal_model_part.submodelpart_solid.chiller"].GetCommunicator().GlobalNumberOfElements(), 46)
        self.assertEqual(self.model["thermal_model_part.submodelpart_solid.chiller"].GetCommunicator().GlobalNumberOfConditions(), 10)
        self.assertEqual(self.model["thermal_model_part.submodelpart_solid"].GetCommunicator().GlobalNumberOfNodes(), 27)
        self.assertEqual(self.model["thermal_model_part.submodelpart_solid"].GetCommunicator().GlobalNumberOfElements(), 46)
        self.assertEqual(self.model["thermal_model_part.submodelpart_solid"].GetCommunicator().GlobalNumberOfConditions(), 10)

        #Liquid
        self.assertEqual(self.model["thermal_model_part.submodelpart_liquid.mainPart.mainPart_57"].GetCommunicator().GlobalNumberOfNodes(), 52)
        self.assertEqual(self.model["thermal_model_part.submodelpart_liquid.mainPart.mainPart_57"].GetCommunicator().GlobalNumberOfElements(), 114)
        self.assertEqual(self.model["thermal_model_part.submodelpart_liquid.mainPart.mainPart_57"].GetCommunicator().GlobalNumberOfConditions(), 98)
        self.assertEqual(self.model["thermal_model_part.submodelpart_liquid.mainPart"].GetCommunicator().GlobalNumberOfNodes(), 52)
        self.assertEqual(self.model["thermal_model_part.submodelpart_liquid.mainPart"].GetCommunicator().GlobalNumberOfElements(), 114)
        self.assertEqual(self.model["thermal_model_part.submodelpart_liquid.mainPart"].GetCommunicator().GlobalNumberOfConditions(), 98)
        self.assertEqual(self.model["thermal_model_part.submodelpart_liquid"].GetCommunicator().GlobalNumberOfNodes(), 52)
        self.assertEqual(self.model["thermal_model_part.submodelpart_liquid"].GetCommunicator().GlobalNumberOfElements(), 114)
        self.assertEqual(self.model["thermal_model_part.submodelpart_liquid"].GetCommunicator().GlobalNumberOfConditions(), 98)

        #Liquid2Solid
        self.assertEqual(self.model["thermal_model_part.submodelpart_liquid.contact_liquid2solid.contact_liquid2chiller_62"].GetCommunicator().GlobalNumberOfNodes(), 11)
        self.assertEqual(self.model["thermal_model_part.submodelpart_liquid.contact_liquid2solid.contact_liquid2chiller_62"].GetCommunicator().GlobalNumberOfElements(), 0)
        self.assertEqual(self.model["thermal_model_part.submodelpart_liquid.contact_liquid2solid.contact_liquid2chiller_62"].GetCommunicator().GlobalNumberOfConditions(), 10)
        self.assertEqual(self.model["thermal_model_part.submodelpart_liquid.contact_liquid2solid"].GetCommunicator().GlobalNumberOfNodes(), 11)
        self.assertEqual(self.model["thermal_model_part.submodelpart_liquid.contact_liquid2solid"].GetCommunicator().GlobalNumberOfElements(), 0)
        self.assertEqual(self.model["thermal_model_part.submodelpart_liquid.contact_liquid2solid"].GetCommunicator().GlobalNumberOfConditions(), 10)

        #Liquid2Air
        self.assertEqual(self.model["thermal_model_part.submodelpart_liquid.contact_liquid2air"].GetCommunicator().GlobalNumberOfNodes(), 50)
        self.assertEqual(self.model["thermal_model_part.submodelpart_liquid.contact_liquid2air"].GetCommunicator().GlobalNumberOfElements(), 0)
        self.assertEqual(self.model["thermal_model_part.submodelpart_liquid.contact_liquid2air"].GetCommunicator().GlobalNumberOfConditions(), 88)
        self.assertEqual(self.model["thermal_model_part.submodelpart_liquid.contact_liquid2air.contact_liquid2mold_158"].GetCommunicator().GlobalNumberOfNodes(), 50)
        self.assertEqual(self.model["thermal_model_part.submodelpart_liquid.contact_liquid2air.contact_liquid2mold_158"].GetCommunicator().GlobalNumberOfElements(), 0)
        self.assertEqual(self.model["thermal_model_part.submodelpart_liquid.contact_liquid2air.contact_liquid2mold_158"].GetCommunicator().GlobalNumberOfConditions(), 88)

        #Solid2Liquid
        self.assertEqual(self.model["thermal_model_part.submodelpart_solid.contact_solid2liquid.contact_chiller_62-2liquid"].GetCommunicator().GlobalNumberOfNodes(), 11)
        self.assertEqual(self.model["thermal_model_part.submodelpart_solid.contact_solid2liquid.contact_chiller_62-2liquid"].GetCommunicator().GlobalNumberOfElements(), 0)
        self.assertEqual(self.model["thermal_model_part.submodelpart_solid.contact_solid2liquid.contact_chiller_62-2liquid"].GetCommunicator().GlobalNumberOfConditions(), 10)
        self.assertEqual(self.model["thermal_model_part.submodelpart_solid.contact_solid2liquid"].GetCommunicator().GlobalNumberOfNodes(), 11)
        self.assertEqual(self.model["thermal_model_part.submodelpart_solid.contact_solid2liquid"].GetCommunicator().GlobalNumberOfElements(), 0)
        self.assertEqual(self.model["thermal_model_part.submodelpart_solid.contact_solid2liquid"].GetCommunicator().GlobalNumberOfConditions(), 10)

        #Solid2Air
        self.assertEqual(self.model["thermal_model_part.submodelpart_solid.contact_solid2air.contact_chiller_62-2air"].GetCommunicator().GlobalNumberOfNodes(), 0)
        self.assertEqual(self.model["thermal_model_part.submodelpart_solid.contact_solid2air.contact_chiller_62-2air"].GetCommunicator().GlobalNumberOfElements(), 0)
        self.assertEqual(self.model["thermal_model_part.submodelpart_solid.contact_solid2air.contact_chiller_62-2air"].GetCommunicator().GlobalNumberOfConditions(), 0)
        self.assertEqual(self.model["thermal_model_part.submodelpart_solid.contact_solid2air"].GetCommunicator().GlobalNumberOfNodes(), 0)
        self.assertEqual(self.model["thermal_model_part.submodelpart_solid.contact_solid2air"].GetCommunicator().GlobalNumberOfElements(), 0)
        self.assertEqual(self.model["thermal_model_part.submodelpart_solid.contact_solid2air"].GetCommunicator().GlobalNumberOfConditions(), 0)

        for elem in self.model["thermal_model_part.submodelpart_liquid"].Elements:
            self.assertEqual(elem.Properties.Id,0)

        for elem in self.model["thermal_model_part.submodelpart_solid"].Elements:
            self.assertEqual(elem.Properties.Id,1)

    def _CreateModelPart(self):
        # Create model part
        self.model = KratosMultiphysics.Model()
        model_part = self.model.CreateModelPart("main_model_part")

        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.TEMPERATURE)

        # Define buffer size
        model_part.SetBufferSize(3)

        # Read mdpa -- Hexa box with X=[-0.25,0.25], Y=[-0.25,0.25], Z=[-0.5,0.5]
        with KratosUnittest.WorkFolderScope(self.work_folder, __file__):
            ReadModelPart(self.file_name, model_part)


        return model_part

if __name__ == '__main__':
    KratosUnittest.main()