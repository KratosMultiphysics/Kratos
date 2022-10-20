# Importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics.modelers.delete_model_parts_modeler import DeleteModelPartsModeler

class TestDeleteModelPartModeler(KratosUnittest.TestCase):
    def testDeleteModelPartModelerOneModelPart(self):
        # Create a fake model part hierarchy to operate with
        model = KratosMultiphysics.Model()
        model.CreateModelPart("ModelPart1")
        model.CreateModelPart("ModelPart2")

        # Set up the import model part modeler
        settings = KratosMultiphysics.Parameters('''{
            "model_part_names" : ["ModelPart2"]
        }''')
        delete_model_parts_modeler = DeleteModelPartsModeler(model, settings)

        # Call the modeler methods
        delete_model_parts_modeler.SetupGeometryModel()
        delete_model_parts_modeler.PrepareGeometryModel()
        delete_model_parts_modeler.SetupModelPart()

        # Check results
        self.assertTrue(model.HasModelPart("ModelPart1"))
        self.assertFalse(model.HasModelPart("ModelPart2"))

    def testDeleteModelPartModelerMultipleModelParts(self):
        # Create a fake model part hierarchy to operate with
        model = KratosMultiphysics.Model()
        model_part_1 = model.CreateModelPart("ModelPart1")
        model_part_1_1 = model_part_1.CreateSubModelPart("SubModelPart11")
        model_part_1_1_2 = model_part_1_1.CreateSubModelPart("SubModelPart112")
        model_part_1_1_2.CreateSubModelPart("SubModelPart1121")
        model_part_1_1_2.CreateSubModelPart("SubModelPart1122")

        model_part_2 = model.CreateModelPart("ModelPart2")
        model_part_2.CreateSubModelPart("SubModelPart21")
        model_part_2_2 = model_part_2.CreateSubModelPart("SubModelPart22")
        model_part_2_2_1 = model_part_2_2.CreateSubModelPart("SubModelPart221")
        model_part_2_2_1.CreateSubModelPart("SubModelPart2211")

        # Set up the import model part modeler
        settings = KratosMultiphysics.Parameters('''{
            "echo_level" : 2,
            "model_part_names" : [
                "ModelPart1.SubModelPart11",
                "ModelPart1.SubModelPart11.SubModelPart112",
                "ModelPart2.SubModelPart22.SubModelPart221.SubModelPart2211",
                "ModelPart2.SubModelPart21",
                "ModelPart1"]
        }''')
        delete_model_parts_modeler = DeleteModelPartsModeler(model, settings)

        # Call the modeler methods
        delete_model_parts_modeler.SetupGeometryModel()
        delete_model_parts_modeler.PrepareGeometryModel()
        delete_model_parts_modeler.SetupModelPart()

        # Check results
        self.assertTrue(model.HasModelPart("ModelPart2"))
        self.assertTrue(model.HasModelPart("ModelPart2.SubModelPart22"))
        self.assertTrue(model.HasModelPart("ModelPart2.SubModelPart22.SubModelPart221"))
        self.assertTrue(model.HasModelPart("ModelPart2.SubModelPart22.SubModelPart221.SubModelPart2211"))
        self.assertFalse(model.HasModelPart("ModelPart1"))

if __name__ == '__main__':
    KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)
    KratosUnittest.main()