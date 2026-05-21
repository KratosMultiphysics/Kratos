# Importing the Kratos Library
import KratosMultiphysics as KM
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.MedApplication as KratosMed
from KratosMultiphysics.MedApplication.modelers.import_med_modeler import ImportMedModeler
from pathlib import Path
from testing_utilities import MedModelPartIOTestCase, GetMedPath, get_num_geometries_by_type


class TestImportMedModeler(MedModelPartIOTestCase):
    def test_import_med_modeler(self):
        # Set up the import model part modeler
        model = KM.Model()
        settings = KM.Parameters(
            """{
            "echo_level" : 0,
            "input_filename" : "",
            "model_part_name" : "Main"
        }"""
        )
        settings["input_filename"].SetString(str(GetMedPath(Path("hexahedral_8N"))))
        import_mdpa_modeler = ImportMedModeler(model, settings)

        # Get the model part created by the modeler
        model_part = model.GetModelPart(settings["model_part_name"].GetString())

        # Call the modeler methods
        import_mdpa_modeler.SetupGeometryModel()
        import_mdpa_modeler.PrepareGeometryModel()
        import_mdpa_modeler.SetupModelPart()

        self._basic_checks(model_part)

        # Check read ModelPart
        self.assertEqual(model_part.NumberOfNodes(), 216)
        self.assertEqual(model_part.NumberOfGeometries(), 371)

        # check how many geoms of each type
        exp_geoms = {KM.Hexahedra3D8: 125, KM.Quadrilateral3D4: 150, KM.Line3D2: 60, KM.Geometry: 36}
        self.assertEqual(sum(exp_geoms.values()), model_part.NumberOfGeometries())
        self.assertDictEqual(exp_geoms, get_num_geometries_by_type(model_part))

        self.assertAlmostEqual(KratosMed.MedTestingUtilities.ComputeLength(model_part), 3200)
        self.assertAlmostEqual(KratosMed.MedTestingUtilities.ComputeArea(model_part), 340000)
        self.assertAlmostEqual(KratosMed.MedTestingUtilities.ComputeVolume(model_part), 10000000)
        self.assertAlmostEqual(KratosMed.MedTestingUtilities.ComputeDomainSize(model_part), 10343200)

        for node in model_part.Nodes:
            self.assertTrue(0.0 <= node.X <= 500.0)
            self.assertTrue(0.0 <= node.X0 <= 500.0)

            self.assertTrue(0.0 <= node.Y <= 100.0)
            self.assertTrue(0.0 <= node.Y0 <= 100.0)

            self.assertTrue(0.0 <= node.Z <= 200.0)
            self.assertTrue(0.0 <= node.Z0 <= 200.0)


if __name__ == "__main__":
    KM.Logger.GetDefaultOutput().SetSeverity(KM.Logger.Severity.WARNING)
    KratosUnittest.main()
