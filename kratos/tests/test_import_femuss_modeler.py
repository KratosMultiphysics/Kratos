# Import necessary modules from KratosMultiphysics
import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest

from KratosMultiphysics.modelers.import_femuss_modeler import ImportFemussModeler

class TestImportFemussModeler(KratosUnittest.TestCase):
    """
    Unit tests for the ImportFemussModeler class.
    """

    def setUp(self):
        """
        Sets up the test environment.
        """


    def test_ImportFemussModeler(self):
        """
        Tests the import of a FEMUSS model geometry.
        """

        # Create the test model part in which the geometry will be imported
        model = KratosMultiphysics.Model()
        model_part = model.CreateModelPart("FemussModelPart")
        model_part.ProcessInfo.SetValue(KratosMultiphysics.DOMAIN_SIZE, 2)

        # Create and execute the modeler to import the FEMUSS geometry
        settings = KratosMultiphysics.Parameters('''{
            "geo_input_filename"  : "auxiliar_files_for_python_unittest/femuss_files/cylinder01.geo.dat",
            "fix_input_filename"  : "auxiliar_files_for_python_unittest/femuss_files/cylinder01.fix",
            "model_part_name" : "FemussModelPart"
        }''')
        modeler = ImportFemussModeler(model, settings)
        modeler.SetupGeometryModel()
        modeler.PrepareGeometryModel()
        modeler.SetupModelPart()

        # Assert the model part content
        self.assertEqual(model_part.NumberOfNodes(), 36)
        self.assertEqual(model_part.NumberOfGeometries(), 53)
        self.assertEqual([model_part.GetGeometry(1)[i].Id for i in range(model_part.GetGeometry(1).PointsNumber())], [18, 20, 17])
        self.assertVectorAlmostEqual([model_part.GetNode(1).X, model_part.GetNode(1).Y, model_part.GetNode(1).Z], [1.600000e+01, 0.0, 0.0])
        self.assertVectorAlmostEqual([model_part.GetNode(2).X, model_part.GetNode(2).Y, model_part.GetNode(2).Z], [1.600000e+01, 3.984547e+00, 0.0])

        # Assert the fixity sub-model parts
        self.assertTrue(model_part.HasSubModelPart("MyBoundary"))
        self.assertTrue(model_part.HasSubModelPart("Fix_11_0_0"))
        self.assertTrue(model_part.HasSubModelPart("Fix_01_0_0"))
        self.assertTrue(model_part.HasSubModelPart("Fix_11_100_0"))
        self.assertEqual(model_part.GetSubModelPart("MyBoundary").NumberOfNodes(), 7)
        self.assertEqual(model_part.GetSubModelPart("Fix_11_0_0").NumberOfNodes(), 3)
        self.assertEqual(model_part.GetSubModelPart("Fix_01_0_0").NumberOfNodes(), 3)
        self.assertEqual(model_part.GetSubModelPart("Fix_11_100_0").NumberOfNodes(), 5)

if __name__ == '__main__':
    # Set the logger severity level and run the KratosUnittest
    KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)
    KratosUnittest.main()