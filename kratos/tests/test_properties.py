import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics as KM

class TestProperties(KratosUnittest.TestCase):

    def test_copy_properties(self):
        current_model = KM.Model()

        model_part= current_model.CreateModelPart("Main")

        model_part.CreateNewProperties(1)
        properties = model_part.GetProperties()[1]

        properties.SetValue(KM.YOUNG_MODULUS, 1.0)
        self.assertEqual(properties.GetValue(KM.YOUNG_MODULUS), 1.0)

        cloned_properties = KM.Properties(properties) #copy constructor
        self.assertEqual(cloned_properties.GetValue(KM.YOUNG_MODULUS), 1.0)

        cloned_properties.SetValue(KM.YOUNG_MODULUS, 10.0)
        self.assertEqual(properties.GetValue(KM.YOUNG_MODULUS), 1.0)
        self.assertEqual(cloned_properties.GetValue(KM.YOUNG_MODULUS), 10.0)

    def test_erase_properties(self):
        current_model = KM.Model()

        model_part= current_model.CreateModelPart("Main")

        model_part.CreateNewProperties(1)
        properties = model_part.GetProperties()[1]

        properties.SetValue(KM.YOUNG_MODULUS, 1.0)
        self.assertEqual(properties.Has(KM.YOUNG_MODULUS), True)
        properties.Erase(KM.YOUNG_MODULUS)
        self.assertEqual(properties.Has(KM.YOUNG_MODULUS), False)
        
    def test_properties_utilities(self):
        current_model = KM.Model()

        model_part= current_model.CreateModelPart("Main")

        properties = model_part.CreateNewProperties(1)
        properties.SetValue(KM.YOUNG_MODULUS, 1.0)
        self.assertEqual(properties.GetValue(KM.YOUNG_MODULUS), 1.0)

        properties_copy = model_part.CreateNewProperties(2)
        properties_copy[KM.TEMPERATURE] = 15.0
        properties_copy[KM.PRESSURE] = 25.0
        KM.PropertiesUtilities.CopyPropertiesValues(properties, properties_copy)
        self.assertEqual(properties_copy.GetValue(KM.YOUNG_MODULUS), 1.0)
        self.assertEqual(properties_copy.Has(KM.TEMPERATURE), False)
        self.assertEqual(properties_copy.Has(KM.PRESSURE), False)

        properties.SetValue(KM.DENSITY, 12.0)
        self.assertEqual(properties.GetValue(KM.DENSITY), 12.0)
        self.assertEqual(properties_copy.Has(KM.DENSITY), False)

    def test_has_table(self):
        current_model = KM.Model()
        model_part= current_model.CreateModelPart("Main")
        properties = model_part.CreateNewProperties(1)
        table = KM.PiecewiseLinearTable()
        table.AddRow(9.0, 1.0)
        table.AddRow(10.0, 1.0)
        properties.SetTable(KM.TEMPERATURE, KM.YOUNG_MODULUS, table)
        self.assertTrue(properties.HasTable(KM.TEMPERATURE, KM.YOUNG_MODULUS))
        self.assertFalse(properties.HasTable(KM.TEMPERATURE, KM.PRESSURE))

if __name__ == '__main__':
    KratosUnittest.main()
