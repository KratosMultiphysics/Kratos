import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics as KM

from unittest.mock import patch
from io import StringIO

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

    def test_accessor(self):
        current_model = KM.Model()
        model_part = current_model.CreateModelPart("Main")
        properties = model_part.CreateNewProperties(1)

        self.assertFalse(properties.HasAccessor(KM.TEMPERATURE))

        properties.SetAccessor(KM.TEMPERATURE, KM.Accessor.Create())

        self.assertTrue(properties.HasAccessor(KM.TEMPERATURE))
        accessor_from_get = properties.GetAccessor(KM.TEMPERATURE)

    def test_accessor_invalid(self):
        current_model = KM.Model()
        model_part= current_model.CreateModelPart("Main")

        properties = model_part.CreateNewProperties(1)

        with self.assertRaisesRegex(Exception, "Error: Trying to retrieve inexisting accessor for 'TEMPERATURE' in properties 1."):
            accessor_from_get = properties.GetAccessor(KM.TEMPERATURE)

    def test_accessor_get_consumed(self):
        current_model = KM.Model()
        model_part= current_model.CreateModelPart("Main")

        properties_1 = model_part.CreateNewProperties(1)
        properties_2 = model_part.CreateNewProperties(2)
        
        properties_1.SetAccessor(KM.TEMPERATURE, KM.Accessor.Create())
        accessor_from_get = properties_1.GetAccessor(KM.TEMPERATURE)
        properties_2.SetAccessor(KM.TEMPERATURE, accessor_from_get)

        with self.assertRaisesRegex(Exception, "Trying to get a consumed or invalid Accessor."):
            accessor_from_get = properties_1.GetAccessor(KM.TEMPERATURE)

    def test_accessor_assign_twice(self):
        current_model = KM.Model()
        model_part= current_model.CreateModelPart("Main")

        properties_1 = model_part.CreateNewProperties(1)
        properties_2 = model_part.CreateNewProperties(2)

        accessor_binding = KM.Accessor.Create()
        
        properties_1.SetAccessor(KM.TEMPERATURE, accessor_binding)

        with self.assertRaisesRegex(Exception, "Trying to set a consumed or invalid Accessor. Accessors are unique. Please create a different one."):
            properties_2.SetAccessor(KM.TEMPERATURE, accessor_binding)

    @patch('sys.stdout', new_callable=StringIO)
    def test_print(self, mock_stdout):
        current_model = KM.Model()
        model_part= current_model.CreateModelPart("Main")
        properties = model_part.CreateNewProperties(1)
        table = KM.PiecewiseLinearTable()
        table.AddRow(9.0, 1.0)
        table.AddRow(10.0, 1.0)
        properties.SetTable(KM.TEMPERATURE, KM.YOUNG_MODULUS, table)
        print(properties)
        output = mock_stdout.getvalue()
        self.assertTrue("Properties" in output)
        self.assertTrue("Id : 1" in output)
        self.assertTrue("This properties contains 1 tables" in output)
        self.assertTrue("\t9\t\t1" in output)
        self.assertTrue("\t10\t\t1" in output)

if __name__ == '__main__':
    KM.Logger.GetDefaultOutput().SetSeverity(KM.Logger.Severity.WARNING)
    KratosUnittest.main()
