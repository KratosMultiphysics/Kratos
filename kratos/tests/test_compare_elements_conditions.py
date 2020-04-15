from __future__ import print_function, absolute_import, division

import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics as KM

class TestCompareElementsAndConditionsUtility(KratosUnittest.TestCase):

    def test_compare_elements(self):
        current_model = KM.Model()

        model_part= current_model.CreateModelPart("Main")

        model_part.CreateNewNode(1, 0.00,0.00,0.00)
        model_part.CreateNewNode(2, 1.00,0.00,0.00)
        model_part.CreateNewNode(3, 1.00,1.00,0.00)
        model_part.AddProperties(KM.Properties(1))
        pelem = model_part.CreateNewElement("Element2D3N", 1, [1,2,3], model_part.GetProperties()[1])

        self.assertEqual(KM.CompareElementsAndConditionsUtility.GetRegisteredName(pelem), "Element2D3N")
        self.assertNotEqual(KM.CompareElementsAndConditionsUtility.GetRegisteredName(pelem), "NotElement2D3N")

    def test_compare_conditions(self):
        current_model = KM.Model()

        model_part= current_model.CreateModelPart("Main")

        model_part.CreateNewNode(1, 0.00,0.00,0.00)
        model_part.CreateNewNode(2, 1.00,0.00,0.00)
        model_part.CreateNewNode(3, 1.00,1.00,0.00)
        model_part.AddProperties(KM.Properties(1))
        pcond = model_part.CreateNewCondition("SurfaceCondition3D3N", 1, [1,2,3], model_part.GetProperties()[1])

        self.assertEqual(KM.CompareElementsAndConditionsUtility.GetRegisteredName(pcond), "SurfaceCondition3D3N")
        self.assertNotEqual(KM.CompareElementsAndConditionsUtility.GetRegisteredName(pcond), "NotSurfaceCondition3D3N")

if __name__ == '__main__':
    KratosUnittest.main()
