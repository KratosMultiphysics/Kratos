from __future__ import print_function, absolute_import, division

import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics as KM

class TestElement(KratosUnittest.TestCase):

    def test_compare_elements(self):
        current_model = KratosMultiphysics.Model()

        model_part= current_model.CreateModelPart("Main")

        model_part.CreateNewNode(1, 0.00,0.00,0.00)
        model_part.CreateNewNode(2, 1.00,0.00,0.00)
        model_part.CreateNewNode(3, 1.00,1.00,0.00)
        model_part.AddProperties(KratosMultiphysics.Properties(1))
        model_part.CreateNewElement("Element3D", 1, [1,2,3], model_part.GetProperties()[1])

        for elem in model_part.Elements:
            geom = elem.GetGeometry()
            self.assertEqual(geom[0].Id, 1)
            self.assertEqual(geom[1].Id, 2)
            self.assertEqual(geom[2].Id, 3)
            for node in geom:


        self.assertEqual(model_part.NumberOfConditions(), 1)
        self.assertEqual(model_part.NumberOfConditions(0), 1)

if __name__ == '__main__':
    KratosUnittest.main()
