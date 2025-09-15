import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest

class TestIntegrationPoints(KratosUnittest.TestCase):

    def test_calculate_on_integration_points(self):
        current_model = KratosMultiphysics.Model()

        mp = current_model.CreateModelPart("Main")

        # Create nodes
        mp.CreateNewNode(1,0.00,3.00,0.00)
        mp.CreateNewNode(2,1.00,2.25,0.00)
        mp.CreateNewNode(3,0.75,1.00,0.00)
        mp.CreateNewNode(4,2.25,2.00,0.00)
        mp.CreateNewNode(5,0.00,0.00,0.00)
        mp.CreateNewNode(6,3.00,3.00,0.00)
        mp.CreateNewNode(7,2.00,0.75,0.00)
        mp.CreateNewNode(8,3.00,0.00,0.00)

        # Create elements
        mp.CreateNewElement("Element2D4N", 1, [8,7,3,5], mp.GetProperties()[1])
        mp.CreateNewElement("Element2D4N", 2, [6,4,7,8], mp.GetProperties()[1])
        mp.CreateNewElement("Element2D4N", 3, [1,2,4,6], mp.GetProperties()[1])
        mp.CreateNewElement("Element2D4N", 4, [4,2,3,7], mp.GetProperties()[1])
        mp.CreateNewElement("Element2D4N", 5, [2,1,5,3], mp.GetProperties()[1])

        for elem in mp.Elements:
            # Boolean
            out = elem.CalculateOnIntegrationPoints(KratosMultiphysics.IS_RESTARTED, mp.ProcessInfo)
            for value in out:
                self.assertFalse(value)

            # Int
            out = elem.CalculateOnIntegrationPoints(KratosMultiphysics.LOAD_RESTART, mp.ProcessInfo)
            for value in out:
                self.assertEqual(value,0)
            # Double
            out = elem.CalculateOnIntegrationPoints(KratosMultiphysics.PRESSURE, mp.ProcessInfo)
            for value in out:
                self.assertAlmostEqual(value,0.0)

            # Array
            out = elem.CalculateOnIntegrationPoints(KratosMultiphysics.DISPLACEMENT, mp.ProcessInfo)
            for value in out:
                self.assertAlmostEqual(value[0],0.0)
                self.assertAlmostEqual(value[1],0.0)
                self.assertAlmostEqual(value[2],0.0)


if __name__ == '__main__':
    KratosUnittest.main()
