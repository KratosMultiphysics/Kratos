from __future__ import print_function, absolute_import, division
import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as UnitTest
import KratosMultiphysics.ConvectionDiffusionApplication as ConvectionDiffusionApplication

import os

class ApplyThermalInterfaceProcessTest(UnitTest.TestCase):
    def runTest(self):
        # Create a model part containing some properties
        self.model = KratosMultiphysics.Model()
        root_model_part = self.model.CreateModelPart("MainModelPart")
        for i in range(5):
            sub_model_part = root_model_part.CreateSubModelPart("SubModelPart" + str(i))
            new_property_i = KratosMultiphysics.Properties(i)
            sub_model_part.AddProperties(new_property_i)

        # Create a fake interface condition
        root_model_part.CreateNewNode(1, 0.0, 0.0, 0.0)
        root_model_part.CreateNewNode(2, 1.0, 0.0, 0.0)
        root_model_part.CreateNewNode(3, 2.0, 0.0, 0.0)
        root_model_part.CreateNewCondition("ThermalFace2D2N", 1, [1,2], root_model_part.GetProperties()[1])
        root_model_part.CreateNewCondition("ThermalFace2D2N", 2, [2,3], root_model_part.GetProperties()[2])

        # Create a fake interface model part
        interface_model_part = root_model_part.CreateSubModelPart("InterfaceModelPart")
        interface_model_part.AddCondition(root_model_part.GetCondition(1))
        interface_model_part.AddCondition(root_model_part.GetCondition(2))

        # Call the apply_thermal_interface_process
        interface_model_part = self.model.GetModelPart("InterfaceModelPart")
        import apply_thermal_interface_process
        settings = KratosMultiphysics.Parameters(r'''{
            "model_part_name": "InterfaceModelPart"
        }''')
        thermal_int_proc = apply_thermal_interface_process.ApplyThermalInterfaceProcess(self.model, settings)

    def checkResults(self):
        # Check the interface properties
        interface_model_part = self.model.GetModelPart("InterfaceModelPart")
        self.assertEqual(interface_model_part.NumberOfProperties(), 1)
        self.assertEqual(interface_model_part.GetRootModelPart().NumberOfProperties(), 6)
        self.assertAlmostEqual(interface_model_part.GetCondition(1).Properties.GetValue(KratosMultiphysics.EMISSIVITY), 0.0, 1e-12)
        self.assertAlmostEqual(interface_model_part.GetCondition(1).Properties.GetValue(KratosMultiphysics.CONVECTION_COEFFICIENT), 0.0, 1e-12)

    def testThermalInterfaceProcess(self):
        self.runTest()
        self.checkResults()

if __name__ == '__main__':
    test = ApplyThermalInterfaceProcessTest()
    test.runTest()
    test.checkResults()
