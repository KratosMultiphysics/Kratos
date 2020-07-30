from __future__ import print_function, absolute_import, division
import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as UnitTest
import KratosMultiphysics.ConvectionDiffusionApplication as ConvectionDiffusionApplication
from KratosMultiphysics.ConvectionDiffusionApplication import apply_thermal_face_process

class ApplyThermalFaceProcessTest(UnitTest.TestCase):
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
        interface_model_part = root_model_part.CreateSubModelPart("FaceModelPart")
        interface_model_part.AddCondition(root_model_part.GetCondition(1))
        interface_model_part.AddCondition(root_model_part.GetCondition(2))

        # Call the apply_thermal_interface_process
        interface_model_part = self.model.GetModelPart("MainModelPart.FaceModelPart")
        settings = KratosMultiphysics.Parameters(r'''{
            "model_part_name": "MainModelPart.FaceModelPart",
            "ambient_temperature": 300.0,
            "add_ambient_radiation": true,
            "emissivity": 0.1,
            "add_ambient_convection": true,
            "convection_coefficient": 0.0
        }''')
        apply_thermal_face_process.ApplyThermalFaceProcess(self.model, settings)

    def checkResults(self):
        # Check the interface properties
        face_model_part = self.model.GetModelPart("MainModelPart.FaceModelPart")
        face_properties = face_model_part.GetCondition(1).Properties
        self.assertEqual(face_model_part.NumberOfProperties(), 1)
        self.assertEqual(face_model_part.GetRootModelPart().NumberOfProperties(), 6)
        self.assertAlmostEqual(face_properties.GetValue(KratosMultiphysics.EMISSIVITY), 0.1, 1e-12)
        self.assertAlmostEqual(face_properties.GetValue(KratosMultiphysics.AMBIENT_TEMPERATURE), 300.0, 1e-12)
        self.assertAlmostEqual(face_properties.GetValue(KratosMultiphysics.CONVECTION_COEFFICIENT), 0.0, 1e-12)

    def testThermalFaceProcess(self):
        self.runTest()
        self.checkResults()

if __name__ == '__main__':
    test = ApplyThermalFaceProcessTest()
    test.runTest()
    test.checkResults()
