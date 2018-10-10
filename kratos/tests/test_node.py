from __future__ import print_function, absolute_import, division
import KratosMultiphysics as KM

import KratosMultiphysics.KratosUnittest as KratosUnittest

class TestNode(KratosUnittest.TestCase):

    def CreateModelPart(self):
        model_part = KM.ModelPart("TestModelPart")
        model_part.AddNodalSolutionStepVariable(KM.VORTICITY)
        model_part.CreateNewNode(1, 1.00,0.00,0.00)
        return model_part

    def test_SetArrayValueFromPython_Array3(self):
        model_part = self.CreateModelPart()
        node = model_part.GetNode(1)

        u = KM.Array3([1,2,3])
        node.SetValue(KM.VORTICITY,u)


    def test_SetArrayValueFromPython_Array3_Implicit(self):
        model_part = self.CreateModelPart()
        node = model_part.GetNode(1)

        node.SetValue(KM.VORTICITY,[1,2,3])
        node.SetValue(KM.VORTICITY,KM.Vector([1,2,3]))

    def test_SetArrayValueFromPython_WrongInput(self):
        model_part = self.CreateModelPart()
        node = model_part.GetNode(1)

        # Implicit conversion from wrong-size vector
        with self.assertRaisesRegex(TypeError, ".*incompatible function arguments.*"):
            node.SetValue(KM.VORTICITY,KM.Vector([1,2,3,4]))

        # Implicit conversion from wrong-size list
        with self.assertRaisesRegex(TypeError, ".*incompatible function arguments.*"):
            node.SetValue(KM.VORTICITY,[1,2,3,4])

    def test_SetArraySolutionStepValueFromPython_Array3(self):
        model_part = self.CreateModelPart()
        node = model_part.GetNode(1)

        u = KM.Array3([1,2,3])
        node.SetSolutionStepValue(KM.VORTICITY,u)
        node.SetSolutionStepValue(KM.VORTICITY,0,u)


    def test_SetArraySolutionStepValueFromPython_Array3_Implicit(self):
        model_part = self.CreateModelPart()
        node = model_part.GetNode(1)

        node.SetSolutionStepValue(KM.VORTICITY,[1,2,3])
        node.SetSolutionStepValue(KM.VORTICITY,KM.Vector([1,2,3]))

        node.SetSolutionStepValue(KM.VORTICITY,0,[1,2,3])
        node.SetSolutionStepValue(KM.VORTICITY,0,KM.Vector([1,2,3]))

    def test_SetArraySolutionStepValueFromPython_WrongInput(self):
        model_part = self.CreateModelPart()
        node = model_part.GetNode(1)

        # Implicit conversion from wrong-size vector
        with self.assertRaisesRegex(TypeError, ".*incompatible function arguments.*"):
            node.SetSolutionStepValue(KM.VORTICITY,KM.Vector([1,2,3,4]))

        # Implicit conversion from wrong-size list
        with self.assertRaisesRegex(TypeError, ".*incompatible function arguments.*"):
            node.SetSolutionStepValue(KM.VORTICITY,[1,2,3,4])

        # Implicit conversion from wrong-size vector
        with self.assertRaisesRegex(TypeError, ".*incompatible function arguments.*"):
            node.SetSolutionStepValue(KM.VORTICITY,0,KM.Vector([1,2,3,4]))

        # Implicit conversion from wrong-size list
        with self.assertRaisesRegex(TypeError, ".*incompatible function arguments.*"):
            node.SetSolutionStepValue(KM.VORTICITY,0,[1,2,3,4])

if __name__ == "__main__":
    KratosUnittest.main()