from __future__ import print_function, absolute_import, division
import KratosMultiphysics as KM

import KratosMultiphysics.KratosUnittest as KratosUnittest

class TestArray1DInterface(KratosUnittest.TestCase):

    def CreateModelPart(self):
        model_part = KM.ModelPart("TestModelPart")
        model_part.AddNodalSolutionStepVariable(KM.VORTICITY)
        model_part.CreateNewNode(1, 0.00,0.00,0.00)
        model_part.CreateNewNode(2, 1.00,0.00,0.00)
        model_part.CreateNewNode(3, 0.00,1.00,0.00)

        model_part.AddProperties(KM.Properties(1))
        model_part.CreateNewElement("Element2D3N", 1, [1,2,3], model_part.GetProperties()[1])
        model_part.CreateNewCondition("Condition3D3N", 1, [1,2,3], model_part.GetProperties()[1])
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

    def test_SetNodalArraySolutionStepValueFromPython_Array3(self):
        model_part = self.CreateModelPart()
        node = model_part.GetNode(1)

        u = KM.Array3([1,2,3])
        node.SetSolutionStepValue(KM.VORTICITY,u)
        node.SetSolutionStepValue(KM.VORTICITY,0,u)

    def test_SetNodalArraySolutionStepValueFromPython_Array3_Implicit(self):
        model_part = self.CreateModelPart()
        node = model_part.GetNode(1)

        node.SetSolutionStepValue(KM.VORTICITY,[1,2,3])
        node.SetSolutionStepValue(KM.VORTICITY,KM.Vector([1,2,3]))

        node.SetSolutionStepValue(KM.VORTICITY,0,[1,2,3])
        node.SetSolutionStepValue(KM.VORTICITY,0,KM.Vector([1,2,3]))

    def test_SetNodalArraySolutionStepValueFromPython_WrongInput(self):
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

    def test_SetElementalArrayValueFromPython_Array3(self):
        model_part = self.CreateModelPart()
        element = model_part.GetElement(1)

        self._ValidAssignments(element)

    def test_SetElementalArrayValueFromPython_WrongInput(self):
        model_part = self.CreateModelPart()
        element = model_part.GetElement(1)

        self._InvalidAssignments(element)

    def test_SetConditionalArrayValueFromPython_Array3(self):
        model_part = self.CreateModelPart()
        condition = model_part.GetCondition(1)

        self._ValidAssignments(condition)

    def test_SetConditionalArrayValueFromPython_WrongInput(self):
        model_part = self.CreateModelPart()
        condition = model_part.GetCondition(1)

        self._InvalidAssignments(condition)

    def test_SetPropertyArrayValueFromPython_Array3(self):
        model_part = self.CreateModelPart()
        prop = model_part.GetProperties(1,0)

        self._ValidAssignments(prop)

    def test_SetPropertyArrayValueFromPython_WrongInput(self):
        model_part = self.CreateModelPart()
        prop = model_part.GetProperties(1,0)

        self._InvalidAssignments(prop)

    def _ValidAssignments(self,tested_object):

        # Set from Array3
        u = KM.Array3([1,2,3])
        tested_object.SetValue(KM.VORTICITY, u)
        value = tested_object.GetValue(KM.VORTICITY)
        for i in range(3):
            self.assertEqual(u[i],value[i])

        # Set from Vector
        v = KM.Vector([4,5,6])
        tested_object.SetValue(KM.VORTICITY, v)
        value = tested_object.GetValue(KM.VORTICITY)
        for i in range(3):
            self.assertEqual(v[i],value[i])

        # Set from list
        l = [7,8,9]
        tested_object.SetValue(KM.VORTICITY, l)
        value = tested_object.GetValue(KM.VORTICITY)
        for i in range(3):
            self.assertEqual(l[i],value[i])

    def _InvalidAssignments(self, tested_object):

        # Forbidden implicit conversion from wrong-size Array
        with self.assertRaisesRegex(TypeError, ".*incompatible function arguments.*"):
            u = KM.Array6([1,2,3,4,5,6])
            tested_object.SetValue(KM.VORTICITY, u)

        # Forbidden implicit conversion from wrong-size Vector
        with self.assertRaisesRegex(ValueError, ".*bad argument*"):
            v = KM.Vector([1,2,3,4])
            tested_object.SetValue(KM.VORTICITY, v)

        # Forbidden implicit conversion from wrong-size list
        with self.assertRaisesRegex(ValueError, ".*bad argument*"):
            l = [1,2,3,4]
            tested_object.SetValue(KM.VORTICITY, l)


if __name__ == "__main__":
    KratosUnittest.main()