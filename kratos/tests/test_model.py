from __future__ import print_function, absolute_import, division

import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest

import sys

class TestModel(KratosUnittest.TestCase):

    def test_model(self):
        current_model = KratosMultiphysics.Model()

        model_part = current_model.CreateModelPart("Main")
        model_part.CreateSubModelPart("Inlets")
        model_part.CreateSubModelPart("Temp")
        outlet = model_part.CreateSubModelPart("Outlet")

        aaa = current_model["Main.Outlet"].CreateSubModelPart("aaa")

        if (sys.version_info < (3, 2)):
            self.assertRaisesRegex = self.assertRaisesRegexp

        self.assertEqual(aaa, current_model["aaa"]) #search by flat name - should be eventually deprecated

        #check that a meaningful error is thrown
        with self.assertRaisesRegex(RuntimeError, "Error: The ModelPart named : \"abc\" was not found either as root-ModelPart or as a flat name. The total input string was \"abc\""):
            current_model["abc"]

        #check that a meaningful error is thrown
        with self.assertRaisesRegex(RuntimeError, "The ModelPart named : \"aaa\" was not found as SubModelPart of : \"Inlets\". The total input string was \"Main.Inlets.aaa\""):
            current_model["Main.Inlets.aaa"]

        #here i create a model part in the lowest level and i check that the other model parts have it
        current_model["Main.Outlet.aaa"].CreateNewNode(1,0.0,0.0,0.0)
        self.assertEqual(len(model_part.Nodes), 1)
        self.assertEqual(len(outlet.Nodes), 1)
        self.assertEqual(len(aaa.Nodes), 1)

        #self.assertEqual(current_model["Main"], model_part )
        #self.assertEqual(current_model["Main.Outlet"].Info(), outlet.Info() )
        #self.assertEqual(current_model["Main.Outlet.aaa"].Info(), aaa.Info() )

        self.assertTrue(current_model.HasModelPart("Main"))
        self.assertTrue(current_model.HasModelPart("Main.Outlet"))
        self.assertFalse(current_model.HasModelPart("Outlet"))


if __name__ == '__main__':
    KratosUnittest.main()
