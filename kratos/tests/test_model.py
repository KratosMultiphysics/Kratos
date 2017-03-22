from __future__ import print_function, absolute_import, division

import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics

import sys

class TestModel(KratosUnittest.TestCase):        

    def test_model(self):
        model_part = KratosMultiphysics.ModelPart("Main")
        model_part.CreateSubModelPart("Inlets")
        model_part.CreateSubModelPart("Temp")
        outlet = model_part.CreateSubModelPart("Outlet")
        
        model = KratosMultiphysics.Model()
        model.AddModelPart(model_part)
                
        aaa = model["Main.Outlet"].CreateSubModelPart("aaa") 
        
        if (sys.version_info < (3, 2)):
            self.assertRaisesRegex = self.assertRaisesRegexp
        
        #check that a meaningful error is thrown
        with self.assertRaisesRegex(RuntimeError, "Error: The model part named : \"aaa\" was not found as root model part.  total input string was \"aaa\""):
            model["aaa"]
        
        #check that a meaningful error is thrown
        with self.assertRaisesRegex(RuntimeError, "The model part named : \"aaa\" was not found as submodelpart of : \"Inlets\" . total input string was \"Main.Inlets.aaa\""):
            model["Main.Inlets.aaa"]
            
        #here i create a model part in the lowest level and i check that the other model parts have it
        model["Main.Outlet.aaa"].CreateNewNode(1,0.0,0.0,0.0)
        self.assertEqual(len(model_part.Nodes), 1)
        self.assertEqual(len(outlet.Nodes), 1)
        self.assertEqual(len(aaa.Nodes), 1)

        #self.assertEqual(model["Main"], model_part )
        #self.assertEqual(model["Main.Outlet"].Info(), outlet.Info() )
        #self.assertEqual(model["Main.Outlet.aaa"].Info(), aaa.Info() )

        
if __name__ == '__main__':
    KratosUnittest.main()
