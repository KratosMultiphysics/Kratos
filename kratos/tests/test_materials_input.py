from __future__ import print_function, absolute_import, division

import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics import *

def GetFilePath(fileName):
    return os.path.dirname(os.path.realpath(__file__)) + "/" + fileName

class TestMaterialsInput(KratosUnittest.TestCase):

    def test_input(self):
        try:
            import KratosMultiphysics.SolidMechanicsApplication
        except:
            self.skipTest("KratosMultiphysics.SolidMechanicsApplication is not available")
        
        model_part = ModelPart("Main")
        model_part.AddNodalSolutionStepVariable(DISPLACEMENT)
        model_part.AddNodalSolutionStepVariable(VISCOSITY)
        model_part_io = ModelPartIO(GetFilePath("test_model_part_io_read")) #reusing the file that is already in the directory
        model_part_io.ReadModelPart(model_part)
    
        #define a Model TODO: replace to use the real Model once available
        Model = {
            "Main" : model_part,
            "Inlets" : model_part.GetSubModelPart("Inlets"),
            "Inlet1" : model_part.GetSubModelPart("Inlets").GetSubModelPart("Inlet1"),
            "Inlet2" : model_part.GetSubModelPart("Inlets").GetSubModelPart("Inlet2"),
            "Outlet" : model_part.GetSubModelPart("Outlet")
            }
        
        test_settings = Parameters("""
            {
                    "Parameters": {
                            "materials_filename": "materials.json"
                    }
            }
            """)
        
        #assign the real path 
        test_settings["Parameters"]["materials_filename"].SetString(GetFilePath("materials.json"))
        
        import read_materials_process
        read_materials_process.Factory(test_settings,Model)
        
        #test if the element properties are assigned correctly to the elements and conditions
        for elem in Model["Inlets"].Elements:
            self.assertTrue(elem.Properties.Id == 1)
        for cond in Model["Inlets"].Conditions:
            self.assertTrue(cond.Properties.Id == 1)
        for elem in Model["Outlet"].Elements:
            self.assertTrue(elem.Properties.Id == 2)
        for cond in Model["Outlet"].Conditions:
            self.assertTrue(cond.Properties.Id == 2)   
                    
        #test that the properties are read correctly
        self.assertTrue(model_part.Properties[1].GetValue(YOUNG_MODULUS) == 200.0)
        self.assertTrue(model_part.Properties[1].GetValue(POISSON_RATIO) == 0.3)
        self.assertTrue(model_part.Properties[1].GetValue(YIELD_STRESS) == 400.0)

        self.assertTrue(model_part.Properties[2].GetValue(YOUNG_MODULUS) == 100.0)
        self.assertTrue(model_part.Properties[2].GetValue(POISSON_RATIO) == 0.1)
        self.assertTrue(model_part.Properties[2].GetValue(YIELD_STRESS) == 800.0)
        self.assertTrue(model_part.Properties[2].GetValue(HTC) == 0.3)
        self.assertTrue(model_part.Properties[2].GetValue(TIME_STEPS) == 159) # int
        self.assertTrue(model_part.Properties[2].GetValue(UPDATE_SENSITIVITIES) == True) # bool
        self.assertTrue(model_part.Properties[2].GetValue(IDENTIFIER) == "MyTestString") # std::string

        mat_vector = model_part.Properties[2].GetValue(CAUCHY_STRESS_VECTOR)
        self.assertAlmostEqual(mat_vector[0],1.5)
        self.assertAlmostEqual(mat_vector[1],0.3)
        self.assertAlmostEqual(mat_vector[2],-2.58)

        mat_matrix = model_part.Properties[2].GetValue(LOCAL_INERTIA_TENSOR)
        self.assertAlmostEqual(mat_matrix[0,0],1.27)
        self.assertAlmostEqual(mat_matrix[0,1],-22.5)
        self.assertAlmostEqual(mat_matrix[1,0],2.01)
        self.assertAlmostEqual(mat_matrix[1,1],0.257)

        table = model_part.Properties[2].GetTable(TEMPERATURE, YOUNG_MODULUS)
        self.assertAlmostEqual(table.GetValue(1.5),11.0)
        self.assertAlmostEqual(table.GetNearestValue(1.1),10.0)
        self.assertAlmostEqual(table.GetDerivative(1.2),2.0)
        

if __name__ == '__main__':
    KratosUnittest.main()
