from __future__ import print_function, absolute_import, division

# Importing the Kratos Library
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics

# Other imports
import os

def GetFilePath(fileName):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), fileName)

class TestMaterialsInput(KratosUnittest.TestCase):

    def test_input(self):
        try:
            import KratosMultiphysics.FluidDynamicsApplication as KratosFluid
        except:
            self.skipTest("KratosMultiphysics.FluidDynamicsApplication is not available")

        try:
            import KratosMultiphysics.StructuralMechanicsApplication
        except:
            self.skipTest("KratosMultiphysics.StructuralMechanicsApplication is not available")

        model_part = KratosMultiphysics.ModelPart("Main")
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VISCOSITY)
        model_part_io = KratosMultiphysics.ModelPartIO(GetFilePath("test_model_part_io_read")) #reusing the file that is already in the directory
        model_part_io.ReadModelPart(model_part)

        # Define a Model
        Model = KratosMultiphysics.Model()
        Model.AddModelPart(model_part)

        test_settings = KratosMultiphysics.Parameters("""
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
        self.assertTrue(model_part.Properties[1].GetValue(KratosMultiphysics.YOUNG_MODULUS) == 200.0)
        self.assertTrue(model_part.Properties[1].GetValue(KratosMultiphysics.POISSON_RATIO) == 0.3)
        self.assertTrue(model_part.Properties[1].GetValue(KratosMultiphysics.YIELD_STRESS) == 400.0)
        self.assertTrue(model_part.Properties[1].GetValue(KratosFluid.SUBSCALE_PRESSURE) == 0.1)
        self.assertTrue(model_part.Properties[1].GetValue(KratosFluid.VORTICITY_MAGNITUDE) == -5.888)

        self.assertTrue(model_part.Properties[2].GetValue(KratosMultiphysics.YOUNG_MODULUS) == 100.0)
        self.assertTrue(model_part.Properties[2].GetValue(KratosMultiphysics.POISSON_RATIO) == 0.1)
        self.assertTrue(model_part.Properties[2].GetValue(KratosMultiphysics.YIELD_STRESS) == 800.0)
        self.assertTrue(model_part.Properties[2].GetValue(KratosMultiphysics.HTC) == 0.3)
        self.assertTrue(model_part.Properties[2].GetValue(KratosMultiphysics.TIME_STEPS) == 159) # int
        self.assertTrue(model_part.Properties[2].GetValue(KratosMultiphysics.UPDATE_SENSITIVITIES) == True) # bool
        self.assertTrue(model_part.Properties[2].GetValue(KratosMultiphysics.IDENTIFIER) == "MyTestString") # std::string

        mat_vector = model_part.Properties[2].GetValue(KratosMultiphysics.CAUCHY_STRESS_VECTOR)
        self.assertAlmostEqual(mat_vector[0],1.5)
        self.assertAlmostEqual(mat_vector[1],0.3)
        self.assertAlmostEqual(mat_vector[2],-2.58)

        mat_matrix = model_part.Properties[2].GetValue(KratosMultiphysics.LOCAL_INERTIA_TENSOR)
        self.assertAlmostEqual(mat_matrix[0,0],1.27)
        self.assertAlmostEqual(mat_matrix[0,1],-22.5)
        self.assertAlmostEqual(mat_matrix[1,0],2.01)
        self.assertAlmostEqual(mat_matrix[1,1],0.257)

        table = model_part.Properties[2].GetTable(KratosMultiphysics.TEMPERATURE, KratosMultiphysics.YOUNG_MODULUS)
        self.assertAlmostEqual(table.GetValue(1.5),11.0)
        self.assertAlmostEqual(table.GetNearestValue(1.1),10.0)
        self.assertAlmostEqual(table.GetDerivative(1.2),2.0)


if __name__ == '__main__':
    KratosUnittest.main()
