from __future__ import print_function, absolute_import, division

# Importing the Kratos Library
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics
try:
    import KratosMultiphysics.FluidDynamicsApplication as KratosFluid
    import KratosMultiphysics.StructuralMechanicsApplication
    missing_external_dependencies = False
    missing_application = ''
except ImportError as e:
    missing_external_dependencies = True
    # extract name of the missing application from the error message
    import re
    missing_application = re.search(r'''.*'KratosMultiphysics\.(.*)'.*''',
                                    '{0}'.format(e)).group(1)

# Other imports
import os

def GetFilePath(fileName):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), fileName)

class TestMaterialsInput(KratosUnittest.TestCase):

    def _prepare_test(self, input_file = "materials.json"):
        self.model_part = KratosMultiphysics.ModelPart("Main")
        self.model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        self.model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VISCOSITY)
        self.model_part_io = KratosMultiphysics.ModelPartIO(GetFilePath("test_model_part_io_read")) #reusing the file that is already in the directory
        self.model_part_io.ReadModelPart(self.model_part)

        # Define a Model
        self.Model = KratosMultiphysics.Model()
        self.Model.AddModelPart(self.model_part)

        self.test_settings = KratosMultiphysics.Parameters("""
        {
            "Parameters": {
                    "materials_filename": "materials.json"
            }
        }
        """)

        #assign the real path
        self.test_settings["Parameters"]["materials_filename"].SetString(GetFilePath(input_file))

    def _check_results(self):
        #test if the element properties are assigned correctly to the elements and conditions
        for elem in self.Model["Inlets"].Elements:
            self.assertEqual(elem.Properties.Id, 1)
        for cond in self.Model["Inlets"].Conditions:
            self.assertEqual(cond.Properties.Id, 1)
        for elem in self.Model["Outlet"].Elements:
            self.assertEqual(elem.Properties.Id, 2)
        for cond in self.Model["Outlet"].Conditions:
            self.assertEqual(cond.Properties.Id, 2)

        #test that the properties are read correctly
        self.assertEqual(self.model_part.Properties[1].GetValue(KratosMultiphysics.YOUNG_MODULUS), 200.0)
        self.assertEqual(self.model_part.Properties[1].GetValue(KratosMultiphysics.POISSON_RATIO), 0.3)
        self.assertEqual(self.model_part.Properties[1].GetValue(KratosMultiphysics.YIELD_STRESS), 400.0)
        self.assertEqual(self.model_part.Properties[1].GetValue(KratosFluid.SUBSCALE_PRESSURE), 0.1)
        self.assertEqual(self.model_part.Properties[1].GetValue(KratosFluid.VORTICITY_MAGNITUDE), -5.888)

        self.assertEqual(self.model_part.Properties[2].GetValue(KratosMultiphysics.YOUNG_MODULUS), 100.0)
        self.assertEqual(self.model_part.Properties[2].GetValue(KratosMultiphysics.POISSON_RATIO), 0.1)
        self.assertEqual(self.model_part.Properties[2].GetValue(KratosMultiphysics.YIELD_STRESS), 800.0)
        self.assertEqual(self.model_part.Properties[2].GetValue(KratosMultiphysics.HTC), 0.3)
        self.assertEqual(self.model_part.Properties[2].GetValue(KratosMultiphysics.TIME_STEPS), 159) # int
        self.assertEqual(self.model_part.Properties[2].GetValue(KratosMultiphysics.UPDATE_SENSITIVITIES), True) # bool
        self.assertEqual(self.model_part.Properties[2].GetValue(KratosMultiphysics.IDENTIFIER), "MyTestString") # std::string

        mat_vector = self.model_part.Properties[2].GetValue(KratosMultiphysics.CAUCHY_STRESS_VECTOR)
        self.assertAlmostEqual(mat_vector[0],1.5)
        self.assertAlmostEqual(mat_vector[1],0.3)
        self.assertAlmostEqual(mat_vector[2],-2.58)

        mat_matrix = self.model_part.Properties[2].GetValue(KratosMultiphysics.LOCAL_INERTIA_TENSOR)
        self.assertAlmostEqual(mat_matrix[0,0],1.27)
        self.assertAlmostEqual(mat_matrix[0,1],-22.5)
        self.assertAlmostEqual(mat_matrix[1,0],2.01)
        self.assertAlmostEqual(mat_matrix[1,1],0.257)

        table = self.model_part.Properties[2].GetTable(KratosMultiphysics.TEMPERATURE, KratosMultiphysics.YOUNG_MODULUS)
        self.assertAlmostEqual(table.GetValue(1.5),11.0)
        self.assertAlmostEqual(table.GetNearestValue(1.1),10.0)
        self.assertAlmostEqual(table.GetDerivative(1.2),2.0)

    def _check_results_with_subproperties(self):
        prop1 = self.model_part.GetProperties()[1]
        self.assertEqual(prop1.NumberOfSubproperties(), 3)

    def test_input_python(self):

        if (missing_external_dependencies is True):
            self.skipTest("{} is not available".format(missing_application))
        self._prepare_test()
        import read_materials_process
        read_materials_process.Factory(self.test_settings,self.Model)
        self._check_results()

    def test_input_cpp(self):

        if (missing_external_dependencies is True):
            self.skipTest("{} is not available".format(missing_application))
        self._prepare_test()
        KratosMultiphysics.ReadMaterialsUtility(self.test_settings, self.Model)
        self._check_results()

    def test_input_with_subproperties_cpp(self):

        if (missing_external_dependencies is True):
            self.skipTest("{} is not available".format(missing_application))
        self._prepare_test("materials_with_subproperties.json")
        KratosMultiphysics.ReadMaterialsUtility(self.test_settings, self.Model)
        self._check_results_with_subproperties()

if __name__ == '__main__':
    KratosUnittest.main()
