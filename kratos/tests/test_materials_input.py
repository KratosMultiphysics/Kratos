from __future__ import print_function, absolute_import, division

import os
import sys

# Importing the Kratos Library
import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as KratosUtils
from KratosMultiphysics import read_materials_process

dependencies_are_available = KratosUtils.CheckIfApplicationsAvailable("StructuralMechanicsApplication", "FluidDynamicsApplication")
if dependencies_are_available:
    import KratosMultiphysics.FluidDynamicsApplication as KratosFluid
    import KratosMultiphysics.StructuralMechanicsApplication


def GetFilePath(fileName):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), fileName)

class TestMaterialsInput(KratosUnittest.TestCase):
    def _prepare_test(self, input_file = "materials.json"):
        # Define a Model
        self.current_model = KratosMultiphysics.Model()

        self.model_part = self.current_model.CreateModelPart("Main")

        self.model_part.AddNodalSolutionStepVariable(KratosMultiphysics.DISPLACEMENT)
        self.model_part.AddNodalSolutionStepVariable(KratosMultiphysics.VISCOSITY)
        self.model_part_io = KratosMultiphysics.ModelPartIO(GetFilePath("auxiliar_files_for_python_unittest/mdpa_files/test_model_part_io_read")) #reusing the file that is already in the directory
        self.model_part_io.ReadModelPart(self.model_part)

        self.test_settings = KratosMultiphysics.Parameters("""
        {
            "Parameters": {
                    "materials_filename": "materials.json"
            }
        }
        """)

        # Assign the real path
        self.test_settings["Parameters"]["materials_filename"].SetString(GetFilePath("auxiliar_files_for_python_unittest/materials_files/" + input_file))

    def _check_results(self):
        #test if the element properties are assigned correctly to the elements and conditions
        for elem in self.current_model["Main.Inlets"].Elements:
            self.assertEqual(elem.Properties.Id, 1)
        for cond in self.current_model["Main.Inlets"].Conditions:
            self.assertEqual(cond.Properties.Id, 1)
        for elem in self.current_model["Main.Outlet"].Elements:
            self.assertEqual(elem.Properties.Id, 2)
        for cond in self.current_model["Main.Outlet"].Conditions:
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

        self.assertEqual(prop1.HasSubProperties(11), True)
        self.assertEqual(prop1.HasSubProperties(120), False)

        sub_prop11 = prop1.GetSubProperties(11)
        self.assertEqual(sub_prop11.GetValue(KratosMultiphysics.YOUNG_MODULUS), 206900000000.0)
        self.assertEqual(sub_prop11.GetValue(KratosMultiphysics.POISSON_RATIO), 0.29)
        self.assertEqual(sub_prop11.GetValue(KratosMultiphysics.THICKNESS), 0.000889)

        self.assertEqual(sub_prop11.NumberOfSubproperties(), 3)

    @KratosUnittest.skipUnless(dependencies_are_available,"StructuralMechanicsApplication or FluidDynamicsApplication are not available")
    def test_input_python(self):
        self._prepare_test()
        read_materials_process.Factory(self.test_settings,self.current_model)
        self._check_results()

    @KratosUnittest.skipUnless(dependencies_are_available,"StructuralMechanicsApplication or FluidDynamicsApplication are not available")
    def test_input_cpp(self):
        self._prepare_test()

        KratosMultiphysics.ReadMaterialsUtility(self.test_settings, self.current_model)
        self._check_results()

    def test_input_with_subproperties_cpp(self):

        if not dependencies_are_available:
            self.skipTest("Dependencies are not available")
        self._prepare_test("materials_with_subproperties.json")
        KratosMultiphysics.ReadMaterialsUtility(self.test_settings, self.current_model)
        self._check_results_with_subproperties()

    @KratosUnittest.expectedFailure
    def test_input_with_subproperties_cpp_expected_failure(self):

        if not dependencies_are_available:
            self.skipTest("Dependencies are not available")
        self._prepare_test("materials_with_subproperties_expected_failure.json")
        KratosMultiphysics.ReadMaterialsUtility(self.test_settings, self.current_model)

    def test_overdefined_materials(self):
        current_model = KratosMultiphysics.Model()
        test_settings = KratosMultiphysics.Parameters(""" { "Parameters": { "materials_filename": ""}} """)


        test_settings["Parameters"]["materials_filename"].SetString(
            GetFilePath(os.path.join("auxiliar_files_for_python_unittest","materials_files","wrong_materials_input","wrong_materials_1.json")))
        expected_error_msg = "Error: Materials for ModelPart \"Main\" are specified multiple times!"
        with self.assertRaisesRegex(RuntimeError, expected_error_msg):
            KratosMultiphysics.ReadMaterialsUtility(test_settings, current_model)

        test_settings["Parameters"]["materials_filename"].SetString(
            GetFilePath(os.path.join("auxiliar_files_for_python_unittest","materials_files","wrong_materials_input","wrong_materials_2.json")))
        expected_error_msg = "Error: Materials for ModelPart \"Main.sub\" are specified multiple times!"
        with self.assertRaisesRegex(RuntimeError, expected_error_msg):
            KratosMultiphysics.ReadMaterialsUtility(test_settings, current_model)

        test_settings["Parameters"]["materials_filename"].SetString(
            GetFilePath(os.path.join("auxiliar_files_for_python_unittest","materials_files","wrong_materials_input","wrong_materials_3.json")))
        expected_error_msg =  "Error: Materials for SubModelPart \"Main.sub\" is being overrided by Parent Model Part \"Main\"!\n"
        with self.assertRaisesRegex(RuntimeError, expected_error_msg):
            KratosMultiphysics.ReadMaterialsUtility(test_settings, current_model)

        test_settings["Parameters"]["materials_filename"].SetString(
            GetFilePath(os.path.join("auxiliar_files_for_python_unittest","materials_files","wrong_materials_input","wrong_materials_4.json")))
        expected_error_msg =  "Error: Materials for SubModelPart \"Main.sub1.subsub\" is being overrided by Parent Model Part \"Main.sub1\"!\n"
        with self.assertRaisesRegex(RuntimeError, expected_error_msg):
            KratosMultiphysics.ReadMaterialsUtility(test_settings, current_model)

        test_settings["Parameters"]["materials_filename"].SetString(
            GetFilePath(os.path.join("auxiliar_files_for_python_unittest","materials_files","wrong_materials_input","wrong_materials_5.json")))
        expected_error_msg =  "Error: Materials for SubModelPart \"Main.sub1.subsub\" is being overrided by Parent Model Part \"Main\"!\n"
        with self.assertRaisesRegex(RuntimeError, expected_error_msg):
            KratosMultiphysics.ReadMaterialsUtility(test_settings, current_model)

    @KratosUnittest.skipUnless(dependencies_are_available,"StructuralMechanicsApplication or FluidDynamicsApplication are not available")
    def test_input_without_tables_and_variables(self):
        self._prepare_test()
        self.test_settings["Parameters"]["materials_filename"].SetString(GetFilePath("auxiliar_files_for_python_unittest/materials_files/material_without_tables_and_variables.json"))

        KratosMultiphysics.ReadMaterialsUtility(self.test_settings, self.current_model)
        for elem in self.current_model["Main.Inlets"].Elements:
            self.assertEqual(elem.Properties.Id, 1)
        for cond in self.current_model["Main.Inlets"].Conditions:
            self.assertEqual(cond.Properties.Id, 1)
        for elem in self.current_model["Main.Outlet"].Elements:
            self.assertEqual(elem.Properties.Id, 2)
        for cond in self.current_model["Main.Outlet"].Conditions:
            self.assertEqual(cond.Properties.Id, 2)

if __name__ == '__main__':
    KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)
    KratosUnittest.main()
