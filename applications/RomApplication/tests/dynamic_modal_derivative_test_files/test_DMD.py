import KratosMultiphysics

try:
    import json
    from .MainKratosDMD import StructuralMechanicsModalDerivativeAnalysisMSConstraints
    json_available = True
except:
    json_available = False

import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.kratos_utilities as kratos_utilities

class TestDynamicModalDerivative(KratosUnittest.TestCase):

    @KratosUnittest.skipIf(json_available == False, "json is required for RomApplication.StructuralMechanicsModalDerivativeAnalysis")
    def test_structural_mechanics_dynamic_modal_derivative(self):

        with KratosUnittest.WorkFolderScope(".", __file__):
            with open("ProjectParameters.json",'r') as parameter_file:
                parameters = KratosMultiphysics.Parameters(parameter_file.read())
            model = KratosMultiphysics.Model()
            simulation = StructuralMechanicsModalDerivativeAnalysisMSConstraints(model,parameters)
            simulation.Run()

            obtained_output = None
            with open("RomParameters.json",'r') as obtained_output_file:
                obtained_output = json.load(obtained_output_file)

            expected_output = None
            with open("ExpectedRomParameters.json",'r') as expected_output_file:
                expected_output = json.load(expected_output_file)

            self.assertEqual(len(expected_output["rom_settings"]["nodal_unknowns"]), len(obtained_output["rom_settings"]["nodal_unknowns"]))
            self.assertEqual(expected_output["rom_settings"]["number_of_rom_dofs"], obtained_output["rom_settings"]["number_of_rom_dofs"])
            self.assertVectorAlmostEqual(expected_output["eigenvalues"], obtained_output["eigenvalues"])
            self.assertEqual(len(expected_output["nodal_modes"]), len(obtained_output["nodal_modes"]))
            
            for (k,v) in expected_output["nodal_modes"].items():
                for i in range(len(v)):
                    self.assertVectorAlmostEqual(v[i], obtained_output["nodal_modes"][k][i])

            # Cleaning
            kratos_utilities.DeleteDirectoryIfExisting("__pycache__")

if __name__ == '__main__':
    KratosMultiphysics.Logger.GetDefaultOutput().SetSeverity(KratosMultiphysics.Logger.Severity.WARNING)
    KratosUnittest.main()