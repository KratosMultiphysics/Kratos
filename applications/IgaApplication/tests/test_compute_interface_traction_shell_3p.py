import KratosMultiphysics as KM
import KratosMultiphysics.IgaApplication as IGA
import KratosMultiphysics.StructuralMechanicsApplication as SMA
import KratosMultiphysics.KratosUnittest as KratosUnittest

from KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_analysis import StructuralMechanicsAnalysis


class ComputeInterfaceTractionShell3pTest(KratosUnittest.TestCase):

    def test_compute_interface_traction_from_input(self):
        model = KM.Model()

        parameters = KM.Parameters("""
        {
            "project_parameters_file_name" : "compute_interface_traction_test/compute_interface_traction_project_parameters.json"
        }
        """)

        with open("compute_interface_traction_test/compute_interface_traction_project_parameters.json", "r") as parameter_file:
            parameters = KM.Parameters(parameter_file.read())

        analysis = StructuralMechanicsAnalysis(model, parameters)
        analysis.Run()

        model_part = model["IgaModelPart"]
        interface_model_part = model["IgaModelPart.Support_3"]

        shell_model_part = model["IgaModelPart.StructuralAnalysis_1"]
        shell_properties = shell_model_part.Elements[1].Properties

        vars_to_copy = [
            KM.THICKNESS,
            KM.YOUNG_MODULUS,
            KM.POISSON_RATIO,
            KM.CONSTITUTIVE_LAW
        ]

        for cond in interface_model_part.Conditions:
            cond_properties = cond.Properties

            for var in vars_to_copy:
                if shell_properties.Has(var) and not cond_properties.Has(var):
                    cond_properties.SetValue(var, shell_properties[var])

        interface_traction_variable = KM.KratosGlobals.GetVariable("INTERFACE_TRACTION")

        IGA.ComputeInterfaceTractionShell3pUtility.ComputeAndSetInterfaceTraction(
            interface_model_part,
            interface_traction_variable
        )

        for cond in interface_model_part.Conditions:
            traction = cond.GetValue(interface_traction_variable)

            self.assertAlmostEqual(traction[0], -1.0, 6)
            self.assertAlmostEqual(traction[1], 0.0, 6)
            self.assertAlmostEqual(traction[2], 0.0, 6)


if __name__ == "__main__":
    KratosUnittest.main()