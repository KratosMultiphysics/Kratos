import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest

from KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_analysis import StructuralMechanicsAnalysis


class TestPuck(KratosUnittest.TestCase):

    def test_PuckDegradation_1(self):
        with KratosUnittest.WorkFolderScope("puck_test", __file__):

            model, simulation = self._RunPuckAnalysis()

            RF_IFF = ...
            RF_FF = ...
            PuckDegradationIFF = ...

            indices = [0, 27, 62, 84, 121]

            expected_RF_FF = [
                1.3761517610369765,
                3.644381497707182,
                1.3000970265279104,
                3.9012832872427974,
                2.3958006885442122
            ]

            expected_RF_IFF = [
                3.952726347844309,
                0.8329307052135733,
                5.006764440715333,
                0.8166329183450516,
                0.5600795767528397
            ]

            expected_PuckDegradationIFF = [
                4.68610351344662,
                0.11722104559856777,
                4.709971846671214,
                0.11722104559856777,
                0.02787511025203488
            ]

            for i, expected in zip(indices, expected_RF_FF):
                self.assertAlmostEqual(RF_FF[i], expected, places=12)

            for i, expected in zip(indices, expected_RF_IFF):
                self.assertAlmostEqual(RF_IFF[i], expected, places=12)

            for i, expected in zip(indices, expected_PuckDegradationIFF):
                self.assertAlmostEqual(PuckDegradationIFF[i], expected, places=12)


    def test_PuckDegradation_2(self):
        with KratosUnittest.WorkFolderScope("puck", __file__):

            model, simulation = self._RunPuckAnalysis()

            PuckDegradationIFF = ...

            indices = [0, 27, 62, 84, 121]

            expected_PuckDegradationIFF = [
                4.68610351344662,
                0.11722104559856776,
                4.709971846671214,
                0.11722104559856777,
                0.027875110252034884
            ]

            for i, expected in zip(indices, expected_PuckDegradationIFF):
                self.assertAlmostEqual(
                    PuckDegradationIFF[i],
                    expected,
                    places=12
                )


    def test_PuckDegradation_3(self):
        with KratosUnittest.WorkFolderScope("puck", __file__):

            model, simulation = self._RunPuckAnalysis()

            PuckDegradationIFF = ...

            indices = [0, 27, 62, 84, 121]

            expected_PuckDegradationIFF = [
                0,
                0,
                0,
                0,
                0
            ]

            for i, expected in zip(indices, expected_PuckDegradationIFF):
                self.assertAlmostEqual(
                    PuckDegradationIFF[i],
                    expected,
                    places=12
                )


    def _RunPuckAnalysis(self):
        with open("ProjectParameters.json", "r") as parameter_file:
            parameters = KratosMultiphysics.Parameters(
                parameter_file.read()
            )

        model = KratosMultiphysics.Model()

        simulation = StructuralMechanicsAnalysis(
            model,
            parameters
        )

        simulation.Run()

        return model, simulation


if __name__ == "__main__":
    KratosUnittest.main()