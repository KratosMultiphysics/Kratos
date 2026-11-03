import KratosMultiphysics
import KratosMultiphysics.KratosUnittest as KratosUnittest

from KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_analysis import StructuralMechanicsAnalysis
from KratosMultiphysics.StructuralMechanicsApplication.structural_verification_process import StructuralVerificationProcess


class TestPuck(KratosUnittest.TestCase):


    def _GetPuckResult(self, simulation):
        for process in simulation._GetListOfProcesses():
            if isinstance(process, StructuralVerificationProcess):
                return process.results[("Shell_gesamt", "puck")]

        self.fail("Could not find StructuralVerificationProcess.")


    def _CreateParametersWithPuckCriterion(self, criterion):
        with open("ProjectParameters.json", "r") as parameter_file:
            parameters = KratosMultiphysics.Parameters(parameter_file.read())

        with open("handbook_config.json", "r") as handbook_file:
            handbook_settings = KratosMultiphysics.Parameters(handbook_file.read())

        handbook_settings["Structural_Components"][0]["metadata"]["criterion"].SetString(criterion)

        process_settings = parameters["processes"]["custom_processes"][0]
        process_settings.RemoveValue("Parameters")
        process_settings.AddValue("Parameters", handbook_settings)

        return parameters

#helper for RF_IFF 
    def _CheckPuckRF_IFF(self, criterion, indices, expected_RF_IFF):
        model, simulation = self._RunPuckAnalysis(criterion)

        puck_result = self._GetPuckResult(simulation)

        RF_IFF = puck_result.metadata["RF_IFF"]

        for index, expected in zip(indices, expected_RF_IFF):
            self.assertAlmostEqual(RF_IFF[index], expected, places=12)

    def _RunPuckAnalysis(self, criterion):
        parameters = self._CreateParametersWithPuckCriterion(criterion)

        model = KratosMultiphysics.Model()
        simulation = StructuralMechanicsAnalysis(model, parameters)
        simulation.Run()

        return model, simulation


#helper for RF_FF 
    def _CheckPuckRF_FF(self, criterion, indices, expected_RF_FF):
        model, simulation = self._RunPuckAnalysis(criterion)

        puck_result = self._GetPuckResult(simulation)

        RF_FF = puck_result.metadata["RF_FF"]

        for index, expected in zip(indices, expected_RF_FF):
            self.assertAlmostEqual(RF_FF[index], expected, places=12)


#helper for PuckDegradation
    def _CheckPuck_DegradationIFF(self, criterion, indices, expected_Puck_DegradationIFF):
        model, simulation = self._RunPuckAnalysis(criterion)

        puck_result = self._GetPuckResult(simulation)

        Puck_DegradationIFF= puck_result.metadata["Puck_DegradationIFF"]

        for index, expected in zip(indices, expected_Puck_DegradationIFF):
            self.assertAlmostEqual(Puck_DegradationIFF[index], expected, places=4)            

    def _CheckPuck_DegradationFF(self, criterion, indices, expected_Puck_DegradationFF):
        model, simulation = self._RunPuckAnalysis(criterion)

        puck_result = self._GetPuckResult(simulation)

        Puck_DegradationFF = puck_result.metadata["Puck_DegradationFF"]

        for index, expected in zip(indices, expected_Puck_DegradationFF):
            self.assertAlmostEqual(Puck_DegradationFF[index], expected, places=12)
    
    def test_PuckDegradation_1(self):
        with KratosUnittest.WorkFolderScope("puck_test", __file__):

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
                4.686103513446605,
                float("inf"),
                4.709971846671217,
                float("inf"),
                float("inf")
            ]

            expected_PuckDegradationFF = [
                1.1898674777487814,
                3.355294616026732,
                1.1168049005964544,
                3.599349043407229,
                1.8172274680886131
            ]

            self._CheckPuckRF_IFF(
                        "last_ply_failure",
                        indices,
                        expected_RF_IFF
                    )
            self._CheckPuckRF_FF(
                        "last_ply_failure",
                        indices,
                        expected_RF_FF
                    )
            self._CheckPuck_DegradationIFF(
                        "last_ply_failure",
                        indices,
                        expected_PuckDegradationIFF
                    )

            self._CheckPuck_DegradationFF(
                        "last_ply_failure",
                        indices,
                        expected_PuckDegradationFF
                    )


    def test_PuckDegradation_2(self):
        with KratosUnittest.WorkFolderScope("puck_test", __file__):

            indices = [0, 27, 62, 84, 121]

            expected_PuckDegradationIFF = [
                4.68610351344662,
                0.11722104559856776,
                4.709971846671214,
                0.11722104559856777,
                0.027875110252034884
            ]

            self._CheckPuck_DegradationIFF(
                        "first_fiber_failure",
                        indices,
                        expected_PuckDegradationIFF
                    )



    def test_PuckDegradation_3(self):
        with KratosUnittest.WorkFolderScope("puck_test", __file__):

            indices = [0, 27, 62, 84, 121]

            expected_PuckDegradationIFF = [
                0,
                0,
                0,
                0,
                0
            ]

            self._CheckPuck_DegradationIFF(
                        "first_ply_failure",
                        indices,
                        expected_PuckDegradationIFF
                    )




if __name__ == "__main__":
    KratosUnittest.main()