import os
import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.GeoMechanicsApplication.context_managers as context_managers
import KratosMultiphysics.GeoMechanicsApplication.geomechanics_analysis as analysis
from KratosMultiphysics.GeoMechanicsApplication.gid_output_file_reader import GiDOutputFileReader

import test_helper

class KratosGeoMechanicsCPhiReductionProcess(KratosUnittest.TestCase):
    """
    This class contains benchmark tests which are checked with the original solution
    """

    @staticmethod
    def _run_stages_with_common_parameters(case_path, common_path, n_stages):
        model = Kratos.Model()
        with context_managers.set_cwd_to(case_path):
            for stage_index in range(1, n_stages + 1):
                project_parameters_file = os.path.join(
                    common_path, f"ProjectParameters_stage{stage_index}.json"
                )
                with open(project_parameters_file, "r") as parameter_file:
                    stage_parameters = Kratos.Parameters(parameter_file.read())
                analysis.GeoMechanicsAnalysis(model, stage_parameters).Run()

    @staticmethod
    def _get_displacement_vectors(case_path):
        common_path = os.path.join(os.path.dirname(case_path), "common")
        KratosGeoMechanicsCPhiReductionProcess._run_stages_with_common_parameters(
            case_path, common_path, 2
        )
        reader = GiDOutputFileReader()
        actual_data = reader.read_output_from(os.path.join(case_path, "stage2.post.res"))
        node_ids = [57, 151, 235]
        return reader.nodal_values_at_time("DISPLACEMENT", 0.2, actual_data, node_ids)

    def _assert_case(self, case_name, expected_displacements_x):
        file_path = test_helper.get_file_path(
            os.path.join("C-Phi_reduction_process", case_name)
        )
        displacement_vectors = self._get_displacement_vectors(file_path)
        for index, expected_x in enumerate(expected_displacements_x):
            self.assertAlmostEqual(expected_x, displacement_vectors[index][0])

    def test_c_phi_reduction_process(self):
        test_cases = {
            "Mohr_Coulomb_model64": [-0.002667, -0.0096777, -0.0115495],
            "Mohr_Coulomb_model": [-0.00268182, -0.00965477,  -0.011628],
        }

        for case_name, expected_displacements_x in test_cases.items():
            with self.subTest(case=case_name):
                self._assert_case(case_name, expected_displacements_x)

if __name__ == '__main__':
    KratosUnittest.main()
