import os
import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics.GeoMechanicsApplication.gid_output_file_reader import GiDOutputFileReader
import KratosMultiphysics.GeoMechanicsApplication.run_multiple_stages as run_multiple_stages

import test_helper

class KratosGeoMechanicsCPhiReductionProcess(KratosUnittest.TestCase):
    """
    This class contains benchmark tests which are checked with the original solution
    """
    @staticmethod
    def _get_displacement_vectors(file_path):
        run_multiple_stages.run_stages(file_path, 2)
        reader = GiDOutputFileReader()
        actual_data = reader.read_output_from(os.path.join(file_path, "stage2.post.res"))
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
