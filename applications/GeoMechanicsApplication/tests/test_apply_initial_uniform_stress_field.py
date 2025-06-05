import unittest
import test_helper
import os.path


class KratosGeoMechanicsApplyInitialUniformStressFieldTests(unittest.TestCase):
    """
    This test applies an initial uniform stress field to a single three dimensional, linear hexahedron (3D8N)
    `UPwSmallStrainElement` using the `ApplyInitialUniformStressField` process.
    The test checks that the initial stress field is applied correctly,
    while all DoF are fixed (i.e. all displacements and water pressures).
    """

    def test_application_of_uniform_stress_field(self):
        output = self.run_simulation("apply_initial_uniform_stress_field")
        stress_vectors = (
            test_helper.GiDOutputFileReader.element_integration_point_values_at_time(
                "CAUCHY_STRESS_VECTOR", 1.0, output, [1]
            )[0]
        )
        expected_stress_vector = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]
        for stress_vector_at_integration_point in stress_vectors:
            for actual_stress_vector_component, expected_stress_vector_component in zip(
                stress_vector_at_integration_point, expected_stress_vector
            ):
                self.assertAlmostEqual(
                    actual_stress_vector_component, expected_stress_vector_component, 6
                )

    @staticmethod
    def run_simulation(test_name):
        file_path = test_helper.get_file_path(test_name)
        test_helper.run_kratos(file_path)
        output_file_path = os.path.join(file_path, "output.post.res")
        output_reader = test_helper.GiDOutputFileReader()
        return output_reader.read_output_from(output_file_path)


if __name__ == "__main__":
    unittest.main()
