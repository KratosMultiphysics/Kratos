import os

import KratosMultiphysics.KratosUnittest as KratosUnittest
import test_helper


class KratosGeoMechanicsSingleElementWithMohrCoulomb(KratosUnittest.TestCase):
    """
    This class contains a test for a single element with a Mohr-Coulomb material model.
    """

    def assert_displacement_field(
        self, time, output_data, expected_displacement_vectors
    ):
        actual_displacement_vectors = (
            test_helper.GiDOutputFileReader.nodal_values_at_time(
                "DISPLACEMENT", time, output_data
            )
        )

        for actual_vector, expected_vector in zip(
            actual_displacement_vectors, expected_displacement_vectors
        ):
            self.assertVectorAlmostEqual(actual_vector[:2], expected_vector)

    def assert_shear_capacity_values(self, time, output_data, expected_values):
        actual_values = (
            test_helper.GiDOutputFileReader.element_integration_point_values_at_time(
                "GEO_SHEAR_CAPACITY", time, output_data
            )[0]
        )
        self.assertVectorAlmostEqual(actual_values, expected_values)

    def test_shear_capacity_calculation_of_UPwSmallStrainElement2D4N(self):
        test_dir_path = test_helper.get_file_path(
            os.path.join(
                "single_element_with_Mohr_Coulomb", "UPwSmallStrainElement2D4N"
            )
        )

        test_helper.run_kratos(test_dir_path)

        output_file_path = os.path.join(test_dir_path, "output.post.res")
        reader = test_helper.GiDOutputFileReader()
        output_data = reader.read_output_from(output_file_path)
        end_time = 1.0

        expected_displacement_vectors = [
            [0.0, 0.0],  # node 1
            [0.015, 0.0],  # node 2
            [0.015, -0.015],  # node 3
            [0.0, -0.015],  # node 4
        ]
        self.assert_displacement_field(
            end_time, output_data, expected_displacement_vectors
        )

        expected_shear_capacity_values = [0.75] * 4
        self.assert_shear_capacity_values(
            end_time, output_data, expected_shear_capacity_values
        )

    def test_shear_capacity_calculattion_of_SmallStrainUPwDiffOrderElement2D8N(self):
        test_dir_path = test_helper.get_file_path(
            os.path.join(
                "single_element_with_Mohr_Coulomb", "SmallStrainUPwDiffOrderElement2D8N"
            )
        )

        test_helper.run_kratos(test_dir_path)

        output_file_path = os.path.join(test_dir_path, "output.post.res")
        reader = test_helper.GiDOutputFileReader()
        output_data = reader.read_output_from(output_file_path)
        end_time = 1.0

        expected_displacements = [
            [0.0, 0.0],  # node 1
            [0.015, 0.0],  # node 2
            [0.015, -0.015],  # node 3
            [0.0, -0.015],  # node 4
            [0.0075, 0.0],  # node 5
            [0.015, -0.0075],  # node 6
            [0.0075, -0.015],  # node 7
            [0.0, -0.0075],  # node 8
        ]
        self.assert_displacement_field(end_time, output_data, expected_displacements)

        expected_shear_capacity_values = [0.75] * 4
        self.assert_shear_capacity_values(
            end_time, output_data, expected_shear_capacity_values
        )


if __name__ == "__main__":
    KratosUnittest.main()
