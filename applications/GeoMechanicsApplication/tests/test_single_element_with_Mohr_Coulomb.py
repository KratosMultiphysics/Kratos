import os

import KratosMultiphysics.KratosUnittest as KratosUnittest
import test_helper


class KratosGeoMechanicsSingleElementWithMohrCoulomb(KratosUnittest.TestCase):
    """
    This class contains a test for a single element with a Mohr-Coulomb material model.
    """

    def assert_displacement_field(self, time, output_data):
        displacements = test_helper.GiDOutputFileReader.nodal_values_at_time(
            "DISPLACEMENT", time, output_data
        )
        self.assertVectorAlmostEqual(displacements[0][:2], [0.0, 0.0])  # node 1
        self.assertVectorAlmostEqual(displacements[1][:2], [0.015, 0.0])  # node 2
        self.assertVectorAlmostEqual(displacements[2][:2], [0.015, -0.015])  # node 3
        self.assertVectorAlmostEqual(displacements[3][:2], [0.0, -0.015])  # node 4

    def assert_shear_capacity_values(self, time, output_data):
        shear_capacity_values_of_element_1 = (
            test_helper.GiDOutputFileReader.element_integration_point_values_at_time(
                "GEO_SHEAR_CAPACITY", time, output_data
            )[0]
        )
        self.assertVectorAlmostEqual(shear_capacity_values_of_element_1, [0.75] * 4)

    def test_shear_capacity_calculation(self):
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
        self.assert_displacement_field(end_time, output_data)
        self.assert_shear_capacity_values(end_time, output_data)


if __name__ == "__main__":
    KratosUnittest.main()
