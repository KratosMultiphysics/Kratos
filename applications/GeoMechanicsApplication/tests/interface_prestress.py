import os
import json
from pathlib import Path

import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics.GeoMechanicsApplication.gid_output_file_reader import (
    GiDOutputFileReader,
)
import KratosMultiphysics.GeoMechanicsApplication.run_multiple_stages as run_multiple_stages
import test_helper


class KratosGeoMechanicsInterfacePreStressTests(KratosUnittest.TestCase):
    def setUp(self):
        super().setUp()
        self.test_path = test_helper.get_file_path("test_interface_prestress")

    def check_vertical_stage_displacements(
        self, output_data, time, expected_vertical_stage_displacement
    ):
        stage_displacements = GiDOutputFileReader.nodal_values_at_time(
            "DISPLACEMENT", time, output_data
        )
        for displacement_vector in stage_displacements:
            self.assertAlmostEqual(
                displacement_vector[1], expected_vertical_stage_displacement, places=6
            )

    def read_interface_output_of_stage_2(self):
        with open(
            Path(self.test_path) / "stage_2_interface_output.json", "r"
        ) as output_file:
            return json.load(output_file)

    def check_output_times(self, output_data, expected_time):
        times = output_data["TIME"]
        self.assertEqual(len(times), 1)
        self.assertAlmostEqual(times[0], expected_time)

    def _check_output_vectors(
        self,
        output_data,
        expected_vector,
        output_item_label,
        name_of_output_item,
        name_of_normal_component,
        name_of_shear_component,
    ):
        element_label = "ELEMENT_3"
        output_vectors_by_integration_point_label = output_data[element_label][
            output_item_label
        ]
        number_of_integration_points = 3
        for integration_point_label in [
            str(i) for i in range(number_of_integration_points)
        ]:
            output_vectors = output_vectors_by_integration_point_label[
                integration_point_label
            ]
            self.assertEqual(
                len(output_vectors),
                1,
                msg=f"number of {name_of_output_item} vectors at integration point {integration_point_label} of {element_label}",
            )
            self.assertEqual(
                len(output_vectors[0]),
                2,
                msg=f"number of {name_of_output_item} components at integration point {integration_point_label} of {element_label}",
            )
            self.assertAlmostEqual(
                output_vectors[0][0],
                expected_vector[0],
                msg=f"{name_of_normal_component} at integration point {integration_point_label} of {element_label}",
            )
            self.assertAlmostEqual(
                output_vectors[0][1],
                expected_vector[1],
                msg=f"{name_of_shear_component} at integration point {integration_point_label} of {element_label}",
            )

    def check_traction_vectors(
        self, output_data, expected_normal_traction, expected_shear_traction
    ):
        expected_traction_vector = [expected_normal_traction, expected_shear_traction]
        self._check_output_vectors(
            output_data,
            expected_traction_vector,
            output_item_label="CAUCHY_STRESS_VECTOR",
            name_of_output_item="traction",
            name_of_normal_component="normal traction",
            name_of_shear_component="shear_traction",
        )

    def check_relative_displacement_vectors(
        self,
        output_data,
        expected_relative_normal_displacement,
        expected_relative_shear_displacement,
    ):
        expected_relative_displacement_vector = [
            expected_relative_normal_displacement,
            expected_relative_shear_displacement,
        ]
        self._check_output_vectors(
            output_data,
            expected_relative_displacement_vector,
            output_item_label="STRAIN",
            name_of_output_item="relative displacement",
            name_of_normal_component="relative normal displacement",
            name_of_shear_component="relative shear displacement",
        )

    def test_interface_pre_stress(self):
        number_of_stages = 2
        run_multiple_stages.run_stages(self.test_path, number_of_stages)

        end_time = 2.0
        reader = GiDOutputFileReader()
        output_data_of_stage_2 = reader.read_output_from(
            os.path.join(
                self.test_path,
                "gid_output",
                "interface_prestress_test_Stage_2.post.res",
            )
        )
        self.check_vertical_stage_displacements(
            output_data_of_stage_2,
            time=end_time,
            expected_vertical_stage_displacement=0.0,
        )

        interface_output_data_of_stage_2 = self.read_interface_output_of_stage_2()
        self.check_output_times(
            interface_output_data_of_stage_2, expected_time=end_time
        )
        self.check_traction_vectors(
            interface_output_data_of_stage_2,
            expected_normal_traction=-1000.0,
            expected_shear_traction=0.0,
        )
        self.check_relative_displacement_vectors(
            interface_output_data_of_stage_2,
            expected_relative_normal_displacement=0.0,
            expected_relative_shear_displacement=0.0,
        )


if __name__ == "__main__":
    KratosUnittest.main()
