import os
import json
from pathlib import Path

import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics.GeoMechanicsApplication.gid_output_file_reader import GiDOutputFileReader
import KratosMultiphysics.GeoMechanicsApplication.run_multiple_stages as run_multiple_stages
import test_helper

class KratosGeoMechanicsInterfacePreStressTests(KratosUnittest.TestCase):
    def check_output_times(self, output_data):
        times = output_data['TIME']
        self.assertEqual(len(times), 1)
        self.assertAlmostEqual(times[0], 2.0)

    def check_traction_vectors(self, output_data, expected_normal_traction, expected_shear_traction):
        element_label = 'ELEMENT_3'
        tractions_at_integration_points = output_data[element_label]['CAUCHY_STRESS_VECTOR']
        integration_point_size_of_interface_element = 3
        for integration_point_label in [str(i) for i in range(integration_point_size_of_interface_element)]:
            output_list = tractions_at_integration_points[integration_point_label]
            self.assertEqual(len(output_list), 1, msg=f"number of traction vectors per integration point of {element_label}")
            self.assertAlmostEqual(output_list[0][0], expected_normal_traction, msg=f'normal traction at integration point {integration_point_label} of {element_label}')
            self.assertAlmostEqual(output_list[0][1], expected_shear_traction, msg=f'shear traction at integration point {integration_point_label} of {element_label}')

    def test_interface_pre_stress(self):
        test_name    = 'test_interface_prestress'
        project_path = test_helper.get_file_path(test_name)
        n_stages     = 2
        run_multiple_stages.run_stages(project_path, n_stages)

        reader = GiDOutputFileReader()
        output_data_stage_2 = reader.read_output_from(os.path.join(project_path, 'gid_output', f'interface_prestress_test_Stage_2.post.res'))

        nodal_displacements = reader.nodal_values_at_time('DISPLACEMENT', 2, output_data_stage_2)

        for nodal_displacement in nodal_displacements:
            y_displacement = nodal_displacement[1]
            self.assertAlmostEqual(y_displacement, 0, places=6)

        interface_output_file_path = Path(project_path) / "stage_2_interface_output.json"
        with open(interface_output_file_path, 'r') as output_file:
            interface_output_data = json.load(output_file)

        self.check_output_times(interface_output_data)
        self.check_traction_vectors(interface_output_data, expected_normal_traction=-1000.0, expected_shear_traction=0.0)



if __name__ == '__main__':
    KratosUnittest.main()
