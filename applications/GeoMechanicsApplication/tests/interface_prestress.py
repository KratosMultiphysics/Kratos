import os

import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics.GeoMechanicsApplication.gid_output_file_reader import GiDOutputFileReader
import KratosMultiphysics.GeoMechanicsApplication.run_multiple_stages as run_multiple_stages
import test_helper

class KratosGeoMechanicsInterfacePreStressTests(KratosUnittest.TestCase):
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



if __name__ == '__main__':
    KratosUnittest.main()
