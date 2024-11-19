import os

import KratosMultiphysics.KratosUnittest as KratosUnittest
import test_helper


class TestThreeDimensionalPipingValidation(KratosUnittest.TestCase):
    def test_three_dimensional_piping(self):
        file_path = test_helper.get_file_path(os.path.join('three_dimensional_piping'))
        simulation = test_helper.run_kratos(file_path)
        result_file_name = os.path.join(file_path, 'Cube_10sS_1.post.res')

        reader = test_helper.GiDOutputFileReader()
        actual_data = reader.read_output_from(result_file_name)

        pipe_elements = [24001, 24002, 24003, 24004, 24005, 24006, 24007, 24008, 24009, 24010, 24011, 24012, 24013,
                         24014, 24015, 24016, 24017, 24018, 24019]
        pipe_active_before_critical_head = \
            reader.element_integration_point_values_at_time("PIPE_ACTIVE", 0.5, actual_data, pipe_elements, [1])
        for active, element in zip(pipe_active_before_critical_head, pipe_elements):
            if element > 24004:
                self.assertEqual(active[0], 1)
            else:
                self.assertEqual(active[0], 0)

        pipe_active_after_critical_head = \
            reader.element_integration_point_values_at_time("PIPE_ACTIVE", 1.0, actual_data, pipe_elements, [1])
        for active in pipe_active_after_critical_head:
            self.assertEqual(active[0], 1)

if __name__ == '__main__':
    KratosUnittest.main()
