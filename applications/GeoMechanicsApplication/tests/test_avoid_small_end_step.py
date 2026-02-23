import os

import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics.GeoMechanicsApplication.gid_output_file_reader import GiDOutputFileReader
import test_helper


class KratosGeoMechanicsAvoidSmallEndStepTests(KratosUnittest.TestCase):
    def test_avoid_small_end_step(self):
        """
        Test calculation that without countermeasure does a very small step to reach the end_time.
        This is triggered by setting the end-time to 1.0 and the time-step to 0.9999995.
        The countermeasure scales the time step to 1.0, such that only one step is taken.
        :return:
        """
        test_name = 'test_avoid_small_end_step'
        file_path = test_helper.get_file_path(os.path.join('.', test_name))

        test_helper.run_kratos(file_path)

        reader = GiDOutputFileReader()
        simulation_output = reader.read_output_from(os.path.join(file_path, "test_avoid_small_end_step.post.res"))

        number_of_steps_taken = len(simulation_output["results"]["DISPLACEMENT"])
        self.assertEqual(number_of_steps_taken, 1)

        end_time = simulation_output["results"]["DISPLACEMENT"][0]["time"]
        self.assertEqual(end_time, 1.0)

if __name__ == '__main__':
    KratosUnittest.main()
