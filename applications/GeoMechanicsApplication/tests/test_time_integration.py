import os

import KratosMultiphysics.KratosUnittest as KratosUnittest
import test_helper


class KratosGeoMechanicsTimeIntegrationTests(KratosUnittest.TestCase):
    """
    This class contains regression tests
    """

    def setUp(self):
        # Code here will be placed BEFORE every test in this TestCase.
        pass

    def tearDown(self):
        # Code here will be placed AFTER every test in this TestCase.
        pass

    def test_avoid_small_last_step(self):
        """
        Test dynamic calculation that without countermeasure does a very small step to reach the end_time.
        :return:
        """
        test_name = 'test_time_integration'
        file_path = test_helper.get_file_path(os.path.join('.', test_name))

        test_helper.run_kratos(file_path)

        reader = test_helper.GiDOutputFileReader()
        simulation_output = reader.read_output_from(os.path.join(file_path, "test_time_integration.post.res"))

        number_of_steps_taken = len(simulation_output["results"]["DISPLACEMENT"])

        self.assertEqual(number_of_steps_taken, 26)

if __name__ == '__main__':
    KratosUnittest.main()
