import os
import shutil
import KratosMultiphysics.LinearSolversApplication
import KratosMultiphysics.StructuralMechanicsApplication
import KratosMultiphysics.GeoMechanicsApplication as GeoMechanicsApplication
import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics.GeoMechanicsApplication.gid_output_file_reader import (
    GiDOutputFileReader,
)
import test_helper


class KratosGeoMechanicsAvoidSmallEndStepTests(KratosUnittest.TestCase):
    def setUp(self):
        super().setUp()

        self.test_root = test_helper.get_file_path("test_avoid_small_end_step")
        self.test_path = os.path.join(self.test_root, self.get_test_dir_name())

        shutil.rmtree(self.test_path, ignore_errors=True)

        os.makedirs(self.test_path)

        self.project_parameters_filenames = ["ProjectParameters.json"]
        input_filenames = self.project_parameters_filenames + [
            "MaterialParameters.json",
            "blocks.mdpa",
        ]

        for filename in input_filenames:
            shutil.copy(
                os.path.join(self.test_root, filename),
                os.path.join(self.test_path, filename),
            )

    def get_test_dir_name(self):
        raise RuntimeError(
            "This base class does not provide a generic test directory name"
        )

    def check_number_of_steps_and_end_time(self):
        reader = GiDOutputFileReader()
        simulation_output = reader.read_output_from(
            os.path.join(self.test_path, "test_avoid_small_end_step.post.res")
        )

        number_of_steps_taken = len(simulation_output["results"]["DISPLACEMENT"])
        self.assertEqual(number_of_steps_taken, 1)

        end_time = simulation_output["results"]["DISPLACEMENT"][0]["time"]
        self.assertEqual(end_time, 1.0)


class KratosGeoMechanicsAvoidSmallTimeStepPythonRoute(
    KratosGeoMechanicsAvoidSmallEndStepTests
):
    """
    This test class checks the settlement workflow using the python route.
    """

    def get_test_dir_name(self):
        return "python"

    def test_avoid_small_end_step(self):
        """
        Test calculation that without countermeasure does a very small step to reach the end_time.
        This is triggered by setting the end-time to 1.0 and the time-step to 0.9999995.
        The countermeasure scales the time step to 1.0, such that only one step is taken.
        :return:
        """
        test_helper.run_kratos(self.test_path)
        self.check_number_of_steps_and_end_time()



class KratosGeoMechanicsAvoidSmallTimeStepCppRoute(
    KratosGeoMechanicsAvoidSmallEndStepTests
):
    """
    This test class checks the settlement workflow using the C++ route.
    """

    def get_test_dir_name(self):
        return "cpp"

    def test_avoid_small_end_step(self):
        """
        Test calculation that without countermeasure does a very small step to reach the end_time.
        This is triggered by setting the end-time to 1.0 and the time-step to 0.9999995.
        The countermeasure scales the time step to 1.0, such that only one step is taken.
        :return:
        """
        import KratosMultiphysics.GeoMechanicsApplication.run_geo_settlement as run_geo_settlement
        status = run_geo_settlement.run_stages(
            self.test_path, self.project_parameters_filenames
        )
        self.assertEqual(status, 0)
        self.check_number_of_steps_and_end_time()


if __name__ == "__main__":
    KratosUnittest.main()
