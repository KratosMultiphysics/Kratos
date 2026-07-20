import os

import KratosMultiphysics.KratosUnittest as KratosUnittest
import test_helper


class KratosGeoMechanicsMuskatTests(KratosUnittest.TestCase):
    """
    To be detailed out
    """

    def _assert_fully_saturated_flow(
        self, parent_name, test_name, node_number, expected_value
    ):
        file_path = test_helper.get_file_path(os.path.join(parent_name, test_name))
        simulation = test_helper.run_kratos(file_path)
        water_pressure = [
            water_pressure
            for water_pressure in test_helper.get_water_pressure(simulation)
        ]
        self.assertAlmostEqual(expected_value, water_pressure[node_number - 1], 6)

    # To be changed, node 100 is a random pick and may be on the prescribed boundary. pick a better result, with computed result
    def test_fully_saturated_hydrostatic(self):
        self._assert_fully_saturated_flow(
            "Muskat", "fully_saturated_hydrostatic", 100, 22800.0
        )

    def test_fully_saturated_hydrostatic_cutoff(self):
        self._assert_fully_saturated_flow(
            "Muskat", "fully_saturated_hydrostatic_cutoff", 100, 0.0
        )


if __name__ == "__main__":
    KratosUnittest.main()
