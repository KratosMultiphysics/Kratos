import os

import KratosMultiphysics.KratosUnittest as KratosUnittest
import test_helper

class KratosGeoMechanicsMovingLoadOnSoil3DTest(KratosUnittest.TestCase):
    """
    This class contains a regression test for a moving load on a 3D soil
    """

    def test_issue_11838(self):
        """
        Regression test for a moving load on the top-middle of a soil in 3D, where dynamic time integration is used.
        No water is added to the soil and the retention law which is used is: SaturatedBelowPhreaticLevel. The issue
        in this specific combination of input was that the density of the soil was set randomly to the saturated or
        unsaturated density, thus resulting in a random mass matrix and results (only in release mode).


        """
        file_path = test_helper.get_file_path(os.path.join('.', 'test_issue_11838'))

        simulation = test_helper.run_kratos(file_path)
        displacements = test_helper.get_displacement(simulation)

        # check if the displacement is consistent
        displacement_node_6 = tuple(displacements[5])
        expected_displacement_node_6 = (0.0, -2.5458817218200053e-07, -4.879563544898403e-08)

        self.assertVectorAlmostEqual(displacement_node_6, expected_displacement_node_6)


if __name__ == '__main__':
    KratosUnittest.main()