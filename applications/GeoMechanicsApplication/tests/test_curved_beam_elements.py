import sys
import os

import KratosMultiphysics.KratosUnittest as KratosUnittest
import test_helper

class KratosGeoMechanicsCurvedBeamElementTests(KratosUnittest.TestCase):
    """
    This class contains benchmark tests which are checked with the analytical solution
    """

    def setUp(self):
        # Code here will be placed BEFORE every test in this TestCase.
        pass

    def tearDown(self):
        # Code here will be placed AFTER every test in this TestCase.
        pass

    def test_inclined_beam_2D3N(self):
        test_name = 'test_inclined_beam_2D3N'
        file_path = test_helper.get_file_path(os.path.join('.', test_name + '.gid'))
        simulation = test_helper.run_kratos(file_path)

        displacements = test_helper.get_displacement(simulation)
        x_displacements = [displacement[0] for displacement in displacements]

        node_number = 15
        expected_value = 32.337500
        self.assertAlmostEqual(expected_value, x_displacements[node_number-1], 6)


if __name__ == '__main__':
    KratosUnittest.main()
