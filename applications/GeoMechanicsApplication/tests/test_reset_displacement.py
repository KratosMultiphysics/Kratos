import sys
import os

sys.path.append(os.path.join('..', '..', '..'))
sys.path.append(os.path.join('..', 'python_scripts'))

import KratosMultiphysics.KratosUnittest as KratosUnittest

import test_helper

class KratosGeoMechanicsResetDisplacementTests(KratosUnittest.TestCase):
    """
    This class contains benchmark tests which are checked with the analytical solution
    """

    def setUp(self):
        # Code here will be placed BEFORE every test in this TestCase.
        pass

    def tearDown(self):
        # Code here will be placed AFTER every test in this TestCase.
        pass

    def test_reset_displacement_truss(self):
        """
        Tests reset displacement in a truss in 4 stages
        stage 1: load is applied / reset displacement is false
        stage 2: load is applied / reset displacement is true
        stage 3: load is applied / reset displacement is false
        stage 4: load is removed / reset displacement is false

        :return:
        """

        # calculate strain
        F = -1e10    # [N]
        E = 2069e8  # [N/m2]
        A = 1       # [m2]

        eps = F/(E*A)

        # get stages
        test_name = 'geo_truss_with_reset_displacemnet'
        project_path = test_helper.get_file_path(os.path.join('.', test_name + '.gid'))
        n_stages = 4

        stages = test_helper.get_stages(project_path, n_stages)

        displacement_stages = [None]*n_stages
        nodal_coordinates_stages = [None]*n_stages

        # run stages and get results
        for idx, stage in enumerate(stages):
            stage.Run()
            displacement_stages[idx] = test_helper.get_displacement(stage)
            nodal_coordinates_stages[idx] = test_helper. get_nodal_coordinates(stage)

        # Assert
        stage_nr = 0
        for idx,node in enumerate(nodal_coordinates_stages[stage_nr]):
            self.assertAlmostEqual(displacement_stages[stage_nr][idx][0], eps*node[0])

        stage_nr = 1
        for idx,node in enumerate(nodal_coordinates_stages[stage_nr]):
            self.assertAlmostEqual(displacement_stages[stage_nr][idx][0], 0)

        stage_nr = 2
        for idx,node in enumerate(nodal_coordinates_stages[stage_nr]):
            self.assertAlmostEqual(displacement_stages[stage_nr][idx][0], 0)

        stage_nr = 3
        for idx,node in enumerate(nodal_coordinates_stages[stage_nr]):
            self.assertAlmostEqual(displacement_stages[stage_nr][idx][0], -eps*node[0])


    def test_reset_displacement_beam(self):
        """
        Tests reset displacement in a beam in 4 stages
        stage 1: load is applied / reset displacement is true
        stage 2: load is applied / reset displacement is true
        stage 3: load is applied / reset displacement is false
        stage 4: load is removed / reset displacement is false

        :return:
        """

        # calculate strain
        F = -1e5    # [N]
        E = 2069e8  # [N/m2]
        I = 1       # [m4]
        L = 1       # [m]

        eps = (F*L**3)/(3*E*I)

        # get stages
        test_name = 'geo_beam_with_reset_displacemnet'
        project_path = test_helper.get_file_path(os.path.join('.', test_name + '.gid'))
        n_stages = 4

        stages = test_helper.get_stages(project_path, n_stages)

        displacement_stages = [None] * n_stages
        nodal_coordinates_stages = [None] * n_stages

        # run stages and get results
        for idx, stage in enumerate(stages):
            stage.Run()
            displacement_stages[idx] = test_helper.get_displacement(stage)
            nodal_coordinates_stages[idx] = test_helper.get_nodal_coordinates(stage)

        # Assert
        stage_nr = 0
        for idx, node in enumerate(nodal_coordinates_stages[stage_nr]):
            self.assertAlmostEqual(displacement_stages[stage_nr][idx][0], eps * node[0], places=5)

        stage_nr = 1
        for idx, node in enumerate(nodal_coordinates_stages[stage_nr]):
            self.assertAlmostEqual(displacement_stages[stage_nr][idx][0], 0, places=5)

        stage_nr = 2
        for idx, node in enumerate(nodal_coordinates_stages[stage_nr]):
            self.assertAlmostEqual(displacement_stages[stage_nr][idx][0], 0, places=5)

        stage_nr = 3
        for idx, node in enumerate(nodal_coordinates_stages[stage_nr]):
            self.assertAlmostEqual(displacement_stages[stage_nr][idx][0], -eps * node[0], places=5)

if __name__ == '__main__':
    suites = KratosUnittest.KratosSuites
    smallSuite = suites['small'] # These tests are executed by the continuous integration tool
    smallSuite.addTests(KratosUnittest.TestLoader().loadTestsFromTestCases([KratosGeoMechanicsResetDisplacementTests]))
    allSuite = suites['all']
    allSuite.addTests(smallSuite)
    KratosUnittest.runTests(suites)
