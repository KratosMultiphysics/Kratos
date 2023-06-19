import os
import json

import KratosMultiphysics.KratosUnittest as KratosUnittest
import test_helper


class KratosGeoMechanicsDynamicsTests(KratosUnittest.TestCase):
    """
    This class contains benchmark tests which are checked with other FE software
    """

    def setUp(self):
        # Code here will be placed BEFORE every test in this TestCase.
        pass

    def tearDown(self):
        # Code here will be placed AFTER every test in this TestCase.
        pass

    def test_wave_through_drained_linear_elastic_soil(self):
        """
        Test dynamic calculation on a drained linear elastic soil column. a line load of -1kN is instantly placed
        on the soil column. The soil parameters are chosen such that after 0.002 seconds, the wave is reflected at the
        bottom of the geometry such that half the stress in the soil column is cancelled out.

        Note that for an accurate results, the timestep size has to be decreased. For regression test purposes, the
        time step size is increase for faster calculation
        :return:
        """
        test_name = 'test_1d_wave_prop_drained_soil.gid'
        file_path = test_helper.get_file_path(os.path.join('.', test_name))

        test_helper.run_kratos(file_path)

        with open(os.path.join(file_path, "calculated_result.json")) as fp:
            calculated_result = json.load(fp)

        with open(os.path.join(file_path, "expected_result.json")) as fp:
            expected_result = json.load(fp)

        where = "NODE_41"
        what = "VELOCITY_Y"
        self.assertVectorAlmostEqual(calculated_result[where][what], expected_result[where][what])


    def test_wave_through_drained_linear_elastic_soil_constant_mass_damping(self):
        """
        Test dynamic calculation on a drained linear elastic soil column. a line load of -1kN is instantly placed
        on the soil column. The soil parameters are chosen such that after 0.002 seconds, the wave is reflected at the
        bottom of the geometry such that half the stress in the soil column is cancelled out. In this test the global
        mass and damping matrix are precalculated in the builder and solver, such that they are not recalculated every
        step.

        Note that for an accurate results, the timestep size has to be decreased. For regression test purposes, the
        time step size is increase for faster calculation
        :return:
        """
        test_name = 'test_1d_wave_prop_drained_soil_constant_mass_damping.gid'
        file_path = test_helper.get_file_path(os.path.join('.', test_name))

        test_helper.run_kratos(file_path)

        with open(os.path.join(file_path, "calculated_result.json")) as fp:
            calculated_result = json.load(fp)

        with open(os.path.join(file_path, "expected_result.json")) as fp:
            expected_result = json.load(fp)

        where = "NODE_41"
        what = "VELOCITY_Y"
        self.assertVectorAlmostEqual(calculated_result[where][what], expected_result[where][what])

    @KratosUnittest.skip("unit test skipped as it is not ready")
    def test_wave_through_undrained_linear_elastic_soil(self):
        test_name = 'test_1d_confined_undrained_wave.gid'
        file_path = test_helper.get_file_path(os.path.join('.', test_name))
        simulation = test_helper.run_kratos(file_path)
        pass

if __name__ == '__main__':
    KratosUnittest.main()
