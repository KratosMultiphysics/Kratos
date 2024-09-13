import os
import json

import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.GeoMechanicsApplication as KratosGeo
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
        time step size is increased for faster calculation
        :return:
        """
        test_name = 'test_1d_wave_prop_drained_soil'
        self.run_wave_through_drained_linear_elastic_soil_test(test_name)

    def test_wave_through_drained_linear_elastic_soil_constant_mass_damping(self):
        """
        Test dynamic calculation on a drained linear elastic soil column. a line load of -1kN is instantly placed
        on the soil column. The soil parameters are chosen such that after 0.002 seconds, the wave is reflected at the
        bottom of the geometry such that half the stress in the soil column is cancelled out. In this test the global
        mass and damping matrix are precalculated in the builder and solver, such that they are not recalculated every
        step.

        Note that for an accurate results, the timestep size has to be decreased. For regression test purposes, the
        time step size is increased for faster calculation
        :return:
        """
        test_name = 'test_1d_wave_prop_drained_soil_constant_mass_damping'

        simulation = self.run_wave_through_drained_linear_elastic_soil_test(test_name)

        # check if the correct builder and solver is used
        self.assertTrue(isinstance(simulation._GetSolver().builder_and_solver,
                                   KratosGeo.ResidualBasedBlockBuilderAndSolverWithMassAndDamping))

    def test_wave_through_drained_linear_elastic_soil_linear_elastic_solver(self):
        """
        Test dynamic calculation on a drained linear elastic soil column. a line load of -1kN is instantly placed
        on the soil column. The soil parameters are chosen such that after 0.002 seconds, the wave is reflected at the
        bottom of the geometry such that half the stress in the soil column is cancelled out. In this test the global
        mass and damping matrix are precalculated in the builder and solver, such that they are not recalculated every
        step.

        Note that for an accurate results, the timestep size has to be decreased. For regression test purposes, the
        time step size is increased for faster calculation
        :return:
        """
        test_name = 'test_1d_wave_prop_drained_soil_linear_elastic_solver'

        simulation = self.run_wave_through_drained_linear_elastic_soil_test(test_name)

        # check if the correct solver is used
        self.assertTrue(isinstance(simulation._GetSolver().solver,
                                   KratosGeo.GeoMechanicNewtonRaphsonStrategyLinearElasticDynamic))

    def test_wave_through_drained_linear_elastic_soil_linear_elastic_solver_initial_acceleration(self):
        """
        Test dynamic calculation on a drained linear elastic soil column. a line load of -1kN is instantly placed
        on the soil column. The soil parameters are chosen such that after 0.002 seconds, the wave is reflected at the
        bottom of the geometry such that half the stress in the soil column is cancelled out. In this a solver is used
        which is designed for linear elastic systems, thus matrices, are not recalculated every step and the linear
        solver is factorized only once. Furthermore, the initial acceleration is calculated such that the first step is
        in equilibrium.

        Note that for an accurate results, the timestep size has to be decreased. For regression test purposes, the
        time step size is increased for faster calculation

        """
        test_name = 'test_1d_wave_prop_drained_soil_linear_elastic_solver_initial_acceleration'

        simulation = self.run_wave_through_drained_linear_elastic_soil_test(test_name)

        # check if the correct solver is used
        self.assertTrue(isinstance(simulation._GetSolver().solver,
                                   KratosGeo.GeoMechanicNewtonRaphsonStrategyLinearElasticDynamic))

    def test_wave_through_drained_linear_elastic_soil_linear_elastic_solver_multi_stage(self):
        """
        Test dynamic calculation on a drained linear elastic soil column. a line load of -1kN is instantly placed
        on the soil column. The soil parameters are chosen such that after 0.002 seconds, the wave is reflected at the
        bottom of the geometry such that half the stress in the soil column is cancelled out. In this a solver is used
        which is designed for linear elastic systems, thus matrices, are not recalculated every step and the linear
        solver is factorized only once. Furthermore, the initial acceleration is calculated such that the first step is
        in equilibrium.

        Note that for an accurate results, the timestep size has to be decreased. For regression test purposes, the
        time step size is increased for faster calculation

        """
        test_name = 'test_1d_wave_prop_drained_soil_linear_elastic_solver_multi_stage'

        file_path = test_helper.get_file_path(os.path.join('.', test_name))
        n_stages = 3

        # run simulation
        test_helper.run_stages(file_path, n_stages)

        where = "NODE_41"
        what = "DISPLACEMENT_Y"
        calculated_displacement = []

        # get calculated results per stage
        for i in range(n_stages):
            with open(os.path.join(file_path, "calculated_result_stage" + str(i+1) + ".json")) as fp:
                calculated_result = json.load(fp)

            calculated_displacement.extend(calculated_result[where][what])

        # get expected results
        with open(os.path.join(file_path, "expected_result.json")) as fp:
            expected_result = json.load(fp)

        self.assertVectorAlmostEqual(calculated_displacement, expected_result[where][what])

    def run_wave_through_drained_linear_elastic_soil_test(self, test_name):
        """
        Test dynamic calculation on a drained linear elastic soil column. a line load of -1kN is instantly placed
        on the soil column. The soil parameters are chosen such that after 0.002 seconds, the wave is reflected at the
        bottom of the geometry such that half the stress in the soil column is cancelled out. In this test the global
        mass and damping matrix are precalculated in the builder and solver, such that they are not recalculated every
        step.

        Note that for an accurate results, the timestep size has to be decreased. For regression test purposes, the
        time step size is increased for faster calculation
        :return:
        """

        file_path = test_helper.get_file_path(os.path.join('.', test_name))

        simulation = test_helper.run_kratos(file_path)

        with open(os.path.join(file_path, "calculated_result.json")) as fp:
            calculated_result = json.load(fp)

        with open(os.path.join(file_path, "expected_result.json")) as fp:
            expected_result = json.load(fp)

        where = "NODE_41"
        what = "VELOCITY_Y"
        self.assertVectorAlmostEqual(calculated_result[where][what], expected_result[where][what])

        return simulation

    def test_load_on_block_2d_no_damping(self):
        """
        Tests a load on a 2d block without damping and a constant mass and stiffness matrix.

        """
        test_name = 'test_load_on_block_2d_no_damping'
        file_path = test_helper.get_file_path(os.path.join('.', test_name))

        # run test
        simulation = test_helper.run_kratos(file_path)
        self.assertTrue(isinstance(simulation._GetSolver().builder_and_solver,
                                   KratosGeo.ResidualBasedBlockBuilderAndSolverWithMassAndDamping))

        # get calculated results
        with open(os.path.join(file_path, "calculated_results.json")) as fp:
            calculated_result = json.load(fp)

        # get expected results
        with open(os.path.join(file_path, "expected_results.json")) as fp:
            expected_result = json.load(fp)

        # check if results are as expected
        nodes = ["NODE_3", "NODE_4"]
        what = "DISPLACEMENT_Y"
        for node in nodes:
            self.assertVectorAlmostEqual(calculated_result[node][what], expected_result[node][what])


if __name__ == '__main__':
    KratosUnittest.main()
