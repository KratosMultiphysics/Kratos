import os
import json
import sys
sys.path.append(r"D:\software_development\Kratos4\bin\Release")
sys.path.append(r"D:\software_development\Kratos4\bin\Release\libs")

import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as KratosUnittest
import KratosMultiphysics.GeoMechanicsApplication as KratosGeo
import test_helper
from strip_load_semi_analytical_solution import StripLoad


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
        test_name = 'test_1d_wave_prop_drained_soil_linear_elastic_solver_initial_acceleration'

        simulation = self.run_wave_through_drained_linear_elastic_soil_test(test_name)

        # check if the correct solver is used
        self.assertTrue(isinstance(simulation._GetSolver().solver,
                                   KratosGeo.GeoMechanicNewtonRaphsonStrategyLinearElasticDynamic))

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
        # test_name = 'test_1d_wave_prop_drained_soil_linear_elastic_solver'
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

    #
    # def test_load_on_block_2d_no_damping(self):
    #     """
    #     Tests a load on a 2d block without damping and a constant mass and stiffness matrix.
    #
    #     """
    #     test_name = 'test_load_on_block_2d_no_damping'
    #     file_path = test_helper.get_file_path(os.path.join('.', test_name))
    #
    #     # run test
    #     simulation = test_helper.run_kratos(file_path)
    #     self.assertTrue(isinstance(simulation._GetSolver().builder_and_solver,
    #                                KratosGeo.ResidualBasedBlockBuilderAndSolverWithMassAndDamping))
    #
    #     # get calculated results
    #     with open(os.path.join(file_path, "calculated_results.json")) as fp:
    #         calculated_result = json.load(fp)
    #
    #     # get expected results
    #     with open(os.path.join(file_path, "expected_results.json")) as fp:
    #         expected_result = json.load(fp)
    #
    #     # check if results are as expected
    #     nodes = ["NODE_3", "NODE_4"]
    #     what = "DISPLACEMENT_Y"
    #     for node in nodes:
    #         self.assertVectorAlmostEqual(calculated_result[node][what], expected_result[node][what])



    # def test_constant_strip_load_2d(self):
    #     """
    #     Tests a constant strip load on a 10m X 10m  block. The block is loaded with a constant strip load of 1kN/m and
    #     a length of 1m.
    #
    #     The solution should roughly follow the semi-analytical solution presented in the following publication:
    #     An Introduction to Soil Dynamics , Verruijt A., 2009, Delft University of Technology, Chapter 12.2
    #
    #
    #     """
    #     # set to true if the results should be checked with the analytical solutio
    #     CHECK_RESULTS = True
    #
    #     test_name = 'test_constant_strip_load_2d'
    #     file_path = test_helper.get_file_path(os.path.join('.', test_name))
    #
    #     simulation = test_helper.run_kratos(file_path)
    #
    #     if CHECK_RESULTS:
    #
    #         # get json output
    #         json_process = simulation._GetListOfProcesses()[4]
    #
    #         elements = json_process.sub_model_part.Elements
    #
    #         calculated_vert_stresses = []
    #         analytical_vert_stresses = []
    #
    #         porosity = 0.0
    #         density_solid = 1020
    #
    #         line_load_length = 1
    #         load_value = -1000
    #         end_time = 1.0
    #
    #         analytical_solution = StripLoad(2.55e3 * 36, 0.25, (1 - porosity) * density_solid, load_value)
    #
    #         x_coords = []
    #         for element in elements:
    #             # get centroid
    #             centroid = element.GetGeometry().Center()
    #             x_coord = centroid.X
    #             y_coord = centroid.Y
    #             depth_check = 10 - y_coord
    #             x_coords.append(centroid.X)
    #
    #             vert_stress_analytic = analytical_solution.calculate_vertical_stress(x_coord, depth_check, end_time,
    #                                                                                  line_load_length, load_value)
    #
    #             analytical_vert_stresses.append(-vert_stress_analytic)
    #
    #             # get element mean vertical cauchy stress
    #             element_stress = element.CalculateOnIntegrationPoints(Kratos.CAUCHY_STRESS_VECTOR,
    #                                                                   json_process.sub_model_part.ProcessInfo)
    #             vert_stress = test_helper.compute_mean_list([gauss_stress[1] for gauss_stress in element_stress])
    #
    #             calculated_vert_stresses.append(vert_stress)
    #
    #         # visualise results if matplotlib is installed
    #         try:
    #             import matplotlib.pyplot as plt
    #             plt.plot(x_coords, analytical_vert_stresses, 'o', color="r")
    #             plt.plot(x_coords, calculated_vert_stresses, 'o', color="b")
    #
    #             plt.legend(["Analytical", "Calculated"])
    #
    #             plt.show()
    #
    #         except ImportError:
    #             print("Matplotlib not installed. Cannot visualise results")


if __name__ == '__main__':
    KratosUnittest.main()
