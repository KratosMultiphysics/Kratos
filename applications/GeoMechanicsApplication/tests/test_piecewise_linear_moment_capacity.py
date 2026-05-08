import os

import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics.GeoMechanicsApplication.gid_output_file_reader import GiDOutputFileReader
import test_helper


def run_piecewise_case(project_subpath, times, expected_moments, expected_displacements, tol=1e-4):
    project_path = test_helper.get_file_path(os.path.join('piecewise_linear_moment_capacity', project_subpath))
    case_project_parameters_file = test_helper.get_file_path(
        os.path.join('piecewise_linear_moment_capacity', project_subpath, 'ProjectParameters.json'))
    if os.path.exists(case_project_parameters_file):
        project_parameters_file = case_project_parameters_file
    else:
        project_parameters_file = test_helper.get_file_path(
            os.path.join('piecewise_linear_moment_capacity', 'common', 'ProjectParameters.json'))

    cwd = os.getcwd()
    try:
        os.chdir(project_path)
        with open(project_parameters_file, 'r') as parameter_file:
            parameters = test_helper.Kratos.Parameters(parameter_file.read())
        simulation = test_helper.analysis.GeoMechanicsAnalysis(test_helper.Kratos.Model(), parameters)
        simulation.Run()
    finally:
        os.chdir(cwd)

    output_file = os.path.join(project_path, 'piecewise_linear_moment_capacity.post.res')
    reader = GiDOutputFileReader()
    output_data = reader.read_output_from(output_file)

    for t, exp_m, exp_d in zip(times, expected_moments, expected_displacements):
        moment = GiDOutputFileReader.element_integration_point_values_at_time("BENDING_MOMENT", t, output_data, [1], [0])[0][0]
        displacement = GiDOutputFileReader.nodal_values_at_time("DISPLACEMENT", t, output_data, [3])[0]
        disp_y = displacement[1]
        abs_tol_m = max(1.0e-12, tol * max(1.0, abs(exp_m)))
        abs_tol_d = max(1.0e-12, tol * max(1.0, abs(exp_d)))
        assert abs(moment - exp_m) <= abs_tol_m, f"Moment mismatch at t={t}: {moment} != {exp_m}"
        assert abs(disp_y - exp_d) <= abs_tol_d, f"Disp mismatch at t={t}: {disp_y} != {exp_d}"


class KratosGeoMechanicsPiecewiseLinearMomentCapacityTests(KratosUnittest.TestCase):
    """
    This class contains tests for the piecewise linear moment capacity material
    """

    def setUp(self):
        # Code here will be placed BEFORE every test in this TestCase.
        pass

    def tearDown(self):
        # Code here will be placed AFTER every test in this TestCase.
        pass


    def no_test_piecewise_linear_moment_capacity_move_up(self):
        """
        move up of the right beam end  for piecewise linear moment capacity material.
        """
        times = [1.0, 2.0, 3.0, 4.0]
        expected_moments_y = [0.120901, 0.168661, 0.205083, 0.240241]
        expected_displacements_y = [0.025, 0.05, 0.075, 0.1]
        run_piecewise_case('move_up', times, expected_moments_y, expected_displacements_y, tol=1e-4)

    def no_test_piecewise_linear_moment_capacity_move_down(self):
        """
        move down of the right beam end for piecewise linear moment capacity material.
        """
        times = [1.0, 2.0, 3.0, 4.0]
        expected_moments_y = [-0.120901, -0.168661, -0.205083, -0.240241]
        expected_displacements_y = [-0.025, -0.05, -0.075, -0.1]
        run_piecewise_case('move_down', times, expected_moments_y, expected_displacements_y, tol=1e-4)

    def no_test_piecewise_linear_moment_capacity_move_up_down(self):
        """
        move up down of the right beam end for piecewise linear moment capacity material.
        """
        times = [1.0, 2.0, 3.0, 4.0]
        expected_moments_y = [0.120901, 0.0, -0.120901, 0.0]
        expected_displacements_y = [0.025, 0.0, -0.025, 0.0]
        run_piecewise_case('move_up_down', times, expected_moments_y, expected_displacements_y, tol=1e-4)

    def test_piecewise_linear_moment_capacity_move_down_up(self):
        """
        move down up of the right beam end for piecewise linear moment capacity material.
        """
        times = [1.0, 2.0, 3.0, 4.0]
        expected_moments_y = [-0.120901, 0.0, 0.120901, 0.0]
        expected_displacements_y = [-0.025, 0.0, 0.025, 0.0]
        run_piecewise_case('move_down_up', times, expected_moments_y, expected_displacements_y, tol=1e-4)


if __name__ == '__main__':
    KratosUnittest.main()
