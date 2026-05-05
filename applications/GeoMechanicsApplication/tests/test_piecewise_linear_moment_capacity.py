import os

import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics.GeoMechanicsApplication.gid_output_file_reader import GiDOutputFileReader
import test_helper


def run_piecewise_case(project_subpath, times, expected_forces, expected_displacements, tol=1e-4):
    project_path = test_helper.get_file_path(os.path.join('piecewise_linear_moment_capacity', project_subpath))
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

    for t, exp_f, exp_d in zip(times, expected_forces, expected_displacements):
        section_force = GiDOutputFileReader.element_integration_point_values_at_time("FORCE", t, output_data, [1], [0])[0][0]
        force_x = section_force[0]
        displacement = GiDOutputFileReader.nodal_values_at_time("DISPLACEMENT", t, output_data, [3])[0]
        disp_x = displacement[0]
        abs_tol_f = max(1.0e-12, tol * max(1.0, abs(exp_f)))
        abs_tol_d = max(1.0e-12, tol * max(1.0, abs(exp_d)))
        assert abs(force_x - exp_f) <= abs_tol_f, f"Force mismatch at t={t}: {force_x} != {exp_f}"
        assert abs(disp_x - exp_d) <= abs_tol_d, f"Disp mismatch at t={t}: {disp_x} != {exp_d}"


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


    def test_piecewise_linear_moment_capacity_tension(self):
        """
        2 element elongation test for piecewise linear moment capacity material.
        """
        times = [1.0, 2.0, 3.0, 4.0]
        expected_forces_x = [0.81055, -1.08659, 1.20015, 1.5]
        expected_displacements_x = [1.0, 0.75, 1.0, 2.0]
        run_piecewise_case('tension', times, expected_forces_x, expected_displacements_x, tol=1e-4)

    def test_piecewise_linear_moment_capacity_compression(self):
        """
        2 element compression test for piecewise linear moment capacity material.
        """
        times = [1.0, 2.0, 3.0, 4.0]
        expected_forces_x = [-0.81055, 1.08659, -1.20015, -1.5]
        expected_displacements_x = [-1.0, -0.75, -1.0, -2.0]
        run_piecewise_case('compression', times, expected_forces_x, expected_displacements_x, tol=1e-4)

    def test_piecewise_linear_moment_capacity_tension_compression(self):
        """
        Elongation-compression loop test for piecewise linear moment capacity material.
        """
        times = [1.0, 2.0, 3.0, 4.0]
        expected_forces_x = [0.81055, -1.38546, -0.81297, 1.38676]
        expected_displacements_x = [1.0, 0.0, -1.0, 0.0]
        run_piecewise_case('tension_compression', times, expected_forces_x, expected_displacements_x, tol=1e-4)

    def test_piecewise_linear_moment_capacity_compression_tension(self):
        """
        Compression-elongation loop test for piecewise linear moment capacity material.
        """
        times = [1.0, 2.0, 3.0, 4.0]
        expected_forces_x = [-0.81055, 1.38546, 0.81297, -1.38676]
        expected_displacements_x = [-1.0, -0.0, 1.0, 0.0]
        run_piecewise_case('compression_tension', times, expected_forces_x, expected_displacements_x, tol=1e-4)


if __name__ == '__main__':
    KratosUnittest.main()
