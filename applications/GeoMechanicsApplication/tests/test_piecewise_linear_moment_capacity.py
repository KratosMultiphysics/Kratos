import os
import json

import KratosMultiphysics.KratosUnittest as KratosUnittest
from KratosMultiphysics.GeoMechanicsApplication.gid_output_file_reader import GiDOutputFileReader
import test_helper


def run_piecewise_case(project_subpath, times, expected_moments, expected_displacements, unreload_modulus=None, tol=1e-4):
    project_path = test_helper.get_file_path(os.path.join('piecewise_linear_moment_capacity', project_subpath))
    project_parameters_file = test_helper.get_file_path(
        os.path.join('piecewise_linear_moment_capacity', project_subpath, 'ProjectParameters.json'))
    materials_file = os.path.join(project_path, 'MaterialParameters.json')

    cwd = os.getcwd()
    try:
        os.chdir(project_path)
        
        # If unreload_modulus is specified, load and modify MaterialParameters
        materials_backup = None
        if unreload_modulus is not None:
            # Read original materials
            with open(materials_file, 'r') as f:
                materials_data = json.load(f)
            materials_backup = json.dumps(materials_data)
            
            # Add UNRELOAD_MODULUS to all properties
            for prop in materials_data.get("properties", []):
                variables = prop["Material"]["Variables"]
                variables["UNRELOAD_MODULUS"] = unreload_modulus
            
            # Write modified materials back to file
            with open(materials_file, 'w') as f:
                json.dump(materials_data, f, indent=4)
        
        try:
            with open(project_parameters_file, 'r') as parameter_file:
                parameters = test_helper.Kratos.Parameters(parameter_file.read())
            simulation = test_helper.analysis.GeoMechanicsAnalysis(test_helper.Kratos.Model(), parameters)
            simulation.Run()
        finally:
            # Restore original materials file if it was modified
            if materials_backup is not None:
                with open(materials_file, 'w') as f:
                    f.write(materials_backup)
    finally:
        os.chdir(cwd)

    output_file = os.path.join(project_path, 'piecewise_linear_moment_capacity.post.res')
    reader = GiDOutputFileReader()
    output_data = reader.read_output_from(output_file)

    for t, exp_m, exp_d in zip(times, expected_moments, expected_displacements):
        moment = GiDOutputFileReader.element_integration_point_values_at_time("BENDING_MOMENT", t, output_data, [1], [0])[0][0]
        displacement = GiDOutputFileReader.nodal_values_at_time("DISPLACEMENT", t, output_data, [2])[0]
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

    def test_piecewise_linear_moment_capacity(self):
        """
        using only backbone moment.
        """
        times = [0.05, 0.1, 0.15, 0.2, 0.25, 0.30, 0.35, 0.5]
        expected_moments_y = [0.39487, 0.260425, 0.494529, -0.327648, -0.641928, -0.494529, -0.82985, -1.1]
        expected_displacements_y = [0.006, 0.004, 0.0075,-0.005, -0.0125, -0.0075, -0.02, -0.0575]
        run_piecewise_case('.', times, expected_moments_y, expected_displacements_y, tol=1e-4)

    def test_piecewise_linear_moment_capacity_with_unreload(self):
        """
        move down up with unload/reload stiffness (UNRELOAD_MODULUS).
        Tests elastic unload/reload behavior during cyclic loading.
        Uses same case files as backbone test but adds UNRELOAD_MODULUS in-flight.
        """
        times = [0.05, 0.1, 0.15, 0.2, 0.25, 0.30, 0.35, 0.5]
        expected_moments_y = [0.39487, 0.128681, 0.494529, -0.641293, -0.830108, -0.164635, -0.95597, -1.1]
        expected_displacements_y = [0.006, 0.004, 0.0075, -0.005, -0.0125, -0.0075, -0.02, -0.0575]
        run_piecewise_case('.', times, expected_moments_y, expected_displacements_y, unreload_modulus=50.0, tol=1e-4)


if __name__ == '__main__':
    KratosUnittest.main()
