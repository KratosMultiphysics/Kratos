from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
import os
os.environ['OMP_NUM_THREADS'] = "1"
#import kratos core and applications
from KratosMultiphysics import *
from KratosMultiphysics.StructuralMechanicsApplication import *
import KratosMultiphysics.KratosUnittest as KratosUnittest
import structural_mechanics_analysis

class TestAdjointSensitivityAnalysisShell3D3NStructure(KratosUnittest.TestCase):

    def setUp(self):
        pass

    def _remove_file(self, file_path):
        if os.path.isfile(file_path):
            os.remove(file_path)

    def _remove_h5_files(self, model_part_name):
        for name in os.listdir():
            if name.find(model_part_name) == 0:
                self._remove_file(name)

    def _solve_primal_problem(self):
        with open("./adjoint_sensitivity_analysis_tests/adjoint_shell_structure_3d3n/linear_shell_test_parameters.json",'r') as parameter_file:
            ProjectParametersPrimal = Parameters( parameter_file.read())
        parameter_file.close()
        primal_analysis = structural_mechanics_analysis.StructuralMechanicsAnalysis(ProjectParametersPrimal)

        primal_analysis.Run()

    def test_local_stress_response(self):
        # Solve primal problem (only here necessary. The other tests corresponding to the same primal problem.)
        self._solve_primal_problem()
        # Create the adjoint solver
        self._solve_primal_problem()
        with open("./adjoint_sensitivity_analysis_tests/adjoint_shell_structure_3d3n/linear_shell_test_local_stress_adjoint_parameters.json",'r') as parameter_file:
            ProjectParametersAdjoint = Parameters( parameter_file.read())
        adjoint_analysis = structural_mechanics_analysis.StructuralMechanicsAnalysis(ProjectParametersAdjoint)
        adjoint_analysis.Run()

        # Check sensitivities for the parameter THICKNESS
        reference_values = [1.7135092490964121, -6.860092387341681, 0.14749301178647778]
        sensitivities_to_check = []
        element_list = [1,2,8]
        for element_id in element_list:
            sensitivities_to_check.append(adjoint_analysis.GetModelPart().Elements[element_id].GetValue(THICKNESS_SENSITIVITY))

        self.assertAlmostEqual(sensitivities_to_check[0], reference_values[0], 5)
        self.assertAlmostEqual(sensitivities_to_check[1], reference_values[1], 5)
        self.assertAlmostEqual(sensitivities_to_check[2], reference_values[2], 5)

    def test_nodal_displacement_response(self):
        # Create the adjoint solver
        with open("./adjoint_sensitivity_analysis_tests/adjoint_shell_structure_3d3n/linear_shell_test_nodal_disp_adjoint_parameters.json",'r') as parameter_file:
            ProjectParametersAdjoint = Parameters( parameter_file.read())
        adjoint_analysis = structural_mechanics_analysis.StructuralMechanicsAnalysis(ProjectParametersAdjoint)
        adjoint_analysis.Run()

        # Check sensitivities for the parameter THICKNESS
        reference_values = [-0.09916013365433643, -0.23348175177098657, -0.04942512089147077]
        sensitivities_to_check = []
        element_list = [1,2,8]
        for element_id in element_list:
            sensitivities_to_check.append(adjoint_analysis.GetModelPart().Elements[element_id].GetValue(THICKNESS_SENSITIVITY))

        self.assertAlmostEqual(sensitivities_to_check[0], reference_values[0], 5)
        self.assertAlmostEqual(sensitivities_to_check[1], reference_values[1], 5)
        self.assertAlmostEqual(sensitivities_to_check[2], reference_values[2], 5)

    def test_strain_energy_response(self):
        # Create the adjoint solver
        with open("./adjoint_sensitivity_analysis_tests/adjoint_shell_structure_3d3n/linear_shell_test_strain_energy_adjoint_parameters.json",'r') as parameter_file:
            ProjectParametersAdjoint = Parameters( parameter_file.read())
        adjoint_analysis = structural_mechanics_analysis.StructuralMechanicsAnalysis(ProjectParametersAdjoint)
        adjoint_analysis.Run()

        # Check sensitivities for the parameter THICKNESS
        reference_values = [-0.4958006682716821, -1.1674087588549331, -0.2471256044520311]
        sensitivities_to_check = []
        element_list = [1,2,8]
        for element_id in element_list:
            sensitivities_to_check.append(adjoint_analysis.GetModelPart().Elements[element_id].GetValue(THICKNESS_SENSITIVITY))

        self.assertAlmostEqual(sensitivities_to_check[0], reference_values[0], 5)
        self.assertAlmostEqual(sensitivities_to_check[1], reference_values[1], 5)
        self.assertAlmostEqual(sensitivities_to_check[2], reference_values[2], 5)

        # Delete *.h5 only after last test case because primal solution is used in each test case
        self._remove_h5_files("Structure")
        self._remove_file("./adjoint_sensitivity_analysis_tests/adjoint_shell_structure_3d3n/rectangular_plate.time")

    def tearDown(self):
        pass
    #TODO: add this test in test_StructualMechanicsApllication.py

if __name__ == '__main__':
    KratosUnittest.main()