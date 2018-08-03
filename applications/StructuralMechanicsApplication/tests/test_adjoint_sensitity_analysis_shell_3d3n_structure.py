from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
#import kratos core and applications
from KratosMultiphysics import *
from KratosMultiphysics.StructuralMechanicsApplication import *
import KratosMultiphysics.KratosUnittest as KratosUnittest
import structural_mechanics_analysis
import KratosMultiphysics.kratos_utilities as kratos_utilities

def solve_primal_problem():
    with open("./adjoint_sensitivity_analysis_tests/adjoint_shell_structure_3d3n/linear_shell_test_parameters.json",'r') as parameter_file:
        ProjectParametersPrimal = Parameters( parameter_file.read())
    model_primal = Model()
    primal_analysis = structural_mechanics_analysis.StructuralMechanicsAnalysis(model_primal, ProjectParametersPrimal)
    primal_analysis.Run()

class TestAdjointSensitivityAnalysisShell3D3NStructure(KratosUnittest.TestCase):

    # called only once for this class, opposed of setUp()
    @classmethod
    def setUpClass(cls):
        solve_primal_problem()

    def test_local_stress_response(self):
        # Create the adjoint solver
        with open("./adjoint_sensitivity_analysis_tests/adjoint_shell_structure_3d3n/linear_shell_test_local_stress_adjoint_parameters.json",'r') as parameter_file:
            ProjectParametersAdjoint = Parameters( parameter_file.read())

        model_part_name = ProjectParametersAdjoint["problem_data"]["model_part_name"].GetString()
        model_adjoint = Model()

        adjoint_analysis = structural_mechanics_analysis.StructuralMechanicsAnalysis(model_adjoint, ProjectParametersAdjoint)
        adjoint_analysis.Run()

        # Check sensitivities for the parameter THICKNESS
        reference_values = [1.7135092490964121, -6.860092387341681, 0.14749301178647778]
        sensitivities_to_check = []
        element_list = [1,2,8]
        for element_id in element_list:
            sensitivities_to_check.append(model_adjoint.GetModelPart(model_part_name).Elements[element_id].GetValue(THICKNESS_SENSITIVITY))

        self.assertAlmostEqual(sensitivities_to_check[0], reference_values[0], 5)
        self.assertAlmostEqual(sensitivities_to_check[1], reference_values[1], 5)
        self.assertAlmostEqual(sensitivities_to_check[2], reference_values[2], 5)

    def test_nodal_displacement_response(self):
        # Create the adjoint solver
        with open("./adjoint_sensitivity_analysis_tests/adjoint_shell_structure_3d3n/linear_shell_test_nodal_disp_adjoint_parameters.json",'r') as parameter_file:
            ProjectParametersAdjoint = Parameters( parameter_file.read())

        model_part_name = ProjectParametersAdjoint["problem_data"]["model_part_name"].GetString()
        model_adjoint = Model()

        adjoint_analysis = structural_mechanics_analysis.StructuralMechanicsAnalysis(model_adjoint, ProjectParametersAdjoint)
        adjoint_analysis.Run()

        # Check sensitivities for the parameter THICKNESS
        reference_values = [-0.09916013365433643, -0.23348175177098657, -0.04942512089147077]
        sensitivities_to_check = []
        element_list = [1,2,8]
        for element_id in element_list:
            sensitivities_to_check.append(model_adjoint.GetModelPart(model_part_name).Elements[element_id].GetValue(THICKNESS_SENSITIVITY))

        self.assertAlmostEqual(sensitivities_to_check[0], reference_values[0], 5)
        self.assertAlmostEqual(sensitivities_to_check[1], reference_values[1], 5)
        self.assertAlmostEqual(sensitivities_to_check[2], reference_values[2], 5)

    def test_strain_energy_response(self):
        # Create the adjoint solver
        with open("./adjoint_sensitivity_analysis_tests/adjoint_shell_structure_3d3n/linear_shell_test_strain_energy_adjoint_parameters.json",'r') as parameter_file:
            ProjectParametersAdjoint = Parameters( parameter_file.read())

        model_part_name = ProjectParametersAdjoint["problem_data"]["model_part_name"].GetString()
        model_adjoint = Model()

        adjoint_analysis = structural_mechanics_analysis.StructuralMechanicsAnalysis(model_adjoint, ProjectParametersAdjoint)
        adjoint_analysis.Run()

        # Check sensitivities for the parameter THICKNESS
        reference_values = [-0.4958006682716821, -1.1674087588549331, -0.2471256044520311]
        sensitivities_to_check = []
        element_list = [1,2,8]
        for element_id in element_list:
            sensitivities_to_check.append(model_adjoint.GetModelPart(model_part_name).Elements[element_id].GetValue(THICKNESS_SENSITIVITY))

        self.assertAlmostEqual(sensitivities_to_check[0], reference_values[0], 5)
        self.assertAlmostEqual(sensitivities_to_check[1], reference_values[1], 5)
        self.assertAlmostEqual(sensitivities_to_check[2], reference_values[2], 5)

    # called only once for this class, opposed of tearDown()
    @classmethod
    def tearDownClass(cls):
        kratos_utilities.DeleteFileIfExisting("./adjoint_sensitivity_analysis_tests/adjoint_shell_structure_3d3n/rectangular_plate.time")
        for file_name in os.listdir():
            if file_name.endswith(".h5"):
                kratos_utilities.DeleteFileIfExisting(file_name)

if __name__ == '__main__':
    KratosUnittest.main()