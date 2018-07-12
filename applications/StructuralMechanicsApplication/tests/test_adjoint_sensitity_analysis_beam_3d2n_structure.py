from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
#import kratos core and applications
from KratosMultiphysics import *
from KratosMultiphysics.StructuralMechanicsApplication import *
import KratosMultiphysics.KratosUnittest as KratosUnittest
import structural_mechanics_analysis

class TestAdjointSensitivityAnalysisBeamStructure(KratosUnittest.TestCase):

    def setUp(self):
        # Solve primal problem (only in one test case necessary)
        self._solve_primal_problem()

    def _remove_file(self, file_path):
        if os.path.isfile(file_path):
            os.remove(file_path)

    def _remove_h5_files(self, model_part_name):
        for name in os.listdir():
            if name.find(model_part_name) == 0:
                self._remove_file(name)

    def _solve_primal_problem(self):
        with open("./adjoint_sensitivity_analysis_tests/adjoint_beam_structure_3d2n/beam_test_parameters.json",'r') as parameter_file:
            ProjectParametersPrimal = Parameters( parameter_file.read())

        model_primal = Model()

        primal_analysis = structural_mechanics_analysis.StructuralMechanicsAnalysis(model_primal, ProjectParametersPrimal)

        primal_analysis.Run()

    def test_local_stress_response(self):
        #Create the adjoint solver
        with open("./adjoint_sensitivity_analysis_tests/adjoint_beam_structure_3d2n/beam_test_local_stress_adjoint_parameters.json",'r') as parameter_file:
            ProjectParametersAdjoint = Parameters( parameter_file.read())

        model_part_name = ProjectParametersAdjoint["problem_data"]["model_part_name"].GetString()
        model_adjoint = Model()

        adjoint_analysis = structural_mechanics_analysis.StructuralMechanicsAnalysis(model_adjoint, ProjectParametersAdjoint)
        adjoint_analysis.Run()

        # Check sensitivities for the parameter I22
        reference_values = [-87622.77099397512, 38125.18144970003, 625.0029074349038]
        sensitivities_to_check = []
        element_list = [1,6,10]
        for element_id in element_list:
            sensitivities_to_check.append(model_adjoint.GetModelPart(model_part_name).Elements[element_id].GetValue(I22_SENSITIVITY))

        self.assertAlmostEqual(sensitivities_to_check[0], reference_values[0], 5)
        self.assertAlmostEqual(sensitivities_to_check[1], reference_values[1], 5)
        self.assertAlmostEqual(sensitivities_to_check[2], reference_values[2], 5)

    def test_nodal_displacement_response(self):
        # Create the adjoint solver
        with open("./adjoint_sensitivity_analysis_tests/adjoint_beam_structure_3d2n/beam_test_nodal_disp_adjoint_parameters.json",'r') as parameter_file:
            ProjectParametersAdjoint = Parameters( parameter_file.read())

        model_part_name = ProjectParametersAdjoint["problem_data"]["model_part_name"].GetString()
        model_adjoint = Model()

        adjoint_analysis = structural_mechanics_analysis.StructuralMechanicsAnalysis(model_adjoint, ProjectParametersAdjoint)

        adjoint_analysis.Run()

        # Check sensitivities for the parameter I22
        reference_values = [-454.1027959305903, -378.2187594016309, -6.200311358415619]
        sensitivities_to_check = []
        element_list = [1,6,10]
        for element_id in element_list:
            sensitivities_to_check.append(model_adjoint.GetModelPart(model_part_name).Elements[element_id].GetValue(I22_SENSITIVITY))

        self.assertAlmostEqual(sensitivities_to_check[0], reference_values[0], 5)
        self.assertAlmostEqual(sensitivities_to_check[1], reference_values[1], 5)
        self.assertAlmostEqual(sensitivities_to_check[2], reference_values[2], 5)

    def test_strain_energy_response(self):
        # Create the adjoint solver
        with open("./adjoint_sensitivity_analysis_tests/adjoint_beam_structure_3d2n/beam_test_strain_energy_adjoint_parameters.json",'r') as parameter_file:
            ProjectParametersAdjoint = Parameters( parameter_file.read())

        model_part_name = ProjectParametersAdjoint["problem_data"]["model_part_name"].GetString()
        model_adjoint = Model()

        adjoint_analysis = structural_mechanics_analysis.StructuralMechanicsAnalysis(model_adjoint, ProjectParametersAdjoint)

        adjoint_analysis.Run()

        # Check sensitivities for the parameter I22
        reference_values = [-9082.055913430777, -7564.375188032661, -124.00616062621255]
        sensitivities_to_check = []
        element_list = [1,6,10]
        for element_id in element_list:
            sensitivities_to_check.append(model_adjoint.GetModelPart(model_part_name).Elements[element_id].GetValue(I22_SENSITIVITY))

        self.assertAlmostEqual(sensitivities_to_check[0], reference_values[0], 5)
        self.assertAlmostEqual(sensitivities_to_check[1], reference_values[1], 5)
        self.assertAlmostEqual(sensitivities_to_check[2], reference_values[2], 5)


        # Delete *.h5 only after last test case because primal solution is used in each test case
        self._remove_h5_files("Structure")
        self._remove_file("./adjoint_sensitivity_analysis_tests/adjoint_beam_structure_3d2n/Beam_structure.time")

    def tearDown(self):
        pass
    #TODO: add this test in test_StructualMechanicsApllication.py

if __name__ == '__main__':
    KratosUnittest.main()
