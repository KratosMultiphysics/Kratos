from __future__ import print_function, absolute_import, division #makes KratosMultiphysics backward compatible with python 2.6 and 2.7
#import kratos core and applications
from KratosMultiphysics import *
from KratosMultiphysics.StructuralMechanicsApplication import *
import KratosMultiphysics.KratosUnittest as KratosUnittest
import structural_mechanics_analysis
import KratosMultiphysics.kratos_utilities as kratos_utilities

def solve_linear_primal_problem():
    with open("./adjoint_sensitivity_analysis_tests/adjoint_truss_stucture_3d2n/linear_truss_test_parameters.json",'r') as parameter_file:
        ProjectParametersPrimal = Parameters( parameter_file.read())

    # To avoid many prints
    if (ProjectParametersPrimal["problem_data"]["echo_level"].GetInt() == 0):
        Logger.GetDefaultOutput().SetSeverity(Logger.Severity.WARNING)

    model_primal = Model()
    primal_analysis = structural_mechanics_analysis.StructuralMechanicsAnalysis(model_primal, ProjectParametersPrimal)
    primal_analysis.Run()

def solve_nonlinear_primal_problem():
    with open("./adjoint_sensitivity_analysis_tests/adjoint_truss_stucture_3d2n/nonlinear_truss_test_parameters.json",'r') as parameter_file:
        ProjectParametersPrimal = Parameters( parameter_file.read())

    # To avoid many prints
    if (ProjectParametersPrimal["problem_data"]["echo_level"].GetInt() == 0):
        Logger.GetDefaultOutput().SetSeverity(Logger.Severity.WARNING)

    model_primal = Model()
    primal_analysis = structural_mechanics_analysis.StructuralMechanicsAnalysis(model_primal, ProjectParametersPrimal)
    primal_analysis.Run()

class TestAdjointSensitivityAnalysisLinearTrussStructure(KratosUnittest.TestCase):

    # called only once for this class, opposed of setUp()
    @classmethod
    def setUpClass(cls):
        solve_linear_primal_problem()

    def test_local_stress_response(self):
        #Create the adjoint solver
        with open("./adjoint_sensitivity_analysis_tests/adjoint_truss_stucture_3d2n/linear_truss_test_local_stress_adjoint_parameters.json",'r') as parameter_file:
            ProjectParametersAdjoint = Parameters( parameter_file.read())

        model_part_name = ProjectParametersAdjoint["problem_data"]["model_part_name"].GetString()
        model_adjoint = Model()

        adjoint_analysis = structural_mechanics_analysis.StructuralMechanicsAnalysis(model_adjoint, ProjectParametersAdjoint)
        adjoint_analysis.Run()

    # called only once for this class, opposed of tearDown()
    @classmethod
    def tearDownClass(cls):
        kratos_utilities.DeleteFileIfExisting("./adjoint_sensitivity_analysis_tests/adjoint_truss_stucture_3d2n/linear_truss_structure.time")
        for file_name in os.listdir():
            if file_name.endswith(".h5"):
                kratos_utilities.DeleteFileIfExisting(file_name)

class TestAdjointSensitivityAnalysisTrussStructure(KratosUnittest.TestCase):

    # called only once for this class, opposed of setUp()
    @classmethod
    def setUpClass(cls):
        solve_nonlinear_primal_problem()

    def test_local_stress_response(self):
        #Create the adjoint solver
        with open("./adjoint_sensitivity_analysis_tests/adjoint_truss_stucture_3d2n/nonlinear_truss_test_local_stress_adjoint_parameters.json",'r') as parameter_file:
            ProjectParametersAdjoint = Parameters( parameter_file.read())

        model_part_name = ProjectParametersAdjoint["problem_data"]["model_part_name"].GetString()
        model_adjoint = Model()

        adjoint_analysis = structural_mechanics_analysis.StructuralMechanicsAnalysis(model_adjoint, ProjectParametersAdjoint)
        adjoint_analysis.Run()

        # Check sensitivities for the parameter I22
        reference_values = [2.018619553, -1.78560606]
        sensitivities_to_check = []
        element_list = [1,2]
        for element_id in element_list:
            sensitivities_to_check.append(model_adjoint.GetModelPart(model_part_name).Elements[element_id].GetValue(CROSS_AREA_SENSITIVITY))

        self.assertAlmostEqual(sensitivities_to_check[0], reference_values[0], 4)
        self.assertAlmostEqual(sensitivities_to_check[1], reference_values[1], 4)

    # called only once for this class, opposed of tearDown()
    @classmethod
    def tearDownClass(cls):
        kratos_utilities.DeleteFileIfExisting("./adjoint_sensitivity_analysis_tests/adjoint_truss_stucture_3d2n/nonlinear_truss_structure.time")
        for file_name in os.listdir():
            if file_name.endswith(".h5"):
                kratos_utilities.DeleteFileIfExisting(file_name)

if __name__ == '__main__':
    KratosUnittest.main()
