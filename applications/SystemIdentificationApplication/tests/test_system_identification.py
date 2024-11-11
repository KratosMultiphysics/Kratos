import numpy
import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as UnitTest
from KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_analysis import StructuralMechanicsAnalysis
from KratosMultiphysics.OptimizationApplication.optimization_analysis import OptimizationAnalysis
from KratosMultiphysics.OptimizationApplication.responses.response_routine import ResponseRoutine
from KratosMultiphysics.kratos_utilities import DeleteFileIfExisting
from KratosMultiphysics.compare_two_files_check_process import CompareTwoFilesCheckProcess

class TestSystemIdentification(UnitTest.TestCase):
    def test_DamagedSystem(self):
        self.addCleanup(DeleteFileIfExisting, "auxiliary_files/damaged_problem/measured_data.csv")

        model = Kratos.Model()
        with open("auxiliary_files/damaged_problem/damaged_project_parameters.json", "r") as file_input:
            params = Kratos.Parameters(file_input.read())
        analysis = StructuralMechanicsAnalysis(model, params)
        analysis.Run()

        data = numpy.loadtxt("auxiliary_files/damaged_problem/measured_data.csv", comments="#", usecols=[0,3,4,5,6], delimiter=",")
        ref_data = numpy.loadtxt("auxiliary_files/damaged_problem/measured_data_ref.csv", comments="#", usecols=[0,3,4,5,6], delimiter=",")
        self.assertTrue(numpy.allclose(data, ref_data, 1e-12, 1e-12))

    def test_SystemIdentification(self):
        model = Kratos.Model()
        with open("auxiliary_files/system_identification/optimization_parameters.json", "r") as file_input:
            params = Kratos.Parameters(file_input.read())

        analysis = OptimizationAnalysis(model, params)
        analysis.Run()

        params = Kratos.Parameters("""{
            "reference_file_name"   : "auxiliary_files/system_identification_summary_ref.csv",
            "output_file_name"      : "auxiliary_files/summary.csv",
            "remove_output_file"    : true,
            "comparison_type"       : "csv_file",
            "tolerance"             : 1e-2,
            "relative_tolerance"    : 1e-3,
            "dimension"             : 3
        }""")
        CompareTwoFilesCheckProcess(params).Execute()

    def test_SystemIdentificationPNorm(self):
        model = Kratos.Model()
        with open("auxiliary_files/system_identification_p_norm/optimization_parameters.json", "r") as file_input:
            params = Kratos.Parameters(file_input.read())

        analysis = OptimizationAnalysis(model, params)
        analysis.Run()

        params = Kratos.Parameters("""{
            "reference_file_name"   : "auxiliary_files/system_identification_p_norm_summary_ref.csv",
            "output_file_name"      : "auxiliary_files/summary_p_norm.csv",
            "remove_output_file"    : true,
            "comparison_type"       : "csv_file",
            "tolerance"             : 1e-2,
            "relative_tolerance"    : 1e-3,
            "dimension"             : 3
        }""")
        CompareTwoFilesCheckProcess(params).Execute()

class TestSystemIdentificationMixedElementsFD(UnitTest.TestCase):
    def __ExecuteDamagedSystemTest(self, case_name: str):
        # self.addCleanup(DeleteFileIfExisting, f"auxiliary_files/{case_name}/measured_data.csv")
        self.addCleanup(DeleteFileIfExisting, f"auxiliary_files/{case_name}/model_part.time")

        model = Kratos.Model()
        with open(f"auxiliary_files/{case_name}/damaged_project_parameters.json", "r") as file_input:
            params = Kratos.Parameters(file_input.read())
        analysis = StructuralMechanicsAnalysis(model, params)
        analysis.Run()

        data = numpy.loadtxt(f"auxiliary_files/{case_name}/measured_data.csv", comments="#", usecols=[0,3,4,5,6], delimiter=",")
        ref_data = numpy.loadtxt(f"auxiliary_files/{case_name}/measured_data_ref.csv", comments="#", usecols=[0,3,4,5,6], delimiter=",")
        self.assertTrue(numpy.allclose(data, ref_data, 1e-12, 1e-12))

    def __ExecuteSystemIdentificationFDTest(self, case_name: str):
        self.addCleanup(DeleteFileIfExisting, f"auxiliary_files/{case_name}/model_part.time")
        with open(f"auxiliary_files/{case_name}/optimization_parameters.json", "r") as file_input:
            parameters = Kratos.Parameters(file_input.read())

        model = Kratos.Model()
        analysis = OptimizationAnalysis(model, parameters)

        analysis.Initialize()
        analysis.Check()

        optimization_problem = analysis.optimization_problem
        objective : ResponseRoutine = optimization_problem.GetComponent("damage_response", ResponseRoutine)
        response = objective.GetReponse()
        ref_value = response.CalculateValue()

        var = objective.GetRequiredPhysicalGradients()
        response.CalculateGradient(var)
        for exp in var[Kratos.YOUNG_MODULUS].GetContainerExpressions():
            adjoint_sensitivity = exp.Evaluate()

        structure_mp = model["AdjointStructure"]

        delta = 1e-6
        for index, elem in enumerate(model["AdjointStructure.all"].Elements):
            Kratos.VariableUtils().SetHistoricalVariableToZero(Kratos.DISPLACEMENT, structure_mp.Nodes)
            elem.Properties[Kratos.YOUNG_MODULUS] += delta * 3e10
            sensitivity = (response.CalculateValue() - ref_value) / delta
            self.assertAlmostEqual(sensitivity, adjoint_sensitivity[index] * 3e10, 12)
            elem.Properties[Kratos.YOUNG_MODULUS] -= delta * 3e10
        analysis.Finalize()

    def test_DamagedSystemSolidTruss(self):
        self.__ExecuteDamagedSystemTest("mixed_element_1")

    def test_GradientsSolidTruss(self):
        self.__ExecuteSystemIdentificationFDTest("mixed_element_1")

    def test_DamagedSystemSolidBeam(self):
        self.__ExecuteDamagedSystemTest("mixed_element_2")

    def test_GradientsSolidBeam(self):
        self.__ExecuteSystemIdentificationFDTest("mixed_element_2")

if __name__ == '__main__':
    UnitTest.main()

