
import KratosMultiphysics as Kratos

# Import KratosUnittest
import KratosMultiphysics.KratosUnittest as kratos_unittest
import KratosMultiphysics.OptimizationApplication as KratosOA
from KratosMultiphysics.OptimizationApplication.utilities.component_data import ComponentData
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.responses.response_function import ResponseFunction
from KratosMultiphysics.OptimizationApplication.utilities.union_utilities import SupportedSensitivityFieldVariableTypes

class TestComponentData(kratos_unittest.TestCase):
    class TestResponse(ResponseFunction):
        def Initialize(self) -> None:
            pass
        def Check(self) -> None:
            pass
        def Finalize(self) -> None:
            pass
        def GetDependentPhysicalKratosVariables(self) -> list[SupportedSensitivityFieldVariableTypes]:
            return None
        def GetModelPart(self) -> Kratos.ModelPart:
            return None
        def CalculateValue(self) -> float:
            return 0.0
        def CalculateGradient(self, _: dict[SupportedSensitivityFieldVariableTypes, KratosOA.ContainerExpression.CollectiveExpressions]) -> None:
            return None

    @classmethod
    def setUpClass(cls) -> None:
        cls.optimization_problem = OptimizationProblem()
        cls.response = TestComponentData.TestResponse()
        cls.optimization_problem.AddResponse("test", cls.response)

    def test_SetDataBuffer(self):
        response_data = ComponentData(self.response, self.optimization_problem)
        response_data.SetDataBuffer(3)

        buffered_data = response_data.GetBufferedData()
        unbuffered_data = response_data.GetUnBufferedData()
        self.assertEqual(buffered_data.GetBufferSize(), 3)
        self.assertEqual(unbuffered_data.GetBufferSize(), 1)

        temp_response_data = ComponentData(self.response, self.optimization_problem)
        self.assertEqual(temp_response_data.GetBufferedData(), buffered_data)
        self.assertEqual(temp_response_data.GetUnBufferedData(), unbuffered_data)

if __name__ == "__main__":
    Kratos.Tester.SetVerbosity(Kratos.Tester.Verbosity.PROGRESS)  # TESTS_OUTPUTS
    kratos_unittest.main()