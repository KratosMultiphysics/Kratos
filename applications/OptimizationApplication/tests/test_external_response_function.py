import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as kratos_unittest

import KratosMultiphysics.OptimizationApplication as KratosOA
from KratosMultiphysics.OptimizationApplication.optimization_analysis import OptimizationAnalysis

class TestExternalResponseFunction(kratos_unittest.TestCase):
    def test_steepest_descent_analysis(self):
        with kratos_unittest.WorkFolderScope(".", __file__):
            with open("auxiliary_files/external_response_function_optimization_parameters.json", "r") as file_input:
                parameters = Kratos.Parameters(file_input.read())

            model = Kratos.Model()
            analysis = OptimizationAnalysis(model, parameters)
            analysis.Run()

            model_part = model["CustomMP"]
            exp = Kratos.Expression.NodalExpression(model_part)
            Kratos.Expression.VariableExpressionIO.Read(exp, KratosOA.CUSTOM_DESIGN_VARIABLE, False)

            self.assertVectorAlmostEqual(exp.Evaluate(), [2,3,4,5])

if __name__ == "__main__":
    kratos_unittest.main()

