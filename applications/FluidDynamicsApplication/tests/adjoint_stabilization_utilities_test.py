from math import exp

import KratosMultiphysics as Kratos
import KratosMultiphysics.KratosUnittest as UnitTest

from KratosMultiphysics.kratos_utilities import DeleteFileIfExisting
from KratosMultiphysics.FluidDynamicsApplication.adjoint_stabilization_utilities import ComputeStabilizationCoefficient

@UnitTest.skipIfApplicationsNotAvailable("StatisticsApplication")
class AdjointStabilizationUtilitiesTests(UnitTest.TestCase):
    def setUp(self):
        self.write_json_output = False

    def testComputeStabilizationCoefficient(self):
        stabilization_parameters = Kratos.Parameters("""{
            "initial_coefficient_bounds"  : [0.0, 1.0],
            "tolerance"                   : 1e-3,
            "plateau_time_range"          : [5.0, 10.0],
            "plateau_max_slope"           : 0.5,
            "max_iterations"              : 20,
            "adjoint_parameters_file_name": "AdjointQSVMSSensitivity2DTest/one_element_test_adjoint_parameters.json"
        }""")

        class DummySolver:
            def __init__(self, model_part):
                self.main_model_part = model_part

        class DummyAnalysisClass:
            def __init__(self, model, parameters):
                self.model = model
                self.model_part = self.model.CreateModelPart("test")
                self.model_part.AddNodalSolutionStepVariable(Kratos.SHAPE_SENSITIVITY)
                self.model_part.CreateNewNode(1, 0, 0, 0)
                self.model_part.CreateNewNode(2, 0, 1, 0)
                self.model_part.CreateNewNode(3, 0, 0, 1)
                self.model_part.CreateNewNode(4, 0, 1, 1)
                self.model_part.CreateNewNode(5, 1, 1, 1)
                self.__solver = DummySolver(self.model_part)
                self.__stabilization_coefficient = parameters["problem_data"]["echo_level"].GetDouble()
                self.time = 1.0

            def OutputSolutionStep(self):
                for node in self.model_part.Nodes:
                    current_value = node.GetSolutionStepValue(Kratos.SHAPE_SENSITIVITY)
                    increment_rate = exp(-self.time * 2.0 * self.__stabilization_coefficient)
                    increment = Kratos.Array3([(node.Id + 1.0) * increment_rate, (node.Id * 2.0) * increment_rate, (node.Id / 3.0) * increment_rate])
                    node.SetSolutionStepValue(Kratos.SHAPE_SENSITIVITY,  current_value + increment)
                self.time += 1.0

            def Finalize(self):
                pass

            def _GetSolver(self):
                return self.__solver

        value = ComputeStabilizationCoefficient(DummyAnalysisClass, stabilization_parameters, execution_method=AdjointStabilizationUtilitiesTests.__ExecuteAnalysis)
        self.assertAlmostEqual(value, 6.860352e-01, 7)
        DeleteFileIfExisting("adjoint_stabilization_data.dat")

    @staticmethod
    def __ExecuteAnalysis(analysis_class_type, adjoint_parameters, stabilization_coefficient, solve_id):
        model = Kratos.Model()
        adjoint_parameters["problem_data"]["echo_level"].SetDouble(stabilization_coefficient)
        analysis = analysis_class_type(model, adjoint_parameters)
        for i in range(0, 10):
            analysis.OutputSolutionStep()
        analysis.Finalize()
        return analysis.time_series_data



if __name__ == '__main__':
    UnitTest.main()