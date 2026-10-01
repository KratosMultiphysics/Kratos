import KratosMultiphysics as Kratos
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.utilities.logger_utilities import time_decorator
from KratosMultiphysics.OptimizationApplication.algorithms.algorithm_steepest_descent import AlgorithmSteepestDescent


def Factory(model: Kratos.Model, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem):
    return AlgorithmAdam(model, parameters, optimization_problem)

class AlgorithmAdam(AlgorithmSteepestDescent):
    """
        Adam (Adaptive Moment Estimation) method to solve unconstrained optimization problems.
        The implementation is based on https://arxiv.org/pdf/1412.6980
    """

    @classmethod
    def GetDefaultParameters(cls):
        return Kratos.Parameters("""{
            "module"            : "KratosMultiphysics.OptimizationApplication.algorithms",
            "type"              : "PLEASE_PROVIDE_AN_ALGORITHM_CLASS_NAME",
            "objective"         : {},
            "controls"          : [],
            "echo_level"        : 0,
            "settings"          : {
                "beta_1"          : 0.9,
                "beta_2"          : 0.999,
                "epsilon"         : 1e-8,
                "echo_level"      : 0,
                "line_search"     : {},
                "conv_settings"   : {}
            }
        }""")

    def __init__(self, model:Kratos.Model, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem):
        super().__init__(model, parameters, optimization_problem)
        self.beta_1 = self.parameters["settings"]["beta_1"].GetDouble()
        self.beta_2 = self.parameters["settings"]["beta_2"].GetDouble()
        self.epsilon = self.parameters["settings"]["epsilon"].GetDouble()
        self.prev_moment_1 = None
        self.prev_moment_2 = None

    @time_decorator()
    def ComputeControlUpdate(self, alpha: Kratos.TensorAdaptors.DoubleCombinedTensorAdaptor):
        search_direction: Kratos.TensorAdaptors.DoubleCombinedTensorAdaptor = self.algorithm_data.GetBufferedData()["search_direction"]
        time_step = self._optimization_problem.GetStep() + 1

        # update the biased first and second raw moment estimates
        moment_1 = Kratos.TensorAdaptors.DoubleCombinedTensorAdaptor(search_direction, perform_collect_data_recursively=False, perform_store_data_recursively=False)
        moment_2 = Kratos.TensorAdaptors.DoubleCombinedTensorAdaptor(search_direction, perform_collect_data_recursively=False, perform_store_data_recursively=False)
        if self.prev_moment_1 and self.prev_moment_2:
            moment_1.data[:] = self.prev_moment_1.data * self.beta_1 + search_direction.data * (1.0 - self.beta_1)
            moment_2.data[:] = self.prev_moment_2.data * self.beta_2 + search_direction.data ** 2 * (1.0 - self.beta_2)
        else:
            moment_1.data[:] = search_direction.data * (1.0 - self.beta_1)
            moment_2.data[:] = search_direction.data ** 2 * (1.0 - self.beta_2)
        moment_1.StoreData()
        moment_2.StoreData()
        self.prev_moment_1 = moment_1
        self.prev_moment_2 = moment_2

        # bias correct the moment estimates and compute the control field update
        moment_1_hat = moment_1.data / (1.0 - self.beta_1 ** time_step)
        moment_2_hat = moment_2.data / (1.0 - self.beta_2 ** time_step)

        control_field_update = Kratos.TensorAdaptors.DoubleCombinedTensorAdaptor(moment_2, perform_store_data_recursively=False)
        control_field_update.data[:] = alpha.data * moment_1_hat / (moment_2_hat ** 0.5 + self.epsilon)
        control_field_update.StoreData()
        self.algorithm_data.GetBufferedData()["control_field_update"] = control_field_update
