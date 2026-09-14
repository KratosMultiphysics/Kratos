import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.utilities.logger_utilities import time_decorator
from KratosMultiphysics.OptimizationApplication.algorithms.algorithm_relaxed_gradient_projection import AlgorithmRelaxedGradientProjection


def Factory(model: Kratos.Model, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem):
    return AlgorithmMomentumRelaxedGradientProjection(model, parameters, optimization_problem)

class AlgorithmMomentumRelaxedGradientProjection(AlgorithmRelaxedGradientProjection):
    """
        Nesterov Accelerated Gradient method to solve unconstrained optimization problems.
        The implementation is based on https://paperswithcode.com/method/nesterov-accelerated-gradient
        In current version to save one response evaluation, the f(x) is only calculated on the momentum points.
    """

    @classmethod
    def GetDefaultParameters(cls):
        return Kratos.Parameters("""{
            "module"            : "KratosMultiphysics.OptimizationApplication.algorithms",
            "type"              : "PLEASE_PROVIDE_AN_ALGORITHM_CLASS_NAME",
            "objective"         : {},
            "constraints"       : [],
            "controls"          : [],
            "echo_level"        : 0,
            "settings"          : {
                "echo_level"                 : 0,
                "line_search"                : {},
                "conv_settings"              : {},
                "linear_solver_settings"     : {},
                "history_size"               : 20,
                "max_inner_iter"             : 20,
                "buffer_coeff_update_factor" : 0.1,
                "eta"                        : 0.9
            }
        }""")

    def __init__(self, model:Kratos.Model, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem):
        super().__init__(model, parameters, optimization_problem)
        self.eta = self.parameters["settings"]["eta"].GetDouble()

    @time_decorator()
    def ComputeControlUpdate(self, alpha: Kratos.TensorAdaptors.DoubleCombinedTensorAdaptor):
        # compute the correction part from momentum point
        algorithm_buffered_data = self.algorithm_data.GetBufferedData()
        search_direction: Kratos.TensorAdaptors.DoubleCombinedTensorAdaptor = algorithm_buffered_data["search_direction"]
        update = Kratos.TensorAdaptors.DoubleCombinedTensorAdaptor(search_direction, perform_store_data_recursively=False)
        update.data[:] *= alpha.data
        update.StoreData()
        # add momentum to the correction update to compute new momentum point. The momentum from
        # the previous call is bookkept in "algorithm"'s own buffered data (step_index 1, i.e. the
        # previous step -- this call writes this step's momentum to step_index 0) rather than a
        # plain Python attribute, so that a restart checkpoint/resume round-trips it correctly;
        # see optimization_problem_restart_input_process.py (mirrors
        # algorithm_nesterov_accelarated_gradient.py's identical fix).
        if algorithm_buffered_data.HasValue("momentum", 1):
            prev_momentum: Kratos.TensorAdaptors.DoubleCombinedTensorAdaptor = algorithm_buffered_data.GetValue("momentum", 1)
            mom_update = Kratos.TensorAdaptors.DoubleCombinedTensorAdaptor(update, perform_store_data_recursively=False)
            mom_update.data[:] += prev_momentum.data * self.eta
            mom_update.StoreData()

            full_update = Kratos.TensorAdaptors.DoubleCombinedTensorAdaptor(update, perform_store_data_recursively=False)
            full_update.data[:] += mom_update.data * self.eta
            full_update.StoreData()
            algorithm_buffered_data["momentum"] = mom_update
        else:
            full_update = Kratos.TensorAdaptors.DoubleCombinedTensorAdaptor(update, perform_store_data_recursively=False)
            full_update.data[:] *= (1 + self.eta)
            full_update.StoreData()
            algorithm_buffered_data["momentum"] = update

        algorithm_buffered_data.SetValue("control_field_update", full_update, overwrite=True)


        return full_update
