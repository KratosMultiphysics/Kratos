import KratosMultiphysics as Kratos
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.utilities.logger_utilities import time_decorator
from KratosMultiphysics.OptimizationApplication.algorithms.algorithm_steepest_descent import AlgorithmSteepestDescent


def Factory(model: Kratos.Model, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem):
    return AlgorithmNesterovAcceleratedGradient(model, parameters, optimization_problem)

class AlgorithmNesterovAcceleratedGradient(AlgorithmSteepestDescent):
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
            "controls"          : [],
            "echo_level"        : 0,
            "settings"          : {
                "eta"             : 0.9,
                "echo_level"      : 0,
                "line_search"     : {},
                "conv_settings"   : {}
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
        update = Kratos.TensorAdaptors.DoubleCombinedTensorAdaptor(search_direction, perform_collect_data_recursively=False, perform_store_data_recursively=False)
        update.data[:] *= alpha.data[:]
        update.StoreData()
        # add momentum to the correction update to compute new momentum point. The momentum from
        # the previous call is bookkept in "algorithm"'s own buffered data (step_index 1, i.e. the
        # previous step -- this call writes this step's momentum to step_index 0) rather than a
        # plain Python attribute, so that a restart checkpoint/resume round-trips it correctly;
        # see optimization_problem_restart_input_process.py.
        if algorithm_buffered_data.HasValue("momentum", 1):
            prev_momentum: Kratos.TensorAdaptors.DoubleCombinedTensorAdaptor = algorithm_buffered_data.GetValue("momentum", 1)
            mom_update = Kratos.TensorAdaptors.DoubleCombinedTensorAdaptor(update, perform_store_data_recursively=False)
            mom_update.data[:] += prev_momentum.data * self.eta
            mom_update.StoreData()
            control_field_update = Kratos.TensorAdaptors.DoubleCombinedTensorAdaptor(mom_update, perform_store_data_recursively=False)
            control_field_update.data[:] = update.data + mom_update.data * self.eta
            control_field_update.StoreData()
            algorithm_buffered_data["control_field_update"] = control_field_update
            algorithm_buffered_data["momentum"] = mom_update
        else:
            control_field_update = Kratos.TensorAdaptors.DoubleCombinedTensorAdaptor(update, perform_store_data_recursively=False)
            control_field_update.data[:] = update.data * (1 + self.eta)
            control_field_update.StoreData()
            algorithm_buffered_data["control_field_update"] = control_field_update
            algorithm_buffered_data["momentum"] = update
