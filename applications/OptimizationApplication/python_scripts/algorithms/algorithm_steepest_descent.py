import KratosMultiphysics as Kratos
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.algorithms.standardized_objective import StandardizedObjective
from KratosMultiphysics.OptimizationApplication.controls.master_control import MasterControl
from KratosMultiphysics.OptimizationApplication.algorithms.algorithm import Algorithm
from KratosMultiphysics.OptimizationApplication.utilities.component_data_view import ComponentDataView
from KratosMultiphysics.OptimizationApplication.utilities.opt_line_search import CreateLineSearch
from KratosMultiphysics.OptimizationApplication.utilities.logger_utilities import time_decorator
from KratosMultiphysics.OptimizationApplication.utilities.logger_utilities import OptimizationAlgorithmTimeLogger
from KratosMultiphysics.OptimizationApplication.convergence_criteria.convergence_criterion import ConvergenceCriterion
from KratosMultiphysics.OptimizationApplication.convergence_criteria.combined_conv_criterion import CombinedConvCriterion
from KratosMultiphysics.OptimizationApplication.convergence_criteria.max_iter_conv_criterion import MaxIterConvCriterion
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem_utilities import OptimizationComponentFactory
from KratosMultiphysics.OptimizationApplication.utilities.logger_utilities import ListLogger


def Factory(model: Kratos.Model, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem):
    return AlgorithmSteepestDescent(model, parameters, optimization_problem)

class AlgorithmSteepestDescent(Algorithm):
    """
        A classical steepest descent algorithm to solve unconstrained optimization problems.
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
                "echo_level"      : 0,
                "line_search"     : {},
                "conv_settings"   : {}
            }
        }""")

    def __init__(self, model:Kratos.Model, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem):
        self.model = model
        self.parameters = parameters
        self._optimization_problem = optimization_problem

        parameters.ValidateAndAssignDefaults(self.GetDefaultParameters())

        self.master_control = MasterControl() # Need to fill it with controls
        self._optimization_problem.AddComponent(self.master_control)

        for control_name in parameters["controls"].GetStringArray():
            control = optimization_problem.GetControl(control_name)
            self.master_control.AddControl(control)


        settings = parameters["settings"]
        settings.ValidateAndAssignDefaults(self.GetDefaultParameters()["settings"])

        self.echo_level = settings["echo_level"].GetInt()

        ComponentDataView("algorithm", self._optimization_problem).SetDataBuffer(self.GetMinimumBufferSize())

        self._convergence_criteria = self.__CreateConvergenceCriteria(settings["conv_settings"])
        self.__line_search_method = CreateLineSearch(settings["line_search"], self._optimization_problem)

        self.__objective = StandardizedObjective(parameters["objective"], self.master_control, self._optimization_problem)
        self._optimization_problem.AddComponent(self.__objective)
        self.__control_field = None
        self.__obj_val = None

    def GetMinimumBufferSize(self) -> int:
        return 2

    def Check(self) -> None:
        self.master_control.Check()
        self.__objective.Check()

    @time_decorator()
    def Initialize(self) -> None:
        self.converged = False
        self.__obj_val = None
        self.master_control.Initialize()
        self.__objective.Initialize()
        self.__control_field = self.master_control.GetControlField()
        self.algorithm_data = ComponentDataView("algorithm", self._optimization_problem)
        self._convergence_criteria.Initialize()

    def Finalize(self) -> None:
        self.master_control.Finalize()
        self.__objective.Finalize()
        self._convergence_criteria.Finalize()

    @time_decorator()
    def ComputeSearchDirection(self, obj_grad: Kratos.TensorAdaptors.DoubleCombinedTensorAdaptor) -> Kratos.TensorAdaptors.DoubleCombinedTensorAdaptor:
        search_direction: Kratos.TensorAdaptors.DoubleCombinedTensorAdaptor = obj_grad.Clone()
        search_direction.data[:] *= -1.0
        Kratos.TensorAdaptors.DoubleCombinedTensorAdaptor(search_direction, perform_store_data_recursively=False, copy=False).StoreData()
        self.algorithm_data.GetBufferedData()["search_direction"] = search_direction
        return search_direction

    @time_decorator()
    def ComputeControlUpdate(self, alpha: Kratos.TensorAdaptors.DoubleCombinedTensorAdaptor) -> None:
        update: Kratos.TensorAdaptors.DoubleCombinedTensorAdaptor = self.algorithm_data.GetBufferedData()["search_direction"].Clone()
        update.data[:] *= alpha.data[:]
        Kratos.TensorAdaptors.DoubleCombinedTensorAdaptor(update, perform_store_data_recursively=False, copy=False).StoreData()
        self.algorithm_data.GetBufferedData()["control_field_update"] = update

    @time_decorator()
    def UpdateControl(self) -> Kratos.TensorAdaptors.DoubleCombinedTensorAdaptor:
        update: Kratos.TensorAdaptors.DoubleCombinedTensorAdaptor = self.algorithm_data.GetBufferedData()["control_field_update"]
        self.__control_field.data[:] += update.data
        Kratos.TensorAdaptors.DoubleCombinedTensorAdaptor(self.__control_field, perform_store_data_recursively=False, copy=False).StoreData()

    @time_decorator()
    def Output(self) -> None:
        self.algorithm_data.GetBufferedData()["control_field"] = Kratos.TensorAdaptors.DoubleCombinedTensorAdaptor(self.__control_field)
        for process in self._optimization_problem.GetListOfProcesses("output_processes"):
            if process.IsOutputStep():
                process.PrintOutput()

    def GetCurrentObjValue(self) -> float:
        return self.__obj_val

    def GetCurrentControlField(self) -> Kratos.TensorAdaptors.DoubleCombinedTensorAdaptor:
        return self.__control_field

    @time_decorator()
    def Solve(self) -> bool:
        while not self.converged:
            with OptimizationAlgorithmTimeLogger("AlgorithmSteepestDescent",self._optimization_problem.GetStep()):
                self._InitializeIteration()

                self.__obj_val = self.__objective.CalculateStandardizedValue(self.__control_field)
                obj_info = self.__objective.GetInfo()
                self.algorithm_data.GetBufferedData()["std_obj_value"] = obj_info["std_value"]
                self.algorithm_data.GetBufferedData()["rel_obj[%]"] = obj_info["rel_change [%]"]
                if "abs_change [%]" in obj_info:
                    self.algorithm_data.GetBufferedData()["abs_obj[%]"] = obj_info["abs_change [%]"]

                obj_grad = self.__objective.CalculateStandardizedGradient()

                self.ComputeSearchDirection(obj_grad)

                alpha = self.__line_search_method.ComputeStep()

                self.ComputeControlUpdate(alpha)

                self._FinalizeIteration()

                self.converged = self._convergence_criteria.IsConverged()

                self.Output()

                self.UpdateControl()

                ListLogger("Convergence info", self._convergence_criteria.GetInfo())

                self._optimization_problem.AdvanceStep()

        return self.converged

    def GetOptimizedObjectiveValue(self) -> float:
        if self.converged:
            return self.__obj_val
        else:
            raise RuntimeError("Optimization problem hasn't been solved.")

    def __CreateConvergenceCriteria(self, settings: Kratos.Parameters) -> ConvergenceCriterion:
        default_settings = Kratos.Parameters("""{
            "max_iter": 0,
            "type"    : "",
            "module"  : "KratosMultiphysics.OptimizationApplication.convergence_criteria",
            "settings": {}
        }""")
        settings.AddMissingParameters(default_settings)

        max_iter_params = Kratos.Parameters("""{
            "max_iter": """ + str(settings["max_iter"].GetInt()) + """
        }""")
        max_iter_conv_criteria = MaxIterConvCriterion(max_iter_params, self._optimization_problem)
        if settings["type"].GetString() == "":
            return max_iter_conv_criteria
        else:
            # create the additional convergence criterions given by the user
            additional_conv_params = Kratos.Parameters("""{
                "type"    : "",
                "module"  : "",
                "settings": {}
            }""")
            additional_conv_params["type"].SetString(settings["type"].GetString())
            additional_conv_params["module"].SetString(settings["module"].GetString())
            additional_conv_params["settings"] = settings["settings"]
            additional_conv: ConvergenceCriterion = OptimizationComponentFactory(self.model, additional_conv_params, self._optimization_problem)

            # now create the combined one and return the combined one
            combined_params = Kratos.Parameters("""{
                "operator": "or"
            }""")
            combined_conv_criteria = CombinedConvCriterion(self.model, combined_params, self._optimization_problem)
            combined_conv_criteria.Add(max_iter_conv_criteria) # add the max iter conv criteria
            combined_conv_criteria.Add(additional_conv) # add the additional conv criteria
            return combined_conv_criteria



