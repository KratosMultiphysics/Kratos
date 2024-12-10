import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.algorithms.standardized_SciPy_objective import StandardizedSciPyObjective
from KratosMultiphysics.OptimizationApplication.algorithms.standardized_SciPy_constraint import StandardizedSciPyConstraint
from KratosMultiphysics.OptimizationApplication.controls.master_control import MasterControl
from KratosMultiphysics.OptimizationApplication.algorithms.algorithm import Algorithm
from KratosMultiphysics.OptimizationApplication.utilities.component_data_view import ComponentDataView
from KratosMultiphysics.OptimizationApplication.utilities.helper_utilities import CallOnAll
from KratosMultiphysics.OptimizationApplication.utilities.logger_utilities import time_decorator

import numpy as np

from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem_utilities import OutputGradientFields

try:
    import scipy
except ImportError:
    raise Exception("SciPy python library is not available")

import scipy.optimize
from scipy.optimize import NonlinearConstraint

def Factory(model: Kratos.Model, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem):
    return SciPyAlgorithms(model, parameters, optimization_problem)

class SciPyAlgorithms(Algorithm):
    """
        A SciPy wrapper to use algorithms from Library. Provide Method and it settings in the options dict.
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
            "SciPy_settings"    : {
                "method"      : "",
                "lower_bound" : -1e16,
                "upper_bound" : 1e16,
                "options"     : {}
            }
        }""")

    def __init__(self, model:Kratos.Model, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem):
        self.model = model
        self.parameters = parameters
        self.parameters.ValidateAndAssignDefaults(self.GetDefaultParameters())
        self._optimization_problem = optimization_problem

        # scipy settings
        self.SciPy_settings = self.parameters["SciPy_settings"]
        self.SciPy_settings.ValidateAndAssignDefaults(self.GetDefaultParameters()["SciPy_settings"])

        # controls
        self.master_control = MasterControl()
        self._optimization_problem.AddComponent(self.master_control)
        for control_name in parameters["controls"].GetStringArray():
            control = optimization_problem.GetControl(control_name)
            self.master_control.AddControl(control)
        self.__control_field = None

        # objective & constraints
        self.__objective = StandardizedSciPyObjective(parameters["objective"], self.master_control, self._optimization_problem)
        self._optimization_problem.AddComponent(self.__objective)
        self.bounds = scipy.optimize.Bounds(lb=self.SciPy_settings["lower_bound"].GetDouble(), ub=self.SciPy_settings["upper_bound"].GetDouble())
        self.__kratos_constraints = []
        self.__scipy_constraints = []
        for constraint_settings in parameters["constraints"]:
            constraint = StandardizedSciPyConstraint(constraint_settings, self.master_control, self._optimization_problem)
            self.__kratos_constraints.append(constraint)
            self._optimization_problem.AddComponent(constraint)

            # create SciPy specific constraints
            scipy_constraint = NonlinearConstraint(constraint.CalculateStandardizedValue, constraint.GetLowerBound(), constraint.GetUpperBound(), constraint.CalculateStandardizedGradient)
            self.__scipy_constraints.append(scipy_constraint)

    def GetMinimumBufferSize(self) -> int:
        return 2

    def Check(self):
        pass

    @time_decorator()
    def Initialize(self):
        self.converged = False
        self.master_control.Initialize()
        self.__objective.Initialize()
        self.__objective.Check()
        CallOnAll(self.__kratos_constraints, StandardizedSciPyConstraint.Initialize)
        self.__control_field = self.master_control.GetControlField()
        self.algorithm_data = ComponentDataView("algorithm", self._optimization_problem)
        self.algorithm_data.SetDataBuffer(self.GetMinimumBufferSize())

        # create SciPy des variables vector
        self.x0 = self.__control_field.Evaluate().reshape(-1)

    @time_decorator()
    def Finalize(self):
        self.__objective.Finalize()
        for constraint in self.__kratos_constraints:
            constraint.Finalize()
        self.master_control.Finalize()

    @time_decorator()
    def Solve(self):

    #     parameters = Kratos.Parameters("""{
    #     "type": "algorithm_steepest_descent",
    #     "settings": {
    #         "echo_level": 0,
    #         "line_search": {
    #             "type": "BB_step",
    #             "init_step": 5e-2,
    #             "max_step": 2e-1,
    #             "gradient_scaling": "l2_norm"
    #         },
    #         "conv_settings": {
    #             "type": "target_value",
    #             "max_iter": 100,
    #             "target_value" : 6e-11
    #         }
    #     },
    #     "controls": [
    #         "material_control"
    #     ],
    #     "objective": {
    #         "response_name": "damage_response",
    #         "type": "minimization",
    #         "scaling": 1.0
    #     }
    # }""")
    #     from KratosMultiphysics.OptimizationApplication.algorithms.algorithm_steepest_descent import AlgorithmSteepestDescent
    #     optimization_problem = self._optimization_problem 
    #     sd_algorithm = AlgorithmSteepestDescent(self.model, parameters, optimization_problem)
    #     sd_algorithm.Initialize()
    #     sd_algorithm.Solve()

    #     from KratosMultiphysics.OptimizationApplication.responses.response_routine import ResponseRoutine

    #     self.__objective.objective = self._optimization_problem.GetComponent("damage_response", ResponseRoutine)
    #     self.algorithm_data.RemoveComponentData("search_direction")

    #     # create SciPy des variables vector
    #     size = self.master_control.GetControlField().GetCollectiveFlattenedDataSize()
    #     self.x0 = sd_algorithm.GetCurrentControlField().Evaluate().reshape(-1)

    #     print("test value", self.__objective.objective.CalculateStandardizedValue(sd_algorithm.GetCurrentControlField(), False) )


        if not self.__scipy_constraints:
            res = scipy.optimize.minimize(self.__objective.CalculateStandardizedValue, self.x0, method=self.SciPy_settings["method"].GetString(), jac=self.__objective.CalculateStandardizedGradient, options=self.__GetOptions(), callback=self.Output, bounds=self.bounds)
        elif self.__scipy_constraints:
            res = scipy.optimize.minimize(self.__objective.CalculateStandardizedValue, self.x0, method=self.SciPy_settings["method"].GetString(), jac=self.__objective.CalculateStandardizedGradient, constraints=self.__scipy_constraints, options=self.__GetOptions(), callback=self.Output, bounds=self.bounds)
        Kratos.Logger.PrintInfo(self.__class__.__name__, f"res::success: {res.success}")
        Kratos.Logger.PrintInfo(self.__class__.__name__, f"res::status: {res.status}")
        Kratos.Logger.PrintInfo(self.__class__.__name__, f"res::message: {res.message}")

    @time_decorator()
    def Output(self, *args) -> KratosOA.CollectiveExpression:
        pass
        # # SciPy calls it at the end of each optimization iterations. Sometime it doesn't call f(x) during iteration, hence to avoid lack of data in the buffer we have a flag "computed" here.
        # if self.__objective.computed:

        #     Kratos.Logger.PrintInfo(self.__class__.__name__, f"Output iteration {self._optimization_problem.GetStep()}")

        #     shape = [c.GetItemShape() for c in self.__control_field.GetContainerExpressions()]
        #     KratosOA.CollectiveExpressionIO.Read(self.__control_field, args[0], shape)

        #     self.algorithm_data.GetBufferedData().SetValue("control_field", self.__control_field.Clone(), overwrite=True)
        #     OutputGradientFields(self.__objective, self._optimization_problem, True)
        #     for constraint in self.__kratos_constraints:
        #         OutputGradientFields(constraint, self._optimization_problem, True)
        #     for process in self._optimization_problem.GetListOfProcesses("output_processes"):
        #         if process.IsOutputStep():
        #             process.PrintOutput()
            
        #     # Advance in Optimization Iteration
        #     self._optimization_problem.AdvanceStep()
        #     self.__objective.computed = False

    def __GetOptions(self):
        options = self.SciPy_settings["options"]
        dict_options = dict()
        for key in options.keys():
            if options[key].IsBool():
                dict_options[key] = options[key].GetBool()
            elif options[key].IsDouble():
                dict_options[key] = options[key].GetDouble()
            elif options[key].IsInt():
                dict_options[key] = options[key].GetInt()
            elif options[key].IsString():
                dict_options[key] = options[key].GetString()
            else:
                raise RuntimeError(f"Unsupported data type with key = \"{key}\" [ value = {options[key]} ].")
        return dict_options