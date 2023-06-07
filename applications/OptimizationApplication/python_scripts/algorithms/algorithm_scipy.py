import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.algorithms.scipy_objective import SciPyObjective
from KratosMultiphysics.OptimizationApplication.algorithms.scipy_constraint import SciPyConstraint
from KratosMultiphysics.OptimizationApplication.controls.master_control import MasterControl
from KratosMultiphysics.OptimizationApplication.utilities.component_data_view import ComponentDataView
from KratosMultiphysics.OptimizationApplication.utilities.opt_convergence import CreateConvergenceCriteria
from KratosMultiphysics.OptimizationApplication.utilities.opt_line_search import CreateLineSearch
from KratosMultiphysics.OptimizationApplication.utilities.logger_utilities import TimeLogger
from KratosMultiphysics.OptimizationApplication.algorithms.algorithm import Algorithm
from scipy.optimize import minimize, show_options
import numpy as np
import json

def GetSciPyDefaultSolversOptions():
        return Kratos.Parameters("""{
                "SLSQP": {"maxiter": 0, "gtol": 1e-3},
                "CG": {"maxiter": 0, "gtol": 1e-3},
                "BFGS": {"maxiter": 0, "gtol": 1e-3, "xrtol":1e-3}
        }""")

def Factory(model: Kratos.Model, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem):
    return AlgorithmSciPy(model, parameters, optimization_problem)

class AlgorithmSciPy(Algorithm):
    """
        Interface to scipy algorithms to solve optimization problems.
    """

    @classmethod
    def GetDefaultParameters(cls):
        return Kratos.Parameters("""{
            "module"            : "KratosMultiphysics.OptimizationApplication.algorithms",
            "type"              : "PLEASE_PROVIDE_AN_ALGORITHM_CLASS_NAME",
            "scipy_method": "",
            "scipy_options": {},
            "control_bounds": [],
            "objective"         : {},
            "constraints"         : [],
            "controls"          : []
        }""")
    
    def __init__(self, model:Kratos.Model, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem):
        self.model = model
        self.parameters = parameters
        self._optimization_problem = optimization_problem

        self.parameters.ValidateAndAssignDefaults(self.GetDefaultParameters())

        default_scipy_settings = GetSciPyDefaultSolversOptions()

        self.scipy_method = self.parameters["scipy_method"].GetString()

        if not self.scipy_method in default_scipy_settings.keys():
            raise RuntimeError("scipy_method "+ self.scipy_method + " is not supported")
        
        self.parameters["scipy_options"].ValidateAndAssignDefaults(default_scipy_settings[self.scipy_method])
    
        self.master_control = MasterControl() # Need to fill it with controls

        for control_name in parameters["controls"].GetStringArray():
            control = optimization_problem.GetControl(control_name)
            self.master_control.AddControl(control)

        ComponentDataView("algorithm", self._optimization_problem).SetDataBuffer(self.GetMinimumBufferSize())

        self.__objective = SciPyObjective(parameters["objective"], self.master_control, self._optimization_problem)

        self.__constraints = []
        self.__sci_py_constraints = []
        if self.parameters["constraints"].size()>0:
            for const_i in range(self.parameters["constraints"].size()):
                constraint_settings = self.parameters["constraints"][const_i]
                constraint = SciPyConstraint(constraint_settings, self.master_control, self._optimization_problem)
                self.__constraints.append(constraint)
                const_type = 'eq'
                if constraint_settings["type"].GetString() != "=":
                    const_type = 'ineq'
                scipy_const_dict = {'type': const_type, 'fun': constraint.CalculateStandardizedValue, 'jac': constraint.CalculateStandardizedGradient}
                self.__sci_py_constraints.append(scipy_const_dict)

        self.__control_field = None

    def GetMinimumBufferSize(self) -> int:
        return 2

    def Check(self):
        pass

    def Initialize(self):
        self.__objective.Initialize()
        self.__objective.Check()
        for constraint in self.__constraints:
            constraint.Initialize()
            constraint.Check()
        self.master_control.Initialize()
        self.__control_field = self.master_control.GetControlField()

    def Finalize(self):
        pass

    def SolveOptimizationProblem(self) -> bool:
        self.Initialize()
        with TimeLogger("Solve Optimization problem", "Start", "End"):
            self.Solve()

    def __ConvertParamsToDict(self,params):
        dict_obj = {}

        for key in params.keys():
            if params[key].IsDouble():
                dict_obj[key] = params[key].GetDouble()
            elif params[key].IsInt():
                dict_obj[key] = params[key].GetInt()
            elif params[key].IsBool():
                dict_obj[key] = params[key].GetBool()
            elif params[key].IsString():
                dict_obj[key] = params[key].GetString()
            elif params[key].IsVector():
                dict_obj[key] = [value for value in params[key].GetVector()]
            elif params[key].IsMatrix():
                dict_obj[key] = [[value for value in row] for row in params[key].GetMatrix()]
            elif params[key].IsSubParameter():
                dict_obj[key] = self.__ConvertParamsToDict(params[key])

        return dict_obj

    def Solve(self):
        # Define the initial guess
        x0 = self.__control_field.Evaluate()

        # Convert kratos params to a Python dictionary.
        scipy_options = self.__ConvertParamsToDict(self.parameters["scipy_options"])

        # Get the bounds if there is any
        bounds_list = self.parameters["control_bounds"].GetVector()
        scipy_bounds = None
        if len(bounds_list)==2:
            scipy_bounds = [(bounds_list[0], bounds_list[1])] * len(x0)

        # Use the minimize function with the "BFGS" algorithm (gradient-based)
        result = minimize(self.__objective.CalculateStandardizedValue, x0, jac=self.__objective.CalculateStandardizedGradient, method=self.scipy_method, constraints=self.__sci_py_constraints, bounds=scipy_bounds, options=scipy_options)

        # Print the optimization result
        print(result)
