try:
    import nlopt
    import numpy
except ImportError:
    raise Exception("NLOPT python library is not available")

import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.algorithms.standardized_NLOPT_objective import StandardizedNLOPTObjective
from KratosMultiphysics.OptimizationApplication.algorithms.standardized_NLOPT_constraint import StandardizedNLOPTConstraint
from KratosMultiphysics.OptimizationApplication.controls.master_control import MasterControl
from KratosMultiphysics.OptimizationApplication.algorithms.algorithm import Algorithm
from KratosMultiphysics.OptimizationApplication.utilities.component_data_view import ComponentDataView


def Factory(model: Kratos.Model, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem):
    return NLOPTAlgorithms(model, parameters, optimization_problem)

class NLOPTAlgorithms(Algorithm):
    """
        A classical steepest descent algorithm to solve unconstrainted optimization problems.
    """

    @classmethod
    def GetDefaultParameters(cls):
        return Kratos.Parameters("""{
            "module"            : "KratosMultiphysics.OptimizationApplication.algorithms",
            "type"              : "PLEASE_PROVIDE_AN_ALGORITHM_CLASS_NAME",
            "objective"         : {},
            "constraints"         : [],
            "controls"          : [],
            "echo_level"        : 0,
            "NLOPT_settings"          : {
                "algorithm_name"      : "MMA",
                "algorithm_settings"     : {},
                "stopping_criteria"   : {
                    "objective_rel_tol": "",
                    "objective_abs_tol": "",
                    "controls_rel_tol": "",
                    "controls_abs_tol": "",
                    "maximum_function_evalualtion": 10
                }
            }
        }""")

    def __init__(self, model:Kratos.Model, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem):
        self.model = model
        self.parameters = parameters
        self._optimization_problem = optimization_problem

        parameters.ValidateAndAssignDefaults(self.GetDefaultParameters())

        # controls
        self.master_control = MasterControl()
        for control_name in parameters["controls"].GetStringArray():
            control = optimization_problem.GetControl(control_name)
            self.master_control.AddControl(control)
        self.__control_field = None

        # objective & constraints
        self.__objective = StandardizedNLOPTObjective(parameters["objective"], self.master_control, self._optimization_problem)
        self.__constraints = []
        for constraint_settings in parameters["constraints"]:
            self.__constraints.append(StandardizedNLOPTConstraint(constraint_settings, self.master_control, self._optimization_problem))

        # nlopt settings
        NLOPT_settings = parameters["NLOPT_settings"]
        NLOPT_settings.ValidateAndAssignDefaults(self.GetDefaultParameters()["NLOPT_settings"])

        # nlopt algorithm settings
        self.algorithm_name = NLOPT_settings["algorithm_name"].GetString()
        self.algorithm_settings = NLOPT_settings["algorithm_settings"]

        # stopping
        self.stopping_criteria = NLOPT_settings["stopping_criteria"]
        self.stopping_criteria.ValidateAndAssignDefaults(self.GetDefaultParameters()["NLOPT_settings"]["stopping_criteria"])

    def GetMinimumBufferSize(self) -> int:
        return 2

    def Check(self):
        pass

    def Initialize(self):
        self.converged = False
        self.__objective.Initialize()
        self.__objective.Check()
        for constraint in self.__constraints:
            constraint.Initialize()
            constraint.Check()
        self.master_control.Initialize()
        self.__control_field = self.master_control.GetControlField()
        self.algorithm_data = ComponentDataView("algorithm", self._optimization_problem)

    def Finalize(self):
        self.__objective.Finalize()
        for constraint in self.__constraints:
            constraint.Finalize()
        self.master_control.Finalize()

    def Solve(self):
        x0 = self.__control_field.Evaluate()
        opt = nlopt.opt(nlopt.LD_MMA, x0.size)
        opt.set_min_objective(self.__objective.CalculateStandardizedValueAndGradients)
        for constraint in self.__constraints:
            opt.add_inequality_constraint(lambda x,grad: constraint.CalculateStandardizedValueAndGradients(x,grad),1e-8)
        opt.set_ftol_rel(1e-2)
        # opt.set_xtol_rel(1e-6)
        opt.set_maxeval(100)
        opt.optimize(x0)
        minf = opt.last_optimum_value()
        print("minimum value = ", minf)
