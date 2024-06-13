import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.utilities.helper_utilities import CallOnAll
from KratosMultiphysics.OptimizationApplication.utilities.logger_utilities import time_decorator
from KratosMultiphysics.OptimizationApplication.controls.master_control import MasterControl
from KratosMultiphysics.OptimizationApplication.algorithms.standardized_PyRol_objective import StandardizedPyRolObjective

from pyrol import getCout, Objective, Problem, Solver
from pyrol.vectors import NumPyVector


def Factory(model: Kratos.Model, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem):
    return PYROLAlgorithms(model, parameters, optimization_problem)

class PYROLAlgorithms():
    """
        A python Wrapper to couple PyRol optimization library. https://pyrol.readthedocs.io/en/latest/intro.html
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
            "settings"          : {}
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
        self.__objective = StandardizedPyRolObjective(parameters["objective"], self.master_control, self._optimization_problem)
        # self.__constraints = []
        # for constraint_settings in parameters["constraints"]:
        #     self.__constraints.append(StandardizedNLOPTConstraint(constraint_settings, self.master_control, self._optimization_problem))

        raise RuntimeError(1)
        # nlopt settings
        NLOPT_settings = parameters["NLOPT_settings"]
        NLOPT_settings.ValidateAndAssignDefaults(self.GetDefaultParameters()["NLOPT_settings"])

        # nlopt verbosity
        self.nlopt_verbosity = NLOPT_settings["verbosity"].GetInt()

        # upper and lower bounds
        self.nlopt_controls_lower_bound = NLOPT_settings["controls_lower_bound"].GetString()
        self.nlopt_controls_upper_bound = NLOPT_settings["controls_upper_bound"].GetString()

        # nlopt algorithm
        self.algorithm_name = NLOPT_settings["algorithm_name"].GetString()
        if len(self.__constraints)==0:
            self.CheckOptimizerSupport(self.algorithm_name,None)
        else:
            for constraint in self.__constraints:
                if constraint.IsEqualityType():
                    self.CheckOptimizerSupport(self.algorithm_name,"equality")
                else:
                    self.CheckOptimizerSupport(self.algorithm_name,"inequality")

        # nlopt subsidiary algorithm
        self.subsidiary_algorithm_name = NLOPT_settings["subsidiary_algorithm_name"].GetString()
        if self.algorithm_name == "augmented_lagrangian":
            if self.subsidiary_algorithm_name == "":
                raise ValueError(f"The algorithm {self.algorithm_name} requires a subsidiary optimizer to be provided.")
            else:
                self.CheckOptimizerSupport(self.subsidiary_algorithm_name,None)

        # stopping
        self.stopping_criteria = NLOPT_settings["stopping_criteria"]
        self.stopping_criteria.ValidateAndAssignDefaults(self.GetDefaultParameters()["NLOPT_settings"]["stopping_criteria"])

        # alg specific settings
        self.opt_algorithm_specific_settings = NLOPT_settings["algorithm_specific_settings"]
    
    def GetMinimumBufferSize(self) -> int:
        return 2

    def Check(self):
        pass

    @time_decorator()
    def Initialize(self):
        self.converged = False