try:
    import scipy
except ImportError:
    raise Exception("SciPy python library is not available")

import KratosMultiphysics as Kratos
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.algorithms.standardized_SciPy_objective import StandardizedSciPyObjective
from KratosMultiphysics.OptimizationApplication.algorithms.standardized_SciPy_constraint import StandardizedSciPyConstraint
from KratosMultiphysics.OptimizationApplication.controls.master_control import MasterControl
from KratosMultiphysics.OptimizationApplication.algorithms.algorithm import Algorithm
from KratosMultiphysics.OptimizationApplication.utilities.component_data_view import ComponentDataView
from KratosMultiphysics.OptimizationApplication.utilities.logger_utilities import DictLogger
from KratosMultiphysics.OptimizationApplication.utilities.helper_utilities import CallOnAll
from KratosMultiphysics.OptimizationApplication.utilities.logger_utilities import time_decorator
import scipy.optimize



def Factory(model: Kratos.Model, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem):
    return SciPyAlgorithms(model, parameters, optimization_problem)

class SciPyAlgorithms(Algorithm):
    """
        A SciPy wrapper to use algorithms from Library
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
            "SciPy_settings"          : {
                "method"      : "",
                "options"     : {}
                },
                "algorithm_specific_settings"   : {}
            }
        }""")

    def __init__(self, model:Kratos.Model, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem):
        self.model = model
        self.parameters = parameters
        self._optimization_problem = optimization_problem

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
        self.__constraints = []
        for constraint_settings in parameters["constraints"]:
            constraint = StandardizedSciPyConstraint(constraint_settings, self.master_control, self._optimization_problem)
            self._optimization_problem.AddComponent(constraint)
            self.__constraints.append(constraint)

        # scipy settings
        self.SciPy_settings = parameters["SciPy_settings"]

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
        CallOnAll(self.__constraints, StandardizedSciPyConstraint.Initialize)
        self.__control_field = self.master_control.GetControlField()
        self.algorithm_data = ComponentDataView("algorithm", self._optimization_problem)

        # create nlopt optimizer
        self.x0 = self.__control_field.Evaluate()

    @time_decorator()
    def Finalize(self):
        self.__objective.Finalize()
        for constraint in self.__constraints:
            constraint.Finalize()
        self.master_control.Finalize()

    @time_decorator()
    def Solve(self):
        res = scipy.optimize.minimize(self.__objective.CalculateStandardizedValue, self.x0, method=self.SciPy_settings["method"].GetString(), jac=self.__objective.CalculateStandardizedGradient, options=self.__GetOptions())
        CallOnAll(self._optimization_problem.GetListOfProcesses("output_processes"), Kratos.OutputProcess.PrintOutput)


    def __GetOptions(self):
        Kratos.Parameters
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
        return dict_options