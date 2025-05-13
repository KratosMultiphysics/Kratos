import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.utilities.helper_utilities import CallOnAll
from KratosMultiphysics.OptimizationApplication.utilities.logger_utilities import time_decorator
from KratosMultiphysics.OptimizationApplication.controls.master_control import MasterControl
from KratosMultiphysics.OptimizationApplication.algorithms.standardized_PyRol_objective import StandardizedPyRolObjective
from KratosMultiphysics.OptimizationApplication.algorithms.standardized_PyRol_constraint import StandardizedPyRolConstraint
from KratosMultiphysics.OptimizationApplication.utilities.component_data_view import ComponentDataView
from KratosMultiphysics.OptimizationApplication.algorithms.algorithm import Algorithm

try:
    import pyrol
except ImportError:
    raise Exception("PyRol python library is not available")


from pyrol import getCout, Objective, Problem, Solver, getParametersFromXmlFile, Bounds
from pyrol.vectors import NumPyVector as myVector
import numpy as np


def Factory(model: Kratos.Model, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem):
    return PyRolAlgorithms(model, parameters, optimization_problem)


class PyRolAlgorithms(Algorithm):
    """
    PyRolAlgorithms is a class that implements optimization algorithms using the PyRol library.
    Provide Method and its settings in the xtml file.
    """

    @classmethod
    def GetDefaultParameters(cls):
        return Kratos.Parameters(cls.create_default_algorithm_config())

    @classmethod
    def create_default_algorithm_config(cls):
        return """{
            "module": "KratosMultiphysics.OptimizationApplication.algorithms",
            "type": "PLEASE_PROVIDE_AN_ALGORITHM_CLASS_NAME",
            "objective": {},
            "constraints": [],
            "controls": [],
            "echo_level": 0,
            "pyrol_settings": {
                "input_filename": "",
                "lower_bound": -1e16,
                "upper_bound": 1e16
            }
        }"""
    
    def __init__(self, model:Kratos.Model, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem):
        """
        Initializes the PyRol algorithm with the given model, parameters, and optimization problem.

        Args:
            model (Kratos.Model): The Kratos model to be used.
            parameters (Kratos.Parameters): The parameters for the optimization algorithm.
            optimization_problem (OptimizationProblem): The optimization problem to be solved.

        Attributes:
            model (Kratos.Model): The Kratos model to be used.
            parameters (Kratos.Parameters): The parameters for the optimization algorithm.
            _optimization_problem (OptimizationProblem): The optimization problem to be solved.
            master_control (MasterControl): The master control object managing all controls.
            __control_field (None): Placeholder for control field, initially set to None.
            settings (Kratos.Parameters): The settings for the PyRol algorithm.
            method_settings (Kratos.Parameters): The method settings loaded from an XML file.
            __objective (StandardizedPyRolObjective): The standardized objective for the PyRol algorithm.
            __constraints (list): A list of standardized constraints for the PyRol algorithm.
        """
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

        # PyRol settings
        self.settings = parameters["pyrol_settings"]
        self.settings.ValidateAndAssignDefaults(self.GetDefaultParameters()["pyrol_settings"])
        self.method_settings = getParametersFromXmlFile(self.settings["input_filename"].GetString())

        # objective & constraints
        self.__objective = StandardizedPyRolObjective(parameters["objective"], self.master_control, self._optimization_problem)
        self.__constraints = []
        for constraint_settings in parameters["constraints"]:
            self.__constraints.append(StandardizedPyRolConstraint(constraint_settings, self.master_control, self._optimization_problem))
    
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
        CallOnAll(self.__constraints, StandardizedPyRolConstraint.Initialize)
        self.__control_field = self.master_control.GetControlField()
        self.algorithm_data = ComponentDataView("algorithm", self._optimization_problem)
        self.algorithm_data.SetDataBuffer(self.GetMinimumBufferSize())
        
        # create SciPy des variables vector
        size = self.master_control.GetControlField().GetCollectiveFlattenedDataSize()
        self.x0 = myVector(self.__control_field.Evaluate().reshape(-1))
        self.grad = self.x0.dual()

        # upper and lower bounds
        self.lower_bounds = myVector(np.full(size, self.settings["lower_bound"].GetDouble()))
        self.upper_bounds = myVector(np.full(size, self.settings["upper_bound"].GetDouble()))
        self.bounds = Bounds(self.lower_bounds, self.upper_bounds)

        # setup the problem
        self.problem = Problem(self.__objective, self.x0, self.grad)
        self.problem.addBoundConstraint(self.bounds)
        for constraint in self.__constraints:
            self.problem.addConstraint(constraint.GetName(), constraint, self.x0.dual())

    @time_decorator()
    def Finalize(self):
        self.__objective.Finalize()
        for constraint in self.__constraints:
            constraint.Finalize()
        self.master_control.Finalize()

    @time_decorator()
    def Solve(self):
        stream = getCout()
        solver = Solver(self.problem, self.method_settings)
        solver.solve(stream)
