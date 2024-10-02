import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.algorithms.standardized_objective import StandardizedObjective
from KratosMultiphysics.OptimizationApplication.controls.master_control import MasterControl
from KratosMultiphysics.OptimizationApplication.utilities.component_data_view import ComponentDataView
from KratosMultiphysics.OptimizationApplication.utilities.logger_utilities import DictLogger
from KratosMultiphysics.OptimizationApplication.utilities.logger_utilities import TimeLogger
from pyrol import Objective
import numpy as np
from pyrol.vectors import NumPyVector as myVector

class StandardizedPyRolObjective(Objective):
    """Standardized objective response function

    This class transformed a user-given optimization problem into the standard format.
    Supported objective types:
        "minimization",
        "maximization"

    """
    def __init__(self, parameters: Kratos.Parameters, master_control: MasterControl, optimization_problem: OptimizationProblem, required_buffer_size: int = 2):
        """
        Initializes the StandardizedPyRolObjective class.

        Args:
            parameters (Kratos.Parameters): The parameters for the objective.
            master_control (MasterControl): The master control object.
            optimization_problem (OptimizationProblem): The optimization problem instance.
            required_buffer_size (int, optional): The required buffer size. Defaults to 2.

        """
        default_parameters = Kratos.Parameters("""{
            "response_name": "",
            "type"         : "",
            "scaling"      : 1.0
        }""")
        parameters.ValidateAndAssignDefaults(default_parameters)

        self.__objective = StandardizedObjective(parameters, master_control, optimization_problem, required_buffer_size)

        self.__optimization_problem = optimization_problem

        super().__init__()

    def value(self, x: myVector, tol: float, save_data: bool = True) -> float:
        """
        Calculate the standardized objective value for a given input vector.
        Args:
            x (myVector): Input vector, expected to be a NumPyVector.
            tol (float): Tolerance value (currently unused in the function).
            save_data (bool, optional): Flag to determine whether to save data. Defaults to True.
        Returns:
            float: The calculated standardized objective value.
        """
        # Convert x into numpy array (x is a NumPyVector)
        list_x = [value for value in x]
        numpy_x = np.array(list_x)  # Direct conversion using np.array() doesn't work, use list 
        
        control_field = self.GetMasterControl().GetEmptyField()
        shape = [c.GetItemShape() for c in control_field.GetContainerExpressions()]
        KratosOA.CollectiveExpressionIO.Read(control_field, numpy_x, shape)
        value = self.__objective.CalculateStandardizedValue(control_field, save_data)

        if save_data:
            for process in self.__optimization_problem.GetListOfProcesses("output_processes"):
                if process.IsOutputStep():
                    process.PrintOutput()

            self.__optimization_problem.AdvanceStep()

        return value

    def gradient(self, g:myVector, x:myVector, tol:float, save_field: bool = True):
        """
        Computes the gradient of the objective function.

        Parameters:
        g (myVector): The vector to store the computed gradient.
        x (myVector): The current point at which the gradient is evaluated.
        tol (float): Tolerance for the gradient computation.
        save_field (bool, optional): Flag to indicate whether to save the field. Defaults to True.

        Raises:
        RuntimeError: If there is an issue with the gradient computation.
        """

        self.value(x, False)  # Compute new primal if x has changed. Does nothing if x the same

        result = self.__objective.CalculateStandardizedGradient(save_field)

        g = myVector(result.Evaluate().reshape(-1) * self.__objective.GetScalingFactor())

    def hessVec(self, hv, v, x, tol):
        raise RuntimeError("Hessian-vector product is not implemented for the pyrol objective response function.")

    def Initialize(self):
        self.__objective.Initialize()

    def Check(self):
        self.__objective.Check()

    def GetMasterControl(self) -> MasterControl:
        return self.__objective.GetMasterControl()
                
    def Finalize(self):
        self.__objective.Finalize()