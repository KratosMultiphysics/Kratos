import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.algorithms.standardized_constraint import StandardizedConstraint
from KratosMultiphysics.OptimizationApplication.controls.master_control import MasterControl
from pyrol import Constraint
import numpy as np
from pyrol.vectors import NumPyVector as myVector
from KratosMultiphysics.kratos_utilities import IssueDeprecationWarning

class StandardizedPyRolConstraint(Constraint):
    """Standardized objective response function

    This class transformed a user-given optimization problem into the standard format.
    Supported objective types:
        "="
    Current implemetation only supports equality constraints and reqires Hessian information to be provided (not implemented).
    """
    def __init__(self, parameters: Kratos.Parameters, master_control: MasterControl, optimization_problem: OptimizationProblem, required_buffer_size: int = 2):
        """
        Initializes the StandardizedPyRolConstraint class.

        Args:
            parameters (Kratos.Parameters): The parameters for the constraint.
            master_control (MasterControl): The master control object.
            optimization_problem (OptimizationProblem): The optimization problem instance.
            required_buffer_size (int, optional): The required buffer size. Defaults to 2.

        """
        # backward compatibility
        if parameters.Has("response_name"):
            IssueDeprecationWarning(self.__class__.__name__, "\"response_name\" is deprecated. Please use \"response_expression\".")
            parameters.AddString("response_expression", parameters["response_name"].GetString())
            parameters.RemoveValue("response_name")
        
        default_parameters = Kratos.Parameters("""{
            "response_expression" : "",
            "type"                : "=",
            "scaled_ref_value"    : "initial_value"
        }""")

        if parameters.Has("scaled_ref_value") and parameters["scaled_ref_value"].IsDouble():
            default_parameters["scaled_ref_value"].SetDouble(0.0)

        if parameters.Has("type") and not parameters["type"].GetString() == "=":
            raise Exception("Only equality constraints are supported at the moment.")

        parameters.ValidateAndAssignDefaults(default_parameters)

        self.__objective = StandardizedConstraint(parameters, master_control, optimization_problem, required_buffer_size)

        self.__optimization_problem = optimization_problem

        super().__init__()

    def value(self, c:myVector, x: myVector, tol: float, save_data: bool = True) -> float:
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

        c[0] = value
        

    def applyJacobian(self, jv:myVector, v:myVector, x:myVector, tol:float, save_field: bool = True):
        '''
        Applies the Jacobian to the given vectors.

        This function computes the Jacobian-vector product and stores the result in `jv`.

            jv (myVector): The vector to store the result of the Jacobian-vector product.
            v (myVector): The vector to be multiplied by the Jacobian.
            x (myVector): The current state vector.
            tol (float): Tolerance for numerical computations. Is not used in this function.
            save_field (bool, optional): Flag indicating whether to save the field. Defaults to True.
        '''


        self.value(x, tol, False)  # Compute new primal if x has changed. Does nothing if x the same
        result = self.__objective.CalculateStandardizedGradient(save_field)
        gAux = result.Evaluate().reshape(-1) * self.__objective.GetScalingFactor()

        jv[0] = 0.0
        for i in range(len(gAux)):
            jv[0] += gAux[i] * v[i]

        def applyAdjointJacobian(self, ajv:myVector, v:myVector, x:myVector, tol:float, save_field: bool = True):
            """
            Applies the adjoint Jacobian to the given vectors.

            This function computes the adjoint Jacobian-vector product and stores the result in `ajv`.

            Args:
                ajv (myVector): The vector to store the result of the adjoint Jacobian-vector product.
                v (myVector): The vector to be multiplied by the adjoint Jacobian.
                x (myVector): The current state vector.
                tol (float): Tolerance for numerical computations.
                save_field (bool, optional): Flag indicating whether to save the field. Defaults to True.

            Returns:
                None
            """
            self.value(x, tol, False)  # Compute new primal if x has changed. Does nothing if x the same
            result = self.__objective.CalculateStandardizedGradient(save_field)
            gAux = result.Evaluate().reshape(-1) * self.__objective.GetScalingFactor()
            ajv[:] = [gAux[i] * v[0] for i in range(len(gAux))] 

    def Initialize(self):
        self.__objective.Initialize()

    def Check(self):
        self.__objective.Check()

    def GetMasterControl(self) -> MasterControl:
        return self.__objective.GetMasterControl()
                
    def Finalize(self):
        self.__objective.Finalize()

    def GetName(self):
        return self.__objective.GetName()