import abc
import KratosMultiphysics as Kratos
from KratosMultiphysics.OptimizationApplication.utilities.union_utilities import SupportedSensitivityFieldVariableTypes
from KratosMultiphysics.SystemIdentificationApplication.utilities.expression_utils import ExpressionUnionType

class ResponseSensitivityInterface(abc.ABC):
    @abc.abstractmethod
    def SetResponseFunction(self, response_function: Kratos.AdjointResponseFunction) -> None:
        """Set the response function to compute the sensitivities using adjoint method.

        Args:
            response_function (Kratos.AdjointResponseFunction): Response function to be used to compute the sensitivities.
        """
        pass

    @abc.abstractmethod
    def GetResponseFunction(self) -> Kratos.AdjointResponseFunction:
        """Returns the response function used to compute the sensitivities.

        Returns:
            Kratos.AdjointResponseFunction: Adjoint response function used in adjoint based sensitivity computation.
        """
        pass

    @abc.abstractmethod
    def GetSensitivityModelPart(self) -> Kratos.ModelPart:
        """Returns the model part used on which the sensitivities are computed.

        This is the model part on which the sensitivities are computed, not the model part on which
        the adjoint solution is computed.

        Returns:
            Kratos.ModelPart: Sensitivity model part on which the sensitivities are computed.
        """
        pass

    @abc.abstractmethod
    def GetSensitivities(self) -> 'dict[SupportedSensitivityFieldVariableTypes, ExpressionUnionType]':
        """Returns sensitivities as a dict having variable and sensitivities as a pair.

        This method returns the computed sensitivities on the sensitivity model part as a pair (variable, expression).
        The variable will be the design variable. Not the variable which is used to store the sensitivities.

        Ex: The design variable may be "SHAPE", and the variable used to store the sensitivities is called "SHAPE_SENSITIVITY".
            In this case, it returns the variable "SHAPE" along with an expression containing the sensitivities.

        Returns:
            dict[SupportedSensitivityFieldVariableTypes, ExpressionUnionType]: Dictionary of design variable and corresponding sensitivities.
        """
        pass

