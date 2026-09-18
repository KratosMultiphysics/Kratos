import abc
import KratosMultiphysics as Kratos

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

