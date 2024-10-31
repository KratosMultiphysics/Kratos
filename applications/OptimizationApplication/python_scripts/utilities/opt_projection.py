# Temporal fix until to avoid list[type] related errors until we move to 3.9 and above.
from __future__ import annotations

from typing import Optional
from abc import ABC, abstractmethod

import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.utilities.union_utilities import ContainerExpressionTypes

class DesignVariableProjection(ABC):
    """Design variable projection methods to convert given x space values to given y space values

    This is the interface class which can be used to implement design variable projection classes
    which project forward and backward given x space values to y space.

    The projection spaces should be scalar spaces.
    """

    @abstractmethod
    def SetProjectionSpaces(self, x_space_values: 'list[float]', y_space_values: 'list[float]') -> None:
        """Space values in x space and y space

        Args:
            x_space_values (list[float]): List of values to be used in x space
            y_space_values (list[float]): List of values to be used in y space
        """
        pass

    @abstractmethod
    def ProjectForward(self, x_values: ContainerExpressionTypes) -> ContainerExpressionTypes:
        """Project x space values to y space

        Args:
            x_values (ContainerExpressionTypes): x space values.

        Returns:
            ContainerExpressionTypes: y space values.
        """
        pass

    @abstractmethod
    def ProjectBackward(self, y_values: ContainerExpressionTypes) -> ContainerExpressionTypes:
        """Project y space values to x space values

        Args:
            y_values (ContainerExpressionTypes): y space values.

        Returns:
            ContainerExpressionTypes: x space values.
        """
        pass

    @abstractmethod
    def ForwardProjectionGradient(self, x_values: ContainerExpressionTypes) -> ContainerExpressionTypes:
        """Computes the forward projected y space values' gradient w.r.t x space.

        Args:
            x_values (ContainerExpressionTypes): x space values.

        Returns:
            ContainerExpressionTypes: gradient of the y space values w.r.t. x space.
        """
        pass

    @abstractmethod
    def Update(self) -> None:
        """Updates the projection method
        """

class IdentityDesignVariableProjection(DesignVariableProjection):
    def __init__(self, parameters: Kratos.Parameters, _):
        default_settings = Kratos.Parameters("""{
            "type"       : "identity_projection",
            "extrapolate": "false"
        }""")
        parameters.ValidateAndAssignDefaults(default_settings)
        pass

    def SetProjectionSpaces(self, _: 'list[float]', __: 'list[float]') -> None:
        # not using any of the spaces provided.
        pass

    def ProjectForward(self, x_values: ContainerExpressionTypes) -> ContainerExpressionTypes:
        # returns the x_values since x and y are the same for IdentityProjection
        return x_values.Clone()

    def ProjectBackward(self, y_values: ContainerExpressionTypes) -> ContainerExpressionTypes:
        # return the y_values since x and y are the same for Identity projection
        return y_values.Clone()

    def ForwardProjectionGradient(self, x_values: ContainerExpressionTypes) -> ContainerExpressionTypes:
        result = x_values.Clone()
        # x = y, hence the gradients are 1.0
        Kratos.Expression.LiteralExpressionIO.SetData(result, 1.0)
        return result

    def Update(self) -> None:
        pass

class SigmoidalDesignVariableProjection(DesignVariableProjection):
    def __init__(self, parameters: Kratos.Parameters, _):
        default_parameters = Kratos.Parameters("""{
            "type"          : "sigmoidal_projection",
            "beta_value"    : 5,
            "penalty_factor": 1
        }""")

        parameters.ValidateAndAssignDefaults(default_parameters)
        self.beta = parameters["beta_value"].GetDouble()
        self.penalty_factor = parameters["penalty_factor"].GetInt()

        self.x_space_values: 'Optional[list[float]]' = None
        self.y_space_values: 'Optional[list[float]]' = None

    def SetProjectionSpaces(self, x_space_values: 'list[float]', y_space_values: 'list[float]') -> None:
        self.x_space_values = x_space_values
        self.y_space_values = y_space_values

    def ProjectForward(self, x_values: ContainerExpressionTypes) -> ContainerExpressionTypes:
        return KratosOA.ControlUtils.SigmoidalProjectionUtils.ProjectForward(x_values, self.x_space_values, self.y_space_values, self.beta, self.penalty_factor)

    def ProjectBackward(self, y_values: ContainerExpressionTypes) -> ContainerExpressionTypes:
        return KratosOA.ControlUtils.SigmoidalProjectionUtils.ProjectBackward(y_values, self.x_space_values, self.y_space_values, self.beta, self.penalty_factor)

    def ForwardProjectionGradient(self, x_values: ContainerExpressionTypes) -> ContainerExpressionTypes:
        return KratosOA.ControlUtils.SigmoidalProjectionUtils.CalculateForwardProjectionGradient(x_values, self.x_space_values, self.y_space_values, self.beta, self.penalty_factor)

    def Update(self) -> None:
        pass

class AdaptiveSigmoidalDesignVariableProjection(DesignVariableProjection):
    def __init__(self, parameters: Kratos.Parameters, optimization_problem: OptimizationProblem):
        default_parameters = Kratos.Parameters("""{
            "type"          : "adaptive_sigmoidal_projection",
            "initial_value" : 5,
            "max_value"     : 30,
            "increase_fac"  : 1.05,
            "update_period" : 50,
            "penalty_factor": 1
        }""")

        parameters.ValidateAndAssignDefaults(default_parameters)
        self.beta = parameters["initial_value"].GetDouble()
        self.max_beta = parameters["max_value"].GetDouble()
        self.increase_fac = parameters["increase_fac"].GetDouble()
        self.update_period = parameters["update_period"].GetInt()
        self.penalty_factor = parameters["penalty_factor"].GetInt()
        self.beta_computed_step = 1

        self.x_space_values: 'Optional[list[float]]' = None
        self.y_space_values: 'Optional[list[float]]' = None

        self.optimization_problem = optimization_problem

    def SetProjectionSpaces(self, x_space_values: 'list[float]', y_space_values: 'list[float]') -> None:
        self.x_space_values = x_space_values
        self.y_space_values = y_space_values

    def ProjectForward(self, x_values: ContainerExpressionTypes) -> ContainerExpressionTypes:
        return KratosOA.ControlUtils.SigmoidalProjectionUtils.ProjectForward(x_values, self.x_space_values, self.y_space_values, self.beta, self.penalty_factor)

    def ProjectBackward(self, y_values: ContainerExpressionTypes) -> ContainerExpressionTypes:
        return KratosOA.ControlUtils.SigmoidalProjectionUtils.ProjectBackward(y_values, self.x_space_values, self.y_space_values, self.beta, self.penalty_factor)

    def ForwardProjectionGradient(self, x_values: ContainerExpressionTypes) -> ContainerExpressionTypes:
        return KratosOA.ControlUtils.SigmoidalProjectionUtils.CalculateForwardProjectionGradient(x_values, self.x_space_values, self.y_space_values, self.beta, self.penalty_factor)

    def Update(self) -> None:
        step = self.optimization_problem.GetStep()
        if step % self.update_period == 0 and self.beta_computed_step != step:
            self.beta_computed_step = step
            self.beta = min(self.beta * self.increase_fac, self.max_beta)
            Kratos.Logger.PrintInfo(self.__class__.__name__, f"Increased beta to {self.beta}.")


def CreateProjection(parameters: Kratos.Parameters, optimization_problem: OptimizationProblem) -> DesignVariableProjection:
    if not parameters.Has("type"):
        raise RuntimeError("DesignVariableProjection \"type\" is not present in following parameters:\n" + str(parameters))

    projection_type = parameters["type"].GetString()

    projection_types_map: 'dict[str, type[DesignVariableProjection]]' = {
        "identity_projection"          : IdentityDesignVariableProjection,
        "sigmoidal_projection"         : SigmoidalDesignVariableProjection,
        "adaptive_sigmoidal_projection": AdaptiveSigmoidalDesignVariableProjection
    }

    if projection_type in projection_types_map.keys():
        return projection_types_map[projection_type](parameters, optimization_problem)
    else:
        raise RuntimeError(f"Unsupported projected type = \"{projection_type}\" requested. Followings are supported:\n\t" + "\n\t".join(list(projection_types_map.keys())))
