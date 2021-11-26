"""
An initial condition process for KratosMultiphysics
license: license.txt
"""

__all__ = ['Factory', 'ImposePerturbedInitialConditionProcessDeprecated']

import warnings

import KratosMultiphysics
import KratosMultiphysics.FluidDynamicsApplication
from KratosMultiphysics.FluidDynamicsApplication.kratos_perturbation_process import Parameters, ImposePerturbedInitialConditionProcess


def Factory(settings, Model):
    return ImposePerturbedInitialConditionProcessDeprecated(Model,Parameters(settings['Parameters']))


class ImposePerturbedInitialConditionProcessDeprecated(ImposePerturbedInitialConditionProcess):
    """
    Process generating divergence-free correlated noise, imposing no-penetrability condition on bluff bodies,
    preserving boundary conditions on boundaries and (if required) adding the contribution of a manually-loaded velocity field.
    The generated initial field reads
        u_0 = u_T + c*P(u_{cn}) + c*\nabla potential ,
    where u_T is the velocity field we load, c a penalty coefficient, and P(u_{cn}) and \nabla potential are discussed next.
    The process first projects the correlated noise u_{cn} into the model part. The projected correlated noise is denoted as P(u_{cn}).
    Then, the associated poisson problem for the unknown TEMPERATURE is solved:
        - \nabla \cdot \nabla TEMPERATURE = \nabla \cdot (c*P(u_{cn}) + u_{T})  on \Omega
        c*P(u_{cn}) \cdot NORMAL = - \nabla TEMPERATURE \cdot NORMAL  on \partial \Omega_{bluff body}
        TEMPERATURE = 0  on \Omega_{outlet}
    Then, the tangential velocity is set to zero on the boundaries \partial \Omega_{bluff body} and \partial \Omega_{external domain} by doing:
        P(u_{cn}) + \nabla TEMPERATURE = 0 on \Gamma
    The correlated noise satisfying the no-penetrability and solenoidal conditions and respecting the boundary conditions is: VELOCITY = u_T + c*P(u_{cn}) + c*\nabla TEMPERATURE.

    Parameters:
    model    : Kratos Model
    settings : Kratos Parameters
    """

    def __init__(self,model,settings):
        super().__init__(model,settings)
        warnings.warn(
                (
                    "\nThe kratos_perturbation_process of KratosMultiphysics.ExaquteSandboxApplication is deprecated "
                    "in favor of the kratos_perturbation_process of KratosMultiphysics.FluidDynamicsApplication. "
                    "\nIt will be removed in the future, tentatively, in July 2022."
                ),
                DeprecationWarning,
            )