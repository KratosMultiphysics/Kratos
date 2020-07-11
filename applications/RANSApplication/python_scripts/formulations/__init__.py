__all__ = ["Factory"]

import KratosMultiphysics as Kratos

# flow solver formulations
from .incompressible_potential_flow import IncompressiblePotentialFlowFormulation
from .monolithic_vms.monolithic_velocity_pressure_formulation import MonolithicVelocityPressureFormulation
from .fractional_step.fractional_step_velocity_pressure_formulation import FractionalStepVelocityPressureFormulation

# turbulence model formulations
### k-epsilon formulations
from .monolithic_vms.monolithic_k_epsilon_formulation import MonolithicKEpsilonFormulation
from .fractional_step.fractional_step_k_epsilon_formulation import FractionalStepKEpsilonFormulation

### k-omega formulations
from .monolithic_vms.monolithic_k_omega_formulation import MonolithicKOmegaFormulation
from .fractional_step.fractional_step_k_omega_formulation import FractionalStepKOmegaFormulation

def Factory(model_part, settings):
    formulation_name = settings["formulation_name"].GetString()
    formulations_list = [
        ["incompressible_potential_flow", IncompressiblePotentialFlowFormulation],
        ["monolithic", MonolithicVelocityPressureFormulation],
        ["monolithic_k_epsilon", MonolithicKEpsilonFormulation],
        ["monolithic_k_omega", MonolithicKOmegaFormulation],
        ["fractional_step", FractionalStepVelocityPressureFormulation],
        ["fractional_step_k_epsilon", FractionalStepKEpsilonFormulation],
        ["fractional_step_k_omega", FractionalStepKOmegaFormulation]
    ]

    formulation_names_list = [
        formulations_list[i][0] for i in range(len(formulations_list))
    ]
    formulation_list = [
        formulations_list[i][1] for i in range(len(formulations_list))
    ]

    if (formulation_name not in formulation_names_list):
        msg = "Unknown formulation_name = \"" + formulation_name + "\". \nFollowing formulations are allowed:\n    "
        msg += "\n    ".join(sorted(formulation_names_list))
        raise Exception(msg + "\n")

    current_formulation = formulation_list[formulation_names_list.index(
        formulation_name)](model_part, settings)

    Kratos.Logger.PrintInfo("RANSFormulationFactory",
                            "Created " + formulation_name + " formulation.")

    return current_formulation
