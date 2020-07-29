from __future__ import print_function, absolute_import, division

# import kratos
import KratosMultiphysics as Kratos

# import formulation interface
from KratosMultiphysics.RANSApplication.formulations.formulation import Formulation

# import formulations
from ..incompressible_potential_flow import IncompressiblePotentialFlowFormulation
from ..k_epsilon import KEpsilonFormulation
from .fractional_step_velocity_pressure_formulation import FractionalStepVelocityPressureFormulation

class FractionalStepKEpsilonFormulation(Formulation):
    def __init__(self, model_part, settings):
        super(FractionalStepKEpsilonFormulation, self).__init__(model_part, settings)

        default_settings = Kratos.Parameters(r'''
        {
            "formulation_name": "fractional_step_k_epsilon",
            "incompressible_potential_flow_initialization_settings": {},
            "fractional_step_flow_solver_settings": {},
            "k_epsilon_solver_settings": {},
            "max_iterations": 1
        }''')
        self.settings.ValidateAndAssignDefaults(default_settings)

        if (not self.settings["incompressible_potential_flow_initialization_settings"].IsEquivalentTo(
                Kratos.Parameters("{}"))):
            self.incompressible_potential_flow_formulation = IncompressiblePotentialFlowFormulation(model_part, settings["incompressible_potential_flow_initialization_settings"])
            self.AddFormulation(self.incompressible_potential_flow_formulation)

        self.fractional_step_formulation = FractionalStepVelocityPressureFormulation(model_part, settings["fractional_step_flow_solver_settings"])
        self.AddFormulation(self.fractional_step_formulation)

        self.k_epsilon_formulation = KEpsilonFormulation(model_part, settings["k_epsilon_solver_settings"])
        self.AddFormulation(self.k_epsilon_formulation)

        self.SetMaxCouplingIterations(self.settings["max_iterations"].GetInt())

    def SetConstants(self, settings):
        self.k_epsilon_formulation.SetConstants(settings)