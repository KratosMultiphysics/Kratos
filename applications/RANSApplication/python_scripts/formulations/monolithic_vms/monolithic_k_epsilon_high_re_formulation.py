from __future__ import print_function, absolute_import, division

# import kratos
import KratosMultiphysics as Kratos

# import RANS
import KratosMultiphysics.RANSApplication as KratosRANS

# import formulation interface
from KratosMultiphysics.RANSApplication.formulations.formulation import Formulation

# import formulations
from ..incompressible_potential_flow import IncompressiblePotentialFlowFormulation
from ..k_epsilon_high_re import KEpsilonHighReFormulation
from .monolithic_velocity_pressure_formulation import MonolithicVelocityPressureFormulation

class MonolithicKEpsilonHighReFormulation(Formulation):
    def __init__(self, model_part, settings):
        super(MonolithicKEpsilonHighReFormulation, self).__init__(model_part, settings)

        default_settings = Kratos.Parameters(r'''
        {
            "formulation_name": "monolithic_k_epsilon_high_re",
            "incompressible_potential_flow_initialization_settings": {},
            "monolithic_flow_solver_settings": {},
            "k_epsilon_high_re_solver_settings": {},
            "max_iterations": 1
        }''')
        self.settings.ValidateAndAssignDefaults(default_settings)

        if (not self.settings["incompressible_potential_flow_initialization_settings"].IsEquivalentTo(
                Kratos.Parameters("{}"))):
            self.incompressible_potential_flow_formulation = IncompressiblePotentialFlowFormulation(model_part, settings["incompressible_potential_flow_initialization_settings"])
            self.AddFormulation(self.incompressible_potential_flow_formulation)

        self.monolithic_formulation = MonolithicVelocityPressureFormulation(model_part, settings["monolithic_flow_solver_settings"])
        self.monolithic_formulation.SetConditionName("RansVMSMonolithicKBasedWallCondition")
        self.AddFormulation(self.monolithic_formulation)

        self.k_epsilon_formulation = KEpsilonHighReFormulation(model_part, settings["k_epsilon_high_re_solver_settings"])
        self.AddFormulation(self.k_epsilon_formulation)

        self.SetMaxCouplingIterations(self.settings["max_iterations"].GetInt())

    def SetConstants(self, settings):
        self.k_epsilon_formulation.SetConstants(settings)
