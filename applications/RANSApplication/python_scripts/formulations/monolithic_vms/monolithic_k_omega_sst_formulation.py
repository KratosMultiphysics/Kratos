from __future__ import print_function, absolute_import, division

# import kratos
import KratosMultiphysics as Kratos

# import formulation interface
from KratosMultiphysics.RANSApplication.formulations.formulation import Formulation

# import formulations
from ..incompressible_potential_flow import IncompressiblePotentialFlowFormulation
from ..k_omega_sst import KOmegaSSTFormulation
from .monolithic_velocity_pressure_formulation import MonolithicVelocityPressureFormulation

class MonolithicKOmegaSSTFormulation(Formulation):
    def __init__(self, model_part, settings):
        super(MonolithicKOmegaSSTFormulation, self).__init__(model_part, settings)

        default_settings = Kratos.Parameters(r'''
        {
            "formulation_name": "monolithic_k_omega_sst",
            "incompressible_potential_flow_initialization_settings": {},
            "monolithic_flow_solver_settings": {},
            "k_omega_sst_solver_settings": {},
            "max_iterations": 1
        }''')
        self.settings.ValidateAndAssignDefaults(default_settings)

        if (not self.settings["incompressible_potential_flow_initialization_settings"].IsEquivalentTo(
                Kratos.Parameters("{}"))):
            self.incompressible_potential_flow_formulation = IncompressiblePotentialFlowFormulation(model_part, settings["incompressible_potential_flow_initialization_settings"])
            self.AddFormulation(self.incompressible_potential_flow_formulation)

        self.monolithic_formulation = MonolithicVelocityPressureFormulation(model_part, settings["monolithic_flow_solver_settings"])
        self.AddFormulation(self.monolithic_formulation)

        self.k_omega_sst_formulation = KOmegaSSTFormulation(model_part, settings["k_omega_sst_solver_settings"])
        self.AddFormulation(self.k_omega_sst_formulation)

        self.SetMaxCouplingIterations(self.settings["max_iterations"].GetInt())

    def SetConstants(self, settings):
        self.k_omega_sst_formulation.SetConstants(settings)
