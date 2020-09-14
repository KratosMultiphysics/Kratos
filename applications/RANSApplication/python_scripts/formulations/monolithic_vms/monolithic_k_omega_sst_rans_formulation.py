from __future__ import print_function, absolute_import, division

# import kratos
import KratosMultiphysics as Kratos

# import formulation interface
from KratosMultiphysics.RANSApplication.formulations.rans_formulation import RansFormulation

# import rans_formulations
from ..incompressible_potential_flow import IncompressiblePotentialFlowRansFormulation
from ..k_omega_sst import KOmegaSSTRansFormulation
from .monolithic_velocity_pressure_rans_formulation import MonolithicVelocityPressureRansFormulation

class MonolithicKOmegaSSTRansFormulation(RansFormulation):
    def __init__(self, model_part, settings):
        super().__init__(model_part, settings)

        default_settings = Kratos.Parameters(r'''
        {
            "formulation_name": "monolithic_k_omega_sst",
            "incompressible_potential_flow_initialization_settings": {},
            "monolithic_flow_solver_settings": {},
            "k_omega_sst_solver_settings": {},
            "max_iterations": 1
        }''')
        settings.ValidateAndAssignDefaults(default_settings)

        if (not settings["incompressible_potential_flow_initialization_settings"].IsEquivalentTo(
                Kratos.Parameters("{}"))):
            self.incompressible_potential_flow_formulation = IncompressiblePotentialFlowRansFormulation(model_part, settings["incompressible_potential_flow_initialization_settings"])
            self.AddRansFormulation(self.incompressible_potential_flow_formulation)

        self.monolithic_formulation = MonolithicVelocityPressureRansFormulation(model_part, settings["monolithic_flow_solver_settings"])
        self.AddRansFormulation(self.monolithic_formulation)

        self.k_omega_sst_formulation = KOmegaSSTRansFormulation(model_part, settings["k_omega_sst_solver_settings"])
        self.AddRansFormulation(self.k_omega_sst_formulation)

        self.SetMaxCouplingIterations(settings["max_iterations"].GetInt())

    def SetConstants(self, settings):
        self.k_omega_sst_formulation.SetConstants(settings)
