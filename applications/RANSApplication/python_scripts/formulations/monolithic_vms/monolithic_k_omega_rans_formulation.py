# import kratos
import KratosMultiphysics as Kratos
import KratosMultiphysics.RANSApplication as KratosRANS

# import formulation interface
from KratosMultiphysics.RANSApplication.formulations.rans_formulation import RansFormulation

# import formulations
from KratosMultiphysics.RANSApplication.formulations.incompressible_potential_flow import IncompressiblePotentialFlowRansFormulation
from KratosMultiphysics.RANSApplication.formulations.turbulence_models.k_omega_rans_formulation import KOmegaRansFormulation
from KratosMultiphysics.RANSApplication.formulations.monolithic_vms.monolithic_velocity_pressure_rans_formulation import MonolithicVelocityPressureRansFormulation

class MonolithicKOmegaRansFormulation(RansFormulation):
    def __init__(self, model_part, settings):
        super().__init__(model_part, settings)

        default_settings = Kratos.Parameters(r'''
        {
            "formulation_name": "monolithic_k_omega",
            "incompressible_potential_flow_initialization_settings": {},
            "monolithic_flow_solver_settings": {},
            "k_omega_solver_settings": {},
            "max_iterations": 1
        }''')

        settings.ValidateAndAssignDefaults(default_settings)

        if (not settings["incompressible_potential_flow_initialization_settings"].IsEquivalentTo(
                Kratos.Parameters("{}"))):
            self.incompressible_potential_flow_formulation = IncompressiblePotentialFlowRansFormulation(model_part, settings["incompressible_potential_flow_initialization_settings"])
            self.AddRansFormulation(self.incompressible_potential_flow_formulation)

        self.monolithic_formulation = MonolithicVelocityPressureRansFormulation(model_part, settings["monolithic_flow_solver_settings"])
        self.AddRansFormulation(self.monolithic_formulation)

        self.k_omega_formulation = KOmegaRansFormulation(model_part, settings["k_omega_solver_settings"])
        self.AddRansFormulation(self.k_omega_formulation)

        self.SetMaxCouplingIterations(settings["max_iterations"].GetInt())

    def SetConstants(self, settings):
        self.k_omega_formulation.SetConstants(settings)

    def Initialize(self):
        super().Initialize()

        if (self.monolithic_formulation.ElementHasNodalProperties()):
            nut_nodal_update_process = KratosRANS.RansNutNodalUpdateProcess(
                                                self.GetBaseModelPart().GetModel(),
                                                self.GetBaseModelPart().Name,
                                                self.k_omega_formulation.echo_level)
            self.k_omega_formulation.AddProcess(nut_nodal_update_process)