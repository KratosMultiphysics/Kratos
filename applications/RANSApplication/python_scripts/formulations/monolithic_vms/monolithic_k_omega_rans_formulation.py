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
    def __init__(self, model_part, settings, deprecated_settings_dict):
        super().__init__(model_part, settings, deprecated_settings_dict)

        settings.ValidateAndAssignDefaults(self.GetDefaultParameters())

        if (not settings["incompressible_potential_flow_initialization_settings"].IsEquivalentTo(
                Kratos.Parameters("{}"))):
            self.incompressible_potential_flow_formulation = IncompressiblePotentialFlowRansFormulation(model_part, settings["incompressible_potential_flow_initialization_settings"], deprecated_settings_dict)
            self.AddRansFormulation(self.incompressible_potential_flow_formulation)

        self.monolithic_formulation = MonolithicVelocityPressureRansFormulation(model_part, settings["monolithic_flow_solver_settings"], deprecated_settings_dict)
        self.AddRansFormulation(self.monolithic_formulation)

        self.k_omega_formulation = KOmegaRansFormulation(model_part, settings["k_omega_solver_settings"], deprecated_settings_dict)
        self.AddRansFormulation(self.k_omega_formulation)

        self.SetMaxCouplingIterations(settings["max_iterations"].GetInt())

    def GetDefaultParameters(self):
        return Kratos.Parameters(r'''
        {
            "formulation_name": "monolithic_k_omega",
            "incompressible_potential_flow_initialization_settings": {},
            "monolithic_flow_solver_settings": {},
            "k_omega_solver_settings": {},
            "max_iterations": 1
        }''')

    def SetConstants(self, settings):
        self.k_omega_formulation.SetConstants(settings)

    def Initialize(self):
        # do not change the order of the initialization. This order is required
        # to add nut_nodal_update_process after nu_t update process. Otherwise
        # nut_nodal_update process (which is responsible for distributing element and condition gauss nut to nodes for old
        # vms and fractional step formulations) is added to the list of processes before its corresponding elemental or condition
        # nut update process.
        super().Initialize()

        if (self.monolithic_formulation.ElementHasNodalProperties()):
            # adding a process to distribute elemental and condition nut to nodes for fractional step
            nut_nodal_update_process = KratosRANS.RansNutNodalUpdateProcess(
                                                self.GetBaseModelPart().GetModel(),
                                                self.GetBaseModelPart().Name,
                                                self.k_omega_formulation.echo_level)
            self.k_omega_formulation.AddProcess(nut_nodal_update_process)

            # calling the execute initialize method here since, it is added after the routine super().Initialize()
            nut_nodal_update_process.ExecuteInitialize()