# import kratos
import KratosMultiphysics as Kratos
import KratosMultiphysics.RANSApplication as KratosRANS

# import formulation interface
from KratosMultiphysics.RANSApplication.formulations.rans_formulation import RansFormulation

# import formulations
from KratosMultiphysics.RANSApplication.formulations.incompressible_potential_flow import IncompressiblePotentialFlowRansFormulation
from KratosMultiphysics.RANSApplication.formulations.turbulence_models.k_omega_sst_rans_formulation import KOmegaSSTRansFormulation
from KratosMultiphysics.RANSApplication.formulations.fractional_step.fractional_step_velocity_pressure_rans_formulation import FractionalStepVelocityPressureRansFormulation

class FractionalStepKOmegaSSTRansFormulation(RansFormulation):
    def __init__(self, model_part, settings):
        super().__init__(model_part, settings)

        default_settings = Kratos.Parameters(r'''
        {
            "formulation_name": "fractional_step_k_epsilon",
            "incompressible_potential_flow_initialization_settings": {},
            "fractional_step_flow_solver_settings": {},
            "k_omega_sst_solver_settings": {},
            "max_iterations": 1
        }''')
        settings.ValidateAndAssignDefaults(default_settings)

        if (not settings["incompressible_potential_flow_initialization_settings"].IsEquivalentTo(
                Kratos.Parameters("{}"))):
            self.incompressible_potential_flow_formulation = IncompressiblePotentialFlowRansFormulation(model_part, settings["incompressible_potential_flow_initialization_settings"])
            self.AddRansFormulation(self.incompressible_potential_flow_formulation)

        self.fractional_step_formulation = FractionalStepVelocityPressureRansFormulation(model_part, settings["fractional_step_flow_solver_settings"])
        self.AddRansFormulation(self.fractional_step_formulation)

        self.k_omega_sst_formulation = KOmegaSSTRansFormulation(model_part, settings["k_omega_sst_solver_settings"])
        self.AddRansFormulation(self.k_omega_sst_formulation)

        self.SetMaxCouplingIterations(settings["max_iterations"].GetInt())

    def SetConstants(self, settings):
        self.k_omega_sst_formulation.SetConstants(settings)

    def Initialize(self):
        super().Initialize()

        nut_nodal_update_process = KratosRANS.RansNutNodalUpdateProcess(
                                            self.GetBaseModelPart().GetModel(),
                                            self.GetBaseModelPart().Name,
                                            self.k_omega_sst_formulation.echo_level)
        self.k_omega_sst_formulation.AddProcess(nut_nodal_update_process)