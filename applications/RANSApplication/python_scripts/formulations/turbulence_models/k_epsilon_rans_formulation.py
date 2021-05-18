# import kratos
import KratosMultiphysics as Kratos

# import RANS
import KratosMultiphysics.RANSApplication as KratosRANS

# import formulation interfaces
from KratosMultiphysics.RANSApplication.formulations.turbulence_models.scalar_turbulence_model_rans_formulation import ScalarTurbulenceModelRansFormulation
from KratosMultiphysics.RANSApplication.formulations.turbulence_models.two_equation_turbulence_model_rans_formulation import TwoEquationTurbulenceModelRansFormulation

class KEpsilonKRansFormulation(ScalarTurbulenceModelRansFormulation):
    def GetSolvingVariable(self):
        return KratosRANS.TURBULENT_KINETIC_ENERGY

    def GetElementNamePrefix(self):
        return "RansKEpsilonK"

    def GetConditionNamePrefix(self):
        return ""


class KEpsilonEpsilonRansFormulation(ScalarTurbulenceModelRansFormulation):
    def GetSolvingVariable(self):
        return KratosRANS.TURBULENT_ENERGY_DISSIPATION_RATE

    def GetElementNamePrefix(self):
        return "RansKEpsilonEpsilon"

    def GetConditionNamePrefix(self):
        return "RansKEpsilonEpsilon"


class KEpsilonRansFormulation(TwoEquationTurbulenceModelRansFormulation):
    def __init__(self, model_part, settings):
        default_settings = Kratos.Parameters(r'''
        {
            "formulation_name": "k_epsilon",
            "stabilization_method": "algebraic_flux_corrected",
            "turbulent_kinetic_energy_solver_settings": {},
            "turbulent_energy_dissipation_rate_solver_settings": {},
            "coupling_settings":
            {
                "relative_tolerance": 1e-3,
                "absolute_tolerance": 1e-5,
                "max_iterations": 10
            },
            "auxiliar_process_list": [],
            "echo_level": 0,
            "minimum_turbulent_viscosity": 1e-12
        }''')

        settings.ValidateAndAssignDefaults(default_settings)

        super().__init__(
            model_part,
            settings,
            KEpsilonKRansFormulation(model_part, settings["turbulent_kinetic_energy_solver_settings"]),
            KEpsilonEpsilonRansFormulation(model_part, settings["turbulent_energy_dissipation_rate_solver_settings"]))

    def AddVariables(self):
        self.GetBaseModelPart().AddNodalSolutionStepVariable(Kratos.DENSITY)
        self.GetBaseModelPart().AddNodalSolutionStepVariable(Kratos.VELOCITY)
        self.GetBaseModelPart().AddNodalSolutionStepVariable(Kratos.MESH_VELOCITY)
        self.GetBaseModelPart().AddNodalSolutionStepVariable(Kratos.NORMAL)
        self.GetBaseModelPart().AddNodalSolutionStepVariable(Kratos.VISCOSITY)
        self.GetBaseModelPart().AddNodalSolutionStepVariable(Kratos.TURBULENT_VISCOSITY)
        self.GetBaseModelPart().AddNodalSolutionStepVariable(KratosRANS.RANS_Y_PLUS)
        self.GetBaseModelPart().AddNodalSolutionStepVariable(KratosRANS.TURBULENT_KINETIC_ENERGY)
        self.GetBaseModelPart().AddNodalSolutionStepVariable(KratosRANS.TURBULENT_KINETIC_ENERGY_RATE)
        self.GetBaseModelPart().AddNodalSolutionStepVariable(KratosRANS.TURBULENT_ENERGY_DISSIPATION_RATE)
        self.GetBaseModelPart().AddNodalSolutionStepVariable(KratosRANS.TURBULENT_ENERGY_DISSIPATION_RATE_2)
        self.GetBaseModelPart().AddNodalSolutionStepVariable(KratosRANS.RANS_AUXILIARY_VARIABLE_1)
        self.GetBaseModelPart().AddNodalSolutionStepVariable(KratosRANS.RANS_AUXILIARY_VARIABLE_2)

        Kratos.Logger.PrintInfo(self.__class__.__name__, "Added solution step variables.")

    def AddDofs(self):
        Kratos.VariableUtils().AddDof(KratosRANS.TURBULENT_KINETIC_ENERGY, self.GetBaseModelPart())
        Kratos.VariableUtils().AddDof(KratosRANS.TURBULENT_ENERGY_DISSIPATION_RATE, self.GetBaseModelPart())

        Kratos.Logger.PrintInfo(self.__class__.__name__, "Added solution step dofs.")

    def Initialize(self):
        model_part = self.GetBaseModelPart()
        model = model_part.GetModel()

        process_info = model_part.ProcessInfo
        wall_model_part_name = process_info[KratosRANS.WALL_MODEL_PART_NAME]
        minimum_nut = self.GetParameters()["minimum_turbulent_viscosity"].GetDouble()

        nut_process = KratosRANS.RansNutKEpsilonUpdateProcess(
                                            model,
                                            self.GetBaseModelPart().Name,
                                            minimum_nut,
                                            self.echo_level)
        self.AddProcess(nut_process)

        nut_wall_process = KratosRANS.RansNutYPlusWallFunctionUpdateProcess(
                                            model,
                                            wall_model_part_name,
                                            minimum_nut,
                                            self.echo_level)
        self.AddProcess(nut_wall_process)

        super().Initialize()

    def SetConstants(self, settings):
        defaults = Kratos.Parameters('''{
            "von_karman"                                    : 0.41,
            "c_mu"                                          : 0.09,
            "c1"                                            : 1.44,
            "c2"                                            : 1.92,
            "sigma_k"                                       : 1.0,
            "sigma_epsilon"                                 : 1.3,
            "stabilizing_upwind_operator_coefficient"       : 1.2,
            "stabilizing_positivity_preserving_coefficient" : 1.2
        }''')

        settings.ValidateAndAssignDefaults(defaults)

        process_info = self.GetBaseModelPart().ProcessInfo
        process_info.SetValue(KratosRANS.VON_KARMAN, settings["von_karman"].GetDouble())
        process_info.SetValue(KratosRANS.TURBULENCE_RANS_C_MU, settings["c_mu"].GetDouble())
        process_info.SetValue(KratosRANS.TURBULENCE_RANS_C1, settings["c1"].GetDouble())
        process_info.SetValue(KratosRANS.TURBULENCE_RANS_C2, settings["c2"].GetDouble())
        process_info.SetValue(KratosRANS.TURBULENT_KINETIC_ENERGY_SIGMA, settings["sigma_k"].GetDouble())
        process_info.SetValue(KratosRANS.TURBULENT_ENERGY_DISSIPATION_RATE_SIGMA, settings["sigma_epsilon"].GetDouble())
        process_info.SetValue(KratosRANS.RANS_STABILIZATION_DISCRETE_UPWIND_OPERATOR_COEFFICIENT, settings["stabilizing_upwind_operator_coefficient"].GetDouble())
        process_info.SetValue(KratosRANS.RANS_STABILIZATION_DIAGONAL_POSITIVITY_PRESERVING_COEFFICIENT, settings["stabilizing_positivity_preserving_coefficient"].GetDouble())
