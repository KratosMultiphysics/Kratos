# import kratos
import KratosMultiphysics as Kratos
from KratosMultiphysics.kratos_utilities import IssueDeprecationWarning

# import RANS
import KratosMultiphysics.RANSApplication as KratosRANS

# import formulation interfaces
from KratosMultiphysics.RANSApplication.formulations.turbulence_models.scalar_turbulence_model_rans_formulation import ScalarTurbulenceModelRansFormulation
from KratosMultiphysics.RANSApplication.formulations.turbulence_models.two_equation_turbulence_model_rans_formulation import TwoEquationTurbulenceModelRansFormulation

# import utilities
from KratosMultiphysics.RANSApplication import RansWallDistanceCalculationProcess

class KOmegaSSTKRansFormulation(ScalarTurbulenceModelRansFormulation):
    def GetSolvingVariable(self):
        return KratosRANS.TURBULENT_KINETIC_ENERGY

    def GetElementNamePrefix(self):
        return "RansKOmegaSSTK"

    def GetConditionNamePrefix(self):
        return ""

    def SetWallFunctionSettings(self, settings=None):
        formulation_settings = self.GetParameters()["wall_function_settings"]
        if settings is not None:
            if not formulation_settings.IsEquivalentTo(Kratos.Parameters("""{}""")):
                Kratos.Logger.PrintWarning(self.__class__.__name__, "Global and specialized \"wall_function_settings\" are defined. Using specialized settings and global settings are discarded for this formulation.")
                settings = formulation_settings
            else:
                IssueDeprecationWarning(self.__class__.__name__, "Using deprecated global \"wall_function_settings\". Please define formulation specialized \"wall_function_settings\" in each leaf formulation.")
        else:
            settings = formulation_settings

        if (settings.Has("wall_function_region_type")):
            wall_function_region_type = settings["wall_function_region_type"].GetString()
        else:
            wall_function_region_type = "logarithmic_region_only"

        if (wall_function_region_type == "logarithmic_region_only"):
            self.condition_name = ""
        elif (wall_function_region_type == "viscous_region_only"):
            self.condition_name = "RansKEpsilonKVisBasedWall"
        else:
            msg = "Unsupported wall function region type provided. [ wall_function_region_type = \"" + wall_function_region_type + "\" ]."
            msg += "Supported wall function region types are:\n"
            msg += "\tlogarithmic_region_only\n"
            msg += "\tviscous_region_only\n"
            raise Exception(msg)

class KOmegaSSTOmegaRansFormulation(ScalarTurbulenceModelRansFormulation):
    def GetSolvingVariable(self):
        return KratosRANS.TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE

    def GetElementNamePrefix(self):
        return "RansKOmegaSSTOmega"

    def GetConditionNamePrefix(self):
        return "RansKOmegaSSTOmega"

    def SetWallFunctionSettings(self, settings=None):
        formulation_settings = self.GetParameters()["wall_function_settings"]
        if settings is not None:
            if not formulation_settings.IsEquivalentTo(Kratos.Parameters("""{}""")):
                Kratos.Logger.PrintWarning(self.__class__.__name__, "Global and specialized \"wall_function_settings\" are defined. Using specialized settings and global settings are discarded for this formulation.")
                settings = formulation_settings
            else:
                IssueDeprecationWarning(self.__class__.__name__, "Using deprecated global \"wall_function_settings\". Please define formulation specialized \"wall_function_settings\" in each leaf formulation.")
        else:
            settings = formulation_settings

        self.condition_name = self.GetConditionNamePrefix()

        if (self.condition_name != ""):
            if (settings.Has("wall_function_region_type")):
                wall_function_region_type = settings["wall_function_region_type"].GetString()
            else:
                wall_function_region_type = "logarithmic_region_only"

            if (settings.Has("wall_friction_velocity_calculation_method")):
                wall_friction_velocity_calculation_method = settings["wall_friction_velocity_calculation_method"].GetString()
            else:
                wall_friction_velocity_calculation_method = "velocity_based"

            if (wall_function_region_type == "logarithmic_region_only"):
                if (wall_friction_velocity_calculation_method == "velocity_based"):
                    self.condition_name = self.condition_name + "UBasedWall"
                elif (wall_friction_velocity_calculation_method ==
                    "turbulent_kinetic_energy_based"):
                    self.condition_name = self.condition_name + "KBasedWall"
                else:
                    msg = "Unsupported wall friction velocity calculation method. [ wall_friction_velocity_calculation_method = \"" + wall_friction_velocity_calculation_method + "\" ].\n"
                    msg += "Supported methods are:\n"
                    msg += "\tvelocity_based\n"
                    msg += "\tturbulent_kinetic_energy_based\n"
                    raise Exception(msg)
            elif (wall_function_region_type == "viscous_region_only"):
                self.condition_name = "RansKOmegaOmegaVisLogBasedWall"
            else:
                msg = "Unsupported wall function region type provided. [ wall_function_region_type = \"" + wall_function_region_type + "\" ]."
                msg += "Supported wall function region types are:\n"
                msg += "\tlogarithmic_region_only\n"
                msg += "\tviscous_region_only\n"
                raise Exception(msg)


class KOmegaSSTRansFormulation(TwoEquationTurbulenceModelRansFormulation):
    def __init__(self, model_part, settings, deprecated_settings_dict):
        settings.ValidateAndAssignDefaults(self.GetDefaultParameters())

        super().__init__(
            model_part,
            settings,
            deprecated_settings_dict,
            KOmegaSSTKRansFormulation(model_part, settings["turbulent_kinetic_energy_solver_settings"], deprecated_settings_dict),
            KOmegaSSTOmegaRansFormulation(model_part, settings["turbulent_specific_energy_dissipation_rate_solver_settings"], deprecated_settings_dict))

    def GetDefaultParameters(self):
        return Kratos.Parameters(r'''
        {
            "formulation_name": "k_omega_sst",
            "stabilization_method": "algebraic_flux_corrected",
            "turbulent_kinetic_energy_solver_settings": {},
            "turbulent_specific_energy_dissipation_rate_solver_settings": {},
            "wall_distance_calculation_settings":
            {
                "max_levels"                           : 100,
                "max_distance"                         : 1e+30,
                "echo_level"                           : 0,
                "distance_variable_name"               : "DISTANCE",
                "nodal_area_variable_name"             : "NODAL_AREA",
                "wall_distance_update_execution_points": ["initialize"]
            },
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

    def AddVariables(self):
        self.GetBaseModelPart().AddNodalSolutionStepVariable(Kratos.VELOCITY)
        self.GetBaseModelPart().AddNodalSolutionStepVariable(Kratos.MESH_VELOCITY)
        self.GetBaseModelPart().AddNodalSolutionStepVariable(Kratos.NORMAL)
        self.GetBaseModelPart().AddNodalSolutionStepVariable(KratosRANS.RANS_Y_PLUS)
        self.GetBaseModelPart().AddNodalSolutionStepVariable(KratosRANS.TURBULENT_KINETIC_ENERGY)
        self.GetBaseModelPart().AddNodalSolutionStepVariable(KratosRANS.TURBULENT_KINETIC_ENERGY_RATE)
        self.GetBaseModelPart().AddNodalSolutionStepVariable(KratosRANS.TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE)
        self.GetBaseModelPart().AddNodalSolutionStepVariable(KratosRANS.TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE_2)
        self.GetBaseModelPart().AddNodalSolutionStepVariable(KratosRANS.RANS_AUXILIARY_VARIABLE_1)
        self.GetBaseModelPart().AddNodalSolutionStepVariable(KratosRANS.RANS_AUXILIARY_VARIABLE_2)

	    # additional variables required for wall distance calculation
        self.GetBaseModelPart().AddNodalSolutionStepVariable(Kratos.DISTANCE)
        self.GetBaseModelPart().AddNodalSolutionStepVariable(Kratos.FLAG_VARIABLE)
        self.GetBaseModelPart().AddNodalSolutionStepVariable(Kratos.NODAL_AREA)

        Kratos.Logger.PrintInfo(self.__class__.__name__, "Added solution step variables.")

    def AddDofs(self):
        Kratos.VariableUtils().AddDof(KratosRANS.TURBULENT_KINETIC_ENERGY, self.GetBaseModelPart())
        Kratos.VariableUtils().AddDof(KratosRANS.TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE, self.GetBaseModelPart())

        Kratos.Logger.PrintInfo(self.__class__.__name__, "Added solution step dofs.")

    def Initialize(self):
        settings = self.GetParameters()
        add_wall_distance_calculation_process = True
        if settings["wall_distance_calculation_settings"].Has("wall_distance_update_execution_points"):
            add_wall_distance_calculation_process = (len(settings["wall_distance_calculation_settings"]["wall_distance_update_execution_points"].GetStringArray()) > 0)

        if add_wall_distance_calculation_process:
            model_part = self.GetBaseModelPart()
            model = model_part.GetModel()
            process_info = model_part.ProcessInfo
            wall_model_part_name = process_info[KratosRANS.WALL_MODEL_PART_NAME]

            wall_distance_calculation_settings = settings["wall_distance_calculation_settings"]
            wall_distance_calculation_settings.AddEmptyValue("main_model_part_name")
            wall_distance_calculation_settings["main_model_part_name"].SetString(self.GetBaseModelPart().Name)
            wall_distance_calculation_settings.AddEmptyValue("wall_model_part_name")
            wall_distance_calculation_settings["wall_model_part_name"].SetString(wall_model_part_name)

            wall_distance_process = RansWallDistanceCalculationProcess(model, wall_distance_calculation_settings)
            self.AddProcess(wall_distance_process)


        super().Initialize()

    def SetConstants(self, settings):
        defaults = Kratos.Parameters('''{
            "wall_law_constants":{
                "kappa"            : 0.41,
                "c_mu"             : 0.09,
                "correction_factor": 1.0
            },
            "k_omega_constants": {
                "sigma_k"    : 0.85,
                "sigma_omega": 0.5,
                "beta"       : 0.0750,
                "a1"         : 0.31
            },
            "k_epsilon_constants": {
                "sigma_k"    : 1.0,
                "sigma_omega": 0.856,
                "beta"       : 0.0828
            },
            "stabilization_constants":{
                "upwind_operator_coefficient"       : 1.2,
                "positivity_preserving_coefficient" : 1.2
            }
        }''')

        settings.RecursivelyValidateAndAssignDefaults(defaults)

        process_info = self.GetBaseModelPart().ProcessInfo
        # wall law constants
        constants = settings["wall_law_constants"]
        process_info.SetValue(KratosRANS.VON_KARMAN, constants["kappa"].GetDouble())
        process_info.SetValue(KratosRANS.TURBULENCE_RANS_C_MU, constants["c_mu"].GetDouble())
        process_info.SetValue(KratosRANS.WALL_CORRECTION_FACTOR, constants["correction_factor"].GetDouble())

        # k-omega constants
        constants = settings["k_omega_constants"]
        process_info.SetValue(KratosRANS.TURBULENT_KINETIC_ENERGY_SIGMA_1, constants["sigma_k"].GetDouble())
        process_info.SetValue(KratosRANS.TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE_SIGMA_1, constants["sigma_omega"].GetDouble())
        process_info.SetValue(KratosRANS.TURBULENCE_RANS_BETA_1, constants["beta"].GetDouble())
        process_info.SetValue(KratosRANS.TURBULENCE_RANS_A1, constants["a1"].GetDouble())

        # k-epsilon constants
        constants = settings["k_epsilon_constants"]
        process_info.SetValue(KratosRANS.TURBULENT_KINETIC_ENERGY_SIGMA_2, constants["sigma_k"].GetDouble())
        process_info.SetValue(KratosRANS.TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE_SIGMA_2, constants["sigma_omega"].GetDouble())
        process_info.SetValue(KratosRANS.TURBULENCE_RANS_BETA_2, constants["beta"].GetDouble())

        # stabilization parameters
        constants = settings["stabilization_constants"]
        process_info.SetValue(KratosRANS.RANS_STABILIZATION_DISCRETE_UPWIND_OPERATOR_COEFFICIENT, constants["upwind_operator_coefficient"].GetDouble())
        process_info.SetValue(KratosRANS.RANS_STABILIZATION_DIAGONAL_POSITIVITY_PRESERVING_COEFFICIENT, constants["positivity_preserving_coefficient"].GetDouble())
