import KratosMultiphysics as Kratos

# import formulation interface
from KratosMultiphysics.RANSApplication.formulations.scalar_rans_formulation import ScalarRansFormulation

# import utilities
from KratosMultiphysics.RANSApplication.formulations.utilities import CalculateNormalsOnConditions
from KratosMultiphysics.RANSApplication.formulations.utilities import InitializeYPlusVariablesInConditions
from KratosMultiphysics.kratos_utilities import IssueDeprecationWarning

class ScalarTurbulenceModelRansFormulation(ScalarRansFormulation):
    def __init__(self, model_part, settings):
        """Scalar turbulence model formulation base class

        This solves the variable given in self.GetSolvingVariable(), using element and conditions
        having prefixes provided by self.GetElementNamePrefix() and self.GetConditionNamePrefix()

        If wall functions are used, then self.GetConditionNamePrefix() should return non-empty prefix
        otherwise it should be empty.

        Args:
            model_part (Kratos.ModelPart): ModelPart to be used in the formulation.
            settings (Kratos.Parameters): Settings to be used in the formulation.
        """

        super().__init__(model_part, settings)

    def Initialize(self):
        InitializeYPlusVariablesInConditions(self.GetBaseModelPart())
        CalculateNormalsOnConditions(self.GetBaseModelPart())

        super().Initialize()

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
            else:
                msg = "Unsupported wall function region type provided. [ wall_function_region_type = \"" + wall_function_region_type + "\" ]."
                msg += "Supported wall function region types are:\n"
                msg += "\tlogarithmic_region_only\n"
                raise Exception(msg)
