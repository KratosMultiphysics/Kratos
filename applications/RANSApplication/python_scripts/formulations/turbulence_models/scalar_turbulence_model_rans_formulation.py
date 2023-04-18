from abc import abstractmethod

# import kratos
import KratosMultiphysics as Kratos
import KratosMultiphysics.FluidDynamicsApplication as KratosCFD

# import formulation interface
from KratosMultiphysics.RANSApplication.formulations.rans_formulation import RansFormulation

# import utilities
from KratosMultiphysics.RANSApplication.formulations.utilities import CreateRansFormulationModelPart
from KratosMultiphysics.RANSApplication.formulations.utilities import CreateBlockBuilderAndSolver
from KratosMultiphysics.RANSApplication.formulations.utilities import InitializePeriodicConditions
from KratosMultiphysics.RANSApplication.formulations.utilities import GetBoundaryFlags
from KratosMultiphysics.RANSApplication.formulations.utilities import GetKratosObjectPrototype
from KratosMultiphysics.RANSApplication.formulations.utilities import CalculateNormalsOnConditions
from KratosMultiphysics.RANSApplication.formulations.utilities import InitializeYPlusVariablesInConditions

class ScalarTurbulenceModelRansFormulation(RansFormulation):
    def __init__(self, model_part, settings, deprecated_settings_dict):
        """Scalar turbulence model formulation base class

        This solves the variable given in self.GetSolvingVariable(), using element and conditions
        having prefixes provided by self.GetElementNamePrefix() and self.GetConditionNamePrefix()

        If wall functions are used, then self.GetConditionNamePrefix() should return non-empty prefix
        otherwise it should be empty.

        Args:
            model_part (Kratos.ModelPart): ModelPart to be used in the formulation.
            settings (Kratos.Parameters): Settings to be used in the formulation.
        """
        self.BackwardCompatibilityHelper(settings, deprecated_settings_dict)
        settings.ValidateAndAssignDefaults(self.GetDefaultParameters())
        self.echo_level = settings["echo_level"].GetInt()

        super().__init__(model_part, settings, deprecated_settings_dict)

    def GetDefaultParameters(self):
        return Kratos.Parameters(r"""{
            "relative_tolerance"    : 1e-3,
            "absolute_tolerance"    : 1e-5,
            "max_iterations"        : 200,
            "relaxation_factor"     : 0.5,
            "echo_level"            : 2,
            "linear_solver_settings": {
                "solver_type"  : "amgcl"
            },
            "boundary_flags": ["INLET", "STRUCTURE"],
            "wall_function_settings": {
                "wall_function_region_type": "logarithmic_region_only",
                "wall_friction_velocity_calculation_method": "velocity_based"
            }
        }""")

    def BackwardCompatibilityHelper(self, settings, deprecated_settings_dict):
        if "wall_function_settings" in deprecated_settings_dict.keys():
            if settings.Has("wall_function_settings"):
                Kratos.Logger.PrintWarning(self.__class__.__name__, "Found \"wall_function_settings\" in deprecated settings as well as in formulation settings. Continuing with formulation based settings.")
            else:
                settings.AddValue("wall_function_settings", deprecated_settings_dict["wall_function_settings"].Clone())

    @abstractmethod
    def GetSolvingVariable(self):
        pass

    @abstractmethod
    def GetElementNamePrefix(self):
        pass

    @abstractmethod
    def GetConditionNamePrefix(self):
        pass

    def PrepareModelPart(self):
        self.turbulence_model_part = CreateRansFormulationModelPart(
            self.GetComputingModelPart(),
            self.__class__.__name__,
            self.GetDomainSize(),
            self.element_name,
            self.condition_name)

        Kratos.Logger.PrintInfo(self.__class__.__name__,
                                "Created formulation model part.")

    def Initialize(self):
        InitializeYPlusVariablesInConditions(self.GetBaseModelPart())
        CalculateNormalsOnConditions(self.GetBaseModelPart())

        settings = self.GetParameters()

        for constraint in self.GetBaseModelPart().MasterSlaveConstraints:
            if (constraint.GetSlaveDofsVector()[0].GetVariable() == self.GetSolvingVariable()):
                self.GetModelPart().AddMasterSlaveConstraint(constraint)

        if (self.IsPeriodic()):
            InitializePeriodicConditions(
                self.GetBaseModelPart(),
                self.GetModelPart(),
                [self.GetSolvingVariable()])

        linear_solver_factory = GetKratosObjectPrototype("LinearSolverFactory")
        linear_solver = linear_solver_factory(settings["linear_solver_settings"])

        builder_and_solver = CreateBlockBuilderAndSolver(
            linear_solver,
            self.IsPeriodic(),
            self.GetCommunicator())

        convergence_criteria_type = GetKratosObjectPrototype("MixedGenericCriteria")
        convergence_criteria = convergence_criteria_type([
            (self.GetSolvingVariable(),
             settings["relative_tolerance"].GetDouble(),
             settings["absolute_tolerance"].GetDouble())])

        if (self.is_steady_simulation):
            scheme = self.scheme_type(settings["relaxation_factor"].GetDouble())
        else:
            scheme_type = GetKratosObjectPrototype("BossakRelaxationScalarScheme")
            scheme = scheme_type(
                self.GetModelPart().ProcessInfo[Kratos.BOSSAK_ALPHA],
                settings["relaxation_factor"].GetDouble(),
                self.GetSolvingVariable())

        solver_type = GetKratosObjectPrototype("ResidualBasedNewtonRaphsonStrategy")
        self.solver = solver_type(
            self.GetModelPart(),
            scheme,
            convergence_criteria,
            builder_and_solver,
            settings["max_iterations"].GetInt(),
            False,
            False,
            False)

        self.solver.SetEchoLevel(self.echo_level)
        convergence_criteria.SetEchoLevel(self.echo_level)

        super().Initialize()
        Kratos.Logger.PrintInfo(self.__class__.__name__, "Initialized formulation")

    def SolveCouplingStep(self):
        if (self.IsBufferInitialized()):
            self.ExecuteBeforeCouplingSolveStep()
            self.solver.SolveSolutionStep()
            self.ExecuteAfterCouplingSolveStep()
            Kratos.Logger.PrintInfo(self.__class__.__name__, "Solved  formulation.")
            return True

        return False

    def GetStrategy(self):
        return self.solver

    def SetTimeSchemeSettings(self, settings):
        if (settings.Has("scheme_type")):
            scheme_type = settings["scheme_type"].GetString()
            if (scheme_type == "steady"):
                self.is_steady_simulation = True
            elif (scheme_type == "bdf2" or scheme_type == "bossak"):
                self.is_steady_simulation = False
            else:
                raise Exception(
                    "Only \"steady\", \"bdf2\" and \"bossak\" scheme types supported. [ scheme_type = \""
                    + scheme_type + "\" ]")
        else:
            raise Exception(
                "\"scheme_type\" is missing in time scheme settings")

    def GetMaxCouplingIterations(self):
        return 0

    def GetModelPart(self):
        return self.turbulence_model_part

    def SetStabilizationMethod(self, stabilization_method):
        self.element_name = self.GetElementNamePrefix()
        if (stabilization_method == "algebraic_flux_corrected"):
            self.element_name = self.element_name + "AFC"
            self.scheme_type = self._CreateAlgebraicFluxCorrectedSteadyScalarScheme
        elif (stabilization_method == "residual_based_flux_corrected"):
            self.element_name = self.element_name + "RFC"
            self.scheme_type = GetKratosObjectPrototype("SteadyScalarScheme")
        elif (stabilization_method == "non_linear_cross_wind_dissipation"):
            self.element_name = self.element_name + "CWD"
            self.scheme_type = GetKratosObjectPrototype("SteadyScalarScheme")
        else:
            raise Exception("Unsupported stabilization method")

    def SetWallFunctionSettings(self):
        wall_function_settings = self.GetParameters()["wall_function_settings"]
        wall_function_settings.ValidateAndAssignDefaults(self.GetDefaultParameters()["wall_function_settings"])
        self.condition_name = self.GetConditionNamePrefix()

        if (self.condition_name != ""):
            wall_function_region_type = wall_function_settings["wall_function_region_type"].GetString()
            wall_friction_velocity_calculation_method = wall_function_settings["wall_friction_velocity_calculation_method"].GetString()

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

    def ElementHasNodalProperties(self):
        return False

    def _CreateAlgebraicFluxCorrectedSteadyScalarScheme(self, relaxation_factor):
        if (self.IsPeriodic()):
            return GetKratosObjectPrototype("AlgebraicFluxCorrectedSteadyScalarScheme")(relaxation_factor, GetBoundaryFlags(self.GetParameters()["boundary_flags"]), KratosCFD.PATCH_INDEX)
        else:
            return GetKratosObjectPrototype("AlgebraicFluxCorrectedSteadyScalarScheme")(relaxation_factor, GetBoundaryFlags(self.GetParameters()["boundary_flags"]))