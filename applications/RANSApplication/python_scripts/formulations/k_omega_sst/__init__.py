from __future__ import print_function, absolute_import, division

# import kratos
import KratosMultiphysics as Kratos
from KratosMultiphysics.process_factory import KratosProcessFactory

# import RANS
import KratosMultiphysics.RANSApplication as KratosRANS
from KratosMultiphysics.RANSApplication import RansVariableUtilities

# import formulation interface
from KratosMultiphysics.RANSApplication.formulations.formulation import Formulation

# import formulations
from .k_omega_sst_k_formulation import KOmegaSSTKFormulation
from .k_omega_sst_omega_formulation import KOmegaSSTOmegaFormulation

# import utilities
from KratosMultiphysics.RANSApplication import RansCalculationUtilities
from KratosMultiphysics.RANSApplication import ScalarVariableDifferenceNormCalculationUtility
from KratosMultiphysics.RANSApplication.formulations.utilities import GetConvergenceInfo
from KratosMultiphysics.RANSApplication.formulations.utilities import CreateWallDistanceCalculationProcess

class KOmegaSSTFormulation(Formulation):
    def __init__(self, model_part, settings):
        super(KOmegaSSTFormulation, self).__init__(model_part, settings)

        default_settings = Kratos.Parameters(r'''
        {
            "formulation_name": "k_omega_sst",
            "stabilization_method": "algebraic_flux_corrected",
            "turbulent_kinetic_energy_solver_settings": {},
            "turbulent_specific_energy_dissipation_rate_solver_settings": {},
            "wall_function_properties":{
                "y_plus_calculation_method": "calculated",
                "y_plus_value": 11.06
            },
            "coupling_settings":
            {
                "relative_tolerance": 1e-3,
                "absolute_tolerance": 1e-5,
                "max_iterations": 10
            },
            "wall_distance_calculation_settings":
            {
                "max_iterations"           : 10,
                "echo_level"               : 0,
                "wall_flag_variable_name"  : "STRUCTURE",
                "wall_flag_variable_value" : true,
                "linear_solver_settings" : {
                    "solver_type"     : "amgcl"
                }
            },
            "auxiliar_process_list": [],
            "echo_level": 0
        }''')
        self.settings.ValidateAndAssignDefaults(default_settings)

        self.stabilization_method = self.settings["stabilization_method"].GetString()

        self.tke_formulation = KOmegaSSTKFormulation(model_part, settings["turbulent_kinetic_energy_solver_settings"])
        self.tke_formulation.SetStabilizationMethod(self.stabilization_method)
        self.AddFormulation(self.tke_formulation)

        self.omega_formulation = KOmegaSSTOmegaFormulation(model_part, settings["turbulent_specific_energy_dissipation_rate_solver_settings"])
        self.omega_formulation.SetStabilizationMethod(self.stabilization_method)
        self.AddFormulation(self.omega_formulation)

        self.echo_level = self.settings["echo_level"].GetInt()
        self.nu_t_convergence_utility = ScalarVariableDifferenceNormCalculationUtility(self.GetBaseModelPart(), Kratos.TURBULENT_VISCOSITY)
        self.SetMaxCouplingIterations(self.settings["coupling_settings"]["max_iterations"].GetInt())

    def AddVariables(self):
        self.GetBaseModelPart().AddNodalSolutionStepVariable(Kratos.DENSITY)
        self.GetBaseModelPart().AddNodalSolutionStepVariable(Kratos.VELOCITY)
        self.GetBaseModelPart().AddNodalSolutionStepVariable(Kratos.MESH_VELOCITY)
        self.GetBaseModelPart().AddNodalSolutionStepVariable(Kratos.NORMAL)
        self.GetBaseModelPart().AddNodalSolutionStepVariable(Kratos.VISCOSITY)
        self.GetBaseModelPart().AddNodalSolutionStepVariable(Kratos.KINEMATIC_VISCOSITY)
        self.GetBaseModelPart().AddNodalSolutionStepVariable(Kratos.TURBULENT_VISCOSITY)
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

        Kratos.Logger.PrintInfo(self.GetName(), "Added solution step variables.")

    def AddDofs(self):
        Kratos.VariableUtils().AddDof(KratosRANS.TURBULENT_KINETIC_ENERGY, self.GetBaseModelPart())
        Kratos.VariableUtils().AddDof(KratosRANS.TURBULENT_SPECIFIC_ENERGY_DISSIPATION_RATE, self.GetBaseModelPart())

        Kratos.Logger.PrintInfo(self.GetName(), "Added solution step dofs.")

    def GetMinimumBufferSize(self):
        if (self.is_steady_simulation):
            return 1
        else:
            return 2

    def Initialize(self):
        model_part = self.GetBaseModelPart()
        model = model_part.GetModel()

        wall_distance_calculation_settings = self.settings["wall_distance_calculation_settings"]
        wall_distance_calculation_settings.AddEmptyValue("model_part_name")
        wall_distance_calculation_settings["model_part_name"].SetString(self.GetBaseModelPart().Name)

        wall_distance_process = CreateWallDistanceCalculationProcess(model, wall_distance_calculation_settings)
        self.AddProcess(wall_distance_process)

        process_info = model_part.ProcessInfo
        a1 = process_info[KratosRANS.TURBULENCE_RANS_A1]
        beta_star = process_info[KratosRANS.TURBULENCE_RANS_C_MU]
        wall_model_part_name = process_info[KratosRANS.WALL_MODEL_PART_NAME]
        kappa = process_info[KratosRANS.WALL_VON_KARMAN]

        nut_process = KratosRANS.RansNutKOmegaSSTUpdateProcess(
                                            model,
                                            self.GetBaseModelPart().Name,
                                            a1,
                                            beta_star,
                                            1e-12,
                                            self.echo_level)
        self.AddProcess(nut_process)

        nut_wall_process = KratosRANS.RansNutYPlusWallFunctionUpdateProcess(
                                            model,
                                            wall_model_part_name,
                                            kappa,
                                            1e-12,
                                            self.echo_level)
        self.AddProcess(nut_wall_process)

        factory = KratosProcessFactory(self.GetBaseModelPart().GetModel())
        self.auxiliar_process_list = factory.ConstructListOfProcesses(
            self.settings["auxiliar_process_list"])
        for process in self.auxiliar_process_list:
            self.AddProcess(process)

        super(KOmegaSSTFormulation, self).Initialize()

    def SetConstants(self, settings):
        defaults = Kratos.Parameters('''{
            "wall_law_constants":{
                "kappa": 0.41,
                "beta" : 5.2,
                "c_mu" : 0.09
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
        process_info.SetValue(KratosRANS.WALL_SMOOTHNESS_BETA, constants["beta"].GetDouble())
        process_info.SetValue(KratosRANS.WALL_VON_KARMAN, constants["kappa"].GetDouble())
        process_info.SetValue(KratosRANS.TURBULENCE_RANS_C_MU, constants["c_mu"].GetDouble())
        process_info.SetValue(KratosRANS.RANS_LINEAR_LOG_LAW_Y_PLUS_LIMIT, RansCalculationUtilities.CalculateLogarithmicYPlusLimit(
                                                                                    process_info[KratosRANS.WALL_VON_KARMAN],
                                                                                    process_info[KratosRANS.WALL_SMOOTHNESS_BETA]
                                                                                    ))
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


    def SetTimeSchemeSettings(self, settings):
        if (settings.Has("scheme_type")):
            scheme_type = settings["scheme_type"].GetString()
            if (scheme_type == "steady"):
                self.is_steady_simulation = True
                self.GetBaseModelPart().ProcessInfo.SetValue(Kratos.BOSSAK_ALPHA, 0.0)
            elif (scheme_type == "transient"):
                self.is_steady_simulation = False
                default_settings = Kratos.Parameters('''{
                    "scheme_type": "transient",
                    "alpha_bossak": -0.3
                }''')
                settings.ValidateAndAssignDefaults(default_settings)
                self.GetBaseModelPart().ProcessInfo.SetValue(Kratos.BOSSAK_ALPHA, settings["alpha_bossak"].GetDouble())
            else:
                raise Exception("Only \"steady\" and \"transient\" scheme types supported. [ scheme_type = \"" + scheme_type  + "\" ]")
        else:
            raise Exception("\"scheme_type\" is missing in time scheme settings")

        super(KOmegaSSTFormulation, self).SetTimeSchemeSettings(settings)

    def SolveCouplingStep(self):
        relative_tolerance = self.settings["coupling_settings"]["relative_tolerance"].GetDouble()
        absolute_tolerance = self.settings["coupling_settings"]["absolute_tolerance"].GetDouble()
        max_iterations = self.GetMaxCouplingIterations()

        for iteration in range(max_iterations):
            self.nu_t_convergence_utility.InitializeCalculation()

            for formulation in self.list_of_formulations:
                if (not formulation.SolveCouplingStep()):
                    return False
            self.ExecuteAfterCouplingSolveStep()

            relative_error, absolute_error = self.nu_t_convergence_utility.CalculateDifferenceNorm()
            info = GetConvergenceInfo(Kratos.TURBULENT_VISCOSITY, relative_error, relative_tolerance, absolute_error, absolute_tolerance)
            Kratos.Logger.PrintInfo(self.GetName() + " CONVERGENCE", info)

            is_converged = relative_error < relative_tolerance or absolute_error < absolute_tolerance
            if (is_converged):
                Kratos.Logger.PrintInfo(self.GetName() + " CONVERGENCE", "TURBULENT_VISCOSITY *** CONVERGENCE ACHIEVED ***")

            Kratos.Logger.PrintInfo(self.GetName(), "Solved coupling itr. " + str(iteration + 1) + "/" + str(max_iterations) + ".")
            if (is_converged):
                return True

        return True