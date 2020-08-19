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
from .constant_energy_spectral_formulation import ConstantEnergySpectralFormulation

# import utilities
from KratosMultiphysics.RANSApplication.generate_velocity_fluctuation import GenerateVelocityFluctuationProcess
from KratosMultiphysics.RANSApplication import ScalarVariableDifferenceNormCalculationUtility
from KratosMultiphysics.RANSApplication.formulations.utilities import GetConvergenceInfo

import math

class RANSKinematicFormulation(Formulation):
    def __init__(self, model_part, settings):
        super(RANSKinematicFormulation, self).__init__(model_part, settings)

        default_settings = Kratos.Parameters(r'''
        {
            "formulation_name": "k_epsilon",
            "constant_energy_spectral_u_solver_settings": {},
            "constant_energy_spectral_v_solver_settings": {},
            "constant_energy_spectral_w_solver_settings": {},
            "effective_wave_number_solver_settings": {},
            "coupling_settings":
            {
                "relative_tolerance": 1e-3,
                "absolute_tolerance": 1e-5,
                "max_iterations": 10
            },
            "auxiliar_process_list": [],
            "echo_level": 0
        }''')
        self.settings.ValidateAndAssignDefaults(default_settings)

        self.spectral_u_formulation = ConstantEnergySpectralFormulation(model_part, settings["constant_energy_spectral_u_solver_settings"])
        self.spectral_u_formulation.element_name = "KinematicSimulationSpectralConstantU"
        self.AddFormulation(self.spectral_u_formulation)

        self.spectral_v_formulation = ConstantEnergySpectralFormulation(model_part, settings["constant_energy_spectral_v_solver_settings"])
        self.spectral_v_formulation.element_name = "KinematicSimulationSpectralConstantV"
        self.AddFormulation(self.spectral_v_formulation)

        self.spectral_w_formulation = ConstantEnergySpectralFormulation(model_part, settings["constant_energy_spectral_w_solver_settings"])
        self.spectral_w_formulation.element_name = "KinematicSimulationSpectralConstantW"
        self.AddFormulation(self.spectral_w_formulation)

        self.effective_wave_number_formulation = ConstantEnergySpectralFormulation(model_part, settings["effective_wave_number_solver_settings"])
        self.effective_wave_number_formulation.element_name = "KinematicSimulationEffectiveWaveNumber"
        self.AddFormulation(self.effective_wave_number_formulation)

        self.echo_level = self.settings["echo_level"].GetInt()
        self.effective_wave_number_convergence_utility = ScalarVariableDifferenceNormCalculationUtility(self.GetBaseModelPart(), KratosRANS.EFFECTIVE_WAVE_NUMBER)
        self.SetMaxCouplingIterations(self.settings["coupling_settings"]["max_iterations"].GetInt())

        self.IsFirstStep = 1 # flag for spectral constant calculation

    def AddVariables(self):
        self.GetBaseModelPart().AddNodalSolutionStepVariable(Kratos.VELOCITY)
        self.GetBaseModelPart().AddNodalSolutionStepVariable(Kratos.KINEMATIC_VISCOSITY)
        self.GetBaseModelPart().AddNodalSolutionStepVariable(KratosRANS.TURBULENT_KINETIC_ENERGY)
        self.GetBaseModelPart().AddNodalSolutionStepVariable(KratosRANS.TURBULENT_ENERGY_DISSIPATION_RATE)
        self.GetBaseModelPart().AddNodalSolutionStepVariable(KratosRANS.SPECTRAL_CONSTANT_U)
        self.GetBaseModelPart().AddNodalSolutionStepVariable(KratosRANS.SPECTRAL_CONSTANT_V)
        self.GetBaseModelPart().AddNodalSolutionStepVariable(KratosRANS.SPECTRAL_CONSTANT_W)
        self.GetBaseModelPart().AddNodalSolutionStepVariable(KratosRANS.TURBULENT_KINETIC_ENERGY_U)
        self.GetBaseModelPart().AddNodalSolutionStepVariable(KratosRANS.TURBULENT_KINETIC_ENERGY_V)
        self.GetBaseModelPart().AddNodalSolutionStepVariable(KratosRANS.TURBULENT_KINETIC_ENERGY_W)
        self.GetBaseModelPart().AddNodalSolutionStepVariable(KratosRANS.EFFECTIVE_WAVE_NUMBER)
        self.GetBaseModelPart().AddNodalSolutionStepVariable(Kratos.LAGRANGE_DISPLACEMENT)

        Kratos.Logger.PrintInfo(self.GetName(), "Added solution step variables.")

    def AddDofs(self):
        Kratos.VariableUtils().AddDof(KratosRANS.SPECTRAL_CONSTANT_U, self.GetBaseModelPart())
        Kratos.VariableUtils().AddDof(KratosRANS.SPECTRAL_CONSTANT_V, self.GetBaseModelPart())
        Kratos.VariableUtils().AddDof(KratosRANS.SPECTRAL_CONSTANT_W, self.GetBaseModelPart())
        Kratos.VariableUtils().AddDof(KratosRANS.EFFECTIVE_WAVE_NUMBER, self.GetBaseModelPart())

        Kratos.Logger.PrintInfo(self.GetName(), "Added solution step dofs.")

    def GetMinimumBufferSize(self):
        return 1

    def Initialize(self):
        model_part = self.GetBaseModelPart()
        model = model_part.GetModel()

        default_settings = Kratos.Parameters('''
            {{
                "model_part_name"      : "{0:s}",
                "ABL_friction_velocity" : 0.375,
                "seed_for_random_samples_generation": 2020,
                "lambda_unsteadiness_parameter" : 1.0
            }}
            '''.format(self.GetBaseModelPart().Name))

        velocity_fluctuations_process = GenerateVelocityFluctuationProcess(
                                            model, default_settings)
        self.AddProcess(velocity_fluctuations_process)

        # calculate components of turbulent kinetic energy
        boundary_layer_height = velocity_fluctuations_process.ABL_friction_velocity*10000/6
        for node in model_part.Nodes:
            beta_v = 1 - 0.22 * pow(math.cos(math.pi*node.Z/(2*boundary_layer_height)),4)
            beta_w = 1 - 0.45 * pow(math.cos(math.pi*node.Z/(2*boundary_layer_height)),4)
            denom = 1+beta_v**2+beta_w**2
            K = node.GetSolutionStepValue(KratosRANS.TURBULENT_KINETIC_ENERGY, 0)
            node.SetSolutionStepValue(KratosRANS.TURBULENT_KINETIC_ENERGY_U, 0, K/denom)
            node.SetSolutionStepValue(KratosRANS.TURBULENT_KINETIC_ENERGY_V, 0, K*(beta_v**2)/denom)
            node.SetSolutionStepValue(KratosRANS.TURBULENT_KINETIC_ENERGY_V, 0, K*(beta_w**2)/denom)

        factory = KratosProcessFactory(self.GetBaseModelPart().GetModel())
        self.auxiliar_process_list = factory.ConstructListOfProcesses(
            self.settings["auxiliar_process_list"])
        for process in self.auxiliar_process_list:
            self.AddProcess(process)

        super(RANSKinematicFormulation, self).Initialize()

    def SetConstants(self, settings):
        defaults = Kratos.Parameters('''{
            "total_wave_number": 10
        }''')

        settings.ValidateAndAssignDefaults(defaults)

        process_info = self.GetBaseModelPart().ProcessInfo
        process_info.SetValue(KratosRANS.TOTAL_WAVE_NUMBER, settings["total_wave_number"].GetInt())

    def SetTimeSchemeSettings(self, settings):
        if (settings.Has("scheme_type")):
            scheme_type = settings["scheme_type"].GetString()
            if (scheme_type == "transient"):
                self.is_steady_simulation = False
            else:
                raise Exception("Only \"transient\" scheme types supported. [ scheme_type = \"" + scheme_type  + "\" ]")
        else:
            raise Exception("\"scheme_type\" is missing in time scheme settings")

        super(RANSKinematicFormulation, self).SetTimeSchemeSettings(settings)

    def SolveCouplingStep(self):

        if self.IsFirstStep == 1:
            self.IsFirstStep = 0 # the spectral constants are calculated only at the first step!
            relative_tolerance = self.settings["coupling_settings"]["relative_tolerance"].GetDouble()
            absolute_tolerance = self.settings["coupling_settings"]["absolute_tolerance"].GetDouble()
            max_iterations = self.GetMaxCouplingIterations()

            for itration in range(max_iterations):
                self.effective_wave_number_convergence_utility.InitializeCalculation()

                for formulation in self.list_of_formulations:
                    if (not formulation.SolveCouplingStep()):
                        return False
                self.ExecuteAfterCouplingSolveStep()

                relative_error, absolute_error = self.effective_wave_number_convergence_utility.CalculateDifferenceNorm()
                info = GetConvergenceInfo(KratosRANS.EFFECTIVE_WAVE_NUMBER, relative_error, relative_tolerance, absolute_error, absolute_tolerance)
                Kratos.Logger.PrintInfo(self.GetName() + " CONVERGENCE", info)

                is_converged = relative_error < relative_tolerance or absolute_error < absolute_tolerance
                if (is_converged):
                    Kratos.Logger.PrintInfo(self.GetName() + " CONVERGENCE", "EFFECTIVE_WAVE_NUMBER *** CONVERGENCE ACHIEVED ***")

                Kratos.Logger.PrintInfo(self.GetName(), "Solved coupling itr. " + str(itration + 1) + "/" + str(max_iterations) + ".")
                if (is_converged):
                    return True

        return True