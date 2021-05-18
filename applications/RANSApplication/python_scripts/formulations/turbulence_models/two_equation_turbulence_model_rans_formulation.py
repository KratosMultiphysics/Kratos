# import kratos
import KratosMultiphysics as Kratos
from KratosMultiphysics.process_factory import KratosProcessFactory

# import formulation interface
from KratosMultiphysics.RANSApplication.formulations.rans_formulation import RansFormulation

# import utilities
from KratosMultiphysics.RANSApplication import ScalarVariableDifferenceNormCalculationUtility
from KratosMultiphysics.RANSApplication.formulations.utilities import GetConvergenceInfo

class TwoEquationTurbulenceModelRansFormulation(RansFormulation):
    def __init__(self, model_part, settings, formulation_1, formulation_2):
        super().__init__(model_part, settings)

        self.stabilization_method = settings["stabilization_method"].GetString()

        self.formulation_1 = formulation_1
        self.formulation_1.SetStabilizationMethod(self.stabilization_method)
        self.AddRansFormulation(self.formulation_1)

        self.formulation_2 = formulation_2
        self.formulation_2.SetStabilizationMethod(self.stabilization_method)
        self.AddRansFormulation(self.formulation_2)

        self.echo_level = settings["echo_level"].GetInt()
        self.nu_t_convergence_utility = ScalarVariableDifferenceNormCalculationUtility(self.GetBaseModelPart(), Kratos.TURBULENT_VISCOSITY)
        self.SetMaxCouplingIterations(settings["coupling_settings"]["max_iterations"].GetInt())

    def GetMinimumBufferSize(self):
        if (self.is_steady_simulation):
            return 1
        else:
            return 2

    def Initialize(self):
        factory = KratosProcessFactory(self.GetBaseModelPart().GetModel())
        self.auxiliar_process_list = factory.ConstructListOfProcesses(
            self.GetParameters()["auxiliar_process_list"])
        for process in self.auxiliar_process_list:
            self.AddProcess(process)

        super().Initialize()

    def SetTimeSchemeSettings(self, settings):
        if (settings.Has("scheme_type")):
            scheme_type = settings["scheme_type"].GetString()
            if (scheme_type == "steady"):
                self.is_steady_simulation = True
                self.GetBaseModelPart().ProcessInfo.SetValue(Kratos.BOSSAK_ALPHA, 0.0)
            elif (scheme_type == "bdf2" or scheme_type == "bossak"):
                self.is_steady_simulation = False
                default_settings = Kratos.Parameters('''{
                    "scheme_type": "PLEASE_SPECIFY_SCHEME_TYPE",
                    "alpha_bossak": -0.3
                }''')
                settings.ValidateAndAssignDefaults(default_settings)
                self.GetBaseModelPart().ProcessInfo.SetValue(Kratos.BOSSAK_ALPHA, settings["alpha_bossak"].GetDouble())
            else:
                raise Exception("Only \"steady\", \"bdf2\" and \"bossak\" scheme types supported. [ scheme_type = \"" + scheme_type  + "\" ]")
        else:
            raise Exception("\"scheme_type\" is missing in time scheme settings")

        super().SetTimeSchemeSettings(settings)

    def SolveCouplingStep(self):
        settings = self.GetParameters()
        relative_tolerance = settings["coupling_settings"]["relative_tolerance"].GetDouble()
        absolute_tolerance = settings["coupling_settings"]["absolute_tolerance"].GetDouble()
        max_iterations = self.GetMaxCouplingIterations()

        for iteration in range(max_iterations):
            self.nu_t_convergence_utility.InitializeCalculation()

            self.ExecuteBeforeCouplingSolveStep()
            for formulation in self.GetRansFormulationsList():
                formulation.SolveCouplingStep()
            self.ExecuteAfterCouplingSolveStep()

            relative_error, absolute_error = self.nu_t_convergence_utility.CalculateDifferenceNorm()
            info = GetConvergenceInfo(Kratos.TURBULENT_VISCOSITY, relative_error, relative_tolerance, absolute_error, absolute_tolerance)
            Kratos.Logger.PrintInfo(self.__class__.__name__ + " CONVERGENCE", info)
            Kratos.Logger.PrintInfo(self.__class__.__name__, "Solved coupling itr. {:d}/{:d}.".format(iteration + 1, max_iterations))

            is_converged = relative_error < relative_tolerance or absolute_error < absolute_tolerance
            if (is_converged):
                Kratos.Logger.PrintInfo(self.__class__.__name__ + " CONVERGENCE", "TURBULENT_VISCOSITY *** CONVERGENCE ACHIEVED ***")
                return True

        return True