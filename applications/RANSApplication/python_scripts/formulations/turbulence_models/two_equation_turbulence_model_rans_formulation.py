# import kratos
import KratosMultiphysics as Kratos
from KratosMultiphysics.process_factory import KratosProcessFactory

# import formulation interface
import KratosMultiphysics.RANSApplication as KratosRANS
from KratosMultiphysics.RANSApplication.formulations.rans_formulation import RansFormulation

# import utilities
from KratosMultiphysics.RANSApplication import RansNutUtility

class TwoEquationTurbulenceModelRansFormulation(RansFormulation):
    def __init__(self, model_part, settings, deprecated_settings_dict, formulation_1, formulation_2):
        super().__init__(model_part, settings)

        self.stabilization_method = settings["stabilization_method"].GetString()

        self.formulation_1 = formulation_1
        self.formulation_1.SetStabilizationMethod(self.stabilization_method)
        self.AddRansFormulation(self.formulation_1)

        self.formulation_2 = formulation_2
        self.formulation_2.SetStabilizationMethod(self.stabilization_method)
        self.AddRansFormulation(self.formulation_2)

        self.echo_level = settings["echo_level"].GetInt()
        self.nu_t_convergence_utility = RansNutUtility(
            self.GetBaseModelPart(),
            settings["coupling_settings"]["relative_tolerance"].GetDouble(),
            settings["coupling_settings"]["absolute_tolerance"].GetDouble(),
            self.echo_level)
        self.SetMaxCouplingIterations(settings["coupling_settings"]["max_iterations"].GetInt())

    def GetDefaultParameters(self):
        return Kratos.Parameters("""{}""")

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

        # adding this process as the last process to copy historical turbulence nodal data to non historical
        # container so user can have custom processes to assign historical nodal turbulence data.
        turbulence_data_copy_process = KratosRANS.RansVariableDataTransferProcess(
            self.GetBaseModelPart().GetModel(),
            self.GetBaseModelPart().Name,
            self.GetBaseModelPart().Name,
            ["initialize", "after_coupling_solve_step"],
            [(self.formulation_1.GetSolvingVariable().Name(), True, 0, self.formulation_1.GetSolvingVariable().Name(), False, 0),
             (self.formulation_2.GetSolvingVariable().Name(), True, 0, self.formulation_2.GetSolvingVariable().Name(), False, 0)],
            self.echo_level
        )
        self.AddProcess(turbulence_data_copy_process)

        super().Initialize()

        self.nu_t_convergence_utility.Initialize()

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
        max_iterations = self.GetMaxCouplingIterations()

        for iteration in range(max_iterations):
            self.ExecuteBeforeCouplingSolveStep()

            self.nu_t_convergence_utility.InitializeCalculation()

            for formulation in self.GetRansFormulationsList():
                formulation.SolveCouplingStep()

            Kratos.Logger.PrintInfo(self.__class__.__name__, "Solved coupling itr. {:d}/{:d}.".format(iteration + 1, max_iterations))
            self.is_converged = self.nu_t_convergence_utility.CheckConvergence()

            self.ExecuteAfterCouplingSolveStep()
            if (self.is_converged):
                return True

        return False

    def IsConverged(self):
        is_converged = super().IsConverged()

        if hasattr(self, "is_converged"):
            return is_converged and self.is_converged
        else:
            return is_converged

    def ExecuteAfterCouplingSolveStep(self):
        super().ExecuteAfterCouplingSolveStep()

        # It is required to update nut here for old FractionalStep and
        # VMS formulations, so this nut can be distributed via rans_nut_nodal_update process
        # otherwise, this can be moved to FinalizeSolvingStep method
        self.nu_t_convergence_utility.UpdateTurbulentViscosity()

