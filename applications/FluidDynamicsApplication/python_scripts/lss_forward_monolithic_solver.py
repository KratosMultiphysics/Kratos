# Importing the Kratos Library
import KratosMultiphysics as Kratos
import KratosMultiphysics.FluidDynamicsApplication as KratosCFD
from KratosMultiphysics.FluidDynamicsApplication.adjoint_monolithic_solver import AdjointMonolithicSolver
from KratosMultiphysics.time_based_ascii_file_writer_utility import TimeBasedAsciiFileWriterUtility

def CreateSolver(model, custom_settings):
    return LSSForwardMonolithicSolver(model, custom_settings)

class LSSForwardMonolithicSolver(AdjointMonolithicSolver):

    def __init__(self, model, settings):
        default_lss_settings = Kratos.Parameters("""
        {
            "delta_time_dialation_alpha" : 1.0,
            "sensitivity_custom_settings": {},
            "convergence_output_settings": {
                "file_name"         : "least_squares_shadowing_delta_time_derivative_difference.dat",
                "output_path"       : "lss_results",
                "write_buffer_size" : -1
            }
        }""")
        if settings.Has("lss_settings"):
            self.lss_settings = settings["lss_settings"].Clone()
            self.lss_settings.ValidateAndAssignDefaults(default_lss_settings)
            settings.RemoveValue("lss_settings")
        else:
            self.lss_settings = default_lss_settings

        super().__init__(model, settings)

        # Overwrite the default buffer size in base AdjointMonolithicSolver
        # TODO: CHECK WHY THE DEFAULT BUFFER SIZE IS 2 IF THE DEFAULT TIME SCHEME IS BOSSAK
        self.min_buffer_size = 2

    def AddVariables(self):
        super().AddVariables()

        self.main_model_part.AddNodalSolutionStepVariable(KratosCFD.LSS_VELOCITY)
        self.main_model_part.AddNodalSolutionStepVariable(KratosCFD.LSS_ACCELERATION)
        self.main_model_part.AddNodalSolutionStepVariable(KratosCFD.LSS_PRESSURE)
        self.main_model_part.AddNodalSolutionStepVariable(Kratos.REACTION)

        Kratos.Logger.PrintInfo(self.__class__.__name__, "LSS forward monolithic solver variables added correctly.")

    def Initialize(self):
        super().Initialize()

        # write time step sensitivity difference to file
        self.time_step_sensitivity_convergence_file = TimeBasedAsciiFileWriterUtility(
                self.GetComputingModelPart(),
                self.lss_settings["convergence_output_settings"],
                "Least squares shadowing delta time total derivative difference square\n Step, Delta time total derivative difference squared")

    def InitializeSolutionStep(self):
        process_info = self.GetComputingModelPart().ProcessInfo
        previous_time_step_sensitivity = process_info[KratosCFD.TIME_STEP_SENSITIVITY]
        self._GetSolutionStrategy().InitializeSolutionStep()
        self.GetResponseFunction().InitializeSolutionStep()
        current_time_step_sensitivity = process_info[KratosCFD.TIME_STEP_SENSITIVITY]
        time_step_sensitivity_difference = (previous_time_step_sensitivity - current_time_step_sensitivity) ** 2
        self.time_step_sensitivity_convergence_file.file.write("{:d}, {:f}\n".format(process_info[Kratos.STEP], time_step_sensitivity_difference))

    def Predict(self):
        self._GetSolutionStrategy().Predict()

    def SolveSolutionStep(self):
        return self._GetSolutionStrategy().SolveSolutionStep()

    def FinalizeSolutionStep(self):
        self._GetSolutionStrategy().FinalizeSolutionStep()
        self.GetResponseFunction().FinalizeSolutionStep()

    def Check(self):
        self._GetSolutionStrategy().Check()

    def _CreateScheme(self):
        response_function = self.GetResponseFunction()
        scheme_type = self.settings["scheme_settings"]["scheme_type"].GetString()
        domain_size = self.main_model_part.ProcessInfo[Kratos.DOMAIN_SIZE]

        # the schemes used in fluid supports SLIP conditions which rotates element/condition
        # matrices based on nodal NORMAL. Hence, the consistent adjoints also need to
        # rotate adjoint element/condition matrices accordingly. Therefore, following
        # schemes are used.
        if scheme_type == "bossak":
            sensitivity_vec = self.settings["sensitivity_settings"]["nodal_solution_step_sensitivity_variables"].GetStringArray()
            if len(sensitivity_vec) > 1:
                raise RuntimeError("LSS only supports sensitivity computation w.r.t. one variable only.")

            sensitivity_custom_settings = self.lss_settings["sensitivity_custom_settings"]
            if sensitivity_vec[0] == "SHAPE_SENSITIVITY":
                self.fluid_lss_sensitivity = KratosCFD.FluidLSSShapeSensitivity(sensitivity_custom_settings, domain_size)
            else:
                raise RuntimeError("LSS only supports SHAPE_SENSITIVITY.")

            def add_variables(domain_variable, other_variable: list):
                result = []
                result.append(Kratos.KratosGlobals.GetVariable(domain_variable.Name() + "_X"))
                result.append(Kratos.KratosGlobals.GetVariable(domain_variable.Name() + "_Y"))
                if domain_size == 3:
                    result.append(Kratos.KratosGlobals.GetVariable(domain_variable.Name() + "_Z"))
                result.append(other_variable)
                return result

            primal_variables_list = add_variables(Kratos.VELOCITY, Kratos.PRESSURE)
            primal_first_derivative_variables_list = add_variables(Kratos.ACCELERATION, None)
            adjoint_variables_list = add_variables(KratosCFD.LSS_VELOCITY, KratosCFD.LSS_PRESSURE)
            adjoint_first_derivative_variables_list = add_variables(KratosCFD.LSS_ACCELERATION, None)
            lss_variables_list = add_variables(KratosCFD.ADJOINT_FLUID_VECTOR_1, KratosCFD.ADJOINT_FLUID_SCALAR_1)
            lss_first_derivative_variables_list = add_variables(KratosCFD.ADJOINT_FLUID_VECTOR_3, None)

            self.lss_variable_utilities = KratosCFD.FluidLSSVariableUtilities(
                                                    primal_variables_list,
                                                    primal_first_derivative_variables_list,
                                                    adjoint_variables_list,
                                                    adjoint_first_derivative_variables_list,
                                                    lss_variables_list,
                                                    lss_first_derivative_variables_list)
            scheme = KratosCFD.LSSBossakForwardScheme(
                                response_function,
                                self.fluid_lss_sensitivity,
                                self.lss_variable_utilities,
                                self.settings["scheme_settings"]["alpha_bossak"].GetDouble(),
                                self.lss_settings["delta_time_dialation_alpha"].GetDouble(),
                                0.0,
                                domain_size,
                                domain_size + 1,
                                self.settings["echo_level"].GetInt())
        else:
            raise Exception("Invalid scheme_type: " + scheme_type)
        return scheme

    def _CreateSensitivityBuilder(self):
        return None


