from os import terminal_size
from pathlib import Path

# import kratos
import KratosMultiphysics as Kratos

# import RANS
import KratosMultiphysics.RANSApplication as KratosRANS

# import utilities
from KratosMultiphysics.RANSApplication.formulations.utilities import SolveProblem

# import formulation interface
from KratosMultiphysics.RANSApplication.formulations.scalar_rans_formulation import ScalarRansFormulation

class CircularConvectionRansFormulation(ScalarRansFormulation):
    def __init__(self, model_part, settings, deprecated_settings):
        default_settings = Kratos.Parameters(r'''
        {
            "formulation_name": "circular_convection",
            "stabilization_method": "residual_based_flux_corrected",
            "circular_convection_solver_settings": {}
        }''')

        settings.ValidateAndAssignDefaults(default_settings)

        self.stabilization_method = settings["stabilization_method"].GetString()
        self.SetStabilizationMethod(self.stabilization_method)

        super().__init__(
            model_part,
            settings["circular_convection_solver_settings"],
            deprecated_settings)

    def AddVariables(self):
        self.GetBaseModelPart().AddNodalSolutionStepVariable(Kratos.NORMAL)

        super().AddVariables()

    def GetSolvingVariable(self):
        return KratosRANS.VELOCITY_POTENTIAL

    def GetElementNamePrefix(self):
        return "RansCircularConvection"

    def GetConditionNamePrefix(self):
        return ""

    def SetConstants(self, settings):
        defaults = Kratos.Parameters('''{
            "circular_convection_constants": {
                "is_clock_wise_rotation": true,
                "rotation_center": [0.0, 0.0, 0.0]
            },
            "stabilization_constants":{
                "dynamic_tau"                       : 0.0,
                "upwind_operator_coefficient"       : 1.2,
                "positivity_preserving_coefficient" : 1.2
            }
        }''')

        settings.RecursivelyValidateAndAssignDefaults(defaults)

        process_info = self.GetBaseModelPart().ProcessInfo

        # stabilization parameters
        constants = settings["circular_convection_constants"]
        process_info.SetValue(KratosRANS.CIRCULAR_CONVECTION_ROTATION_CLOCKWISE, constants["is_clock_wise_rotation"].GetBool())
        process_info.SetValue(KratosRANS.CIRCULAR_CONVECTION_ROTATION_CENTER, constants["rotation_center"].GetVector())

        constants = settings["stabilization_constants"]
        if (self.is_steady_simulation):
            self.GetBaseModelPart().ProcessInfo.SetValue(Kratos.DYNAMIC_TAU, 0.0)
        else:
            self.GetBaseModelPart().ProcessInfo.SetValue(Kratos.DYNAMIC_TAU, constants["dynamic_tau"].GetDouble())

        # stabilization parameters
        process_info.SetValue(KratosRANS.RANS_STABILIZATION_DISCRETE_UPWIND_OPERATOR_COEFFICIENT, constants["upwind_operator_coefficient"].GetDouble())
        process_info.SetValue(KratosRANS.RANS_STABILIZATION_DIAGONAL_POSITIVITY_PRESERVING_COEFFICIENT, constants["positivity_preserving_coefficient"].GetDouble())

    def ComputeTransientResponseFunctionInterpolationError(self, settings):
        # in here we need to calculate response function interpolation error
        default_settings = Kratos.Parameters("""{
            "primal_transient_project_parameters_file_name": "PLEASE_SPECIFY_TRANSIENT_PRIMAL_PROJECT_PARAMETERS_FILE",
            "adjoint_steady_project_parameters_file_name"  : "PLEASE_SPECIFY_STEADY_ADJOINT_PROJECT_PARAMETERS_FILE",
            "logs_path"                                    : "logs",
            "number_of_transient_steps_to_consider"        : 50
        }""")

        settings.ValidateAndAssignDefaults(default_settings)

        current_model_part = self.GetBaseModelPart()
        current_model = current_model_part.GetModel()
        process_info = current_model_part.ProcessInfo
        current_time = process_info[Kratos.TIME]
        current_step = process_info[Kratos.STEP]
        current_dt   = process_info[Kratos.DELTA_TIME]

        log_path = Path(settings["logs_path"].GetString())
        log_path.mkdir(exist_ok=True, parents=True)
        execution_prefix = str(log_path / "response_function_interpolation_error_evaluation_step_{:d}".format(current_step))

        # run the forward primal solution for user defined steps in the transient problem
        # while calculating time averaged quantities
        time_steps = settings["number_of_transient_steps_to_consider"].GetInt()
        start_time = current_time
        end_time = current_time + current_dt * time_steps

        Kratos.Logger.PrintInfo(self.__class__.__name__, "Solving forward primal solution for transient response function interpolation error calculation...")

        # open primal parameters json file
        with open(settings["primal_transient_project_parameters_file_name"].GetString(), "r") as file_input:
            primal_parameters = Kratos.Parameters(file_input.read().replace("<error_computation_step>", str(current_step)))

        # set start time and end time
        primal_parameters["problem_data"]["start_time"].SetDouble(start_time)
        primal_parameters["problem_data"]["end_time"].SetDouble(end_time)

        # solve primal problem
        from KratosMultiphysics.RANSApplication.rans_analysis import RANSAnalysis
        primal_model, primal_simulation = SolveProblem(RANSAnalysis, primal_parameters, execution_prefix + "_primal")

        # copy time averaged quantities from the primal_simulation
        time_averaged_variable_data = [
            ("TIME_AVERAGED_{:s}".format(var.Name()), False, "TIME_AVERAGED_{:s}".format(var.Name()), False) for var in self.GetSolvingVariables()
        ]
        KratosRANS.RansVariableDataTransferProcess(
            primal_model,
            current_model,
            primal_simulation._GetSolver().GetComputingModelPart().FullName(),
            current_model_part.FullName(),
            ["execute"],
            time_averaged_variable_data,
            self.echo_level).Execute()

        # free the memory consumed by primal analysis
        del primal_simulation
        del primal_model

        # now run the backward adjoint steady problem
        Kratos.Logger.PrintInfo(self.__class__.__name__, "Solving adjoint solution for transient response function interpolation error calculation...")

        # open adjoint parameters json file
        with open(settings["adjoint_steady_project_parameters_file_name"].GetString(), "r") as file_input:
            adjoint_parameters = Kratos.Parameters(file_input.read().replace("<error_computation_step>", str(current_step)))

        # set start time and end time
        adjoint_parameters["problem_data"]["start_time"].SetDouble(0.0)
        adjoint_parameters["problem_data"]["end_time"].SetDouble(1.0)

        # solve adjoint problem
        from KratosMultiphysics.RANSApplication.adjoint_rans_analysis import AdjointRANSAnalysis
        adjoint_model, adjoint_simulation = SolveProblem(AdjointRANSAnalysis, adjoint_parameters, execution_prefix + "_adjoint")

        # now transfer RESPONSE_FUNCTION_INTERPOLATION_ERROR data to the current model part
        KratosRANS.RansVariableDataTransferProcess(
            adjoint_model,
            current_model,
            adjoint_simulation._GetSolver().GetComputingModelPart().FullName(),
            current_model_part.FullName(),
            ["execute"],
            [("RESPONSE_FUNCTION_INTERPOLATION_ERROR", False, "RESPONSE_FUNCTION_INTERPOLATION_ERROR", False)],
            self.echo_level).Execute()

        # free the memory consumed by primal analysis
        del adjoint_simulation
        del adjoint_model

        Kratos.Logger.PrintInfo(self.__class__.__name__, "Computed response based interpolation errors for adaptive mesh refinement.")
