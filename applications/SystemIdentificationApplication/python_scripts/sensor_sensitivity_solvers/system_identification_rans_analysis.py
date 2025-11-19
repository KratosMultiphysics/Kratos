import KratosMultiphysics as Kratos
import KratosMultiphysics.SystemIdentificationApplication as KratosSI

from KratosMultiphysics.RANSApplication.adjoint_rans_analysis import AdjointRANSAnalysis
from KratosMultiphysics.SystemIdentificationApplication.sensor_sensitivity_solvers.response_sensitivity_analysis import ResponseSensitivityAnalysis
from KratosMultiphysics.SystemIdentificationApplication.sensor_sensitivity_solvers.sensor_sensitivity_adjoint_rans_solver import SensorSensitivityAdjointRANSSolver

def Factory(model: Kratos.Model, parameters: Kratos.Parameters, _):
    return SystemIdentificationRANSAnalysis(model, parameters["settings"])

class SystemIdentificationRANSAnalysis(AdjointRANSAnalysis, ResponseSensitivityAnalysis):
    def _CreateSolver(self) -> SensorSensitivityAdjointRANSSolver:
        return SensorSensitivityAdjointRANSSolver(self.model, self.project_parameters["solver_settings"])

    def _GetSimulationName(self) -> str:
        return "::[SystemIdentificationRANSAnalysis]:: "

    def CalculateGradient(self, response_function: Kratos.AdjointResponseFunction) -> None:
        self._GetSolver().SetResponseFunction(response_function)
        self.time = self._AdvanceTime()
        self.InitializeSolutionStep()
        response_function.InitializeSolutionStep()
        self._GetSolver().Predict()
        self._GetSolver().SolveSolutionStep()
        response_function.FinalizeSolutionStep()
        self.FinalizeSolutionStep()
        self.OutputSolutionStep()

    def PrintAnalysisStageProgressInformation(self):
        process_info = self._GetSolver().GetComputingModelPart().ProcessInfo
        Kratos.Logger.PrintInfo(self._GetSimulationName(), f"Step {process_info[Kratos.STEP]}: Computed sensitivities for response \"{process_info[KratosSI.SENSOR_NAME]}\" using \"{process_info[KratosSI.TEST_ANALYSIS_NAME]}\" analysis.")

if __name__ == "__main__":
    from sys import argv

    if len(argv) > 2:
        err_msg =  'Too many input arguments!\n'
        err_msg += 'Use this script in the following way:\n'
        err_msg += '- With default ProjectParameters (read from "ProjectParameters.json"):\n'
        err_msg += '    "python3 structural_mechanics_analysis.py"\n'
        err_msg += '- With custom ProjectParameters:\n'
        err_msg += '    "python3 structural_mechanics_analysis.py CustomProjectParameters.json"\n'
        raise Exception(err_msg)

    if len(argv) == 2: # ProjectParameters is being passed from outside
        project_parameters_file_name = argv[1]
    else: # using default name
        project_parameters_file_name = "ProjectParameters.json"

    with open(project_parameters_file_name,'r') as parameter_file:
        parameters = Kratos.Parameters(parameter_file.read())

    model = Kratos.Model()
    simulation = SystemIdentificationRANSAnalysis(model, parameters)
    simulation.Run()
