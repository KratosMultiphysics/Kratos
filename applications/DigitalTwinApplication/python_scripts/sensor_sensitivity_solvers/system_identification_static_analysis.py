import typing
import KratosMultiphysics as Kratos
import KratosMultiphysics.DigitalTwinApplication as KratosDT

from KratosMultiphysics.analysis_stage import AnalysisStage
from KratosMultiphysics.DigitalTwinApplication.sensor_sensitivity_solvers.sensor_sensitivity_adjoint_static_solver import SensorSensitivityAdjointStaticSolver
from KratosMultiphysics.DigitalTwinApplication.utilities.sensor_utils import GetSensors
from KratosMultiphysics.DigitalTwinApplication.utilities.expression_utils import GetContainerExpression
from KratosMultiphysics.DigitalTwinApplication.utilities.expression_utils import ExpressionUnionType
from KratosMultiphysics.DigitalTwinApplication.utilities.expression_utils import ExpressionDataLocation

class SystemIdentificationStaticAnalysis(AnalysisStage):
    def __init__(self, model: Kratos.Model, project_parameters: Kratos.Parameters):
        super().__init__(model, project_parameters)
        self.listof_sensors: 'list[KratosDT.Sensors.Sensor]' = {}

    def Initialize(self):
        super().Initialize()

        default_sensor_settings = Kratos.Parameters("""{
            "perturbation_size"            : 1e-8,
            "adapt_perturbation_size"      : true,
            "p_coefficient"                : 2.0,
            "list_of_sensors"              : []
        }""")

        sensor_settings = self.project_parameters["sensor_settings"]
        sensor_settings.ValidateAndAssignDefaults(default_sensor_settings)

        model_part: Kratos.ModelPart = self._GetSolver().GetComputingModelPart()
        model_part.ProcessInfo[KratosDT.PERTURBATION_SIZE] = sensor_settings["perturbation_size"].GetDouble()
        model_part.ProcessInfo[KratosDT.ADAPT_PERTURBATION_SIZE] = sensor_settings["adapt_perturbation_size"].GetBool()
        self.listof_sensors = GetSensors(model_part, sensor_settings["list_of_sensors"].values())

        self.measurement_residual_response_function = KratosDT.Sensors.MeasurementResidualPNormResponseFunction(sensor_settings["p_coefficient"].GetDouble())
        for sensor in self.listof_sensors:
            sensor.SetValue(KratosDT.SENSOR_MEASURED_VALUE, 0.0)
            self.measurement_residual_response_function.AddSensor(sensor)

        self.measurement_residual_response_function.Initialize()

        # since we have to run adjoint per test scenario, we have
        # to rebuild the LHS and RHS for every adjoint solve.
        self._GetSolver()._GetSolutionStrategy().SetRebuildLevel(1)

        self._GetSolver().SetSensor(self.measurement_residual_response_function)

    def _CreateSolver(self) -> SensorSensitivityAdjointStaticSolver:
        return SensorSensitivityAdjointStaticSolver(self.model, self.project_parameters["solver_settings"])

    def _GetSimulationName(self) -> str:
        return "::[SystemIdentificationAnalysis]:: "

    def RunSolutionLoop(self):
        self.InitializeSolutionStep()
        self.measurement_residual_response_function.InitializeSolutionStep()
        self._GetSolver().SolveSolutionStep()
        self.measurement_residual_response_function.FinalizeSolutionStep()
        self.FinalizeSolutionStep()
        self.OutputSolutionStep()

    def GetListOfSensors(self) -> 'list[KratosDT.Sensors.Sensor]':
        return self.listof_sensors

    def GetResponseFunction(self) -> Kratos.AdjointResponseFunction:
        return self.measurement_residual_response_function

    def PrintAnalysisStageProgressInformation(self):
        process_info = self._GetSolver().GetComputingModelPart().ProcessInfo
        Kratos.Logger.PrintInfo(self._GetSimulationName(), f"Computed sensitivities for \"{process_info[KratosDT.TEST_ANALYSIS_NAME]}\".")

    def GetSensitivities(self) -> 'dict[typing.Union[Kratos.DoubleVariable, Kratos.Array1DVariable3], ExpressionUnionType]':
        sensitivity_model_part: Kratos.ModelPart = self._GetSolver().GetSensitivityModelPart()
        sensitivity_variables: 'dict[ExpressionDataLocation, list[typing.Union[Kratos.DoubleVariable, Kratos.Array1DVariable3]]]' = self._GetSolver().GetSensitivtyVariables()

        result: 'dict[typing.Union[Kratos.DoubleVariable, Kratos.Array1DVariable3], ExpressionUnionType]' = {}
        for data_location, variables in sensitivity_variables.items():
            for variable in variables:
                result[variable] = GetContainerExpression(sensitivity_model_part, data_location, variable)
        return result

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
    simulation = SystemIdentificationStaticAnalysis(model, parameters)
    simulation.Run()
