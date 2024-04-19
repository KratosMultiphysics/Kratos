from pathlib import Path
import typing
import KratosMultiphysics as Kratos
import KratosMultiphysics.SystemIdentificationApplication as KratosSI

from KratosMultiphysics.analysis_stage import AnalysisStage
from KratosMultiphysics.SystemIdentificationApplication.sensor_sensitivity_solvers.sensor_sensitivity_adjoint_static_solver import SensorSensitivityAdjointStaticSolver
from KratosMultiphysics.SystemIdentificationApplication.utilities.sensor_utils import GetSensors
from KratosMultiphysics.SystemIdentificationApplication.utilities.expression_utils import GetContainerExpression
from KratosMultiphysics.SystemIdentificationApplication.utilities.expression_utils import ExpressionUnionType
from KratosMultiphysics.SystemIdentificationApplication.utilities.expression_utils import ExpressionDataLocation

class SystemIdentificationStaticAnalysis(AnalysisStage):
    def __init__(self, model: Kratos.Model, project_parameters: Kratos.Parameters):
        super().__init__(model, project_parameters)
        self.listof_sensors: 'list[KratosSI.Sensors.Sensor]' = {}

    def Initialize(self):
        super().Initialize()

        default_sensor_settings = Kratos.Parameters("""{
            "perturbation_size"            : 1e-8,
            "adapt_perturbation_size"      : true,
            "list_of_sensors"              : [],
            "output_settings"              : {
                "output_sensor_sensitivity_fields": false,
                "output_folder"                   : "Optimization_Results/sensor_sensitivity_fields"
            }
        }""")

        sensor_settings = self.project_parameters["sensor_settings"]
        sensor_settings.ValidateAndAssignDefaults(default_sensor_settings)

        model_part: Kratos.ModelPart = self._GetSolver().GetComputingModelPart()
        model_part.ProcessInfo[KratosSI.PERTURBATION_SIZE] = sensor_settings["perturbation_size"].GetDouble()
        model_part.ProcessInfo[KratosSI.ADAPT_PERTURBATION_SIZE] = sensor_settings["adapt_perturbation_size"].GetBool()
        self.listof_sensors = GetSensors(model_part, sensor_settings["list_of_sensors"].values())

        self.measurement_residual_response_function = KratosSI.Sensors.MeasurementResidualResponseFunction()
        for sensor in self.listof_sensors:
            sensor.SetValue(KratosSI.SENSOR_MEASURED_VALUE, 0.0)
            self.measurement_residual_response_function.AddSensor(sensor)

        self.measurement_residual_response_function.Initialize()

        # since we have to run adjoint per test scenario, we have
        # to rebuild the LHS and RHS for every adjoint solve.
        self._GetSolver()._GetSolutionStrategy().SetRebuildLevel(1)

        self._GetSolver().SetSensor(self.measurement_residual_response_function)

        output_settings = sensor_settings["output_settings"]
        output_settings.ValidateAndAssignDefaults(default_sensor_settings["output_settings"])
        self.output_sensor_sensitivity_fields = output_settings["output_sensor_sensitivity_fields"].GetBool()
        self.output_sensor_sensitivity_path = output_settings["output_folder"].GetString()

    def _CreateSolver(self) -> SensorSensitivityAdjointStaticSolver:
        return SensorSensitivityAdjointStaticSolver(self.model, self.project_parameters["solver_settings"])

    def _GetSimulationName(self) -> str:
        return "::[SystemIdentificationAnalysis]:: "

    def RunSolutionLoop(self):
        process_info = self._GetSolver().GetComputingModelPart().ProcessInfo

        # first we calculate sensor specific distribution if required because, at the end
        # of this function we should return the residual sensitivities
        if self.output_sensor_sensitivity_fields:
            import KratosMultiphysics.HDF5Application as KratosHDF5
            from KratosMultiphysics.HDF5Application.core.file_io import OpenHDF5File

            output_path = Path(self.output_sensor_sensitivity_path) / f"iteration_{process_info[Kratos.STEP]:05d}.h5"
            output_path.parent.mkdir(parents=True, exist_ok=True)

            hdf5_params = Kratos.Parameters("""{
                "file_name"       : "",
                "file_access_mode": "truncate"
            }""")
            hdf5_params["file_name"].SetString(str(output_path))
            with OpenHDF5File(hdf5_params, self._GetSolver().GetSensitivityModelPart()) as h5_file:
                exp_io = KratosHDF5.ExpressionIO(Kratos.Parameters("""{"prefix": "/SensitivityFieldData/"}"""), h5_file)

                for sensor in self.listof_sensors:
                    process_info[KratosSI.SENSOR_NAME] = f"sensor \"{sensor.GetName()}\""
                    self._GetSolver().SetSensor(sensor)
                    self.InitializeSolutionStep()
                    sensor.InitializeSolutionStep()
                    self._GetSolver().SolveSolutionStep()
                    sensor.FinalizeSolutionStep()
                    self.FinalizeSolutionStep()

                    sensitivities = self.GetSensitivities()
                    for var, cexp in sensitivities.items():
                        exp_io.Write(f"{sensor.GetName()}_{var.Name()}", cexp)

            self._GetSolver().SetSensor(self.measurement_residual_response_function)

        # now we calculate the residual sensitivities to be used in the optimization algorithm
        process_info[KratosSI.SENSOR_NAME] = "for aggregated MeasurementResidualResponseFunction"
        self.InitializeSolutionStep()
        self.measurement_residual_response_function.InitializeSolutionStep()
        self._GetSolver().SolveSolutionStep()
        self.measurement_residual_response_function.FinalizeSolutionStep()
        self.FinalizeSolutionStep()
        self.OutputSolutionStep()

    def GetListOfSensors(self) -> 'list[KratosSI.Sensors.Sensor]':
        return self.listof_sensors

    def GetResponseFunction(self) -> Kratos.AdjointResponseFunction:
        return self.measurement_residual_response_function

    def PrintAnalysisStageProgressInformation(self):
        process_info = self._GetSolver().GetComputingModelPart().ProcessInfo
        Kratos.Logger.PrintInfo(self._GetSimulationName(), f"Computed sensitivities for {process_info[KratosSI.SENSOR_NAME]} using \"{process_info[KratosSI.TEST_ANALYSIS_NAME]}\" analysis.")

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
