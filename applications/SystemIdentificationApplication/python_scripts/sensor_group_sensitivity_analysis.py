import KratosMultiphysics as Kratos
import KratosMultiphysics.SystemIdentificationApplication as KratosSI

from KratosMultiphysics.OptimizationApplication.optimization_analysis import OptimizationAnalysis
from KratosMultiphysics.OptimizationApplication.utilities.component_data_view import ComponentDataView
from KratosMultiphysics.OptimizationApplication.utilities.helper_utilities import CallOnAll
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem_utilities import OptimizationComponentFactory
from KratosMultiphysics.OptimizationApplication.utilities.logger_utilities import OptimizationAnalysisTimeLogger
from KratosMultiphysics.SystemIdentificationApplication.utilities.sensor_utils import GetSensors
from KratosMultiphysics.SystemIdentificationApplication.sensor_sensitivity_solvers.response_sensitivity_analysis import ResponseSensitivityAnalysis

class SensorGroupSensitivityAnalysis(OptimizationAnalysis):
    @classmethod
    def GetDefaultParameters(cls) -> Kratos.Parameters:
        return Kratos.Parameters("""{
            "problem_data"      : {},
            "model_parts"       : [],
            "analyses"          : [],
            "algorithm_settings": {},
            "processes"  : {
                "kratos_processes"           : {},
                "optimization_data_processes": {}
            }
        }""")

    def Run(self):
        self.Initialize()
        self.Check()
        self.__Solve()
        self.Finalize()

    def _CreateResponses(self):
        pass

    def _CreateControls(self):
        pass

    def _CreateAlgorithm(self):
        default_settings = Kratos.Parameters("""{
            "module"           : "KratosMultiphysics.SystemIdentificationApplication.sensor_sensitivity_solvers",
            "sensor_group_name": ""
        }""")
        algorithm_settings = self.project_parameters["algorithm_settings"]
        algorithm_settings.AddMissingParameters(default_settings)
        self.sensor_group_name = algorithm_settings["sensor_group_name"].GetString()
        self._algorithm = OptimizationComponentFactory(self.model, algorithm_settings, self.optimization_problem)

    def __Solve(self) -> None:
        # first get the sensors list
        sensor_group_data = ComponentDataView(self.sensor_group_name, self.optimization_problem)
        list_of_sensors = GetSensors(sensor_group_data)

        sensor_model_part = self.model[self.sensor_group_name]
        computing_model_part: Kratos.ModelPart = self.GetAlgorithm()._GetSolver().GetComputingModelPart()

        # clear sensor data containers
        for sensor in list_of_sensors:
            sensor.ClearContainerExpressions()

        for process_type in self.GetAlgorithm().GetProcessesOrder():
            CallOnAll(self.optimization_problem.GetListOfProcesses(process_type), Kratos.Process.ExecuteBeforeSolutionLoop)

        for execution_policy in self.optimization_problem.GetListOfExecutionPolicies():
            sensor_model_part.ProcessInfo[KratosSI.TEST_ANALYSIS_NAME] = execution_policy.GetName()
            computing_model_part.ProcessInfo[KratosSI.TEST_ANALYSIS_NAME] = execution_policy.GetName()

            # run the test analysis corresponding to a given load case
            execution_policy.Execute()

            for sensor in list_of_sensors:
                self.optimization_problem.AdvanceStep()

                sensor_model_part.ProcessInfo[Kratos.STEP] = self.optimization_problem.GetStep()
                computing_model_part.ProcessInfo[Kratos.STEP] = self.optimization_problem.GetStep()

                sensor_model_part.ProcessInfo[KratosSI.SENSOR_NAME] = sensor.GetName()
                computing_model_part.ProcessInfo[KratosSI.SENSOR_NAME] = sensor.GetName()

                for process_type in self.GetAlgorithm().GetProcessesOrder():
                    CallOnAll(self.optimization_problem.GetListOfProcesses(process_type), Kratos.Process.ExecuteInitializeSolutionStep)

                # now run the sensor analysis to get sensor sensitivities.
                sensitivities = self.GetAlgorithm().CalculateGradient(sensor)

                for variable, expression in sensitivities.items():
                    sensor_group_data.GetUnBufferedData().SetValue(f"{variable.Name()}", expression.Clone(), overwrite=True)

                for process_type in self.GetAlgorithm().GetProcessesOrder():
                    CallOnAll(self.optimization_problem.GetListOfProcesses(process_type), Kratos.Process.ExecuteFinalizeSolutionStep)

                for process_type in self.GetAlgorithm().GetProcessesOrder():
                    CallOnAll(self.optimization_problem.GetListOfProcesses(process_type), Kratos.Process.ExecuteBeforeOutputStep)

                for process in self.optimization_problem.GetListOfProcesses("output_processes"):
                    if process.IsOutputStep():
                        process.PrintOutput()

            Kratos.Logger.PrintInfo(self.__class__.__name__, f"Finished sensor analysis for \"{execution_policy.GetName()}\"." )
