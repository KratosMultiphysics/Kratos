import importlib

import KratosMultiphysics as Kratos
import KratosMultiphysics.SystemIdentificationApplication as KratosSI

from KratosMultiphysics.analysis_stage import AnalysisStage
from KratosMultiphysics.process_factory import KratosProcessFactory
from KratosMultiphysics.OptimizationApplication.model_part_controllers.model_part_controller import ModelPartController
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.utilities.component_data_view import ComponentDataView
from KratosMultiphysics.OptimizationApplication.utilities.helper_utilities import CallOnAll
from KratosMultiphysics.OptimizationApplication.utilities.helper_utilities import OptimizationComponentFactory
from KratosMultiphysics.SystemIdentificationApplication.utilities.sensor_utils import GetSensors
from KratosMultiphysics.SystemIdentificationApplication.utilities.sensor_utils import GetSensorModelPart
from KratosMultiphysics.SystemIdentificationApplication.sensor_sensitivity_solvers.sensor_sensitivity_analysis import SensorSensitivityAnalysis

class SensorListSensitivityAnalysis:
    @classmethod
    def GetDefaultParameters(cls) -> Kratos.Parameters:
        return Kratos.Parameters("""{
            "problem_data"            : {},
            "model_parts"             : [],
            "test_analyses"           : [],
            "sensor_analysis_settings": {},
            "processes"               : {
                "kratos_processes"           : {},
                "optimization_data_processes": {}
            }
        }""")

    def __init__(self, model: Kratos.Model, project_parameters: Kratos.Parameters) -> None:
        self.model = model
        self.project_parameters = project_parameters
        self.project_parameters.ValidateAndAssignDefaults(self.GetDefaultParameters())
        self.optimization_problem = OptimizationProblem(0)

        self.__list_of_model_part_controllers: 'list[ModelPartController]' = []

        self._CreateModelPartControllers()

        self.__list_of_test_analyses: 'list[tuple[str, AnalysisStage]]' = []
        for analyses_settings in self.project_parameters["test_analyses"]:
            self.__list_of_test_analyses.append(self._CreateAnalysis(analyses_settings))

        self.__sensor_analysis: SensorSensitivityAnalysis = self._CreateAnalysis(self.project_parameters["sensor_analysis_settings"])[1]

        self._CreateProcesses()

    def Initialize(self) -> None:
        CallOnAll(self.__list_of_model_part_controllers, ModelPartController.ImportModelPart)
        CallOnAll(self.__list_of_model_part_controllers, ModelPartController.Initialize)

        for process_type in self.__sensor_analysis.GetProcessesOrder():
            CallOnAll(self.optimization_problem.GetListOfProcesses(process_type), Kratos.Process.ExecuteInitialize)

        self.__sensor_analysis.Initialize()

    def Check(self) -> None:
        for process_type in self.__sensor_analysis.GetProcessesOrder():
            CallOnAll(self.optimization_problem.GetListOfProcesses(process_type), Kratos.Process.Check)

    def Finalize(self) -> None:
        for process_type in self.__sensor_analysis.GetProcessesOrder():
            CallOnAll(self.optimization_problem.GetListOfProcesses(process_type), Kratos.Process.ExecuteFinalize)

    def Run(self):
        self.Initialize()
        self.Check()

        sensor_model_part = GetSensorModelPart(self.optimization_problem)
        sensors_list = GetSensors(self.optimization_problem)

        # clear sensor data containers
        for sensor in sensors_list:
            sensor.Clear()

        for process_type in self.__sensor_analysis.GetProcessesOrder():
            CallOnAll(self.optimization_problem.GetListOfProcesses(process_type), Kratos.Process.ExecuteBeforeSolutionLoop)

        for test_analysis_name, test_analysis in self.__list_of_test_analyses:
            self.optimization_problem.AdvanceStep()

            sensor_model_part.ProcessInfo[KratosSI.TEST_ANALYSIS_NAME] = test_analysis_name

            # run the test analysis corresponding to a given load case
            test_analysis.Run()

            for sensor in sensors_list:
                sensor_model_part.ProcessInfo[KratosSI.SENSOR_NAME] = sensor.GetName()

                for process_type in self.__sensor_analysis.GetProcessesOrder():
                    CallOnAll(self.optimization_problem.GetListOfProcesses(process_type), Kratos.Process.ExecuteInitializeSolutionStep)

                # now run the sensor analysis to get sensor sensitivities.
                sensitivities = self.__sensor_analysis.CalculateGradient(sensor)

                for var_name, expression in sensitivities.items():
                    # add to sensor
                    sensor.AddContainerExpression(f"{test_analysis_name}_{var_name}", expression.Clone())

                    # add to optimization problem for output and other uses
                    ComponentDataView("sensors", self.optimization_problem).GetUnBufferedData().SetValue(f"{sensor.GetName()}/{sensor.GetName()}_{var_name}", expression.Clone())

                for process_type in self.__sensor_analysis.GetProcessesOrder():
                    CallOnAll(self.optimization_problem.GetListOfProcesses(process_type), Kratos.Process.ExecuteFinalizeSolutionStep)

            for process_type in self.__sensor_analysis.GetProcessesOrder():
                CallOnAll(self.optimization_problem.GetListOfProcesses(process_type), Kratos.Process.ExecuteBeforeOutputStep)

            for process in self.optimization_problem.GetListOfProcesses("output_processes"):
                if process.IsOutputStep():
                    process.PrintOutput()

            Kratos.Logger.PrintInfo(self.__class__.__name__, f"Finished sensor analysis for \"{test_analysis_name}\"." )

        self.Finalize()

    def _CreateModelPartControllers(self):
        default_settings = Kratos.Parameters("""{
            "type": "mdpa_model_part_controller",
            "module": "KratosMultiphysics.OptimizationApplication.model_part_controllers"
        }""")
        for model_part_controller_settings in self.project_parameters["model_parts"].values():
            model_part_controller_settings.AddMissingParameters(default_settings)
            model_part_controller: ModelPartController = OptimizationComponentFactory(self.model, model_part_controller_settings, self.optimization_problem)
            self.__list_of_model_part_controllers.append(model_part_controller)

    def _CreateAnalysis(self, analyses_settings: Kratos.Parameters) -> 'tuple[str, AnalysisStage]':
        problem_name = analyses_settings["problem_data"]["problem_name"].GetString()
        analysis_stage_module_name = analyses_settings["analysis_stage"].GetString()
        analysis_stage_class_name = analysis_stage_module_name.split('.')[-1]
        analysis_stage_class_name = ''.join(x.title() for x in analysis_stage_class_name.split('_'))
        analysis_stage_module = importlib.import_module(analysis_stage_module_name)

        return (problem_name, getattr(analysis_stage_module, analysis_stage_class_name)(self.model, analyses_settings))

    def _CreateProcesses(self):
        process_settings = self.project_parameters["processes"]
        process_settings.ValidateAndAssignDefaults(self.GetDefaultParameters()["processes"])

        kratos_processes = process_settings["kratos_processes"]
        optimization_data_processes = process_settings["optimization_data_processes"]

        factory = KratosProcessFactory(self.model)

        optimization_data_process_default_settings = Kratos.Parameters("""{
            "module" : "KratosMultiphysics.OptimizationApplication.optimization_data_processes"
        }""")

        for process_type in self.__sensor_analysis.GetProcessesOrder():
            self.optimization_problem.AddProcessType(process_type)
            if kratos_processes.Has(process_type):
                for process in factory.ConstructListOfProcesses(kratos_processes[process_type]):
                    self.optimization_problem.AddProcess(process_type, process)
            if optimization_data_processes.Has(process_type):
                for process_settings in optimization_data_processes[process_type].values():
                    process_settings.AddMissingParameters(optimization_data_process_default_settings)
                    process = OptimizationComponentFactory(self.model, process_settings, self.optimization_problem)
                    self.optimization_problem.AddProcess(process_type, process)
