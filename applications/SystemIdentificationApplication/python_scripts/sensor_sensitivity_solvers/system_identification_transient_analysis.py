import csv
import typing
import KratosMultiphysics as Kratos
import KratosMultiphysics.StructuralMechanicsApplication as KratosSA
import KratosMultiphysics.SystemIdentificationApplication as KratosSI
import KratosMultiphysics.OptimizationApplication as KratosOA

from KratosMultiphysics.analysis_stage import AnalysisStage
from KratosMultiphysics.SystemIdentificationApplication.sensor_sensitivity_solvers.sensor_sensitivity_adjoint_transient_solver import SensorSensitivityAdjointTransientSolver
from KratosMultiphysics.SystemIdentificationApplication.utilities.sensor_utils import GetSensors
from KratosMultiphysics.SystemIdentificationApplication.utilities.expression_utils import GetContainerExpression
from KratosMultiphysics.SystemIdentificationApplication.utilities.expression_utils import ExpressionUnionType
from KratosMultiphysics.SystemIdentificationApplication.utilities.expression_utils import ExpressionDataLocation

class SystemIdentificationTransientAnalysis(AnalysisStage):
    def __init__(self, model: Kratos.Model, project_parameters: Kratos.Parameters):
        solver_settings = project_parameters["solver_settings"]

        # Making sure that time step is negative
        if solver_settings["time_stepping"].Has("time_step"):
            if not solver_settings["time_stepping"]["time_step"].GetDouble() < 0:
                raise Exception("StructuralMechanicsAdjointDynamicAnalysis: " + '"time_step" in adjoint problem has to be negative!')
        elif solver_settings["time_stepping"].Has("time_step_table"):
            raise Exception("StructuralMechanicsAdjointDynamicAnalysis: there is currently no variable time stepping possible!")

        # Setting start and end time
        if not project_parameters["problem_data"].Has("start_time"):
            project_parameters["problem_data"].AddEmptyValue("start_time")
            project_parameters["problem_data"]["start_time"].SetDouble( \
                            project_parameters["problem_data"]["start_step"].GetDouble() \
                                )    

        if not project_parameters["problem_data"].Has("end_time"):
            project_parameters["problem_data"].AddEmptyValue("end_time")
            project_parameters["problem_data"]["end_time"].SetDouble( \
                            project_parameters["problem_data"]["start_step"].GetDouble() + \
                            (project_parameters["problem_data"]["nsteps"].GetInt()-0.5)*solver_settings["time_stepping"]["time_step"].GetDouble()
                        )
        self.nsteps = project_parameters["problem_data"]["nsteps"].GetInt()
        super().__init__(model, project_parameters)
        self.listof_sensors: 'list[KratosSI.Sensors.Sensor]' = []

    def Initialize(self):
        super().Initialize()

        print('--------------------------------------')
        print("primal model part time:", self.model["Structure"].ProcessInfo[Kratos.TIME])
        print("adjoint model part time", self.model["AdjointStructure"].ProcessInfo[Kratos.TIME])
        print('--------------------------------------')

        default_sensor_settings = Kratos.Parameters("""{
            "perturbation_size"            : 1e-8,
            "adapt_perturbation_size"      : true,
            "list_of_sensors"              : [],
            "p_coefficient"                : 1
        }""")
            # "output_settings"              : {
            #     "output_sensor_sensitivity_fields": false,
            #     "output_folder"                   : "Optimization_Results/sensor_sensitivity_fields"
            # "output_hdf5_file_name": "sensor_data_<step>.h5"

        sensor_settings = self.project_parameters["sensor_settings"]
        sensor_settings.ValidateAndAssignDefaults(default_sensor_settings)

        model_part: Kratos.ModelPart = self._GetSolver().GetComputingModelPart()
        model_part.ProcessInfo[KratosSI.PERTURBATION_SIZE] = sensor_settings["perturbation_size"].GetDouble()
        model_part.ProcessInfo[KratosSI.ADAPT_PERTURBATION_SIZE] = sensor_settings["adapt_perturbation_size"].GetBool()
        self.list_of_sensors = GetSensors(model_part, sensor_settings["list_of_sensors"].values())
        
        self.sensor_measurement_data_file_name: str = ""

        p_coefficient = sensor_settings["p_coefficient"].GetDouble()
        self.least_squares_response_function = KratosSI.Sensors.LeastSquaresResponseFunction()

        for sensor in self.list_of_sensors:
            sensor.SetValue(KratosSI.SENSOR_MEASURED_VALUE, 0.0)
            self.least_squares_response_function.AddSensor(sensor)

        self.sensor_name_dict: 'dict[str, KratosSI.Sensors.Sensor]' = {}
        for sensor in self.list_of_sensors:
            self.sensor_name_dict[sensor.GetName()] = sensor

        self.least_squares_response_function.Initialize()

        # since we have to run adjoint per test scenario, we have
        # to rebuild the LHS and RHS for every adjoint solve.
        self._GetSolver()._GetSolutionStrategy().SetRebuildLevel(1)

        self._GetSolver().SetSensor(self.least_squares_response_function)

        # dummy time step to correctly calculate the first DELTA_TIME
        self._GetSolver().main_model_part.CloneTimeStep(self.time)

        # output_settings = sensor_settings["output_settings"]
        # output_settings.ValidateAndAssignDefaults(default_sensor_settings["output_settings"])
        # self.output_sensor_sensitivity_fields = output_settings["output_sensor_sensitivity_fields"].GetBool()
        # self.output_sensor_sensitivity_path = output_settings["output_folder"].GetString()
    

    def _CreateSolver(self) -> SensorSensitivityAdjointTransientSolver:
        return SensorSensitivityAdjointTransientSolver(self.model, self.project_parameters["solver_settings"])
    
    def _SetSensorMeasurementDataFileName(self, data_file_name: str) -> None:
        self.sensor_measurement_data_file_name = data_file_name

    def _GetSimulationName(self) -> str:
        return "::[SystemIdentificationAnalysis]:: "

    def KeepAdvancingSolutionLoop(self):
        # Note that the adjoint problem is solved in reverse time
        return self.time > self.end_time
    
    def InitializeSolutionStep(self):
        super().InitializeSolutionStep()
        # set sensor measured value for current time step
        adjoint_step = self._GetSolver().GetComputingModelPart().ProcessInfo[Kratos.STEP]
        step = self.nsteps - adjoint_step + 1
        self.__SetSensorMeasuredValue(self.sensor_measurement_data_file_name.replace("<step>", f"{step:06d}"))

    def GetListOfSensors(self) -> 'list[KratosSI.Sensors.Sensor]':
        return self.list_of_sensors
    
    def GetSensorNameDictionary(self) -> 'dict[str, KratosSI.Sensors.Sensor]':
        return self.sensor_name_dict

    def GetResponseFunction(self) -> Kratos.AdjointResponseFunction:
        return self.least_squares_response_function

    def PrintAnalysisStageProgressInformation(self):
        process_info = self._GetSolver().GetComputingModelPart().ProcessInfo
        Kratos.Logger.PrintInfo(self._GetSimulationName(), f"Computed sensitivities for {process_info[KratosSI.SENSOR_NAME]} using \"{process_info[KratosSI.TEST_ANALYSIS_NAME]}\" analysis.")

    def __GetSensor(self, sensor_name: str) -> KratosSI.Sensors.Sensor:
        return self.sensor_name_dict[sensor_name]

    def __GetHeaderIndices(self, csv_stream: csv.reader) -> 'tuple[int, int]':
        headers = [s.strip() for s in next(csv_stream)]
        name_index = headers.index("name")
        value_index = headers.index("value")
        return name_index, value_index

    def __SetSensorMeasuredValue(self, sensor_measurement_data_file_name: str) -> None:
        with open(sensor_measurement_data_file_name, "r") as csv_measurement_file:
            csv_measurement_stream = csv.reader(csv_measurement_file, delimiter=",")
            measured_name_index, measured_value_index = self.__GetHeaderIndices(csv_measurement_stream)

            for measured_row in csv_measurement_stream:
                measured_sensor_name = measured_row[measured_name_index].strip()
                measured_value = float(measured_row[measured_value_index])
                self.__GetSensor(measured_sensor_name).SetValue(KratosSI.SENSOR_MEASURED_VALUE, measured_value)

    def GetSensitivities(self, model_part: Kratos.ModelPart) -> 'dict[typing.Union[Kratos.DoubleVariable, Kratos.Array1DVariable3], ExpressionUnionType]':
        sensitivity_variables: 'dict[ExpressionDataLocation, list[typing.Union[Kratos.DoubleVariable, Kratos.Array1DVariable3]]]' = self._GetSolver().GetSensitivtyVariables()

        result: 'dict[typing.Union[Kratos.DoubleVariable, Kratos.Array1DVariable3], ExpressionUnionType]' = {}
        for data_location, variables in sensitivity_variables.items():
            for variable in variables:
                result[variable] = GetContainerExpression(model_part, data_location, variable)
        return result
    
#     def ReadMeasrumentData(self) -> None:
#         process_info = self._GetSolver().GetComputingModelPart().ProcessInfo
#         sensor_measurement_data_file_name = self.measurement_data_file_name.replace("<step>", str(process_info[Kratos.STEP]))
#         with open(sensor_measurement_data_file_name, "r") as csv_measurement_file:
#             csv_measurement_stream = csv.reader(csv_measurement_file, delimiter=",")
#             measured_name_index, measured_value_index = self.__GetHeaderIndices(csv_measurement_stream)

#             for measured_row in csv_measurement_stream:
#                 measured_sensor_name = measured_row[measured_name_index].strip()
#                 measured_value = float(measured_row[measured_value_index])
#                 self.__GetSensor(measured_sensor_name).SetValue(KratosSI.SENSOR_MEASURED_VALUE, measured_value)
    
#     def OutputSensorData(self) -> None:
#         # read the following from the json params
#         sensitiviy_variable = KratosSA.YOUNG_MODULUS_SENSITIVITY
#         elem_exp = Kratos.Expression.ElementExpression(self.model["all_nodes_elements_model_part"])
#         KratosOA.PropertiesVariablexpressionIO.Read(elem_exp, sensitiviy_variable)

#         import h5py
#         with h5py.File("sensor_data.h5", "w") as h5_output:
#             h5_output["/data"] = elem_exp.Evaluate()