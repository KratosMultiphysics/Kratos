from pathlib import Path
import typing
import KratosMultiphysics as Kratos
import KratosMultiphysics.SystemIdentificationApplication as KratosSI
import KratosMultiphysics.HDF5Application as KratosHDF5
import KratosMultiphysics.OptimizationApplication as KratosOA

from KratosMultiphysics.HDF5Application.core.file_io import OpenHDF5File
from KratosMultiphysics.HDF5Application.single_mesh_temporal_output_process import SingleMeshTemporalOutputProcessFactory

from KratosMultiphysics.SystemIdentificationApplication.utilities.sensor_utils import CreateSensors, GetSensors
from KratosMultiphysics.SystemIdentificationApplication.utilities.sensor_utils import PrintSensorListToCSV
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.utilities.component_data_view import ComponentDataView
from KratosMultiphysics.OptimizationApplication.responses.response_function import ResponseFunction


def Factory(arg1: Kratos.Parameters, arg2: Kratos.Model, optimization_problem: 'typing.Optional[OptimizationProblem]' = None):
    if isinstance(arg1, Kratos.Parameters) and isinstance(arg2, Kratos.Model):
        # this is to use SensorHDF5OutputProcess as a normal kratos process
        return SensorHDF5OutputProcess(arg2, arg1["Parameters"], optimization_problem)
    elif isinstance(arg1, Kratos.Model) and isinstance(arg2, Kratos.Parameters):
        # this is to use SensorHDF5OutputProcess as a process for optimization application where optimization_problem cannot be none
        return SensorHDF5OutputProcess(arg1, arg2["settings"], optimization_problem)
    else:
        raise RuntimeError("Argument mismatch.")

class SensorHDF5OutputProcess(Kratos.OutputProcess):
    def __init__(self, model: Kratos.Model, settings: Kratos.Parameters, optimization_problem: 'typing.Optional[OptimizationProblem]' = None):
        Kratos.OutputProcess.__init__(self)

        default_settings = Kratos.Parameters("""{
            "model_part_name": "PLEASE_PROVIDE_A_MODEL_PART_NAME",
            "sensor_model_part_name": "SensorModelPart",
            "sensor_value_variable_name": "PLEASE_ENTER_THE_SENSOR_VARIABLE_NAME",
            "output_file_name": "hdf5_output/<model_part_name>_T_<step>.h5",
            "list_of_sensors": []
        }""")

        if optimization_problem is not None:
            default_settings.RemoveValue("list_of_sensors")
            default_settings.AddEmptyValue("response_name")
            default_settings["response_name"].SetString("")
        settings.ValidateAndAssignDefaults(default_settings)

        self.model_part = model[settings["model_part_name"].GetString()]
        #self.properties_list: 'list[str]' = settings["properties_list"].GetStringArray()
        self.output_file_name = Path(settings["output_file_name"].GetString())
        self.optimization_problem = optimization_problem
        print("SensorHDF5OutputProcess opt problem is ", self.optimization_problem)

        self.list_of_sensors: 'list[KratosSI.Sensors.Sensor]' = []
        if optimization_problem is None:
            self.initialized_sensors = True
            self.compute_sensor_value = True
            sensor_model_part_name = settings["sensor_model_part_name"].GetString()
            if not model.HasModelPart(sensor_model_part_name):
                model.CreateModelPart(sensor_model_part_name)
            #print("sensor_mp is ", model[sensor_model_part_name])
            
            self.list_of_sensors = CreateSensors(model[sensor_model_part_name], self.model_part, settings["list_of_sensors"].values())
            # else: 
            #     self.response_name = self.model_part.GetValue(KratosOA.EXECUTION_POLICY_NAME)
            #     self.list_of_sensors = ComponentDataView(self.optimization_problem.GetComponent(self.response_name, ResponseFunction), self.optimization_problem).GetUnBufferedData().GetValue("sensors")
        else:
            self.initialized_sensors = False
            self.compute_sensor_value = False
            #self.response_name = settings["response_name"].GetString()
            self.response_name = self.model_part.GetValue(KratosOA.EXECUTION_POLICY_NAME)
            print("SensorHDF5OutputProcess response name is ", self.response_name )

        self.sensor_model_part = model[sensor_model_part_name]
        self.sensor_value_variable = Kratos.KratosGlobals.GetVariable(settings["sensor_value_variable_name"].GetString())

    def ExecuteInitializeSolutionStep(self):
        self.sensor_model_part.ProcessInfo[Kratos.STEP] = self.model_part.ProcessInfo[Kratos.STEP]
        self.sensor_model_part.ProcessInfo[Kratos.TIME] = self.model_part.ProcessInfo[Kratos.TIME]

    def ExecuteFinalizeSolutionStep(self):
        if not self.initialized_sensors:
            self.list_of_sensors = ComponentDataView(self.optimization_problem.GetComponent(self.response_name, ResponseFunction), self.optimization_problem).GetUnBufferedData().GetValue("sensors")
            self.initialized_sensors = True

        # calculate and assign values for each sensor
        if self.compute_sensor_value:
            for sensor in self.list_of_sensors:
                sensor.SetSensorValue(sensor.CalculateValue(self.model_part))
                sensor.GetNode().SetValue(self.sensor_value_variable, sensor.GetSensorValue())


    def PrintOutput(self) -> None:
        s_name = str(self.output_file_name)
        s_name = s_name.replace("<model_part_name>", f"{self.sensor_model_part.Name}")
        s_name_orig = s_name
        if not self.initialized_sensors:
            self.list_of_sensors = ComponentDataView(self.optimization_problem.GetComponent(self.response_name, ResponseFunction), self.optimization_problem).GetUnBufferedData().GetValue("sensors")
            self.initialized_sensors = True
            s_name = s_name.replace("<step>", f"{self.optimization_problem.GetStep()}")
        else:
            s_name = s_name.replace("<step>", f"{self.model_part.ProcessInfo[Kratos.STEP]}")
       
        #print("s_name is ", s_name)
        params = Kratos.Parameters("""{
            "file_access_mode": "truncate"
        }""")
        params.AddString("file_name", s_name)

        prefix_settings = Kratos.Parameters("""{                         
        }""")
        if not prefix_settings.Has("list_of_variables"):
            prefix_settings.AddStringArray("list_of_variables", [f"{self.sensor_value_variable.Name()}"])
        if not prefix_settings.Has("prefix"):
            prefix = "/ResultsData"
            prefix_settings.AddString("prefix", prefix)

        with OpenHDF5File(params, self.sensor_model_part) as h5_file:
            expio = KratosHDF5.HDF5NodalDataValueIO(prefix_settings, h5_file)
            expio.Write(self.sensor_model_part)

        exec_policy_name = self.model_part.GetValue(KratosOA.EXECUTION_POLICY_NAME)
        if "MEASURED" in self.sensor_value_variable.Name():
            KratosOA.OptAppModelPartUtils.LogModelPartStatus(self.model_part, f"{exec_policy_name}_MEASURED_SENSOR_DATA_written")
            KratosOA.OptAppModelPartUtils.LogModelPartStatus(self.model_part, f"{exec_policy_name}_MEASURED_SENSOR_DATA_Path:{s_name_orig}:{prefix}")
        elif "COMPUTED" in self.sensor_value_variable.Name():
            KratosOA.OptAppModelPartUtils.LogModelPartStatus(self.model_part, f"{exec_policy_name}_COMPUTED_SENSOR_DATA_written")
            KratosOA.OptAppModelPartUtils.LogModelPartStatus(self.model_part, f"{exec_policy_name}_COMPUTED_SENSOR_DATA_Path:{s_name_orig}:{prefix}")


        




