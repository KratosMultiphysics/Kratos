import typing
from pathlib import Path

import KratosMultiphysics as Kratos
import KratosMultiphysics.SystemIdentificationApplication as KratosSI
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem
from KratosMultiphysics.OptimizationApplication.utilities.component_data_view import ComponentDataView
from KratosMultiphysics.OptimizationApplication.utilities.buffered_dict import BufferedDict
from KratosMultiphysics.SystemIdentificationApplication.utilities.data_utils import SupportedValueUnionType
from KratosMultiphysics.SystemIdentificationApplication.utilities.data_utils import SupportedVariableUnionType
from KratosMultiphysics.SystemIdentificationApplication.utilities.data_utils import GetParameterToKratosValuesConverter
from KratosMultiphysics.SystemIdentificationApplication.utilities.data_utils import GetKratosValueToCSVStringConverter
from KratosMultiphysics.SystemIdentificationApplication.utilities.data_utils import GetNameToCSVString

def CreateSensors(sensor_model_part: Kratos.ModelPart, analysis_model_part: Kratos.ModelPart, list_of_parameters: 'list[Kratos.Parameters]') -> 'list[KratosSI.Sensors.Sensor]':
    """Get list of sensors from given parameters.

    This generates list of sensors from given list of parameters. Each parameter
    corresponds to one settings set for one sensor.

    Args:
        sensor_model_part (Kratos.ModelPart): Model part to store sensor nodes.
        analysis_model_part (Kratos.ModelPart): Analysis model part on which the sensors will be generated.
        list_of_parameters (list[Kratos.Parameters]): List of parameters.

    Returns:
        list[KratosSI.Sensors.Sensor]: List of sensors generated.
    """
    point_locator = Kratos.BruteForcePointLocator(analysis_model_part)

    dict_of_sensor_types: 'dict[str, typing.Type[KratosSI.Sensors.Sensor]]' = {}
    for sensor_type_name in dir(KratosSI.Sensors):
        sensor_type = getattr(KratosSI.Sensors, sensor_type_name)
        if sensor_type_name.endswith("Sensor") and issubclass(sensor_type, KratosSI.Sensors.Sensor):
            dict_of_sensor_types[Kratos.StringUtilities.ConvertCamelCaseToSnakeCase(sensor_type_name)] = sensor_type

    list_of_sensors: 'list[KratosSI.Sensors.Sensor]' = []
    shape_funcs = Kratos.Vector()
    for sensor_index, parameters in enumerate(list_of_parameters):
        if not parameters.Has("type"):
            raise RuntimeError(f"The sensor parameters does not contain \"type\".")
        sensor_type_name = parameters["type"].GetString()

        if not sensor_type_name in dict_of_sensor_types.keys():
            raise RuntimeError(f"Unsupported sensor type = \"{sensor_type_name}\" requested. Followings are supported:\n\t" + "\n\t".join(dict_of_sensor_types.keys()))
        parameters.ValidateAndAssignDefaults(dict_of_sensor_types[sensor_type_name].GetDefaultParameters())

        if sensor_type_name == "displacement_sensor":
            name = parameters["name"].GetString()
            loc = parameters["location"].GetVector()
            loc = Kratos.Point(loc[0], loc[1], loc[2])
            weight = parameters["weight"].GetDouble()
            direction = parameters["direction"].GetVector()
            elem_id = point_locator.FindElement(loc, shape_funcs, Kratos.Configuration.Initial, 1e-8)
            node = sensor_model_part.CreateNewNode(sensor_index + 1, loc[0], loc[1], loc[2])
            sensor = KratosSI.Sensors.DisplacementSensor(name, node, direction, analysis_model_part.GetElement(elem_id), weight)
            AddSensorVariableData(sensor, parameters["variable_data"])
            list_of_sensors.append(sensor)
        elif sensor_type_name == "strain_sensor":
            name = parameters["name"].GetString()
            loc = parameters["location"].GetVector()
            loc = Kratos.Point(loc[0], loc[1], loc[2])
            weight = parameters["weight"].GetDouble()
            strain_variable: Kratos.MatrixVariable = Kratos.KratosGlobals.GetVariable(parameters["strain_variable"].GetString())
            strain_type = parameters["strain_type"].GetString()
            if strain_type == "strain_xx":
                strain_type_value = KratosSI.Sensors.StrainSensor.STRAIN_XX
            elif strain_type == "strain_yy":
                strain_type_value = KratosSI.Sensors.StrainSensor.STRAIN_YY
            elif strain_type == "strain_zz":
                strain_type_value = KratosSI.Sensors.StrainSensor.STRAIN_ZZ
            elif strain_type == "strain_xy":
                strain_type_value = KratosSI.Sensors.StrainSensor.STRAIN_XY
            elif strain_type == "strain_xz":
                strain_type_value = KratosSI.Sensors.StrainSensor.STRAIN_XZ
            elif strain_type == "strain_yz":
                strain_type_value = KratosSI.Sensors.StrainSensor.STRAIN_YZ
            elem_id = point_locator.FindElement(loc, shape_funcs, Kratos.Configuration.Initial, 1e-8)
            node = sensor_model_part.CreateNewNode(sensor_index + 1, loc[0], loc[1], loc[2])
            sensor = KratosSI.Sensors.StrainSensor(name, node, strain_variable, strain_type_value, analysis_model_part.GetElement(elem_id), weight)
            AddSensorVariableData(sensor, parameters["variable_data"])
            list_of_sensors.append(sensor)
    return list_of_sensors

def PrintSensorListToCSV(output_file_name: Path, list_of_sensors: 'list[KratosSI.Sensors.Sensor]', list_of_sensor_properties: 'list[str]') -> None:
    """Writes data in the list_of_sensors to CSV file.

    This method writes data from all the list_of_sensors to the CSV file. The columns of the
    CSV is determined by list_of_sensor_properties.

    Args:
        output_file_name (Path): CSV file name (including the .csv extension).
        list_of_sensors (list[Sensor]): List of sensors.
        list_of_sensor_properties (list[str]): List of columns to write.
    """
    number_of_sensor_views = len(list_of_sensors)

    # do nothing if number of clusters is zero
    if number_of_sensor_views == 0:
        return

    # create output path
    output_file_name.parent.mkdir(exist_ok=True, parents=True)

    with open(str(output_file_name), "w") as file_output:
        # print the headers
        first_sensor = list_of_sensors[0]

        # write the csv id
        file_output.write("#")

        # store converters to CSV in dicts
        params_value_dict: 'dict[str, typing.Callable[[Kratos.Parameters], str]]' = {}
        params_str_dict: 'dict[str, typing.Callable[[SupportedValueUnionType], str]]' = {}
        data_value_dict: 'dict[SupportedVariableUnionType, typing.Callable[[SupportedValueUnionType], str]]' = {}

        dummy_params = first_sensor.GetSensorParameters()
        list_of_params_keys: 'list[str]' = []
        list_of_data_keys: 'list[SupportedVariableUnionType]' = []
        param_property_headers = ""
        value_property_headers = ""
        for sensor_property in list_of_sensor_properties:
            if dummy_params.Has(sensor_property):
                list_of_params_keys.append(sensor_property)
                value_func = GetParameterToKratosValuesConverter(dummy_params[sensor_property])
                value = value_func(dummy_params[sensor_property])
                params_value_dict[sensor_property] = value_func
                params_str_dict[sensor_property] = GetKratosValueToCSVStringConverter(value)
                param_property_headers += f",{GetNameToCSVString(sensor_property, value)}"
            else:
                var = Kratos.KratosGlobals.GetVariable(sensor_property)
                if first_sensor.Has(var):
                    list_of_data_keys.append(var)
                    value = first_sensor.GetNode().GetValue(var)
                    data_value_dict[var] = GetKratosValueToCSVStringConverter(value)
                    value_property_headers += f",{GetNameToCSVString(sensor_property, value)}"
                else:
                    raise RuntimeError(f"Sensor property = \"{sensor_property}\" not found in either sensor parameters nor in data value container. [ Sensor = {first_sensor} ].")

        file_output.write(f"{param_property_headers}{value_property_headers}\n")

        for i, sensor in enumerate(list_of_sensors):
            # write the id
            file_output.write(f"{i+1}")

            # first write the values from sensor params
            sensor_params = sensor.GetSensorParameters()
            for sensor_property in list_of_params_keys:
                file_output.write(f",{params_str_dict[sensor_property](params_value_dict[sensor_property](sensor_params[sensor_property]))}")

            # now write the data values from sensor
            for var in list_of_data_keys:
                file_output.write(f",{data_value_dict[var](sensor.GetNode().GetValue(var))}")

            file_output.write("\n")

def AddSensorVariableData(sensor: KratosSI.Sensors.Sensor, variable_data: Kratos.Parameters) -> None:
    """Adds sensor variable data.

    Args:
        sensor (KratosSI.Sensors.Sensor): Sensor to add data.
        variable_data (Kratos.Parameters): Parameters containing variable data.
    """
    for var_name, var_value in variable_data.items():
        var = Kratos.KratosGlobals.GetVariable(var_name)
        value_func =  GetParameterToKratosValuesConverter(var_value)
        sensor.GetNode().SetValue(var, value_func(var_value))

def SetSensorModelPart(sensor_model_part: Kratos.ModelPart, optimization_problem: OptimizationProblem) -> None:
    ComponentDataView("sensors", optimization_problem).GetUnBufferedData().SetValue("sensor_model_part", sensor_model_part)

def GetSensorModelPart(optimization_problem: OptimizationProblem) -> Kratos.ModelPart:
    return ComponentDataView("sensors", optimization_problem).GetUnBufferedData().GetValue("sensor_model_part")

def SetSensors(list_of_sensors: 'list[KratosSI.Sensors.Sensor]', optimization_problem: OptimizationProblem) -> None:
    data = ComponentDataView("sensors", optimization_problem).GetUnBufferedData()
    for sensor in list_of_sensors:
        data.SetValue(f"{sensor.GetName()}/sensor", sensor)

def GetSensors(optimization_problem: OptimizationProblem) -> 'list[KratosSI.Sensors.Sensor]':
    data = ComponentDataView("sensors", optimization_problem).GetUnBufferedData()

    list_of_sensors: 'list[KratosSI.Sensors.Sensor]' = []
    for _, sensor_data in data.GetSubItems().items():
        list_of_sensors.append(sensor_data.GetValue("sensor"))

    return list_of_sensors

def AddSensorStatusControlUpdater(name: str, sensor_status_control_updater: typing.Any, optimization_problem: OptimizationProblem) -> None:
    if not hasattr(sensor_status_control_updater, "Update"):
        raise RuntimeError(f"The sensor status control update {sensor_status_control_updater} does not have the Update method.")

    data = ComponentDataView("sensors", optimization_problem).GetUnBufferedData()
    if not data.HasValue("sensor_status_control_updaters"):
        data.SetValue("sensor_status_control_updaters", [])
    list_of_sensor_status_control_updaters: 'list[tuple[str, typing.Any]]' = data.GetValue("sensor_status_control_updaters")
    list_of_sensor_status_control_updaters.append((name, sensor_status_control_updater))

def GetSensorStatusControlUpdater(query_name: str, optimization_problem: OptimizationProblem) -> typing.Any:
    for updater_name, updater in GetListOfSensorStatusControlUpdaters(optimization_problem):
        if updater_name == query_name:
            return updater
    raise RuntimeError(f"The queried sensor status control updater with name = \"{query_name}\" not found.")

def GetListOfSensorStatusControlUpdaters(optimization_problem: OptimizationProblem) -> 'list[tuple[str, typing.Any]]':
    data = ComponentDataView("sensors", optimization_problem).GetUnBufferedData()
    if data.HasValue("sensor_status_control_updaters"):
        return data.GetValue("sensor_status_control_updaters")
    else:
        return []

def HasSensorStatusControlUpdater(query_name: str, optimization_problem: OptimizationProblem) -> bool:
    for updater_name, _ in GetListOfSensorStatusControlUpdaters(optimization_problem):
        if updater_name == query_name:
            return True
    return False

def UpdateSensorStatusControlUpdaters(optimization_problem: OptimizationProblem) -> None:
    for _, updater in GetListOfSensorStatusControlUpdaters(optimization_problem):
        updater.Update()

