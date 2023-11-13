import typing
from pathlib import Path

import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
import KratosMultiphysics.DigitalTwinApplication as KratosDT
from KratosMultiphysics.DigitalTwinApplication.utilities.data_utils import SensorViewUnionType
from KratosMultiphysics.DigitalTwinApplication.utilities.data_utils import SupportedValueUnionType
from KratosMultiphysics.DigitalTwinApplication.utilities.data_utils import SupportedVariableUnionType
from KratosMultiphysics.DigitalTwinApplication.utilities.data_utils import GetParameterToKratosValuesConverter
from KratosMultiphysics.DigitalTwinApplication.utilities.data_utils import GetKratosValueToCSVStringConverter
from KratosMultiphysics.DigitalTwinApplication.utilities.data_utils import GetNameToCSVString

def GetSensors(model_part: Kratos.ModelPart, list_of_parameters: 'list[Kratos.Parameters]') -> 'list[KratosDT.Sensors.Sensor]':
    """Get list of sensors from given parameters.

    This generates list of sensors from given list of parameters. Each parameter
    corresponds to one settings set for one sensor.

    Args:
        model_part (Kratos.ModelPart): Model part on which the sensors will be generated.
        list_of_parameters (list[Kratos.Parameters]): List of parameters.

    Returns:
        list[KratosDT.Sensors.Sensor]: List of sensors generated.
    """
    point_locator = Kratos.BruteForcePointLocator(model_part)

    defaults_map:'dict[str, Kratos.Parameters]' = {
            "displacement_sensor": Kratos.Parameters(
                                        """{
                                                "type"         : "displacement_sensor",
                                                "name"         : "",
                                                "value"        : 0,
                                                "location"     : [0.0, 0.0, 0.0],
                                                "direction"    : [0.0, 0.0, 0.0],
                                                "weight"       : 1.0,
                                                "variable_data": {}
                                        }""")
    }

    list_of_sensors: 'list[KratosDT.Sensors.Sensor]' = []
    shape_funcs = Kratos.Vector()
    for parameters in list_of_parameters:
        if not parameters.Has("type"):
            raise RuntimeError(f"The sensor parameters does not contain \"type\".")
        sensor_type = parameters["type"].GetString()

        if not sensor_type in defaults_map.keys():
            raise RuntimeError(f"Unsupported sensor type = \"{sensor_type}\" requested. Followings are supported:\n\t" + "\n\t".join(defaults_map.keys()))
        parameters.ValidateAndAssignDefaults(defaults_map[sensor_type])

        if sensor_type == "displacement_sensor":
            name = parameters["name"].GetString()
            loc = parameters["location"].GetVector()
            loc = Kratos.Point(loc[0], loc[1], loc[2])
            weight = parameters["weight"].GetDouble()
            direction = parameters["direction"].GetVector()
            elem_id = point_locator.FindElement(loc, shape_funcs, Kratos.Configuration.Initial, 1e-8)
            sensor = KratosDT.Sensors.DisplacementSensor(name, loc, direction, model_part.GetElement(elem_id), weight)
            AddSensorVariableData(sensor, parameters["variable_data"])
            list_of_sensors.append(sensor)
    return list_of_sensors

def PrintSensorViewsListToCSV(output_file_name: Path, list_of_sensors: 'list[KratosDT.Sensors.Sensor]', list_of_sensor_properties: 'list[str]') -> None:
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
                    value = first_sensor.GetValue(var)
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
                file_output.write(f",{data_value_dict[var](sensor.GetValue(var))}")

            file_output.write("\n")

def GetNormalizedSensorViews(list_of_sensor_views: 'list[SensorViewUnionType]') -> 'list[SensorViewUnionType]':
    """Get the normalized sensor views.

    This method returns sensor views which will contain normalized container expressions.
    The sensor views with zero norm will not be included in the return list

    Args:
        list_of_sensor_views (list[SensorViewUnionType]): List of sensors to normalize.

    Returns:
        list[SensorViewUnionType]: Normalized list of sensor views.
    """
    result: 'list[SensorViewUnionType]' = []
    for sensor_view in list_of_sensor_views:
        cexp = sensor_view.GetContainerExpression()
        norm_l2 = KratosOA.ExpressionUtils.NormL2(cexp)
        if norm_l2 > 0.0:
            cexp.SetExpression(cexp.GetExpression() / norm_l2)
            cexp.SetExpression(cexp.Flatten().GetExpression())
            result.append(sensor_view)
    return result

def GetCosineDistances(list_of_sensor_views: 'list[SensorViewUnionType]') -> 'list[float]':
    """Get the cosine distances in a flat array.

    The distances are computed as follows

         d(u, v) = 1.0 - u\cdot v

    The vectors in u, v must be unit vectors.

    Args:
        list_of_sensor_views (list[SensorViewUnionType]): List of sensor views

    Returns:
        list[float]: d(u, v) calculated for all u and v
    """
    results: 'list[float]' = []
    for i, sensor_view_i in enumerate(list_of_sensor_views):
        for sensor_view_j in list_of_sensor_views[i+1:]:
            results.append(1.0 - KratosOA.ExpressionUtils.InnerProduct(sensor_view_i.GetContainerExpression(), sensor_view_j.GetContainerExpression()))
    return results

def AddSensorVariableData(sensor: KratosDT.Sensors.Sensor, variable_data: Kratos.Parameters) -> None:
    """Adds sensor variable data.

    Args:
        sensor (KratosDT.Sensors.Sensor): Sensor to add data.
        variable_data (Kratos.Parameters): Parameters containing variable data.
    """
    for var_name, var_value in variable_data.items():
        var = Kratos.KratosGlobals.GetVariable(var_name)
        value_func =  GetParameterToKratosValuesConverter(var_value)
        sensor.SetValue(var, value_func(var_value))

def GetSensorTypeDict(list_of_sensors: 'list[SensorViewUnionType]') -> 'dict[str, list[KratosDT.Sensors.Sensor]]':
    """Get the sensor type dict from list of sensors with keys being type of the sensor.

    Args:
        list_of_sensors (list[KratosDT.Sensors.Sensor]): List of sensors

    Returns:
        dict[str, list[KratosDT.Sensors.Sensor]]: Sensor type dict.
    """
    result: 'dict[str, list[KratosDT.Sensors.Sensor]]' = {}
    for sensor in list_of_sensors:
        sensor_type = type(sensor.GetSensor()).__name__
        if sensor_type not in result.keys():
            result[sensor_type] = []
        result[sensor_type].append(sensor)
    return result

