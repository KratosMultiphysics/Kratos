import typing
from pathlib import Path

import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
import KratosMultiphysics.SystemIdentificationApplication as KratosSI
from KratosMultiphysics.SystemIdentificationApplication.utilities.data_utils import SupportedValueUnionType
from KratosMultiphysics.SystemIdentificationApplication.utilities.data_utils import SupportedVariableUnionType
from KratosMultiphysics.SystemIdentificationApplication.utilities.data_utils import GetParameterToKratosValuesConverter
from KratosMultiphysics.SystemIdentificationApplication.utilities.data_utils import GetKratosValueToCSVStringConverter
from KratosMultiphysics.SystemIdentificationApplication.utilities.data_utils import GetNameToCSVString
from KratosMultiphysics.SystemIdentificationApplication.utilities.expression_utils import GetContainerExpressionType
from KratosMultiphysics.OptimizationApplication.utilities.component_data_view import ComponentDataView

def CreateSensors(sensor_model_part: Kratos.ModelPart, domain_model_part: Kratos.ModelPart, list_of_parameters: 'list[Kratos.Parameters]') -> 'list[KratosSI.Sensors.Sensor]':
    """Create list of sensors from given parameters.

    This creates a list of sensors from given list of parameters. Each parameter
    corresponds to one settings set for one sensor.

    Args:
        sensor_model_part (Kratos.ModelPart): Model part to which the sensors will be stored.
        domain_model_part (Kratos.ModelPart): Model part on which the sensors will be located.
        list_of_parameters (list[Kratos.Parameters]): List of parameters.

    Returns:
        list[KratosSI.Sensors.Sensor]: List of sensors generated.
    """
    point_locator = Kratos.BruteForcePointLocator(domain_model_part)

    if sensor_model_part.NumberOfNodes() != 0:
        raise RuntimeError(f"The sensor model part \"{sensor_model_part.FullName()}\" is not empty.")

    dict_of_sensor_types: 'dict[str, typing.Type[KratosSI.Sensors.Sensor]]' = {}
    for sensor_type_name in dir(KratosSI.Sensors):
        sensor_type = getattr(KratosSI.Sensors, sensor_type_name)
        if sensor_type_name.endswith("Sensor") and issubclass(sensor_type, KratosSI.Sensors.Sensor):
            dict_of_sensor_types[Kratos.StringUtilities.ConvertCamelCaseToSnakeCase(sensor_type_name)] = sensor_type

    list_of_sensors: 'list[KratosSI.Sensors.Sensor]' = []
    shape_funcs = Kratos.Vector()
    sensor_id = 0
    for parameters in list_of_parameters:
        sensor_id += 1

        if not parameters.Has("type"):
            raise RuntimeError(f"The sensor parameters does not contain \"type\".")
        sensor_type_name = parameters["type"].GetString()

        if not sensor_type_name in dict_of_sensor_types.keys():
            raise RuntimeError(f"Unsupported sensor type = \"{sensor_type_name}\" requested. Followings are supported:\n\t" + "\n\t".join(dict_of_sensor_types.keys()))
        parameters.ValidateAndAssignDefaults(dict_of_sensor_types[sensor_type_name].GetDefaultParameters())

        name = parameters["name"].GetString()
        loc = Kratos.Point(parameters["location"].GetVector())
        node = sensor_model_part.CreateNewNode(sensor_id, loc.X, loc.Y, loc.Z)
        weight = parameters["weight"].GetDouble()

        if sensor_type_name == "displacement_sensor":
            direction = parameters["direction"].GetVector()
            elem_id = point_locator.FindElement(loc, shape_funcs, Kratos.Configuration.Initial, 1e-8)
            sensor = KratosSI.Sensors.DisplacementSensor(name, node, direction, domain_model_part.GetElement(elem_id), weight)
            AddSensorVariableData(sensor, parameters["variable_data"])
            list_of_sensors.append(sensor)
        elif sensor_type_name == "strain_sensor":
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
            sensor = KratosSI.Sensors.StrainSensor(name, node, strain_variable, strain_type_value, domain_model_part.GetElement(elem_id), weight)
            AddSensorVariableData(sensor, parameters["variable_data"])
            list_of_sensors.append(sensor)
        else:
            raise RuntimeError(f"Unsupported sensor type name = \"{sensor_type_name}\".")
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

def SetSensors(sensor_group_data: ComponentDataView, list_of_sensors: 'list[KratosSI.Sensors.Sensor]') -> None:
    for sensor in list_of_sensors:
        sensor_group_data.GetUnBufferedData().SetValue(f"{sensor.GetName()}/sensor", sensor)

def GetSensors(sensor_group_data: ComponentDataView) -> 'list[KratosSI.Sensors.Sensor]':
    list_of_sensors: 'list[KratosSI.Sensors.Sensor]' = []
    for _, sensor_data in sensor_group_data.GetUnBufferedData().GetSubItems().items():
        list_of_sensors.append(sensor_data.GetValue("sensor"))
    return list_of_sensors
