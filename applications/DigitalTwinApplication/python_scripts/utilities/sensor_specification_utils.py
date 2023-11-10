import typing
from pathlib import Path

import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA
import KratosMultiphysics.DigitalTwinApplication as KratosDT

def GetSpecifications(model_part: Kratos.ModelPart, list_of_parameters: 'list[Kratos.Parameters]') -> 'list[KratosDT.Sensors.SensorSpecification]':
    point_locator = Kratos.BruteForcePointLocator(model_part)

    list_of_specifications: 'list[KratosDT.Sensors.SensorSpecification]' = []
    shape_funcs = Kratos.Vector()
    for parameters in list_of_parameters:
        sensor_id = parameters["id"].GetInt()
        sensor_name = parameters["name"].GetString()
        location = Kratos.Point(parameters["location"].GetVector())

        spec = KratosDT.Sensors.SensorSpecification(sensor_name, sensor_id)
        spec.SetLocation(location)
        spec.SetValue(KratosDT.SENSOR_ELEMENT_ID, point_locator.FindElement(location, shape_funcs, Kratos.Configuration.Initial, 1e-8))

        if parameters.Has("variable_data"):
            var_value: Kratos.Parameters
            for var_name, var_value in parameters["variable_data"].items():
                var = Kratos.KratosGlobals.GetVariable(var_name)
                if isinstance(var, Kratos.BoolVariable):
                    spec.SetValue(var, var_value.GetBool())
                elif isinstance(var, Kratos.IntegerVariable):
                    spec.SetValue(var, var_value.GetInt())
                elif isinstance(var, Kratos.DoubleVariable):
                    spec.SetValue(var, var_value.GetDouble())
                elif isinstance(var, Kratos.Array1DVariable3):
                    vec = var_value.GetVector()
                    spec.SetValue(var, Kratos.Array3([vec[0], vec[1], vec[2]]))
                elif isinstance(var, Kratos.Array1DVariable4):
                    vec = var_value.GetVector()
                    spec.SetValue(var, Kratos.Array4([vec[0], vec[1], vec[2], vec[3]]))
                elif isinstance(var, Kratos.Array1DVariable6):
                    vec = var_value.GetVector()
                    spec.SetValue(var, Kratos.Array6([vec[0], vec[1], vec[2], vec[3], vec[4], vec[6]]))
                elif isinstance(var, Kratos.Array1DVariable9):
                    vec = var_value.GetVector()
                    spec.SetValue(var, Kratos.Array9([vec[0], vec[1], vec[2], vec[3], vec[4], vec[6], vec[7], vec[8]]))
                elif isinstance(var, Kratos.VectorVariable):
                    spec.SetValue(var, var_value.GetVector())

        list_of_specifications.append(spec)

    return list_of_specifications

def PrintSpecificationDataToCSV(list_of_spec_views: 'list[typing.Union[KratosDT.Sensors.NodalSensorSpecificationView, KratosDT.Sensors.ConditionSensorSpecificationView, KratosDT.Sensors.ElementSensorSpecificationView]]', output_file_name: Path) -> None:
    number_of_clusters = len(list_of_spec_views)

    # do nothing if number of clusters is zero
    if number_of_clusters == 0:
        return

    output_file_name.parent.mkdir(exist_ok=True, parents=True)

    with open(str(output_file_name), "w") as csv_output:
        # write the headers of spec
        csv_output.write("#; name; value; location_x; location_y; location_z")
        # now write the headers of data containers
        list_of_vars = []
        for var_name in list_of_spec_views[0].GetSensorSpecification().GetDataVariableNames():
            csv_output.write(f"; {var_name}")
            list_of_vars.append(Kratos.KratosGlobals.GetVariable(var_name))
        csv_output.write("\n")

        # now write the data
        sensor_id = 0
        for spec_view in list_of_spec_views:
            spec = spec_view.GetSensorSpecification()
            sensor_id += 1
            loc = spec.GetLocation()
            csv_output.write(f"{sensor_id}; {spec.GetName()}; {spec.GetSensorValue()}; {loc[0]}; {loc[1]}; {loc[2]}")
            for var in list_of_vars:
                csv_output.write(f"; {spec.GetValue(var)}")
            csv_output.write("\n")

def GetNormalizedSpecifications(list_of_spec_views: 'list[typing.Union[KratosDT.Sensors.NodalSensorSpecificationView, KratosDT.Sensors.ConditionSensorSpecificationView, KratosDT.Sensors.ElementSensorSpecificationView]]') -> 'list[typing.Union[KratosDT.Sensors.NodalSensorSpecificationView, KratosDT.Sensors.ConditionSensorSpecificationView, KratosDT.Sensors.ElementSensorSpecificationView]]':
    result: 'list[typing.Union[KratosDT.Sensors.NodalSensorSpecificationView, KratosDT.Sensors.ConditionSensorSpecificationView, KratosDT.Sensors.ElementSensorSpecificationView]]' = []
    for spec_view in list_of_spec_views:
        cexp = spec_view.GetContainerExpression()
        norm_l2 = KratosOA.ExpressionUtils.NormL2(cexp)
        if norm_l2 > 0.0:
            cexp.SetExpression(cexp.GetExpression() / norm_l2)
            cexp.SetExpression(cexp.Flatten().GetExpression())
            result.append(spec_view)
    return result

def GetCosineDistances(list_of_spec_views: 'list[typing.Union[KratosDT.Sensors.NodalSensorSpecificationView, KratosDT.Sensors.ConditionSensorSpecificationView, KratosDT.Sensors.ElementSensorSpecificationView]]') -> 'list[float]':
    """Get the cosine distances in a flat array.

    The distances are computed as follows

         d(u, v) = 1.0 - u\cdot v

    The vectors in u, v must be unit vectors

    Args:
        list_of_specification_views (list[KratosDT.Sensors.SensorSpecificationView]): List of specification views

    Returns:
        list[float]: d(u, v) calculated for all u and v
    """
    results: 'list[float]' = []
    for i, spec_i in enumerate(list_of_spec_views):
        for spec_j in list_of_spec_views[i+1:]:
            results.append(1.0 - KratosOA.ExpressionUtils.InnerProduct(spec_i.GetContainerExpression(), spec_j.GetContainerExpression()))
    return results

def AddSensorSpecificationVariableData(specification: KratosDT.Sensors.SensorSpecification, variable_data: Kratos.Parameters) -> None:
    for var_name, var_value in variable_data:
        var = Kratos.KratosGlobals.GetVariable(var_name)
        if isinstance(var, Kratos.BoolVariable):
            specification.SetValue(var, var_value.GetBool())
        elif isinstance(var, Kratos.IntegerVariable):
            specification.SetValue(var, var_value.GetInt())
        elif isinstance(var, Kratos.DoubleVariable):
            specification.SetValue(var, var_value.GetDouble())
        elif isinstance(var, Kratos.StringVariable):
            specification.SetValue(var, var_value.GetString())
        elif isinstance(var, Kratos.Array1DVariable3):
            value = var_value.GetVector()
            specification.SetValue(var, [value[0], value[1], value[2]])
        elif isinstance(var, Kratos.Array1DVariable4):
            value = var_value.GetVector()
            specification.SetValue(var, [value[0], value[1], value[2], value[3]])
        elif isinstance(var, Kratos.Array1DVariable6):
            value = var_value.GetVector()
            specification.SetValue(var, [value[0], value[1], value[2], value[3], value[4], value[5]])
        elif isinstance(var, Kratos.Array1DVariable9):
            value = var_value.GetVector()
            specification.SetValue(var, [value[0], value[1], value[2], value[3], value[4], value[5], value[6], value[7], value[8]])
        elif isinstance(var, Kratos.VectorVariable):
            specification.SetValue(var, var_value.GetVector())
