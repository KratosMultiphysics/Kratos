import typing

import KratosMultiphysics as Kratos
import KratosMultiphysics.SystemIdentificationApplication as KratosSI

SensorViewUnionType = typing.Union[
                            KratosSI.Sensors.NodalSensorView,
                            KratosSI.Sensors.ConditionSensorView,
                            KratosSI.Sensors.ElementSensorView
                        ]

SupportedVariableUnionType = typing.Union[
                                    Kratos.BoolVariable,
                                    Kratos.IntegerVariable,
                                    Kratos.DoubleVariable,
                                    Kratos.StringVariable,
                                    Kratos.Array1DVariable3,
                                    Kratos.Array1DVariable4,
                                    Kratos.Array1DVariable6,
                                    Kratos.Array1DVariable9,
                                    Kratos.VectorVariable
                                ]

SupportedValueUnionType = typing.Union[
                                    bool,
                                    int,
                                    float,
                                    str,
                                    Kratos.Array3,
                                    Kratos.Array4,
                                    Kratos.Array6,
                                    Kratos.Array9,
                                    Kratos.Vector
                                ]


def GetNameToCSVString(name: str, kratos_value: SupportedValueUnionType) -> str:
    if type(kratos_value) in [bool, int, float, str]:
        return name
    elif type(kratos_value) in [Kratos.Array3, Kratos.Array4, Kratos.Array6, Kratos.Array9, Kratos.Vector]:
        return ",".join([f"{name}_{i}" for i in range(kratos_value.Size())])
    else:
        raise RuntimeError(f"Unsupported kratos_value type = {type(kratos_value)} with name = {name} and kratos_value = {kratos_value}.")

def GetKratosValueToPythonValueConverter(kratos_value: SupportedValueUnionType) -> 'typing.Callable[[SupportedValueUnionType], typing.Union[bool, int, float, str, list[float]]]':
    if type(kratos_value) in [bool, int, float, str]:
        return lambda x: x
    elif type(kratos_value) in [Kratos.Array3, Kratos.Array4, Kratos.Array6, Kratos.Array9, Kratos.Vector]:
        return lambda x: [v for v in x]
    else:
        raise RuntimeError(f"Unsupported kratos_value type = {type(kratos_value)} with kratos_value = {kratos_value}.")

def GetKratosValueToCSVStringConverter(kratos_value: SupportedValueUnionType) -> 'typing.Callable[[SupportedValueUnionType], str]':
    if type(kratos_value) in [bool, int, float, str]:
        return lambda x: str(x)
    elif type(kratos_value) in [Kratos.Array3, Kratos.Array4, Kratos.Array6, Kratos.Array9, Kratos.Vector]:
        return lambda x: ",".join([str(v) for v in x])
    else:
        raise RuntimeError(f"Unsupported kratos_value type = {type(kratos_value)} with kratos_value = {kratos_value}.")

def GetParameterToKratosValuesConverter(params: Kratos.Parameters) -> 'typing.Callable[[Kratos.Parameters], SupportedValueUnionType]':
    if params.IsVector():
        return lambda x: x.GetVector()
    elif params.IsBool():
        return lambda x: x.GetBool()
    elif params.IsInt():
        return lambda x: x.GetInt()
    elif params.IsDouble():
        return lambda x: x.GetDouble()
    elif params.IsString():
        return lambda x: x.GetString()
    else:
        raise RuntimeError(f"Unsupported parameters type = {params}.")