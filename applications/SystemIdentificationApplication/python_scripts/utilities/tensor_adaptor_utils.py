import typing
from enum import Enum
import KratosMultiphysics as Kratos
import KratosMultiphysics.OptimizationApplication as KratosOA

class PropertiesDataLocation(Enum):
    ElementProperties = 100
    ConditionProperties = 200

class TensorAdaptorDataLocation(Kratos.Globals.DataLocation):
    ElementProperties = PropertiesDataLocation.ElementProperties,
    ConditionProperties = PropertiesDataLocation.ConditionProperties

def GetTensorAdaptor(model_part: Kratos.ModelPart, data_location: TensorAdaptorDataLocation, variable: typing.Any) -> Kratos.TensorAdaptors.DoubleTensorAdaptor:
    if data_location == TensorAdaptorDataLocation.NodeHistorical:
        ta = Kratos.TensorAdaptors.HistoricalVariableTensorAdaptor(model_part.Nodes, variable)
    elif data_location == TensorAdaptorDataLocation.NodeNonHistorical:
        ta = Kratos.TensorAdaptors.VariableTensorAdaptor(model_part.Nodes, variable)
    elif data_location == TensorAdaptorDataLocation.Condition:
        ta = Kratos.TensorAdaptors.VariableTensorAdaptor(model_part.Conditions, variable)
    elif data_location == TensorAdaptorDataLocation.Element:
        ta = Kratos.TensorAdaptors.VariableTensorAdaptor(model_part.Elements, variable)
    elif data_location == TensorAdaptorDataLocation.ConditionProperties:
        ta = KratosOA.TensorAdaptors.PropertiesVariableTensorAdaptor(model_part.Conditions, variable)
    elif data_location == TensorAdaptorDataLocation.ElementProperties:
        ta = KratosOA.TensorAdaptors.PropertiesVariableTensorAdaptor(model_part.Elements, variable)
    else:
        raise RuntimeError(f"Unsupported {data_location}.")

    ta.CollectData()
    return ta
