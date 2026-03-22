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

class TensorAdaptorBoundingManager:
    def __init__(self, bounds: 'list[float]') -> None:
        """
        This class is used to bound values of an tensor adaptor to specified interval.
            If a value in an tensor adaptor is above the max value of bounds, then it will have a value larger than 1.0
            if a value in an tensor adaptor is below the min value of bounds, then it will have a value smaller than 0.0
            Values in an tensor adaptor will be linearly interpolated and mapped from interval [min(bounds), max(bounds)]
            to [0, 1].

        Args:
            bounds (list[float]): Bounds of the given tensor adaptor values.

        Raises:
            RuntimeError: If bounds does not contain two values.
        """
        if len(bounds) != 2:
            raise RuntimeError(f"The bounds should be of size 2. [bounds = {bounds}]")
        self.bounds = sorted(bounds)

    def GetBoundGap(self) -> float:
        return self.bounds[1] - self.bounds[0]

    def GetBoundedTensorAdaptor(self, unbounded_tensor_adaptor: Kratos.TensorAdaptors.DoubleTensorAdaptor) -> Kratos.TensorAdaptors.DoubleTensorAdaptor:
        result = unbounded_tensor_adaptor.Clone()
        result.data[:] = (result.data[:] - self.bounds[0]) / self.GetBoundGap()
        return result

    def GetUnboundedTensorAdaptor(self, bounded_tensor_adaptor: Kratos.TensorAdaptors.DoubleTensorAdaptor) -> Kratos.TensorAdaptors.DoubleTensorAdaptor:
        result = bounded_tensor_adaptor.Clone()
        result.data[:] = (result.data[:] * self.GetBoundGap() + self.bounds[0])
        return result
