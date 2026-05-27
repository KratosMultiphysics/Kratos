import typing
import KratosMultiphysics as Kratos
from KratosMultiphysics.OptimizationApplication.utilities.buffered_dict import BufferedDict
from KratosMultiphysics.OptimizationApplication.utilities.component_data_view import ComponentDataView
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem

class AlgorithmDataManager:
    def __init__(self, optimization_problem: OptimizationProblem):
        self.__optimization_problem = optimization_problem
        self.__list_of_control_names: 'list[str]' = ComponentDataView("algorithm", self.__optimization_problem).GetUnBufferedData().GetValue("controls").split(",")

    def HasValue(self, key: str, is_historical = True, step_index = 0) -> bool:
        data_view = self.__GetBufferedDataDict(is_historical)
        return data_view.HasValue(key, step_index)

    def SetValue(self, key: str, value: 'typing.Union[bool, int, float, str, Kratos.TensorAdaptors.DoubleCombinedTensorAdaptor]', is_historical = True, step_index = 0, overwrite = False) -> None:
        data_view = self.__GetBufferedDataDict(is_historical)

        if isinstance(value, Kratos.TensorAdaptors.DoubleCombinedTensorAdaptor):
            for control_name, tensor_adaptor in zip(self.__list_of_control_names, value.GetTensorAdaptors()):
                data_view.SetValue(f"{key}/{control_name}_{key}", tensor_adaptor, step_index, overwrite)
        else:
            data_view.SetValue(key, value, step_index, overwrite)

    def GetValue(self, key: str, is_historical = True, step_index = 0) -> 'typing.Union[bool, int, float, str, Kratos.TensorAdaptors.DoubleCombinedTensorAdaptor]':
        data_view = self.__GetBufferedDataDict(is_historical)
        value = data_view.GetValue(key, step_index)

        if isinstance(value, BufferedDict):
            list_of_tensor_adaptors = [data_view.GetValue(f"{key}/{control_name}_{key}", step_index) for control_name in self.__list_of_control_names]
            ta = Kratos.TensorAdaptors.DoubleCombinedTensorAdaptor(list_of_tensor_adaptors, perform_collect_data_recursively=False, perform_store_data_recursively=False, copy=False)
            ta.CollectData()
            return ta
        else:
            return value

    def __GetBufferedDataDict(self, is_historical: bool) -> BufferedDict:
        if is_historical:
            return ComponentDataView("algorithm", self.__optimization_problem).GetBufferedData()
        else:
            return ComponentDataView("algorithm", self.__optimization_problem).GetUnBufferedData()

    def __str__(self) -> str:
        return f"Buffered dict:\n{ComponentDataView("algorithm", self.__optimization_problem).GetBufferedData()}\nUnbuffered dict:\n{ComponentDataView("algorithm", self.__optimization_problem).GetUnBufferedData()}"
