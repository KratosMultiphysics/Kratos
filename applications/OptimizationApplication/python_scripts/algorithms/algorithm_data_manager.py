import typing
import KratosMultiphysics as Kratos
from KratosMultiphysics.OptimizationApplication.controls.control import Control
from KratosMultiphysics.OptimizationApplication.utilities.component_data_view import ComponentDataView
from KratosMultiphysics.OptimizationApplication.utilities.optimization_problem import OptimizationProblem

class AlgorithmDataManager:
    def __init__(self, optimization_problem: OptimizationProblem):
        self.__optimization_problem = optimization_problem
        self.__list_of_controls = [optimization_problem.GetControl(control_name) for control_name in ComponentDataView("algorithm").GetUnBufferedData().GetValue("controls").split(",")]

    def HasValue(self, key: str, is_historical = True, step_index = 0) -> bool:
        if is_historical:
            data_view = ComponentDataView("algorithm", self.__optimization_problem).GetBufferedData()
        else:
            data_view = ComponentDataView("algorithm", self.__optimization_problem).GetUnBufferedData()

        if data_view.HasValue(key, step_index):
            return True
        else:
            if is_historical:
                list_of_tensor_adaptors = [ComponentDataView(control, self.__optimization_problem).GetBufferedData().HasValue(f"{control.GetName()}_{key}", step_index) for control in self.__list_of_controls]
            else:
                list_of_tensor_adaptors = [ComponentDataView(control, self.__optimization_problem).GetUnBufferedData().HasValue(f"{control.GetName()}_{key}", step_index) for control in self.__list_of_controls]

            return all(list_of_tensor_adaptors)

    def SetValue(self, key: str, value: 'typing.Union[bool, int, float, str, Kratos.TensorAdaptors.DoubleCombinedTensorAdaptor]', is_historical = True, step_index = 0, overwrite = False) -> None:
        if type(value) in [bool, int, float, str]:
            if is_historical:
                ComponentDataView("algorithm", self.__optimization_problem).GetBufferedData().SetValue(key, value, step_index, overwrite)
            else:
                ComponentDataView("algorithm", self.__optimization_problem).GetUnBufferedData().SetValue(key, value, step_index, overwrite)
        elif isinstance(value, Kratos.TensorAdaptors.DoubleCombinedTensorAdaptor):
            for control, tensor_adaptor in zip(self.__list_of_controls, value.GetTensorAdaptors()):
                if is_historical:
                    ComponentDataView(control, self.__optimization_problem).GetBufferedData().SetValue(f"{control.GetName()}_{key}", tensor_adaptor, step_index, overwrite)
                else:
                    ComponentDataView(control, self.__optimization_problem).GetUnBufferedData().SetValue(f"{control.GetName()}_{key}", tensor_adaptor, step_index, overwrite)
        else:
            raise RuntimeError(f"Unsupported value type [ key = {key}, value = {value} ].")

    def GetValue(self, key: str, is_historical = True, step_index = 0) -> 'typing.Union[bool, int, float, str, Kratos.TensorAdaptors.DoubleCombinedTensorAdaptor]':
        if is_historical:
            data_view = ComponentDataView("algorithm", self.__optimization_problem).GetBufferedData()
        else:
            data_view = ComponentDataView("algorithm", self.__optimization_problem).GetUnBufferedData()

        if data_view.HasValue(key, step_index):
            return data_view.GetValue(key, step_index)
        else:
            if is_historical:
                list_of_tensor_adaptors = [ComponentDataView(control, self.__optimization_problem).GetBufferedData().GetValue(f"{control.GetName()}_{key}", step_index) for control in self.__list_of_controls]
            else:
                list_of_tensor_adaptors = [ComponentDataView(control, self.__optimization_problem).GetUnBufferedData().GetValue(f"{control.GetName()}_{key}", step_index) for control in self.__list_of_controls]

            return Kratos.TensorAdaptors.DoubleCombinedTensorAdaptor(list_of_tensor_adaptors, perform_collect_data_recursively=False, perform_store_data_recursively=False, copy=False)

    def __str__(self):
        pass
