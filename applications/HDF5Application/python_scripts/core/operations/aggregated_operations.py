import typing

import KratosMultiphysics as Kratos
import KratosMultiphysics.HDF5Application as KratosHDF5

from KratosMultiphysics.HDF5Application.core.file_io import OpenHDF5File

class ControlledOperation:
    def __init__(self,
                 operation_type: 'typing.Type[Kratos.Operation]',
                 parameters: Kratos.Parameters,
                 controller: Kratos.Controller,
                 *args) -> None:
        self.__parameters = parameters
        self.__operation_type = operation_type
        self.__controller = controller
        self.__other_args = args

        if self.__parameters is not None:
            self.__parameters.AddMissingParameters(self.__operation_type.GetDefaultParameters())

    def Check(self) -> None:
        self.__controller.Check()

    def Evaluate(self) -> bool:
        return self.__controller.Evaluate()

    def Execute(self, model_part: Kratos.ModelPart, hdf5_file: KratosHDF5.HDF5File) -> None:
        if self.Evaluate():
            self.__operation_type(model_part, self.__parameters, hdf5_file, *self.__other_args).Execute()

    def Update(self) -> None:
        self.__controller.Update()

    def __str__(self) -> str:
        return f"{self.__class__.__name__}: operation type = {self.__operation_type}, controller = {self.__controller}, parameters = \n{self.__parameters}"

class AggregatedControlledOperations:
    def __init__(self,
                 model_part: Kratos.ModelPart,
                 hdf5_file_settings: Kratos.Parameters) -> None:
        self.__model_part = model_part
        self.__hdf5_file_parameters = hdf5_file_settings
        self.__list_of_controlled_operations: 'list[ControlledOperation]' = []

    def Check(self):
        list(map(lambda x: x.Check(), self.__list_of_controlled_operations))

    def Evaluate(self) -> bool:
        return any([controlled_operation.Evaluate() for controlled_operation in self.__list_of_controlled_operations])

    def Execute(self) -> None:
        with OpenHDF5File(self.__hdf5_file_parameters, self.__model_part) as h5_file:
            list(map(lambda x: x.Execute(self.__model_part, h5_file) , self.__list_of_controlled_operations))

    def Update(self) -> None:
        list(map(lambda x: x.Update() , self.__list_of_controlled_operations))

    def AddControlledOperation(self, controlled_operation: ControlledOperation) -> None:
        self.__list_of_controlled_operations.append(controlled_operation)

    def GetListOfControlledOperations(self) -> 'list[ControlledOperation]':
        return self.__list_of_controlled_operations

    def __str__(self) -> str:
        msg = f"{self.__class__.__name__}: model part name = {self.__model_part.FullName()}, hdf5_file_params =\n{self.__hdf5_file_parameters}"
        msg += "List of aggregated operations:" + "\n\t --- ".join(str(operation) for operation in self.__list_of_controlled_operations)
        return msg