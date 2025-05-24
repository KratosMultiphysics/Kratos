'''HDF5 model part operations.

license: HDF5Application/license.txt
'''


# --- Code Imports ---
import KratosMultiphysics

# --- HDF5 Imports ---
import KratosMultiphysics.HDF5Application as KratosHDF5
from ..utils import ParametersWrapper
from ..file_io import OpenHDF5File

# --- STD Imports ---
from importlib import import_module
import sys
import typing
#import abc


def Prefix(pattern, model_part, time_format=''):
    if hasattr(model_part, 'ProcessInfo'):
        time = model_part.ProcessInfo[KratosMultiphysics.TIME]
        prefix = format(time, time_format).join(pattern.split('<time>'))
        if KratosMultiphysics.STEP in model_part.ProcessInfo:
            prefix = prefix.replace('<step>', str(model_part.ProcessInfo[KratosMultiphysics.STEP]))
        else:
            # to be removed once analysis stage sets the STEP variable.
            prefix = prefix.replace('<step>', "0")
    else:
        prefix = pattern
    if hasattr(model_part, 'Name'):
        prefix = prefix.replace('<model_part_name>', model_part.Name)
    return prefix


class IOOperation(KratosMultiphysics.Operation): # IOOperation(KratosMultiphysics.Operation, metaclass = abc.ABCMeta) # <== conflicting metaclasses
    """ @brief Base class for HDF5 IO operations."""

    def __init__(self,
                 model_part: KratosMultiphysics.ModelPart,
                 parameters: KratosMultiphysics.Parameters,
                 file: KratosHDF5.HDF5File):
        super().__init__()
        self.__model_part = model_part
        self.__file = file
        self.__time_format = parameters["time_format"].GetString()

        # Pick necessary entries from the input parameters
        # (assuming they're all present)
        self.__parameters = KratosMultiphysics.Parameters()
        self.__parameters.AddString("prefix", Prefix(parameters["prefix"].GetString(), self.model_part, self.time_format))

    #@abc.abstractmethod # <== missing the ABCMeta metaclass
    def Execute(self) -> None:
        raise NotImplementedError("Call to a pure abstract method")

    #@abc.abstractclassmethod # <== missing the ABCMeta metaclass
    @classmethod
    def GetDefaultParameters(cls) -> KratosMultiphysics.Parameters:
        return KratosMultiphysics.Parameters("""{
            "operation_type" : "",
            "prefix" : "/",
            "time_format" : "0.4f"
        }""")

    @property
    def model_part(self) -> KratosMultiphysics.ModelPart:
        return self.__model_part

    @property
    def parameters(self) -> KratosMultiphysics.Parameters:
        return self.__parameters

    @property
    def file(self) -> KratosHDF5.HDF5File:
        return self.__file

    @property
    def time_format(self) -> str:
        return self.__time_format

    @property
    def prefix(self) -> str:
        return Prefix(self.parameters["prefix"].GetString(),
                      self.model_part,
                      self.time_format)


class ModelPartIOOperation(IOOperation):
    """ @brief Base class for HDF5 IO operations on @ref ModelPart s."""

    @classmethod
    def GetDefaultParameters(cls: "typing.Type[ModelPartIOOperation]") -> KratosMultiphysics.Parameters:
        parameters = super().GetDefaultParameters()
        parameters["prefix"].SetString("/ModelData")
        return parameters


class ModelPartInput(ModelPartIOOperation):
    '''Reads a @ref ModelPart from a file.'''

    def Execute(self) -> None:
        KratosHDF5.HDF5ModelPartIO(self.file, self.prefix).ReadModelPart(self.model_part)


class ModelPartOutput(ModelPartIOOperation):
    '''Writes a @ref ModelPart to a file.'''

    def Execute(self) -> None:
        KratosHDF5.HDF5ModelPartIO(self.file, self.prefix).WriteModelPart(self.model_part)


class PartitionedModelPartOutput(ModelPartIOOperation):
    '''Writes a partitioned model part to a file.'''

    def Execute(self) -> None:
        KratosHDF5.HDF5PartitionedModelPartIO(self.file, self.prefix).WriteModelPart(self.model_part)


class ProcessInfoIOOperation(IOOperation):

    @classmethod
    def GetDefaultParameters(cls: "typing.Type[ProcessInfoIOOperation]") -> KratosMultiphysics.Parameters:
        parameters = super().GetDefaultParameters()
        parameters["prefix"].SetString("/ProcessInfo")
        return parameters


class ProcessInfoOutput(ProcessInfoIOOperation):
    '''Writes a @ref ProcessInfo to a file.'''

    def Execute(self) -> None:
        KratosHDF5.WriteDataValueContainer(self.file, self.prefix, self.model_part.ProcessInfo)


class ProcessInfoInput(ProcessInfoIOOperation):
    '''Reads a @ref ProcessInfo from a file.'''

    def Execute(self) -> None:
        KratosHDF5.ReadDataValueContainer(self.file, self.prefix, self.model_part.ProcessInfo)


class VariableIOOperation(IOOperation):
    '''Generates json settings for variable data IO.'''

    def __init__(self,
                 model_part: KratosMultiphysics.ModelPart,
                 parameters: KratosMultiphysics.Parameters,
                 file: KratosHDF5.HDF5File):
        super().__init__(model_part, parameters, file)
        self.parameters.AddValue("list_of_variables", parameters["list_of_variables"])

    @classmethod
    def GetDefaultParameters(cls: "typing.Type[VariableIOOperation]") -> KratosMultiphysics.Parameters:
        parameters = super().GetDefaultParameters()
        parameters["prefix"].SetString("/ResultsData")
        parameters.AddEmptyArray("list_of_variables")
        return parameters


class ElementDataValueOutput(VariableIOOperation):
    '''Writes non-historical element data values to a file.'''

    def Execute(self) -> None:
        KratosHDF5.HDF5ElementDataValueIO(
            self.parameters,
            self.file).WriteElementResults(self.model_part.Elements)


class ElementDataValueInput(VariableIOOperation):
    '''Reads non-historical element data values from a file.'''

    def Execute(self) -> None:
        KratosHDF5.HDF5ElementDataValueIO(
            self.parameters,
            self.file).ReadElementResults(self.model_part.Elements,
                                          self.model_part.GetCommunicator())


class ElementFlagValueOutput(VariableIOOperation):
    '''Writes non-historical element flag values to a file.'''

    def Execute(self) -> None:
        KratosHDF5.HDF5ElementFlagValueIO(
            self.parameters,
            self.file).WriteElementFlags(self.model_part.Elements)


class ElementFlagValueInput(VariableIOOperation):
    '''Reads non-historical element flag values from a file.'''

    def Execute(self) -> None:
        KratosHDF5.HDF5ElementFlagValueIO(
            self.parameters,
            self.file).ReadElementFlags(self.model_part.Elements,
                                        self.model_part.GetCommunicator())


class ElementGaussPointOutput(VariableIOOperation):
    '''Write element integration point values to a file.'''

    def Execute(self) -> None:
        KratosHDF5.HDF5ElementGaussPointIO(
            self.parameters,
            self.file).WriteElementGaussPointValues(self.model_part.Elements,
                                                    self.model_part.GetCommunicator().GetDataCommunicator(),
                                                    self.model_part.ProcessInfo)


class ConditionDataValueOutput(VariableIOOperation):
    '''Writes non-historical element data values to a file.'''

    def Execute(self) -> None:
        KratosHDF5.HDF5ConditionDataValueIO(
            self.parameters,
            self.file).WriteConditionResults(self.model_part.Conditions)


class ConditionDataValueInput(VariableIOOperation):
    '''Reads non-historical element data values from a file.'''

    def Execute(self) -> None:
        KratosHDF5.HDF5ConditionDataValueIO(
            self.parameters,
            self.file).ReadConditionResults(self.model_part.Conditions,
                                            self.model_part.GetCommunicator())


class ConditionFlagValueOutput(VariableIOOperation):
    '''Writes non-historical element flag values to a file.'''

    def Execute(self) -> None:
        KratosHDF5.HDF5ConditionFlagValueIO(
            self.parameters,
            self.file).WriteConditionFlags(self.model_part.Conditions)


class ConditionFlagValueInput(VariableIOOperation):
    '''Reads non-historical element flag values from a file.'''

    def Execute(self) -> None:
        KratosHDF5.HDF5ConditionFlagValueIO(
            self.parameters,
            self.file).ReadConditionFlags(self.model_part.Conditions,
                                          self.model_part.GetCommunicator())


class ConditionGaussPointOutput(VariableIOOperation):
    '''Write condition integration point values to a file.'''

    def Execute(self) -> None:
        KratosHDF5.HDF5ConditionGaussPointIO(
            self.parameters,
            self.file).WriteConditionGaussPointValues(self.model_part.Conditions,
                                                      self.model_part.GetCommunicator().GetDataCommunicator(),
                                                      self.model_part.ProcessInfo)


class NodalSolutionStepDataOutput(VariableIOOperation):
    '''Writes nodal solution step data to a file.'''

    def Execute(self) -> None:
        KratosHDF5.HDF5NodalSolutionStepDataIO(
            self.parameters,
            self.file).WriteNodalResults(self.model_part, 0)


class NodalSolutionStepDataInput(VariableIOOperation):
    '''Reads nodal solution step data from a file.'''

    def Execute(self) -> None:
        KratosHDF5.HDF5NodalSolutionStepDataIO(
            self.parameters,
            self.file).ReadNodalResults(self.model_part, 0)


class NodalDataValueOutput(VariableIOOperation):
    '''Writes non-historical nodal data values to a file.'''

    def Execute(self) -> None:
        KratosHDF5.HDF5NodalDataValueIO(
            self.parameters,
            self.file).WriteNodalResults(self.model_part.Nodes)


class NodalDataValueInput(VariableIOOperation):
    '''Reads non-historical nodal data values from a file.'''

    def Execute(self) -> None:
        KratosHDF5.HDF5NodalDataValueIO(
            self.parameters,
            self.file).ReadNodalResults(self.model_part.Nodes,
                                        self.model_part.GetCommunicator())


class NodalFlagValueOutput(VariableIOOperation):
    '''Writes non-historical nodal flag values to a file.'''

    def Execute(self) -> None:
        KratosHDF5.HDF5NodalFlagValueIO(
            self.parameters,
            self.file).WriteNodalFlags(self.model_part.Nodes)


class NodalFlagValueInput(VariableIOOperation):
    '''Reads non-historical nodal flag values from a file.'''

    def Execute(self) -> None:
        KratosHDF5.HDF5NodalFlagValueIO(
            self.parameters,
            self.file).ReadNodalFlags(self.model_part.Nodes,
                                      self.model_part.GetCommunicator())


class PrimalBossakOutput(VariableIOOperation):
    '''Writes nodal solution step data to a file for Bossak time schemes.

    Behaves the same as NodalSolutionStepDataOutput except for ACCELERATION,
    which is computed as
    (1 - alpha_bossak) * node.GetSolutionStepValue(ACCELERATION, 0) +
          alpha_bossak * node.GetSolutionStepValue(ACCELERATION, 1)
    and written to the file for the current time step. This is used by the
    transient adjoint solvers.
    '''

    def __init__(self,
                 model_part: KratosMultiphysics.ModelPart,
                 parameters: KratosMultiphysics.Parameters,
                 file: KratosHDF5.HDF5File):
        self.__alpha_bossak = parameters["alpha_bossak"].GetDouble()
        super().__init__(model_part, parameters, file)

    @property
    def alpha_bossak(self) -> float:
        return self.__alpha_bossak

    def Execute(self) -> None:
        primal_io = KratosHDF5.HDF5NodalSolutionStepBossakIO(self.parameters, self.file)
        primal_io.SetAlphaBossak(self.__alpha_bossak)
        primal_io.WriteNodalResults(self.model_part)

    @classmethod
    def GetDefaultParameters(cls: "typing.Type[PrimalBossakOutput]") -> KratosMultiphysics.Parameters:
        parameters = super().GetDefaultParameters()
        parameters.AddDouble("alpha_bossak", -0.3)
        return parameters


class PrimalBossakInput(VariableIOOperation):
    '''Reads nodal solution step data from a file.

    This is used by the transient adjoint solvers.
    '''

    def Execute(self) -> None:
        KratosHDF5.HDF5NodalSolutionStepBossakIO(
            self.parameters,
            self.file).ReadNodalResults(self.model_part)


class MoveMesh(ModelPartIOOperation):
    '''Perform a mesh move operation on a model part.

    The primary use case is to set the mesh to the current configuration after
    reading the model part.
    '''

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        print("HDF5 MoveMesh operation is deprecated and will be removed", file = sys.stderr)

    def Execute(self) -> None:
        KratosMultiphysics.ImplicitSolvingStrategy(self.model_part, True).MoveMesh()


def GetSubclasses(base_class: type) -> "list[type]":
    """Recursively find all subclasses of a base class"""
    subclasses = base_class.__subclasses__()
    for subclass in base_class.__subclasses__():
        subclasses += GetSubclasses(subclass)
    return subclasses


class AggregateOperation(KratosMultiphysics.Operation):
    """ @brief Class for aggregating HDF5 IO operations on the same file."""

    def __init__(self, model: KratosMultiphysics.Model, parameters: typing.Union[KratosMultiphysics.Parameters, ParametersWrapper]):
        super().__init__()
        if isinstance(parameters, ParametersWrapper):
            parameters = parameters.Get()
        parameters.AddMissingParameters(self.GetDefaultParameters())

        self.__model_part = model.GetModelPart(parameters["model_part_name"].GetString())
        self.__io_parameters = parameters["io_settings"]

        # {operation_type, operation_parameters}
        self.__operations: "list[tuple[type, KratosMultiphysics.Parameters]]" = []
        for i in range(parameters["list_of_operations"].size()):
            self.__Add(parameters["list_of_operations"][i])

    def Execute(self) -> None:
        with OpenHDF5File(self.__io_parameters, self.__model_part) as file:
            for operation, operation_parameters in self.__operations:
                operation(self.__model_part, operation_parameters, file).Execute()

    def __Add(self, operation_parameters: KratosMultiphysics.Parameters) -> None:
        # Convert input snake case name to the internal camel case name
        operation_type = KratosMultiphysics.StringUtilities.ConvertSnakeCaseToCamelCase(operation_parameters["operation_type"].GetString())
        operation: typing.Type[IOOperation] = next((op for op in GetSubclasses(IOOperation) if op.__name__ == operation_type), None)

        # Found an operation with a matching name
        if not (operation is None):
            operation_parameters.AddMissingParameters(operation.GetDefaultParameters())
            self.__operations.append((operation, operation_parameters))

        # No operation with the specified name in this scope
        else:
            if operation_parameters.Has("module_name"):
                module_name = operation_parameters["module_name"].GetString()
                module = import_module(f"KratosMultiphysics.HDF5Application.core.{module_name}")
                operation = module.Create(operation_parameters)

            if operation is None:
                raise ValueError(f"Invalid operation type '{operation_parameters['operation_type'].GetString()}'. Available options: {[KratosMultiphysics.StringUtilities.ConvertCamelCaseToSnakeCase(op.__name__) for op in GetSubclasses(IOOperation)]}")
            else:
                self.__operations.append((operation, operation_parameters))

    @classmethod
    def GetDefaultParameters(cls: "typing.Type[AggregateOperation]") -> KratosMultiphysics.Parameters:
        return KratosMultiphysics.Parameters("""{
            "model_part_name" : "",
            "list_of_operations" : [],
            "io_settings" : {}
        }""")


def Create(model: KratosMultiphysics.Model, settings: typing.Union[KratosMultiphysics.Parameters, ParametersWrapper]) -> AggregateOperation:
    '''Return the operation factory specified by the setting 'operation_type'.

    If the 'operation_type' is not found and the settings have a 'module_name',
    the module is imported and used to create the operation. If 'module_name'
    is not found, an exception is raised. Empty settings will contain default
    values after returning from the function call.
    '''
    return AggregateOperation(model, settings)
