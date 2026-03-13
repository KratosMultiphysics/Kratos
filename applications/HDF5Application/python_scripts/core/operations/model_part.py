'''HDF5 model part operations.

license: HDF5Application/license.txt
'''


# --- Code Imports ---
import KratosMultiphysics

# --- HDF5 Imports ---
import KratosMultiphysics.HDF5Application as KratosHDF5
from ..utils import ParametersWrapper
from ..utils import EvaluatePattern
from ..file_io import OpenHDF5File

# --- STD Imports ---
from importlib import import_module
import sys
import typing
#import abc

class IOOperation(KratosMultiphysics.Operation): # IOOperation(KratosMultiphysics.Operation, metaclass = abc.ABCMeta) # <== conflicting metaclasses
    """ @brief Base class for HDF5 IO operations."""

    def __init__(self,
                 model_part: KratosMultiphysics.ModelPart,
                 parameters: KratosMultiphysics.Parameters,
                 file: KratosHDF5.HDF5File):
        super().__init__()
        self.__model_part = model_part
        self.__file = file
        self.__parameters = parameters

        self.__parameters.ValidateAndAssignDefaults(self.GetDefaultParameters())

        self.__time_format = ""
        if self.__parameters.Has("time_format"):
            self.__time_format = self.__parameters["time_format"].GetString()

        if not self.__parameters.Has("prefix"):
            raise RuntimeError(f"The parameters does not have prefix. Parameters: {self.__parameters}")

    #@abc.abstractmethod # <== missing the ABCMeta metaclass
    def Execute(self) -> None:
        raise NotImplementedError("Call to a pure abstract method")

    #@abc.abstractclassmethod # <== missing the ABCMeta metaclass
    @classmethod
    def GetDefaultParameters(cls) -> KratosMultiphysics.Parameters:
        raise NotImplementedError("Call to a pure abstract method")

    @property
    def model_part(self) -> KratosMultiphysics.ModelPart:
        return self.__model_part

    @property
    def parameters(self) -> KratosMultiphysics.Parameters:
        current_parameters = self.__parameters.Clone()
        current_parameters["prefix"].SetString(self.prefix)
        return current_parameters

    @property
    def file(self) -> KratosHDF5.HDF5File:
        return self.__file

    @property
    def prefix(self) -> str:
        prefix = self.__parameters["prefix"].GetString()
        evaluated_prefix = EvaluatePattern(prefix, self.model_part, self.time_format)
        return evaluated_prefix

    @property
    def time_format(self) -> str:
        return self.__time_format

class ModelPartIOOperation(IOOperation):
    """ @brief Base class for HDF5 IO operations on @ref ModelPart s."""

    @classmethod
    def GetDefaultParameters(cls) -> KratosMultiphysics.Parameters:
        return KratosMultiphysics.Parameters("""{
            "prefix"        : "/ModelData",
            "time_format"   : "0.4f",
            "operation_type": ""
        }""")


class ModelPartInput(ModelPartIOOperation):
    '''Reads a @ref ModelPart from a file.'''

    def Execute(self) -> None:
        KratosHDF5.HDF5ModelPartIO(self.parameters, self.file).ReadModelPart(self.model_part)


class ModelPartOutput(ModelPartIOOperation):
    '''Writes a @ref ModelPart to a file.'''

    def Execute(self) -> None:
        KratosHDF5.HDF5ModelPartIO(self.parameters, self.file).WriteModelPart(self.model_part)


class PartitionedModelPartOutput(ModelPartIOOperation):
    '''Writes a partitioned model part to a file.'''

    def Execute(self) -> None:
        KratosHDF5.HDF5PartitionedModelPartIO(self.parameters, self.file).WriteModelPart(self.model_part)


class ProcessInfoIOOperation(IOOperation):
    @classmethod
    def GetDefaultParameters(cls) -> KratosMultiphysics.Parameters:
        return KratosMultiphysics.Parameters("""{
            "prefix"     : "/ProcessInfo",
            "time_format": "0.4f"
        }""")

class ProcessInfoOutput(ProcessInfoIOOperation):
    '''Writes a @ref ProcessInfo to a file.'''

    def Execute(self) -> None:
        KratosHDF5.WriteDataValueContainer(self.file, self.prefix, self.model_part.ProcessInfo)


class ProcessInfoInput(ProcessInfoIOOperation):
    '''Reads a @ref ProcessInfo from a file.'''

    def Execute(self) -> None:
        KratosHDF5.ReadDataValueContainer(self.file, self.prefix, self.model_part.ProcessInfo)

class VariableOutputOperation(IOOperation):
    def _MakeIO(self) -> typing.Any:
        class_name = self.__class__.__name__
        if class_name.endswith("Output"):
            # This maps the output classes to their corresponding IO classes
            # Eg. ElementDataValueOutput is mapped to ElementDataValueIO
            class_type_name = f"HDF5{class_name[:-6]}IO"
            if hasattr(KratosHDF5, class_type_name):
                return getattr(KratosHDF5, class_type_name)(self.parameters, self.file)
            else:
                raise RuntimeError(f"The {class_type_name} not found in HDF5Application.")
        else:
            raise RuntimeError(f"Unsupported operation type {class_name}")

    def Execute(self) -> None:
        self._MakeIO().Write(self.model_part, self.parameters["custom_attributes"])

    '''Generates json settings for variable data output.'''
    @classmethod
    def GetDefaultParameters(cls) -> KratosMultiphysics.Parameters:
        return KratosMultiphysics.Parameters("""{
            "prefix"           : "/ResultsData",
            "list_of_variables": [],
            "custom_attributes": {},
            "time_format"      : "0.4f",
            "operation_type"   : ""
        }""")

class VariableInputOperation(IOOperation):
    def _MakeIO(self) -> typing.Any:
        class_name = self.__class__.__name__
        if class_name.endswith("Input"):
            class_type_name = f"HDF5{class_name[:-5]}IO"
            # This maps the input classes to their corresponding IO classes
            # Eg. ElementDataValueInput is mapped to ElementDataValueIO
            if hasattr(KratosHDF5, class_type_name):
                return getattr(KratosHDF5, class_type_name)(self.parameters, self.file)
            else:
                raise RuntimeError(f"The {class_type_name} not found in HDF5Application.")
        else:
            raise RuntimeError(f"Unsupported operation type {class_name}")

    def Execute(self) -> None:
        self._attributes = self._MakeIO().Read(self.model_part)

    '''Generates json settings for variable data output.'''
    @classmethod
    def GetDefaultParameters(cls) -> KratosMultiphysics.Parameters:
        return KratosMultiphysics.Parameters("""{
            "prefix"           : "/ResultsData",
            "list_of_variables": [],
            "time_format"      : "0.4f",
            "operation_type"   : ""
        }""")

    @property
    def attributes(self) -> 'dict[str, KratosMultiphysics.Parameters]':
        if hasattr("_attributes", self):
            return self._attributes
        else:
            raise RuntimeError("Please run Execute method before retireving attributes.")


class ElementDataValueOutput(VariableOutputOperation): pass
class ElementDataValueInput(VariableInputOperation): pass
class ElementFlagValueOutput(VariableOutputOperation): pass
class ElementFlagValueInput(VariableInputOperation): pass
class ElementGaussPointOutput(VariableOutputOperation): pass
class ConditionDataValueOutput(VariableOutputOperation): pass
class ConditionDataValueInput(VariableInputOperation): pass
class ConditionFlagValueOutput(VariableOutputOperation): pass
class ConditionFlagValueInput(VariableInputOperation): pass
class ConditionGaussPointOutput(VariableOutputOperation): pass
class NodalSolutionStepDataOutput(VariableOutputOperation): pass
class NodalSolutionStepDataInput(VariableInputOperation): pass
class NodalDataValueOutput(VariableOutputOperation): pass
class NodalDataValueInput(VariableInputOperation): pass
class NodalFlagValueOutput(VariableOutputOperation): pass
class NodalFlagValueInput(VariableInputOperation): pass

class PrimalBossakOutput(VariableOutputOperation):
    def __init__(self,
                 model_part: KratosMultiphysics.ModelPart,
                 parameters: KratosMultiphysics.Parameters,
                 file: KratosHDF5.HDF5File):
        super().__init__(model_part, parameters, file)
        self.__alpha_bossak = parameters["alpha_bossak"].GetDouble()

    def _MakeIO(self) -> typing.Any:
        io = KratosHDF5.HDF5NodalSolutionStepBossakIO(self.parameters, self.file)
        io.SetAlphaBossak(self.__alpha_bossak)
        return io

    @property
    def alpha_bossak(self) -> float:
        return self.__alpha_bossak

    @classmethod
    def GetDefaultParameters(cls: "typing.Type[PrimalBossakOutput]") -> KratosMultiphysics.Parameters:
        return KratosMultiphysics.Parameters("""{
            "prefix"           : "/ResultsData",
            "list_of_variables": [],
            "custom_attributes": {},
            "time_format"      : "0.4f",
            "alpha_bossak"     : -0.3,
            "operation_type"   : ""
        }""")


class PrimalBossakInput(VariableInputOperation):
    '''Reads nodal solution step data from a file.

    This is used by the transient adjoint solvers.
    '''
    def _MakeIO(self) -> typing.Any:
        return KratosHDF5.HDF5NodalSolutionStepBossakIO(self.parameters, self.file)


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
