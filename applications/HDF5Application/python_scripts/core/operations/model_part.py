'''HDF5 model part operations.

license: HDF5Application/license.txt
'''


# --- Code Imports ---
import KratosMultiphysics

# --- HDF5 Imports ---
import KratosMultiphysics.HDF5Application as KratosHDF5
from ..pattern import EvaluatePattern

# --- STD Imports ---
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

        self.__parameters.AddMissingParameters(self.GetDefaultParameters())

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
            "prefix"           : "/ModelData",
            "time_format"      : "0.4f",
            "custom_attributes": {}
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
            "time_format"      : "0.4f",
            "custom_attributes": {}
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
            "time_format"      : "0.4f"
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
            "alpha_bossak"     : -0.3
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


class VertexCoordinateOutput(IOOperation):
    @classmethod
    def GetDefaultParameters(cls) -> KratosMultiphysics.Parameters:
        return KratosMultiphysics.Parameters("""{
            "prefix"           : "/ResultsData",
            "custom_attributes": {},
            "time_format"      : "0.4f"
        }""")

    def __init__(self,
                 model_part: KratosMultiphysics.ModelPart,
                 parameters: KratosMultiphysics.Parameters,
                 file: KratosHDF5.HDF5File,
                 vertices: KratosHDF5.VertexContainer):
        super().__init__(model_part, parameters, file)
        self.__vertices = vertices

    def Execute(self) -> None:
        KratosHDF5.VertexContainerCoordinateIO(self.parameters, self.file).Write(self.__vertices, self.parameters["custom_attributes"])


class VertexValueOutput(IOOperation):
    @classmethod
    def GetDefaultParameters(cls) -> KratosMultiphysics.Parameters:
        return KratosMultiphysics.Parameters("""{
            "prefix"           : "/VertexData",
            "list_of_variables": [],
            "custom_attributes": {},
            "time_format"      : "0.4f"
        }""")

    def _MakeIO(self) -> typing.Any:
        raise NotImplementedError("Please implement _MakeIO method in derrived class.")

    def __init__(self,
                 model_part: KratosMultiphysics.ModelPart,
                 parameters: KratosMultiphysics.Parameters,
                 file: KratosHDF5.HDF5File,
                 vertices: KratosHDF5.VertexContainer):
        super().__init__(model_part, parameters, file)
        self.__vertices = vertices

    def Execute(self) -> None:
        self._MakeIO().Write(self.__vertices, self.parameters["custom_attributes"])

class VertexHistoricalValueOutput(VertexValueOutput):
    def _MakeIO(self) -> typing.Any:
        return KratosHDF5.VertexContainerHistoricalVariableIO(self.parameters, self.file)

class VertexNonHistoricalValueOutput(VertexValueOutput):
    def _MakeIO(self) -> typing.Any:
        return KratosHDF5.VertexContainerNonHistoricalVariableIO(self.parameters, self.file)
