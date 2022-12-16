'''HDF5 model part operations.

license: HDF5Application/license.txt
'''


# --- Code Imports ---
import KratosMultiphysics

# --- HDF5 Imports ---
import KratosMultiphysics.HDF5Application as KratosHDF5
from ..utils import ParametersWrapper

# --- STD Imports ---
from importlib import import_module
import abc


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

    # @abc.abstractmethod # <== missing the ABCMeta metaclass
    def Execute(self) -> None:
        pass


class IOFactory(metaclass = abc.ABCMeta):
    """ @brief Factory for creating IO operations.
        @details This intermediate factory class takes care of parsing input parameters,
                 that has two purposes:
                 - prevent reparsing the same parameters on each call to an operation
                 - fail fast on invalid input parameters, well before the first operation\n
                   is constructed
    """

    def __init__(self, parameters: KratosMultiphysics.Parameters):
        parameters.AddMissingParameters(self.GetDefaultParameters())
        if isinstance(parameters, KratosMultiphysics.Parameters):
            self.__prefix_pattern = parameters["prefix"].GetString()
        elif isinstance(parameters, ParametersWrapper):
            self.__prefix_pattern = parameters["prefix"]
        else:
            raise TypeError(f"Expecting KratosMultiphysics::Parameters as input parameters, but got {type(parameters)}")

    @abc.abstractmethod
    def __call__(self, model_part: KratosMultiphysics.ModelPart, file: KratosHDF5.HDF5File) -> IOOperation:
        pass

    @abc.abstractstaticmethod
    def GetDefaultParameters() -> KratosMultiphysics.Parameters:
        pass

    @property
    def prefix_pattern(self) -> str:
        return self.__prefix_pattern


class ModelPartIOOperation(IOOperation):
    """ @brief Base class for HDF5 IO operations on @ref ModelPart s."""

    def __init__(self,
                 model_part: KratosMultiphysics.ModelPart,
                 file: KratosHDF5.HDF5File,
                 prefix_pattern: str,
                 time_format: str):
        super().__init__()
        self.__model_part = model_part
        self.__file = file
        self.__prefix = Prefix(prefix_pattern, model_part, time_format)

    @property
    def model_part(self) -> KratosMultiphysics.ModelPart:
        return self.__model_part

    @property
    def file(self) -> KratosHDF5.HDF5File:
        return self.__file

    @property
    def prefix(self) -> str:
        return self.__prefix


class ModelPartInput(ModelPartIOOperation):
    '''Reads a @ref ModelPart from a file.'''

    def Execute(self) -> None:
        KratosHDF5.HDF5ModelPartIO(self.file, self.prefix).ReadModelPart(self.model_part)

    class Factory(IOFactory):

        def __init__(self, parameters: KratosMultiphysics.Parameters):
            super().__init__(parameters)
            self.__time_format = parameters["time_format"]

        def __call__(self, model_part: KratosMultiphysics.ModelPart, file: KratosHDF5.HDF5File) -> "ModelPartInput":
            return ModelPartInput(model_part, file, self.prefix_pattern, self.__time_format)

        @staticmethod
        def GetDefaultParameters() -> KratosMultiphysics.Parameters:
            return KratosMultiphysics.Parameters("""{
                "prefix" : "/ModelData",
                "time_format" : "0.4f"
            }""")


class ModelPartOutput(ModelPartIOOperation):
    '''Writes a @ref ModelPart to a file.'''

    def Execute(self) -> None:
        KratosHDF5.HDF5ModelPartIO(self.file, self.prefix).WriteModelPart(self.model_part)

    class Factory(IOFactory):

        def __init__(self, parameters: KratosMultiphysics.Parameters):
            super().__init__(parameters)
            self.__time_format = parameters["time_format"]

        def __call__(self, model_part: KratosMultiphysics.ModelPart, file: KratosHDF5.HDF5File) -> "ModelPartOutput":
            return ModelPartOutput(model_part, file, self.prefix_pattern, self.__time_format)

        @staticmethod
        def GetDefaultParameters() -> KratosMultiphysics.Parameters:
            return KratosMultiphysics.Parameters("""{
                "prefix" : "/ModelData",
                "time_format" : "0.4f"
            }""")


class PartitionedModelPartOutput(ModelPartIOOperation):
    '''Writes a partitioned model part to a file.'''

    def Execute(self) -> None:
        KratosHDF5.HDF5PartitionedModelPartIO(self.file, self.prefix).WriteModelPart(self.model_part)

    class Factory(IOFactory):

        def __init__(self, parameters: KratosMultiphysics.Parameters):
            super().__init__(parameters)
            self.__time_format = parameters["time_format"]

        def __call__(self, model_part: KratosMultiphysics.ModelPart, file: KratosHDF5.HDF5File) -> "PartitionedModelPartOutput":
            return PartitionedModelPartOutput(model_part, file, self.prefix_pattern, self.__time_format)

        @staticmethod
        def GetDefaultParameters() -> KratosMultiphysics.Parameters:
            return KratosMultiphysics.Parameters("""{
                "prefix" : "/ModelData",
                "time_format" : "0.4f"
            }""")


class ProcessInfoOutput(ModelPartIOOperation):
    '''Writes a @ref ProcessInfo to a file.'''

    def Execute(self) -> None:
        KratosHDF5.WriteDataValueContainer(self.file, self.prefix, self.model_part.ProcessInfo)

    class Factory(IOFactory):

        def __init__(self, parameters: KratosMultiphysics.Parameters):
            super().__init__(parameters)
            self.__time_format = parameters["time_format"]

        def __call__(self, model_part: KratosMultiphysics.ModelPart, file: KratosHDF5.HDF5File) -> "ProcessInfoOutput":
            return ProcessInfoOutput(model_part, file, self.prefix_pattern, self.__time_format)

        @staticmethod
        def GetDefaultParameters() -> KratosMultiphysics.Parameters:
            return KratosMultiphysics.Parameters("""{
                "prefix" : "/ProcessInfo",
                "time_format" : "0.4f"
            }""")


class ProcessInfoInput(ModelPartIOOperation):
    '''Reads a @ref ProcessInfo from a file.'''

    def Execute(self) -> None:
        KratosHDF5.ReadDataValueContainer(self.file, self.prefix, self.model_part.ProcessInfo)

    class Factory(IOFactory):

        def __init__(self, parameters: KratosMultiphysics.Parameters):
            super().__init__(parameters)
            self.__time_format = parameters["time_format"]

        def __call__(self, model_part: KratosMultiphysics.ModelPart, file: KratosHDF5.HDF5File) -> "ProcessInfoInput":
            return ProcessInfoInput(model_part, file, self.prefix_pattern, self.__time_format)

        @staticmethod
        def GetDefaultParameters() -> KratosMultiphysics.Parameters:
            return KratosMultiphysics.Parameters("""{
                "prefix" : "/ProcessInfo",
                "time_format" : "0.4f"
            }""")


class VariableIOOperation(IOOperation):
    '''Generates json settings for variable data IO.'''

    def __init__(self,
                 model_part: KratosMultiphysics.ModelPart,
                 file: KratosHDF5.HDF5File,
                 io_parameters: ParametersWrapper,
                 time_format: str):
        super().__init__()
        self.__model_part = model_part
        self.__file = file
        self.__io_parameters: KratosMultiphysics.Parameters = io_parameters.Get()
        self.__io_parameters["prefix"].SetString(Prefix(
            self.__io_parameters["prefix"].GetString(),
            self.__model_part,
            time_format))

    @staticmethod
    def GetDefaultParameters() -> KratosMultiphysics.Parameters:
        return KratosMultiphysics.Parameters("""{
            "prefix" : "/ResultsData",
            "list_of_variables" : [],
            "time_format" : "0.4f"
        }""")

    @property
    def model_part(self) -> KratosMultiphysics.ModelPart:
        return self.__model_part

    @property
    def file(self) -> KratosHDF5.HDF5File:
        return self.__file

    @property
    def io_parameters(self) -> KratosMultiphysics.Parameters:
        return self.__io_parameters

    class Factory(IOFactory):

        def __init__(self, parameters: KratosMultiphysics.Parameters):
            super().__init__(parameters)
            parameters.AddMissingParameters(self.GetDefaultParameters())
            self.__io_parameters = ParametersWrapper()
            self.__io_parameters.Get().AddValue("prefix", parameters.Get()["prefix"])
            self.__io_parameters.Get().AddValue("list_of_variables", parameters.Get()["list_of_variables"])
            self.__time_format = parameters["time_format"]

        @property
        def io_parameters(self) -> KratosMultiphysics.Parameters:
            return self.__io_parameters

        @property
        def time_format(self) -> str:
            return self.__time_format

        @staticmethod
        def GetDefaultParameters() -> KratosMultiphysics.Parameters:
            return VariableIOOperation.GetDefaultParameters()


class ElementDataValueOutput(VariableIOOperation):
    '''Writes non-historical element data values to a file.'''

    def Execute(self) -> None:
        KratosHDF5.HDF5ElementDataValueIO(
            self.io_parameters,
            self.file).WriteElementResults(self.model_part.Elements)

    class Factory(VariableIOOperation.Factory):

        def __call__(self, model_part: KratosMultiphysics.ModelPart, file: KratosHDF5.HDF5File) -> VariableIOOperation:
            return ElementDataValueOutput(model_part, file, self.io_parameters, self.time_format)


class ElementDataValueInput(VariableIOOperation):
    '''Reads non-historical element data values from a file.'''

    def Execute(self) -> None:
        KratosHDF5.HDF5ElementDataValueIO(
            self.io_parameters,
            self.file).ReadElementResults(self.model_part.Elements,
                                          self.model_part.GetCommunicator())

    class Factory(VariableIOOperation.Factory):

        def __call__(self, model_part: KratosMultiphysics.ModelPart, file: KratosHDF5.HDF5File) -> VariableIOOperation:
            return ElementDataValueInput(model_part, file, self.io_parameters, self.time_format)

class ElementFlagValueOutput(VariableIOOperation):
    '''Writes non-historical element flag values to a file.'''

    def Execute(self) -> None:
        KratosHDF5.HDF5ElementFlagValueIO(
            self.io_parameters,
            self.file).WriteElementFlags(self.model_part.Elements)

    class Factory(VariableIOOperation.Factory):

        def __call__(self, model_part: KratosMultiphysics.ModelPart, file: KratosHDF5.HDF5File) -> VariableIOOperation:
            return ElementFlagValueOutput(model_part, file, self.io_parameters, self.time_format)


class ElementFlagValueInput(VariableIOOperation):
    '''Reads non-historical element flag values from a file.'''

    def Execute(self) -> None:
        KratosHDF5.HDF5ElementFlagValueIO(
            self.io_parameters,
            self.file).ReadElementFlags(self.model_part.Elements,
                                        self.model_part.GetCommunicator())

    class Factory(VariableIOOperation.Factory):

        def __call__(self, model_part: KratosMultiphysics.ModelPart, file: KratosHDF5.HDF5File) -> VariableIOOperation:
            return ElementFlagValueInput(model_part, file, self.io_parameters, self.time_format)

class ElementGaussPointOutput(VariableIOOperation):
    '''Write element integration point values to a file.'''

    def Execute(self) -> None:
        KratosHDF5.HDF5ElementGaussPointOutput(
            self.io_parameters,
            self.file).WriteElementGaussPointValues(self.model_part.Elements,
                                                    self.model_part.GetCommunicator().GetDataCommunicator(),
                                                    self.model_part.ProcessInfo)

    class Factory(VariableIOOperation.Factory):

        def __call__(self, model_part: KratosMultiphysics.ModelPart, file: KratosHDF5.HDF5File) -> VariableIOOperation:
            return ElementGaussPointOutput(model_part, file, self.io_parameters, self.time_format)

class ConditionDataValueOutput(VariableIOOperation):
    '''Writes non-historical element data values to a file.'''

    def Execute(self) -> None:
        KratosHDF5.HDF5ConditionDataValueIO(
            self.io_parameters,
            self.file).WriteConditionResults(self.model_part.Conditions)

    class Factory(VariableIOOperation.Factory):

        def __call__(self, model_part: KratosMultiphysics.ModelPart, file: KratosHDF5.HDF5File) -> VariableIOOperation:
            return ConditionDataValueOutput(model_part, file, self.io_parameters, self.time_format)


class ConditionDataValueInput(VariableIOOperation):
    '''Reads non-historical element data values from a file.'''

    def Execute(self) -> None:
        KratosHDF5.HDF5ConditionDataValueIO(
            self.io_parameters,
            self.file).ReadConditionResults(self.model_part.Conditions,
                                            self.model_part.GetCommunicator())

    class Factory(VariableIOOperation.Factory):

        def __call__(self, model_part: KratosMultiphysics.ModelPart, file: KratosHDF5.HDF5File) -> VariableIOOperation:
            return ConditionDataValueInput(model_part, file, self.io_parameters, self.time_format)

class ConditionFlagValueOutput(VariableIOOperation):
    '''Writes non-historical element flag values to a file.'''

    def Execute(self) -> None:
        KratosHDF5.HDF5ConditionFlagValueIO(
            self.io_parameters,
            self.file).WriteConditionFlags(self.model_part.Conditions)

    class Factory(VariableIOOperation.Factory):

        def __call__(self, model_part: KratosMultiphysics.ModelPart, file: KratosHDF5.HDF5File) -> VariableIOOperation:
            return ConditionFlagValueOutput(model_part, file, self.io_parameters, self.time_format)


class ConditionFlagValueInput(VariableIOOperation):
    '''Reads non-historical element flag values from a file.'''

    def Execute(self) -> None:
        KratosHDF5.HDF5ConditionFlagValueIO(
            self.io_parameters,
            self.file).ReadConditionFlags(self.model_part.Conditions,
                                          self.model_part.GetCommunicator())

    class Factory(VariableIOOperation.Factory):

        def __call__(self, model_part: KratosMultiphysics.ModelPart, file: KratosHDF5.HDF5File) -> VariableIOOperation:
            return ConditionFlagValueInput(model_part, file, self.io_parameters, self.time_format)

class ConditionGaussPointOutput(VariableIOOperation):
    '''Write condition integration point values to a file.'''

    def Execute(self) -> None:
        KratosHDF5.HDF5ConditionGaussPointOutput(
            self.io_parameters,
            self.file).WriteConditionGaussPointValues(self.model_part.Conditions,
                                                      self.model_part.GetCommunicator().GetDataCommunicator(),
                                                      self.model_part.ProcessInfo)

    class Factory(VariableIOOperation.Factory):

        def __call__(self, model_part: KratosMultiphysics.ModelPart, file: KratosHDF5.HDF5File) -> VariableIOOperation:
            return ConditionGaussPointOutput(model_part, file, self.io_parameters, self.time_format)


class NodalSolutionStepDataOutput(VariableIOOperation):
    '''Writes nodal solution step data to a file.'''

    def Execute(self) -> None:
        KratosHDF5.HDF5NodalSolutionStepDataIO(
            self.io_parameters,
            self.file).WriteNodalResults(self.model_part, 0)

    class Factory(VariableIOOperation.Factory):

        def __call__(self, model_part: KratosMultiphysics.ModelPart, file: KratosHDF5.HDF5File) -> VariableIOOperation:
            return NodalSolutionStepDataOutput(model_part, file, self.io_parameters, self.time_format)


class NodalSolutionStepDataInput(VariableIOOperation):
    '''Reads nodal solution step data from a file.'''

    def Execute(self) -> None:
        KratosHDF5.HDF5NodalSolutionStepDataIO(
            self.io_parameters,
            self.file).ReadNodalResults(self.model_part, 0)

    class Factory(VariableIOOperation.Factory):

        def __call__(self, model_part: KratosMultiphysics.ModelPart, file: KratosHDF5.HDF5File) -> VariableIOOperation:
            return NodalSolutionStepDataInput(model_part, file, self.io_parameters, self.time_format)


class NodalDataValueOutput(VariableIOOperation):
    '''Writes non-historical nodal data values to a file.'''

    def Execute(self) -> None:
        KratosHDF5.HDF5NodalDataValueIO(
            self.io_parameters,
            self.file).WriteNodalResults(self.model_part.Nodes)

    class Factory(VariableIOOperation.Factory):

        def __call__(self, model_part: KratosMultiphysics.ModelPart, file: KratosHDF5.HDF5File) -> VariableIOOperation:
            return NodalDataValueOutput(model_part, file, self.io_parameters, self.time_format)


class NodalDataValueInput(VariableIOOperation):
    '''Reads non-historical nodal data values from a file.'''

    def Execute(self) -> None:
        KratosHDF5.HDF5NodalDataValueIO(
            self.io_parameters,
            self.file).ReadNodalResults(self.model_part.Nodes,
                                        self.model_part.GetCommunicator())

    class Factory(VariableIOOperation.Factory):

        def __call__(self, model_part: KratosMultiphysics.ModelPart, file: KratosHDF5.HDF5File) -> VariableIOOperation:
            return NodalDataValueInput(model_part, file, self.io_parameters, self.time_format)

class NodalFlagValueOutput(VariableIOOperation):
    '''Writes non-historical nodal flag values to a file.'''

    def Execute(self) -> None:
        KratosHDF5.HDF5NodalFlagValueIO(
            self.io_parameters,
            self.file).WriteNodalFlags(self.model_part.Nodes)

    class Factory(VariableIOOperation.Factory):

        def __call__(self, model_part: KratosMultiphysics.ModelPart, file: KratosHDF5.HDF5File) -> VariableIOOperation:
            return NodalFlagValueOutput(model_part, file, self.io_parameters, self.time_format)


class NodalFlagValueInput(VariableIOOperation):
    '''Reads non-historical nodal flag values from a file.'''

    def Execute(self) -> None:
        KratosHDF5.HDF5NodalFlagValueIO(
            self.io_parameters,
            self.file).ReadNodalFlags(self.model_part.Nodes,
                                      self.model_part.GetCommunicator())

    class Factory(VariableIOOperation.Factory):

        def __call__(self, model_part: KratosMultiphysics.ModelPart, file: KratosHDF5.HDF5File) -> VariableIOOperation:
            return NodalFlagValueInput(model_part, file, self.io_parameters, self.time_format)


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
                 file: KratosHDF5.HDF5File,
                 parameters: KratosMultiphysics.Parameters,
                 io_parameters: KratosMultiphysics.Parameters,
                 time_format: str):
        super().__init__(model_part, file, io_parameters, time_format)
        self.__parameters = parameters

    def Execute(self) -> None:
        primal_io = KratosHDF5.HDF5NodalSolutionStepBossakIO(self.io_parameters, self.file)
        primal_io.SetAlphaBossak(self.__parameters["alpha_bossak"].GetDouble())
        primal_io.WriteNodalResults(self.model_part)

    @staticmethod
    def GetDefaultParameters() -> KratosMultiphysics.Parameters:
        return KratosMultiphysics.Parameters("""{
            "prefix" : "/ResultsData",
            "list_of_variables" : [],
            "time_format" : "0.4f",
            "alpha_bossak" : -0.3
        }""")

    class Factory(IOFactory):

        def __init__(self, parameters: KratosMultiphysics.Parameters):
            parameters.AddMissingParameters(self.GetDefaultParameters())
            super().__init__(parameters)
            self.__parameters = parameters.Get()
            self.__io_parameters = ParametersWrapper()
            self.__io_parameters.Get().AddValue("prefix", self.__parameters["prefix"].Clone())
            self.__io_parameters.Get().AddValue("list_of_variables", self.__parameters["list_of_variables"].Clone())

            if self.__parameters.Has("time_format"):
                self.__time_format = self.__parameters["time_format"].GetString()
            else:
                self.__time_format = PrimalBossakOutput.GetDefaultParameters()["time_format"].GetString()

        def __call__(self, model_part: KratosMultiphysics.ModelPart, file: KratosHDF5.HDF5File) -> VariableIOOperation:
            return PrimalBossakOutput(model_part, file, self.parameters, self.__io_parameters, self.__time_format)

        @staticmethod
        def GetDefaultParameters() -> KratosMultiphysics.Parameters:
            return PrimalBossakOutput.GetDefaultParameters()

        @property
        def parameters(self) -> KratosMultiphysics.Parameters:
            return self.__parameters


class PrimalBossakInput(VariableIOOperation):
    '''Reads nodal solution step data from a file.

    This is used by the transient adjoint solvers.
    '''

    def Execute(self) -> None:
        KratosHDF5.HDF5NodalSolutionStepBossakIO(
            self.io_parameters,
            self.file).ReadNodalResults(self.model_part)

    class Factory(VariableIOOperation.Factory):

        def __call__(self, model_part: KratosMultiphysics.ModelPart, file: KratosHDF5.HDF5File) -> VariableIOOperation:
            return PrimalBossakInput(model_part, file, self.io_parameters, self.time_format)


class MoveMesh(ModelPartIOOperation):
    '''Perform a mesh move operation on a model part.

    The primary use case is to set the mesh to the current configuration after
    reading the model part.
    '''

    def Execute(self) -> None:
        KratosMultiphysics.ImplicitSolvingStrategy(self.model_part, True).MoveMesh()

    class Factory(IOFactory):

        def __call__(self, model_part: KratosMultiphysics.ModelPart, file: KratosHDF5.HDF5File) -> "ProcessInfoInput":
            return MoveMesh(model_part, file, self.prefix_pattern, "")


def GetSubclasses(base_class: type) -> "list[type]":
    """Recursively find all subclasses of a base class"""
    subclasses = base_class.__subclasses__()
    for subclass in base_class.__subclasses__():
        subclasses += GetSubclasses(subclass)
    return subclasses


def Create(settings: ParametersWrapper) -> IOFactory:
    '''Return the operation factory specified by the setting 'operation_type'.

    If the 'operation_type' is not found and the settings have a 'module_name',
    the module is imported and used to create the operation. If 'module_name'
    is not found, an exception is raised. Empty settings will contain default
    values after returning from the function call.
    '''
    settings.SetDefault('operation_type', 'model_part_output')
    operation_type = settings['operation_type']

    # Find operation in the local definitions
    snake_to_camel = lambda string: "".join(part.title() for part in string.split('_'))
    operation_type_camel = snake_to_camel(operation_type)
    factory = next((io.Factory for io in GetSubclasses(IOOperation) if io.__name__ == operation_type_camel), None)

    if factory is None: # the requested operation was not defined in this script
        if settings.Has('module_name'):
            module_name = settings['module_name']
            module = import_module(
                'KratosMultiphysics.HDF5Application.core.' + module_name)
            instance =  module.Create(settings)
        else:
            raise ValueError(
                '"operation_type" has invalid value "' + operation_type + '"')
    else:
        instance = factory(settings)

    if instance is None:
        raise RuntimeError(f"no instance created for operation: {operation_type_camel}")
    return instance