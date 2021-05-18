# Importing the Kratos Library
import KratosMultiphysics as KM

# Other imports
from abc import ABCMeta, abstractmethod

class PythonMapper(metaclass=ABCMeta):
    """Baseclass for python based mappers in Kratos
    The inteface matches the C++ version ("custom_mappers/mapper.h")
    The py-mappers are intentionally NOT derived from the c++ version.
    Reasons:
    - Doing so would require some special treatment of the pure virtual functions exposed to python
    - They are more or less temporary until Kratos has more Mappers
    """
    def __init__(self, model_part_origin, model_part_destination, mapper_settings):
        self.model_part_origin = model_part_origin
        self.model_part_destination = model_part_destination

        self.mapper_settings = mapper_settings
        self.mapper_settings.ValidateAndAssignDefaults(self._GetDefaultParameters())

        self.echo_level = self.mapper_settings["echo_level"].GetInt()

    # public methods, same as in "custom_mappers/mapper.h"
    def Map(self, variable_origin, variable_destination, mapper_flags=KM.Flags()):
        CheckVariables(variable_origin, variable_destination)
        self._MapInternal(variable_origin, variable_destination, mapper_flags)

    def InverseMap(self, variable_origin, variable_destination, mapper_flags=KM.Flags()):
        CheckVariables(variable_origin, variable_destination)
        self._InverseMapInternal(variable_origin, variable_destination, mapper_flags)

    @abstractmethod
    def UpdateInterface(self): pass

    # protected methods
    @abstractmethod
    def _MapInternal(self, variable_origin, variable_destination, mapper_flags): pass

    @abstractmethod
    def _InverseMapInternal(self, variable_origin, variable_destination, mapper_flags): pass

    @classmethod
    def _GetDefaultParameters(cls):
        return KM.Parameters("""{
            "mapper_type" : "",
            "echo_level"  : 0
        }""")

    @classmethod
    def _ClassName(cls):
        return cls.__name__

def CheckVariables(variable_origin, variable_destination):
    pass