import KratosMultiphysics as KM
from KratosMultiphysics.MappingApplication.python_mapper import PythonMapper

def CreateMapper(model_part_origin, model_part_destination, mapper_settings):
    return EmpireMapperWrapper(model_part_origin, model_part_destination, mapper_settings)

class EmpireMapperWrapper(PythonMapper):
    def __init__(self, model_part_origin, model_part_destination, mapper_settings):
        super(EmpireMapperWrapper, self).__init__(model_part_origin, model_part_destination, mapper_settings)

        default_parameters = KM.Parameters("""{
            "mapper_type" : "PLEASE_SPECIFY_MAPPER_TYPE"
        }""")

        self.mapper_settings.ValidateAndAssignDefaults(default_parameters)

        raise NotImplementedError("This mapper was not yet implemented!")
