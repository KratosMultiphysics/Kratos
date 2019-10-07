import KratosMultiphysics as KM
from KratosMultiphysics.MappingApplication.python_mapper import PythonMapper

class EmpireMapperWrapper(PythonMapper):
    def __init__(model_part_origin, model_part_destination, mapper_settings):
        super(EmpireMapperWrapper, self).__init__(model_part_origin, model_part_destination, mapper_settings)

        default_parameters = KM.Parameters("""{
            "mapper_type" : "PLEASE_SPECIFY_MAPPER_TYPE"
        }""")

        self.settings.ValidateAndAssignDefault(default_parameters)
