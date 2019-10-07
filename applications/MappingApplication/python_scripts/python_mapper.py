import KratosMultiphysics as KM

class PythonMapper(object):
    """Baseclass for python based mappers in Kratos
    The inteface matches the C++ version ("custom_mappers/mapper.h")
    """
    def __init__(model_part_origin, model_part_destination, mapper_settings):
        self.model_part_origin = model_part_origin
        self.model_part_destination = model_part_destination
        self.mapper_settings = mapper_settings # Note: no validation done here

    def Map(Variable_origin, variable_destination, mapper_flags=KM.Flags()):
        raise NotImplementedError('"Map" was not implemented for "{}"'.format(self.__class__.__name__))

    def InverseMap(Variable_origin, variable_destination, mapper_flags=KM.Flags()):
        raise NotImplementedError('"InverseMap" was not implemented for "{}"'.format(self.__class__.__name__))

    def UpdateInterface():
        raise NotImplementedError('"UpdateInterface" was not implemented for "{}"'.format(self.__class__.__name__))
