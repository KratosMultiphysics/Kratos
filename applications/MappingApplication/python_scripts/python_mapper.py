import KratosMultiphysics as KM

class PythonMapper(object):
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
        self.mapper_settings = mapper_settings # Note: no validation done here, should be done in derived class

    def Map(self, variable_origin, variable_destination, mapper_flags=KM.Flags()):
        raise NotImplementedError('"Map" was not implemented for "{}"'.format(self.__class__.__name__))

    def InverseMap(self, variable_origin, variable_destination, mapper_flags=KM.Flags()):
        raise NotImplementedError('"InverseMap" was not implemented for "{}"'.format(self.__class__.__name__))

    def UpdateInterface(self):
        raise NotImplementedError('"UpdateInterface" was not implemented for "{}"'.format(self.__class__.__name__))
