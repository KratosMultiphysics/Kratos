# Application dependent names and paths
from KratosMultiphysics import _ImportApplication
from KratosMultiphysics import Mapper, MapperFactory # for backward compatibility
from KratosMappingApplication import *
application = KratosMappingApplication()
application_name = "KratosMappingApplication"

_ImportApplication(application, application_name)


# hack for backward compatibility
from KratosMultiphysics.kratos_utilities import IssueDeprecationWarning
def CreateMPIMapper(*args):
    IssueDeprecationWarning("MappingApplication", 'CreateMPIMapper is deprecated, please use "MappingApplication.MPIExtension.MPIMapperFactory.CreateMapper" instead')
    from KratosMultiphysics.MappingApplication.MPIExtension import MPIMapperFactory
    return MPIMapperFactory.CreateMapper(*args)

def HasMPIMapper(*args):
    IssueDeprecationWarning("MappingApplication", 'HasMPIMapper is deprecated, please use "MappingApplication.MPIExtension.MPIMapperFactory.HasMapper" instead')
    from KratosMultiphysics.MappingApplication.MPIExtension import MPIMapperFactory
    return MPIMapperFactory.HasMapper(*args)

def GetRegisteredMPIMapperNames(*args):
    IssueDeprecationWarning("MappingApplication", 'GetRegisteredMPIMapperNames is deprecated, please use "MappingApplication.MPIExtension.MPIMapperFactory.GetRegisteredMapperNames" instead')
    from KratosMultiphysics.MappingApplication.MPIExtension import MPIMapperFactory
    return MPIMapperFactory.GetRegisteredMapperNames(*args)

setattr(MapperFactory, 'CreateMPIMapper', CreateMPIMapper)
setattr(MapperFactory, 'HasMPIMapper', HasMPIMapper)
setattr(MapperFactory, 'GetRegisteredMPIMapperNames', GetRegisteredMPIMapperNames)

'''
TODO:
    - Test for Serialization
    - Test for local-search?
    - Cleanup how the MapperParams are used
    - Further cleanup Trilinos and try some things (read up on opt-stuff)
    - use std::unordered_set for row & column indices-vectors in trilinos => does the map need sorted indices?
    - For Trilinos: What happens if a rank does not have local nodes???
    - Function-Documentation
    - Delete copy and assignment-constructors?
    - testing => do some logical tests with USE_TRANSPOSE
'''
