# Application dependent names and paths
from KratosMultiphysics import _ImportApplication
from KratosMultiphysics import Mapper as _CoreMapper
from KratosMultiphysics import MapperFactory as _CoreMapperFactory
from KratosMappingApplication import *
application = KratosMappingApplication()
application_name = "KratosMappingApplication"

_ImportApplication(application, application_name)

# hacks for backward compatibility
from KratosMultiphysics.kratos_utilities import IssueDeprecationWarning

class _DeprecatedMapper:
    def __getattribute__(self, method_name):
        IssueDeprecationWarning("MappingApplication", 'The "Mapper" was moved to the Core! (used for "{}")'.format(method_name))
        return getattr(_CoreMapper, method_name)

class _DeprecatedMapperFactory:
    def __getattr__(self, method_name):
        IssueDeprecationWarning("MappingApplication", 'The "MapperFactory" was moved to the Core! (used for "{}")'.format(method_name))
        return getattr(_CoreMapperFactory, method_name)

    def CreateMPIMapper(self, *args):
        IssueDeprecationWarning("MappingApplication", 'CreateMPIMapper is deprecated, please use "MappingApplication.MPIExtension.MPIMapperFactory.CreateMapper" instead')
        from KratosMultiphysics.MappingApplication.MPIExtension import MPIMapperFactory
        return MPIMapperFactory.CreateMapper(*args)

    def HasMPIMapper(self, *args):
        IssueDeprecationWarning("MappingApplication", 'HasMPIMapper is deprecated, please use "MappingApplication.MPIExtension.MPIMapperFactory.HasMapper" instead')
        from KratosMultiphysics.MappingApplication.MPIExtension import MPIMapperFactory
        return MPIMapperFactory.HasMapper(*args)

    def GetRegisteredMPIMapperNames(self, *args):
        IssueDeprecationWarning("MappingApplication", 'GetRegisteredMPIMapperNames is deprecated, please use "MappingApplication.MPIExtension.MPIMapperFactory.GetRegisteredMapperNames" instead')
        from KratosMultiphysics.MappingApplication.MPIExtension import MPIMapperFactory
        return MPIMapperFactory.GetRegisteredMapperNames(*args)

Mapper = _DeprecatedMapper()
MapperFactory = _DeprecatedMapperFactory()

'''
TODO:
    - Test for Serialization
    - Test for local-search?
    - Cleanup how the MapperParams are used
    - Further cleanup Trilinos and try some things (read up on opt-stuff)
    - use std::unordered_set for row & column indices-vectors in trilinos => does the map need sorted indices?
    - Function-Documentation
    - Delete copy and assignment-constructors?
'''
