# makes KratosMultiphysics backward compatible with python 2.6 and 2.7
from __future__ import print_function, absolute_import, division

# Application dependent names and paths
import KratosMultiphysics as KM
from KratosMappingApplication import *
application = KratosMappingApplication()
application_name = "KratosMappingApplication"
application_folder = "MappingApplication"

KM._ImportApplicationAsModule(application, application_name, application_folder, __path__)

'''
TODO:
    - Test for Serialization
    - Test for local-search?
    - Finish new tests => requires first implementation in core
    - Cleanup & remove old tests
    - Cleanup how the MapperParams are used
    - Further cleanup Trilinos and try some things (read up on opt-stuff)
    - use std::unordered_set for row & column indices-vectors in trilinos => does the map need sorted indices?
    - For Trilinos: What happens if a rank does not have local nodes???
    - Function-Documentation
    - Delete copy and assignment-constructors?
    - testing => do some logical tests with USE_TRANSPOSE
    - MapperFlags: Check that they are used correctly and all of them are used in tests (USE_TRANSPOSE)
'''