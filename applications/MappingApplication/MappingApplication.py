# makes KratosMultiphysics backward compatible with python 2.6 and 2.7
from __future__ import print_function, absolute_import, division

# Application dependent names and paths
from KratosMultiphysics import _ImportApplication
from KratosMappingApplication import *
application = KratosMappingApplication()
application_name = "KratosMappingApplication"

_ImportApplication(application, application_name)

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