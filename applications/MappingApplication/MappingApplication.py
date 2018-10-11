# makes KratosMultiphysics backward compatible with python 2.6 and 2.7
from __future__ import print_function, absolute_import, division

# Application dependent names and paths
from KratosMappingApplication import *
application = KratosMappingApplication()
application_name = "KratosMappingApplication"
application_folder = "MappingApplication"

# The following lines are common for all applications
from . import application_importer
import inspect
caller = inspect.stack()[1] # Information about the file that imported this, to check for unexpected imports
application_importer.ImportApplication(application,application_name,application_folder,caller)

'''
TODO:
    - Nearest Element with Volume
    - Nearest Element tests
    - Test for Serialization
    - Test for local-search?
    - Serial: Construct Matrix Structure
    - Finish new tests => requires first implementation in core
    - Cleanup & remove old tests
    - Cleanup how the MapperParams are used
    - Cleanup how the AssemblyUtility is used (in Mapper-BaseClass)
    - Further cleanup Trilinos and try some things (read up on opt-stuff)
    - For Trilinos: What happens if a rank does not have local nodes???
    - Function-Documentation
    - Delete copy and assignment-constructors
    - use explicit
    - testing => do some logical tests with CONSERVATIVE
    - MapperFlags: Check that they are used correctly and all of them are used in tests (CONSERVATIVE & USE_TRANSPOSE)
'''