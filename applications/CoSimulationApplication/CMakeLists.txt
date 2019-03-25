set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosCoSimulationApplication ****")

if(${CO_SIMULATION_APPLICATION} MATCHES ON)

    include_directories( ${CMAKE_SOURCE_DIR}/kratos )

    # generate variables with the sources
    set( KRATOS_CO_SIMULATION_APPLICATION_SOURCES
        ${CMAKE_CURRENT_SOURCE_DIR}/co_simulation_application.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/co_simulation_python_application.cpp
    )

    # define library Kratos which defines the basic python interface
    pybind11_add_module(KratosCoSimulationApplication MODULE THIN_LTO ${KRATOS_CO_SIMULATION_APPLICATION_SOURCES})
    target_link_libraries(KratosCoSimulationApplication PRIVATE KratosCore )
    set_target_properties(KratosCoSimulationApplication PROPERTIES PREFIX "")
    install(TARGETS KratosCoSimulationApplication DESTINATION libs )

    # changing the .dll suffix to .pyd (Windows)
    if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
        set_target_properties(KratosCoSimulationApplication PROPERTIES SUFFIX .pyd)
    endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

    # changing the .dylib suffix to .so (OS X)
    if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
        set_target_properties(KratosCoSimulationApplication PROPERTIES SUFFIX .so)
    endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

    # install the python files
    if(${INSTALL_PYTHON_FILES} MATCHES ON)
        get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
        install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/python_scripts DESTINATION applications/${CURRENT_DIR_NAME}  FILES_MATCHING PATTERN "*.py"  PATTERN ".svn" EXCLUDE)
    endif(${INSTALL_PYTHON_FILES} MATCHES ON)

    # Add to the KratosMultiphisics Python module
    install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/CoSimulationApplication.py" DESTINATION "KratosMultiphysics/CoSimulationApplication" RENAME "__init__.py")

elseif(NOT ${CO_SIMULATION_APPLICATION_PYTHON} MATCHES ON)

    project (KratosMultiphysics)
    cmake_minimum_required (VERSION 2.8.6)

    MESSAGE("Configuring the python only CoSimulationApplication ... ")
    get_filename_component(PARENT_DIR ${CMAKE_CURRENT_SOURCE_DIR} DIRECTORY)
    get_filename_component(KRATOS_DIR ${PARENT_DIR} DIRECTORY)

    ## Write KratosMultiphisics __init__.py 
    #take care: do not indent the lines below
    file(WRITE ${KRATOS_DIR}/KratosMultiphysics/__init__.py "import sys
import os.path
kratos_scripts=os.path.abspath(os.path.join(os.path.dirname(__file__),'../kratos/python_scripts'))
kratos_tests=os.path.abspath(os.path.join(os.path.dirname(__file__),'../kratos/tests'))
sys.path.append(kratos_scripts)
sys.path.append(kratos_tests)

def _ImportApplicationAsModule(application, application_name, application_folder, mod_path):
    print('Importing python   ' + application_name)

    # adding the scripts in APP_NAME/python_scripts such that they are treated as a regular python-module
    application_path = os.path.join(\"${KRATOS_DIR}/applications\", application_folder)
    python_path = os.path.join(application_path, 'python_scripts')
    mod_path.append(python_path)

")

    ## Write CoSimulationApplication __init__.py
    file(WRITE ${KRATOS_DIR}/KratosMultiphysics/CoSimulationApplication/__init__.py "# makes KratosMultiphysics backward compatible with python 2.6 and 2.7
from __future__ import print_function, absolute_import, division
import os
import sys

# Application dependent names and paths
import KratosMultiphysics as KM
application = 'KratosCoSimulationApplication'
application_name = 'KratosCoSimulationApplication'
application_folder = 'CoSimulationApplication'

KM._ImportApplicationAsModule(application, application_name, application_folder, __path__)
## adding pyKratos to path
application_path = os.path.join(\"${KRATOS_DIR}/applications\", application_folder)
python_path = os.path.join(application_path, 'custom_data_structure')
sys.path.append(python_path)
")

endif(NOT ${CO_SIMULATION_APPLICATION_PYTHON} MATCHES ON)
