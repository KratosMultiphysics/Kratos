set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosFSIApplication ****")

include_directories( ${KRATOS_SOURCE_DIR}/kratos )

## define application core sources (pure C++ parts)
set( KRATOS_FSI_APPLICATION_CORE_SOURCES
	# application class
	${CMAKE_CURRENT_SOURCE_DIR}/fsi_application.cpp
)

## define python interface sources
set( KRATOS_FSI_APPLICATION_PYTHON_INTERFACE_SOURCES
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/fsi_python_application.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_mappers_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_utilities_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_convergence_accelerators_to_python.cpp
)

## define fsi test sources
if(${KRATOS_BUILD_TESTING} MATCHES ON)
 	file(GLOB_RECURSE KRATOS_FSI_TESTING_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/tests/cpp_tests/*.cpp)
    kratos_add_gtests(TARGET KratosFSICore SOURCES "${KRATOS_FSI_TESTING_SOURCES}" USE_CUSTOM_MAIN)
endif(${KRATOS_BUILD_TESTING} MATCHES ON)

###############################################################################
## FSIApplication core library (C++ parts)
add_library( KratosFSICore SHARED ${KRATOS_FSI_APPLICATION_CORE_SOURCES})
target_link_libraries( KratosFSICore PUBLIC KratosCore)
set_target_properties( KratosFSICore PROPERTIES COMPILE_DEFINITIONS "FSI_APPLICATION=EXPORT,API")

## FSIApplication python module
pybind11_add_module( KratosFSIApplication MODULE THIN_LTO ${KRATOS_FSI_APPLICATION_PYTHON_INTERFACE_SOURCES} )
target_link_libraries( KratosFSIApplication PRIVATE KratosFSICore)
set_target_properties( KratosFSIApplication PROPERTIES PREFIX "")

# Set batch size in the unity build
IF(CMAKE_UNITY_BUILD MATCHES ON)
    set_target_properties(KratosFSICore PROPERTIES UNITY_BUILD_BATCH_SIZE ${KRATOS_UNITY_BUILD_BATCH_SIZE})
    set_target_properties(KratosFSIApplication PROPERTIES UNITY_BUILD_BATCH_SIZE ${KRATOS_UNITY_BUILD_BATCH_SIZE})
ENDIF(CMAKE_UNITY_BUILD MATCHES ON)

###############################################################################
# changing the .dll suffix to .pyd
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
	set_target_properties(KratosFSIApplication PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

# changing the .dylib suffix to .so
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
	set_target_properties(KratosFSIApplication PROPERTIES SUFFIX .so)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

# Add to the KratosMultiphisics Python module
kratos_python_install(${INSTALL_PYTHON_USING_LINKS} ${CMAKE_CURRENT_SOURCE_DIR}/FSIApplication.py KratosMultiphysics/FSIApplication/__init__.py )

# Install python files
get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
kratos_python_install_directory(${INSTALL_PYTHON_USING_LINKS} ${CMAKE_CURRENT_SOURCE_DIR}/python_scripts KratosMultiphysics/${CURRENT_DIR_NAME} )

# Kratos Testing. Install everything except sources to ensure that reference and configuration files are copied.
if(${INSTALL_TESTING_FILES} MATCHES ON )
    get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
    install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/tests DESTINATION applications/${CURRENT_DIR_NAME}
        PATTERN "*.git" EXCLUDE
        PATTERN "*.c" EXCLUDE
        PATTERN "*.h" EXCLUDE
        PATTERN "*.cpp" EXCLUDE
        PATTERN "*.hpp" EXCLUDE
  )
endif(${INSTALL_TESTING_FILES} MATCHES ON)

# Install targets
install(TARGETS KratosFSICore DESTINATION libs )
install(TARGETS KratosFSIApplication DESTINATION libs )

# Define custom targets
set(KRATOS_KERNEL "${KRATOS_KERNEL};KratosFSICore" PARENT_SCOPE)
set(KRATOS_PYTHON_INTERFACE "${KRATOS_PYTHON_INTERFACE};KratosFSIApplication" PARENT_SCOPE)