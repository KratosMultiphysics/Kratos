set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosFreeSurfaceApplication ****")

include(pybind11Tools)
include_directories( ${KRATOS_SOURCE_DIR}/kratos )

# generate variables with the sources
set( KRATOS_FREE_SURFACE_APPLICATION_CORE_SOURCES
	${CMAKE_CURRENT_SOURCE_DIR}/free_surface_application.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/free_surface_application_variables.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/edgebased_gid_io.cpp
)

set( KRATOS_FREE_SURFACE_PYTHON_INTERFACE_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_utilities_to_python.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_io_to_python.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_edgebased_levelset_solver_to_python.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/free_surface_python_application.cpp
)

# FreeSurface core library (C++ parts)
add_library(KratosFreeSurfaceCore SHARED ${KRATOS_FREE_SURFACE_APPLICATION_CORE_SOURCES})
target_link_libraries(KratosFreeSurfaceCore PUBLIC KratosCore)
set_target_properties(KratosFreeSurfaceCore PROPERTIES COMPILE_DEFINITIONS "FREE_SURFACE_APPLICATION=EXPORT,API")

# FreeSurface python module
pybind11_add_module(KratosFreeSurfaceApplication MODULE THIN_LTO ${KRATOS_FREE_SURFACE_PYTHON_INTERFACE_SOURCES})
target_link_libraries(KratosFreeSurfaceApplication PRIVATE KratosFreeSurfaceCore )
set_target_properties(KratosFreeSurfaceApplication PROPERTIES PREFIX "")

# Set batch size in the unity build
IF(CMAKE_UNITY_BUILD MATCHES ON)
    set_target_properties(KratosFreeSurfaceCore PROPERTIES UNITY_BUILD_BATCH_SIZE ${KRATOS_UNITY_BUILD_BATCH_SIZE})
    set_target_properties(KratosFreeSurfaceApplication PROPERTIES UNITY_BUILD_BATCH_SIZE ${KRATOS_UNITY_BUILD_BATCH_SIZE})
ENDIF(CMAKE_UNITY_BUILD MATCHES ON)

# changing the .dll suffix to .pyd (Windows)
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
	set_target_properties(KratosFreeSurfaceApplication PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

# changing the .dylib suffix to .so (OS X)
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
	set_target_properties(KratosFreeSurfaceApplication PROPERTIES SUFFIX .so)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

# Add to the KratosMultiphisics Python module
kratos_python_install(${INSTALL_PYTHON_USING_LINKS} ${CMAKE_CURRENT_SOURCE_DIR}/FreeSurfaceApplication.py KratosMultiphysics/FreeSurfaceApplication/__init__.py )

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
install(TARGETS KratosFreeSurfaceCore DESTINATION libs )
install(TARGETS KratosFreeSurfaceApplication DESTINATION libs )

# Define custom targets
set(KRATOS_KERNEL "${KRATOS_KERNEL};KratosFreeSurfaceCore" PARENT_SCOPE)
set(KRATOS_PYTHON_INTERFACE "${KRATOS_PYTHON_INTERFACE};KratosFreeSurfaceApplication" PARENT_SCOPE)