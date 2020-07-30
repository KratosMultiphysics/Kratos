set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosMetisApplication ****")

include(pybind11Tools)

include_directories( ${KRATOS_SOURCE_DIR}/kratos )

# Use Metis 5 by default
if(NOT DEFINED USE_METIS_5)
	message("USE_METIS_5 not defined. Setting to ON")
	SET (USE_METIS_5 ON)
endif(NOT DEFINED USE_METIS_5)

if(${USE_METIS_5} MATCHES ON)
    message("using metis (and not Parmetis) version 5")
    add_definitions(-DKRATOS_USE_METIS_5)
    find_package(METIS REQUIRED)
    include_directories(${METIS_INCLUDE_DIR})
    link_directories(${METIS_LIBRARIES})
    install(FILES ${METIS_LIBRARIES} DESTINATION libs)
else(${USE_METIS_5} MATCHES ON)
    message("using PARMETIS < 4 as default. Set USE_METIS_5=ON to compile with METIS v5 instead- note that current compile option is deprecated and will be removed soon")
    find_package(PARMETIS REQUIRED)
    include_directories(${PARMETIS_INCLUDE_DIR})
    link_directories(${PARMETIS_LIBRARIES})
endif(${USE_METIS_5} MATCHES ON)

## generate variables with the sources
set( KRATOS_METIS_APPLICATION_SOURCES
	${CMAKE_CURRENT_SOURCE_DIR}/metis_application.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_processes_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/kratos_metis_python_application.cpp
)

# ###############################################################
pybind11_add_module(KratosMetisApplication MODULE THIN_LTO ${KRATOS_METIS_APPLICATION_SOURCES})
target_link_libraries(KratosMetisApplication PRIVATE KratosCore ${METIS_LIBRARY} ${MPI_LIBRARIES} )
set_target_properties(KratosMetisApplication PROPERTIES PREFIX "")

# changing the .dll suffix to .pyd
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
	set_target_properties(KratosMetisApplication PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

# changing the .dylib suffix to .so
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
	set_target_properties(KratosMetisApplication PROPERTIES SUFFIX .so)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

# Cotire
if(USE_COTIRE MATCHES ON)
    cotire(KratosMetisApplication)
endif(USE_COTIRE MATCHES ON)

# Add to the KratosMultiphisics Python module
kratos_python_install(${INSTALL_PYTHON_USING_LINKS} ${CMAKE_CURRENT_SOURCE_DIR}/MetisApplication.py KratosMultiphysics/MetisApplication/__init__.py )

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
install(TARGETS KratosMetisApplication DESTINATION libs )