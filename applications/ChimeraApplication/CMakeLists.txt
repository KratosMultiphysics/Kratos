set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** Configuring KratosChimeraApplication ****")

if(${FLUID_DYNAMICS_APPLICATION} MATCHES OFF )
    MESSAGE( FATAL_ERROR "ChimeraApplication needs FluidDynamicsApplicaiton to work. Please set FLUID_DYNAMICS_APPLICATION=ON in the configuration file." )
endif(${FLUID_DYNAMICS_APPLICATION} MATCHES OFF)

include(pybind11Tools)

include_directories( ${CMAKE_SOURCE_DIR}/kratos )
#dependecies on other applications
include_directories( ${CMAKE_SOURCE_DIR}/applications/FluidDynamicsApplication )
include_directories( ${CMAKE_SOURCE_DIR}/applications/ChimeraApplication )


# generate variables with the sources
set( KRATOS_CHIMERA_APPLICATION_SOURCES
	${CMAKE_CURRENT_SOURCE_DIR}/chimera_application.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/chimera_application_variables.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/chimera_python_application.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_processes_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_strategies_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_utilities_to_python.cpp
)

set(KRATOS_CHIMERA_APPLICATION_EXTRA_LIBS "")
## If we have mpi support, use it
if(${MPI_NEEDED} MATCHES ON )
	set(KRATOS_CHIMERA_APPLICATION_EXTRA_LIBS ${KRATOS_CHIMERA_APPLICATION_EXTRA_LIBS} ${MPI_LIBRARIES})
endif(${MPI_NEEDED} MATCHES ON )


# define library Kratos which defines the basic python interface
pybind11_add_module(KratosChimeraApplication MODULE THIN_LTO ${KRATOS_CHIMERA_APPLICATION_SOURCES})
target_link_libraries(KratosChimeraApplication PRIVATE KratosCore KratosFluidDynamicsCore ${KRATOS_CHIMERA_APPLICATION_EXTRA_LIBS})
set_target_properties(KratosChimeraApplication PROPERTIES PREFIX "")
install(TARGETS KratosChimeraApplication DESTINATION libs )

# changing the .dll suffix to .pyd (Windows)
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
	set_target_properties(KratosChimeraApplication PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

# changing the .dylib suffix to .so (OS X)
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
	set_target_properties(KratosChimeraApplication PROPERTIES SUFFIX .so)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

# install the python files
if(${INSTALL_PYTHON_FILES} MATCHES ON)
  get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
  install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/python_scripts DESTINATION applications/${CURRENT_DIR_NAME}  FILES_MATCHING PATTERN "*.py"  PATTERN ".svn" EXCLUDE)
endif(${INSTALL_PYTHON_FILES} MATCHES ON)

# Add to the KratosMultiphisics Python module
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/ChimeraApplication.py" DESTINATION "KratosMultiphysics/ChimeraApplication" RENAME "__init__.py")
