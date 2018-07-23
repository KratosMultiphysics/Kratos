set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosStabilizedCFDApplication ****")

include_directories( ${CMAKE_SOURCE_DIR}/kratos )
include_directories( ${CMAKE_SOURCE_DIR}/external_libraries )

# generate variables with the sources
set( KRATOS_STABILIZED_CFD_APPLICATION_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/stabilized_cfd_application.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/stabilized_cfd_application_variables.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/stabilized_cfd_python_application.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_utilities_to_python.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/turbulence_statistics_container.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/dss.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/dss_fic.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/dss_fic_limited.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/dss_gls.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/dynss.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/dss_notau2.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/dynss_notau2.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/dss_ps.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/dss_face.cpp
)
## If we have mpi support, use it
set(KRATOS_STABILIZED_CDF_APPLICATION_EXTRA_LIBS "")
if(${MPI_NEEDED} MATCHES ON )
	set(KRATOS_STABILIZED_CDF_APPLICATION_EXTRA_LIBS ${MPI_LIBRARIES})
endif(${MPI_NEEDED} MATCHES ON )

## Define PyBind module for the application
pybind11_add_module(KratosStabilizedCFDApplication MODULE THIN_LTO ${KRATOS_STABILIZED_CFD_APPLICATION_SOURCES})

target_link_libraries(KratosStabilizedCFDApplication PRIVATE KratosCore KratosFluidDynamicsCore ${KRATOS_STABILIZED_CDF_APPLICATION_EXTRA_LIBS})
set_target_properties(KratosStabilizedCFDApplication PROPERTIES PREFIX "")

# changing the .dll suffix to .pyd (Windows)
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
	set_target_properties(KratosStabilizedCFDApplication PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

# changing the .dylib suffix to .so (OS X)
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
	set_target_properties(KratosStabilizedCFDApplication PROPERTIES SUFFIX .so)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

install(TARGETS KratosStabilizedCFDApplication DESTINATION libs )

# install the python files
if(${INSTALL_PYTHON_FILES} MATCHES ON)
  get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
  install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/python_scripts DESTINATION applications/${CURRENT_DIR_NAME}  FILES_MATCHING PATTERN "*.py"  PATTERN ".svn" EXCLUDE)
endif(${INSTALL_PYTHON_FILES} MATCHES ON)

# Add to the KratosMultiphisics Python module
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/StabilizedCFDApplication.py" DESTINATION KratosMultiphysics )
