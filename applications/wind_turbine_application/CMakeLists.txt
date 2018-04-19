set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosWindTurbineApplication ****")

add_definitions( -DTRILIBRARY )


include_directories( ${CMAKE_SOURCE_DIR}/kratos )
include_directories( ${CMAKE_SOURCE_DIR}/external_libraries/triangle )

if(${USE_TETGEN_NONFREE_TPL} MATCHES ON )
    add_definitions( -DTETGEN143 )
    add_definitions( -DTETLIBRARY )
    add_definitions( -DUSE_TETGEN_NONFREE_TPL )
    include_directories(${TETGEN_INCLUDES})
else(${USE_TETGEN_NONFREE_TPL} MATCHES ON )
    message("*********************************************************************************")
    message("********                  NEEDS TETGEN TO COMPILE                 ***************")
    message("*********************************************************************************")
endif(${USE_TETGEN_NONFREE_TPL} MATCHES ON )

## generate variables with the sources
set( KRATOS_WIND_TURBINE_APPLICATION_SOURCES     
	${CMAKE_CURRENT_SOURCE_DIR}//wind_turbine_application.cpp 
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/wind_turbine_python_application.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_utilities_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/wind_turbine_rotation_utilities.cpp

)



###############################################################
## define library Kratos which defines the basic python interface
add_library(KratosWindTurbineApplication SHARED ${KRATOS_WIND_TURBINE_APPLICATION_SOURCES})
if(${WIND_TURBINE_USE_PARALLEL_EXTENSION} MATCHES ON)
    add_definitions(-DWIND_TURBINE_USE_PARALLEL_EXTENSION)
    target_link_libraries(KratosWindTurbineApplication KratosCore triangle ${TETGEN_LIBRARIES} ${MPI_LIBRARIES})
else(${WIND_TURBINE_USE_PARALLEL_EXTENSION} MATCHES ON)
    target_link_libraries(KratosWindTurbineApplication KratosCore triangle ${TETGEN_LIBRARIES})
endif(${WIND_TURBINE_USE_PARALLEL_EXTENSION} MATCHES ON)
set_target_properties(KratosWindTurbineApplication PROPERTIES PREFIX "")
install(TARGETS KratosWindTurbineApplication DESTINATION libs )

# changing the .dll suffix to .pyd 
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
	set_target_properties(KratosWindTurbineApplication PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

# changing the .dylib suffix to .so 
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
	set_target_properties(KratosWindTurbineApplication PROPERTIES SUFFIX .so)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

if(${INSTALL_PYTHON_FILES} MATCHES ON)
  get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
  install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/python_scripts DESTINATION applications/${CURRENT_DIR_NAME}  FILES_MATCHING PATTERN "*.py"  PATTERN ".svn" EXCLUDE)
endif(${INSTALL_PYTHON_FILES} MATCHES ON)


# Add to the KratosMultiphisics Python module
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/WindTurbineApplication.py" DESTINATION KratosMultiphysics )
