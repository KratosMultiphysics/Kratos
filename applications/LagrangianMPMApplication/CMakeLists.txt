set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosLagrangianMPMApplication ****")

include_directories( ${CMAKE_SOURCE_DIR}/kratos )


# generate variables with the sources
set( KRATOS_LAGRANGIAN_MPM_APPLICATION_SOURCES
     ${CMAKE_CURRENT_SOURCE_DIR}/lagrangian_mpm_application.cpp
     ${CMAKE_CURRENT_SOURCE_DIR}/lagrangian_mpm_application_variables.cpp
     ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/lagrangian_mpm_python_application.cpp
     ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_utilities_to_python.cpp
     ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_processes_to_python.cpp
     ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/meshless_base_element.cpp
     ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/updated_lagrangian_mpm_element.cpp
     #${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/lagrange_multiplier_condition_y.cpp
     #${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/lagrange_multiplier_condition_x.cpp
     ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/lagrange_multiplier_condition.cpp
     ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/lagrange_multiplier_condition_0.cpp
     ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/pointforce2D.cpp
   )

# define library Kratos which defines the basic python interface
add_library(KratosLagrangianMPMApplication SHARED ${KRATOS_LAGRANGIAN_MPM_APPLICATION_SOURCES})
target_link_libraries(KratosLagrangianMPMApplication KratosCore )
set_target_properties(KratosLagrangianMPMApplication PROPERTIES PREFIX "")

# changing the .dll suffix to .pyd (Windows)
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    set_target_properties(KratosLagrangianMPMApplication PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

# changing the .dylib suffix to .so (OS X)
    if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
        set_target_properties(KratosLagrangianMPMApplication PROPERTIES SUFFIX .so)
        endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")



if(USE_COTIRE MATCHES ON)
    cotire(KratosLagrangianMPMApplication)
endif(USE_COTIRE MATCHES ON)

install(TARGETS KratosLagrangianMPMApplication DESTINATION libs )

# install the python files
if(${INSTALL_PYTHON_FILES} MATCHES ON)
    get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
    install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/python_scripts DESTINATION applications/${CURRENT_DIR_NAME}  FILES_MATCHING PATTERN "*.py"  PATTERN ".svn" EXCLUDE)
endif(${INSTALL_PYTHON_FILES} MATCHES ON)

# Add to the KratosMultiphisics Python module
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/LagrangianMPMApplication.py" DESTINATION KratosMultiphysics )
