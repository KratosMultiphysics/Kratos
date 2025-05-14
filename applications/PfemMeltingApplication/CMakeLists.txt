set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosPfemMeltingApplication ****")

kratos_add_dependency(${KRATOS_SOURCE_DIR}/applications/FluidDynamicsApplication)
kratos_add_dependency(${KRATOS_SOURCE_DIR}/applications/ConvectionDiffusionApplication)
kratos_add_dependency(${KRATOS_SOURCE_DIR}/applications/MeshingApplication)

include_directories( ${KRATOS_SOURCE_DIR}/kratos )
include_directories( ${CMAKE_SOURCE_DIR}/applications/FluidDynamicsApplication )

## generate variables with the sources
set( KRATOS_PFEM_MELTING_APPLICATION_CORE_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/pfem_melting_application.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/pfem_melting_application_variables.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/hypo.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/qfluid.cpp
    #${CMAKE_SOURCE_DIR}/applications/FluidDynamicsApplication/custom_elements/hypo.cpp
    #${CMAKE_SOURCE_DIR}/applications/FluidDynamicsApplication/fluid_dynamics_application_variables.cpp  #TODO: this should REALLY
    #${CMAKE_SOURCE_DIR}/applications/FluidDynamicsApplication/custom_elements/vms.cpp  #TODO: this should REALLY
)



## define python interface sources
set( KRATOS_PFEM_MELTING_PYTHON_INTERFACE_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/kratos_pfem_melting_python_application.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_utilities_to_python.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_processes_to_python.cpp
##    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_strategies_to_python.cpp
)

#    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_strategies_to_python.cpp
#    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_response_functions_to_python.cpp

## define test sources
#if(${KRATOS_BUILD_TESTING} MATCHES ON)
#     file(GLOB_RECURSE KRATOS_PFEM_MELTING_TESTING_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/tests/cpp_tests/*.cpp)
#endif(${KRATOS_BUILD_TESTING} MATCHES ON)

###############################################################
## ConvectionDiffusionApplication core library (C++ parts)
add_library( KratosPfemMeltingCore SHARED ${KRATOS_PFEM_MELTING_APPLICATION_CORE_SOURCES} ${KRATOS_PFEM_MELTING_TESTING_SOURCES} )
target_link_libraries( KratosPfemMeltingCore PUBLIC KratosFluidDynamicsCore KratosConvectionDiffusionCore KratosCore )
set_target_properties( KratosPfemMeltingCore PROPERTIES COMPILE_DEFINITIONS "PFEM_MELTING_APPLICATION=EXPORT,API")

## FluidDynamicsApplication python module
pybind11_add_module( KratosPfemMeltingApplication MODULE THIN_LTO ${KRATOS_PFEM_MELTING_PYTHON_INTERFACE_SOURCES} )
target_link_libraries( KratosPfemMeltingApplication PRIVATE KratosPfemMeltingCore)
set_target_properties( KratosPfemMeltingApplication PROPERTIES PREFIX "")

# changing the .dll suffix to .pyd
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    set_target_properties(KratosPfemMeltingApplication PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

# changing the .dylib suffix to .so
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    set_target_properties(KratosPfemMeltingApplication PROPERTIES SUFFIX .so)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

# Cotire
if(USE_COTIRE MATCHES ON)
    cotire(KratosPfemMeltingCore)
    cotire(KratosPfemMeltingApplication)
endif(USE_COTIRE MATCHES ON)

# Add to the KratosMultiphisics Python module
kratos_python_install(${INSTALL_PYTHON_USING_LINKS} ${CMAKE_CURRENT_SOURCE_DIR}/PfemMeltingApplication.py KratosMultiphysics/PfemMeltingApplication/__init__.py )

# Install python files
get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
kratos_python_install_directory(${INSTALL_PYTHON_USING_LINKS} ${CMAKE_CURRENT_SOURCE_DIR}/python_scripts KratosMultiphysics/${CURRENT_DIR_NAME} )

# Kratos Testing. Install everything except sources to ensure that reference and configuration files are copied.
#if(${INSTALL_TESTING_FILES} MATCHES ON )
#    get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
#    install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/tests DESTINATION applications/${CURRENT_DIR_NAME}
#        PATTERN "*.git" EXCLUDE
#        PATTERN "*.c" EXCLUDE
#        PATTERN "*.h" EXCLUDE
#        PATTERN "*.cpp" EXCLUDE
#        PATTERN "*.hpp" EXCLUDE
#  )
#endif(${INSTALL_TESTING_FILES} MATCHES ON)

# Install Targets
install(TARGETS KratosPfemMeltingCore DESTINATION libs )
install(TARGETS KratosPfemMeltingApplication DESTINATION libs )
