set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosConvectionDiffusionApplication ****")

include_directories( ${CMAKE_SOURCE_DIR}/kratos )

## generate variables with the sources
set( KRATOS_CONVECTION_DIFFUSION_APPLICATION_CORE_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/convection_diffusion_application.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/convection_diffusion_application_variables.cpp

    ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/eulerian_conv_diff.cpp

    ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/conv_diff_2d.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/conv_diff_3d.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/laplacian_element.cpp

    ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/thermal_face2D.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/thermal_face3D.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/flux_condition.cpp
)

## define python interface sources
set( KRATOS_CONVECTION_DIFFUSION_PYTHON_INTERFACE_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/kratos_convection_diffusion_python_application.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_strategies_to_python.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_utilities_to_python.cpp
)

## define fluid dynamics test sources
if(${KRATOS_BUILD_TESTING} MATCHES ON)
     file(GLOB_RECURSE KRATOS_CONVECTION_DIFFUSION_TESTING_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/tests/cpp_tests/*.cpp)
endif(${KRATOS_BUILD_TESTING} MATCHES ON)

###############################################################
## ConvectionDiffusionApplication core library (C++ parts)
add_library( KratosConvectionDiffusionCore SHARED ${KRATOS_CONVECTION_DIFFUSION_APPLICATION_CORE_SOURCES} ${KRATOS_CONVECTION_DIFFUSION_TESTING_SOURCES} )
target_link_libraries( KratosConvectionDiffusionCore PUBLIC KratosCore)
set_target_properties( KratosConvectionDiffusionCore PROPERTIES COMPILE_DEFINITIONS "CONVECTION_DIFFUSION_APPLICATION=EXPORT,API")

## FluidDynamicsApplication python module
pybind11_add_module( KratosConvectionDiffusionApplication MODULE THIN_LTO ${KRATOS_CONVECTION_DIFFUSION_PYTHON_INTERFACE_SOURCES} )
target_link_libraries( KratosConvectionDiffusionApplication PRIVATE KratosConvectionDiffusionCore)
set_target_properties( KratosConvectionDiffusionApplication PROPERTIES PREFIX "")

# changing the .dll suffix to .pyd
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    set_target_properties(KratosConvectionDiffusionApplication PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

# changing the .dylib suffix to .so
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    set_target_properties(KratosConvectionDiffusionApplication PROPERTIES SUFFIX .so)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

if(${INSTALL_PYTHON_FILES} MATCHES ON)
  get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
  install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/python_scripts DESTINATION applications/${CURRENT_DIR_NAME}  FILES_MATCHING PATTERN "*.py"  PATTERN ".svn" EXCLUDE)
endif(${INSTALL_PYTHON_FILES} MATCHES ON)

###############################################################################
if(USE_COTIRE MATCHES ON)
    cotire(KratosConvectionDiffusionCore)
    cotire(KratosConvectionDiffusionApplication)
endif(USE_COTIRE MATCHES ON)

###############################################################################
## installing the resulting libraries

install(TARGETS KratosConvectionDiffusionCore DESTINATION libs )
install(TARGETS KratosConvectionDiffusionApplication DESTINATION libs )

# get_property(inc_dirs DIRECTORY PROPERTY INCLUDE_DIRECTORIES)
# message("KratosIncompressibleFluidApplication subdir inc_dirs = ${inc_dirs}")

# Add to the KratosMultiphisics Python module
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/ConvectionDiffusionApplication.py" DESTINATION KratosMultiphysics )
