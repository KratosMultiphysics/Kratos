message("FOUND TRILINOS EXTENSIONS")

include_directories(
    ${CMAKE_SOURCE_DIR}/kratos
    ${CMAKE_SOURCE_DIR}/applications
    ${CMAKE_SOURCE_DIR}/applications/FluidDynamicsApplication
)
message("CMAKE_BUILD_DIR: ${CMAKE_BUILD_DIR}")
link_directories(
    ${CMAKE_SOURCE_DIR}/applications/FluidDynamicsApplication
#    ${CMAKE_BUILD_DIR}/applications/TrilinosApplication
)

file(
    GLOB_RECURSE
    KRATOS_FLUID_DYNAMICS_TRILINOS_EXTENSION_PYTHON_INTERFACE_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/*.cpp
)

## Python module
pybind11_add_module( KratosFluidDynamicsTrilinosExtension MODULE THIN_LTO ${KRATOS_FLUID_DYNAMICS_TRILINOS_EXTENSION_PYTHON_INTERFACE_SOURCES} )
target_link_libraries( KratosFluidDynamicsTrilinosExtension PRIVATE KratosFluidDynamicsCore, KratosTrilinosApplication, KratosMPICore)
set_target_properties( KratosFluidDynamicsTrilinosExtension PROPERTIES PREFIX "")
