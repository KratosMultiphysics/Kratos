set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosRANSModellingApplication ****")

# check for fluid dynamics application
if(NOT (${FLUID_DYNAMICS_APPLICATION} MATCHES ON))
  message(FATAL_ERROR "RANSModellingApplication requires FluidDynamicsApplication, Please re-compile with -DFLUID_DYNAMICS_APPLICATION=ON flag.")
endif(NOT (${FLUID_DYNAMICS_APPLICATION} MATCHES ON))

################### PYBIND11
include(pybind11Tools)

include_directories( ${CMAKE_SOURCE_DIR}/kratos )
include_directories( ${CMAKE_SOURCE_DIR}/applications/FluidDynamicsApplication )

# generate variables with the sources
set( KRATOS_RANS_MODELLING_APPLICATION_CORE
    ${CMAKE_CURRENT_SOURCE_DIR}/rans_modelling_application.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/rans_modelling_application_variables.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_strategies/relaxed_dof_updater.cpp
)

# add element sources
set( KRATOS_RANS_MODELLING_APPLICATION_CORE ${KRATOS_RANS_MODELLING_APPLICATION_CORE}
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/evm_k_epsilon/evm_k_epsilon_utilities.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/evm_k_epsilon/rans_evm_k_epsilon_k.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/evm_k_epsilon/rans_evm_k_epsilon_low_re_k.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/evm_k_epsilon/rans_evm_k_epsilon_epsilon.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/evm_k_epsilon/rans_evm_k_epsilon_low_re_epsilon.cpp
)

# add condition sources
set( KRATOS_RANS_MODELLING_APPLICATION_CORE ${KRATOS_RANS_MODELLING_APPLICATION_CORE}
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/evm_k_epsilon/rans_evm_k_epsilon_epsilon_wall.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/evm_k_epsilon/rans_evm_k_epsilon_vms_monolithic_wall.cpp
)
# add utilities sources
set( KRATOS_RANS_MODELLING_APPLICATION_CORE ${KRATOS_RANS_MODELLING_APPLICATION_CORE}
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/rans_calculation_utilities.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/rans_variable_utilities.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_elements/evm_k_epsilon/evm_k_epsilon_utilities.cpp
)

# add auxiliary process sources
set( KRATOS_RANS_MODELLING_APPLICATION_CORE ${KRATOS_RANS_MODELLING_APPLICATION_CORE}
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_processes/auxiliary_processes/rans_apply_exact_nodal_periodic_condition_process.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_processes/auxiliary_processes/rans_apply_flag_process.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_processes/auxiliary_processes/rans_check_scalar_bounds_process.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_processes/auxiliary_processes/rans_check_vector_bounds_process.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_processes/auxiliary_processes/rans_clip_scalar_variable_process.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_processes/auxiliary_processes/rans_epsilon_turbulent_mixing_inlet_process.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_processes/auxiliary_processes/rans_find_condition_parent_process.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_processes/auxiliary_processes/rans_k_turbulent_intensity_inlet_process.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_processes/auxiliary_processes/rans_line_output_process.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_processes/auxiliary_processes/rans_logarithmic_y_plus_calculation_process.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_processes/auxiliary_processes/rans_nut_k_epsilon_high_re_calculation_process.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_processes/auxiliary_processes/rans_nut_low_re_calculation_process.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_processes/auxiliary_processes/rans_nut_y_plus_wall_function_process.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_processes/auxiliary_processes/rans_scalar_cell_center_averaging_process.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_processes/auxiliary_processes/rans_vector_align_process.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_processes/auxiliary_processes/rans_vector_cell_center_averaging_process.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_processes/auxiliary_processes/rans_wall_distance_calculation_process.cpp
)

set( KRATOS_RANS_MODELLING_APPLICATION_PYTHON_INTERFACE
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/rans_modelling_python_application.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_utilities_to_python.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_strategies_to_python.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_strategies_to_python.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_processes_to_python.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_auxiliary_processes_to_python.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_solving_processes_to_python.cpp
)

set(KRATOS_RANS_MODELLING_APPLICATION_SOURCES_MPI "")
set(KRATOS_RANS_MODELLING_APPLICATION_EXTRA_LIBS "")
if(${MPI_NEEDED} MATCHES ON )
    include_directories( ${CMAKE_SOURCE_DIR}/applications/TrilinosApplication )
    find_package(TRILINOS REQUIRED)
    include_directories(SYSTEM ${TRILINOS_INCLUDE_DIR}) # Included as "SYSTEM" to suppress warnings from it

    set(KRATOS_RANS_MODELLING_APPLICATION_EXTRA_LIBS ${KRATOS_RANS_MODELLING_APPLICATION_EXTRA_LIBS} ${MPI_LIBRARIES})
    set(KRATOS_RANS_MODELLING_APPLICATION_EXTRA_LIBS ${KRATOS_RANS_MODELLING_APPLICATION_EXTRA_LIBS} ${TRILINOS_LIBRARIES})
endif(${MPI_NEEDED} MATCHES ON )

if(${KRATOS_BUILD_TESTING} MATCHES ON)
  file(GLOB_RECURSE
    KRATOS_RANS_MODELLING_APPLICATION_TESTING
    ${CMAKE_CURRENT_SOURCE_DIR}/tests/cpp/*.cpp
  )
endif(${KRATOS_BUILD_TESTING} MATCHES ON)

add_library(KratosRANSModellingCore SHARED ${KRATOS_RANS_MODELLING_APPLICATION_CORE} ${KRATOS_RANS_MODELLING_APPLICATION_TESTING})
target_link_libraries(KratosRANSModellingCore PUBLIC KratosCore KratosFluidDynamicsCore ${KRATOS_RANS_MODELLING_APPLICATION_EXTRA_LIBS})
set_target_properties(KratosRANSModellingCore PROPERTIES COMPILE_DEFINITIONS "RANS_MODELLING_APPLICATION=EXPORT,API")

pybind11_add_module(KratosRANSModellingApplication MODULE THIN_LTO ${KRATOS_RANS_MODELLING_APPLICATION_PYTHON_INTERFACE})
target_link_libraries(KratosRANSModellingApplication PRIVATE KratosRANSModellingCore)
set_target_properties(KratosRANSModellingApplication PROPERTIES PREFIX "")

if(USE_COTIRE MATCHES ON)
    cotire(KratosRANSModellingApplication)
endif(USE_COTIRE MATCHES ON)

install(TARGETS KratosRANSModellingCore DESTINATION libs )
install(TARGETS KratosRANSModellingApplication DESTINATION libs )

# changing the .dll suffix to .pyd (Windows)
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    set_target_properties(KratosRANSModellingApplication PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

# changing the .dylib suffix to .so (OS X)
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    set_target_properties(KratosRANSModellingApplication PROPERTIES SUFFIX .so)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

# install the python files
if(${INSTALL_PYTHON_FILES} MATCHES ON)
  get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
  install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/python_scripts DESTINATION applications/${CURRENT_DIR_NAME}  FILES_MATCHING PATTERN "*.py"  PATTERN ".svn" EXCLUDE)
endif(${INSTALL_PYTHON_FILES} MATCHES ON)

# Add to the KratosMultiphisics Python module
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/RANSModellingApplication.py" DESTINATION "KratosMultiphysics/RANSModellingApplication" RENAME "__init__.py")
