set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosCoSimulationApplication ****")

################### PYBIND11
include(pybind11Tools)

include_directories( ${KRATOS_SOURCE_DIR}/kratos )

## CoSimulation Core sources
file(GLOB_RECURSE KRATOS_CO_SIMULATION_APPLICATION_CORE
    ${CMAKE_CURRENT_SOURCE_DIR}/co_simulation_application.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/co_simulation_application_variables.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_io/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_processes/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/*.cpp
)

## CoSimulation testing sources
if(${KRATOS_BUILD_TESTING} MATCHES ON)
    file(GLOB_RECURSE KRATOS_CO_SIMULATION_APPLICATION_TESTING_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/tests/cpp_tests/*.cpp)
endif(${KRATOS_BUILD_TESTING} MATCHES ON)

## CoSimulation python interface sources
file(GLOB_RECURSE KRATOS_CO_SIMULATION_APPLICATION_PYTHON_INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/*.cpp)

set(KRATOS_CO_SIMULATION_TESTING_CPP_PING_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/tests/cpp_ping_pong_solvers/cpp_test_solver_ping/cpp_test_solver_ping.cpp)
set(KRATOS_CO_SIMULATION_TESTING_CPP_PONG_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/tests/cpp_ping_pong_solvers/cpp_test_solver_pong/cpp_test_solver_pong.cpp)


add_library(KratosCoSimulationCore SHARED ${KRATOS_CO_SIMULATION_APPLICATION_CORE} ${KRATOS_CO_SIMULATION_APPLICATION_TESTING_SOURCES})
target_link_libraries(KratosCoSimulationCore PUBLIC KratosCore)
set_target_properties(KratosCoSimulationCore PROPERTIES COMPILE_DEFINITIONS "CO_SIMULATION_APPLICATION=EXPORT,API")

###############################################################
## define library Kratos which defines the basic python interface
pybind11_add_module(KratosCoSimulationApplication MODULE THIN_LTO ${KRATOS_CO_SIMULATION_APPLICATION_PYTHON_INTERFACE})
target_link_libraries(KratosCoSimulationApplication PRIVATE KratosCoSimulationCore)
set_target_properties(KratosCoSimulationApplication PROPERTIES PREFIX "")

# Set batch size in the unity build
IF(CMAKE_UNITY_BUILD MATCHES ON)
    set_target_properties(KratosCoSimulationCore PROPERTIES UNITY_BUILD_BATCH_SIZE ${KRATOS_UNITY_BUILD_BATCH_SIZE})
    set_target_properties(KratosCoSimulationApplication PROPERTIES UNITY_BUILD_BATCH_SIZE ${KRATOS_UNITY_BUILD_BATCH_SIZE})
ENDIF(CMAKE_UNITY_BUILD MATCHES ON)

# Making the cpp_test_solver_one and cpp_test_solver_two executables
add_executable(ping ${KRATOS_CO_SIMULATION_TESTING_CPP_PING_SOURCES})
add_executable(pong ${KRATOS_CO_SIMULATION_TESTING_CPP_PONG_SOURCES})

# changing the .dll suffix to .pyd (Windows)
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    set_target_properties(KratosCoSimulationApplication PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

# changing the .dylib suffix to .so (OS X)
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    set_target_properties(KratosCoSimulationApplication PROPERTIES SUFFIX .so)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

# Add to the KratosMultiphisics Python module
kratos_python_install(${INSTALL_PYTHON_USING_LINKS} ${CMAKE_CURRENT_SOURCE_DIR}/CoSimulationApplication.py KratosMultiphysics/CoSimulationApplication/__init__.py )

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
install(TARGETS KratosCoSimulationCore DESTINATION libs )
install(TARGETS KratosCoSimulationApplication DESTINATION libs )

# Define custom targets
set(KRATOS_KERNEL "${KRATOS_KERNEL};KratosCoSimulationCore" PARENT_SCOPE)
set(KRATOS_PYTHON_INTERFACE "${KRATOS_PYTHON_INTERFACE};KratosCoSimulationApplication" PARENT_SCOPE)