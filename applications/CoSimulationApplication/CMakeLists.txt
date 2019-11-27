set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosCoSimulationApplication ****")

if(${CO_SIMULATION_APPLICATION_PYTHON} MATCHES ON)
    message("**** configuring KratosCoSimulationApplication for PYTHON-ONLY ****")
    project (KratosMultiphysics)
    cmake_minimum_required (VERSION 2.8.6)

    get_filename_component(PARENT_DIR ${CMAKE_CURRENT_SOURCE_DIR} DIRECTORY)
    get_filename_component(KRATOS_DIR ${PARENT_DIR} DIRECTORY)

    message(STATUS "Deleting: ${KRATOS_DIR}/KratosMultiphysics")
    file(REMOVE_RECURSE ${KRATOS_DIR}/KratosMultiphysics)
    message(STATUS "Deleting: ${KRATOS_DIR}/libs")
    file(REMOVE_RECURSE ${KRATOS_DIR}/libs)

    install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/custom_data_structure/KratosMultiphysics.py" DESTINATION ${KRATOS_DIR}/KratosMultiphysics RENAME "__init__.py")
    install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/custom_data_structure/CoSimulationApplication.py" DESTINATION ${KRATOS_DIR}/KratosMultiphysics/CoSimulationApplication RENAME "__init__.py")

else(${CO_SIMULATION_APPLICATION_PYTHON} MATCHES ON)

    ################### PYBIND11
    include(pybind11Tools)

    include_directories( ${CMAKE_SOURCE_DIR}/kratos )

    # generate variables with the sources
    set( KRATOS_CO_SIMULATION_APPLICATION_SOURCES
        ${CMAKE_CURRENT_SOURCE_DIR}/co_simulation_application.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/co_simulation_application_variables.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/co_simulation_python_application.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_io_to_python.cpp
    )

    ## define test sources (cpp tests)
    if(${KRATOS_BUILD_TESTING} MATCHES ON)
         file(GLOB_RECURSE KRATOS_CO_SIMULATION_TESTING_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/tests/cpp_tests/*.cpp)
    endif(${KRATOS_BUILD_TESTING} MATCHES ON)
         
    set(KRATOS_CO_SIMULATION_TESTING_CPP_PING_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/tests/cpp_ping_pong_solvers/cpp_test_solver_ping/cpp_test_solver_ping.cpp)
    set(KRATOS_CO_SIMULATION_TESTING_CPP_PONG_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/tests/cpp_ping_pong_solvers/cpp_test_solver_pong/cpp_test_solver_pong.cpp)
    # define library Kratos which defines the basic python interface
    pybind11_add_module(KratosCoSimulationApplication MODULE THIN_LTO ${KRATOS_CO_SIMULATION_APPLICATION_SOURCES} ${KRATOS_CO_SIMULATION_TESTING_SOURCES})
    target_link_libraries(KratosCoSimulationApplication PRIVATE KratosCore )
    set_target_properties(KratosCoSimulationApplication PROPERTIES PREFIX "")
    install(TARGETS KratosCoSimulationApplication DESTINATION libs )

    # making the cpp_test_solver_one and cpp_test_solver_two executables
    MESSAGE(${KRATOS_CO_SIMULATION_TESTING_SOURCES})
    add_executable (ping ${KRATOS_CO_SIMULATION_TESTING_CPP_PING_SOURCES})
    add_executable (pong ${KRATOS_CO_SIMULATION_TESTING_CPP_PONG_SOURCES})

    # changing the .dll suffix to .pyd (Windows)
    if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
        set_target_properties(KratosCoSimulationApplication PROPERTIES SUFFIX .pyd)
    endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

    # changing the .dylib suffix to .so (OS X)
    if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
        set_target_properties(KratosCoSimulationApplication PROPERTIES SUFFIX .so)
    endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

    # install the python files
    if(${INSTALL_PYTHON_FILES} MATCHES ON)
    get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
    install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/python_scripts DESTINATION applications/${CURRENT_DIR_NAME}  FILES_MATCHING PATTERN "*.py"  PATTERN ".svn" EXCLUDE)
    endif(${INSTALL_PYTHON_FILES} MATCHES ON)

    # Add to the KratosMultiphisics Python module
    install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/CoSimulationApplication.py" DESTINATION "KratosMultiphysics/CoSimulationApplication" RENAME "__init__.py")

endif(${CO_SIMULATION_APPLICATION_PYTHON} MATCHES ON)









