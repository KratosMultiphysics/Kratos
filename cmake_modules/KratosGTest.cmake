# This function automatically configures a given application to build its tests
macro(kratos_add_gtests)
    set(options USE_MPI USE_CUSTOM_MAIN)
    set(oneValueArgs TARGET WORKING_DIRECTORY)
    set(multiValueArgs SOURCES)

    cmake_parse_arguments(KRATOS_ADD_GTEST "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

    if(KRATOS_ADD_GTEST_SOURCES)
        include(GoogleTest)

        # Need to add the custom main initializing the verbosity and other customizations
        if(NOT KRATOS_ADD_GTEST_USE_CUSTOM_MAIN)
            set(KRATOS_GTEST_MAIN_SOURCE ${KRATOS_BASE_FOLDER}/kratos/testing/testing.cpp)
        endif()

        if(KRATOS_ADD_GTEST_USE_MPI)
            set(TESTING_MPI_UTILITIES "KratosMPICoreTestUtilities" ${MPI_LIBRARIES})
            # Need to add the custom main initializing the kernel with MPI
            if(NOT KRATOS_ADD_GTEST_USE_CUSTOM_MAIN)
                set(KRATOS_GTEST_MAIN_SOURCE ${KRATOS_BASE_FOLDER}/kratos/mpi/testing/mpi_testing.cpp)
            endif()
        endif()

        add_executable("${KRATOS_ADD_GTEST_TARGET}Test" ${KRATOS_ADD_GTEST_SOURCES} ${KRATOS_GTEST_MAIN_SOURCE})
        target_link_libraries("${KRATOS_ADD_GTEST_TARGET}Test" ${KRATOS_ADD_GTEST_TARGET} KratosCoreTestUtilities "${TESTING_MPI_UTILITIES}" GTest::gmock_main)
        set_target_properties("${KRATOS_ADD_GTEST_TARGET}Test" PROPERTIES COMPILE_DEFINITIONS "KRATOS_TEST_CORE=IMPORT,API")

        install(TARGETS ${KRATOS_ADD_GTEST_TARGET}Test DESTINATION test)

        if(DEFINED KRATOS_ADD_GTEST_WORKING_DIRECTORY)
            gtest_discover_tests(${KRATOS_ADD_GTEST_TARGET}Test DISCOVERY_MODE PRE_TEST WORKING_DIRECTORY "${KRATOS_ADD_GTEST_WORKING_DIRECTORY}")
        else()
            gtest_discover_tests(${KRATOS_ADD_GTEST_TARGET}Test DISCOVERY_MODE PRE_TEST )
        endif()
    endif()

endmacro(kratos_add_gtests)