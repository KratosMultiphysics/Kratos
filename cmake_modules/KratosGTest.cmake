# This function automatically configures a given application to build its tests
macro(kratos_add_gtests)
	set(options USE_MPI)
	set(oneValueArgs TARGET WORKING_DIRECTORY)
	set(multiValueArgs SOURCES)
	
	cmake_parse_arguments(KRATOS_ADD_GTEST "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

	include(GoogleTest)

	if(KRATOS_ADD_GTEST_USE_MPI)
		set(TESTING_MPI_UTILITIES "KratosMPICoreTestUtilities")
	endif()
	
	add_executable("${KRATOS_ADD_GTEST_TARGET}Test" ${KRATOS_ADD_GTEST_SOURCES})
	target_link_libraries("${KRATOS_ADD_GTEST_TARGET}Test" ${KRATOS_ADD_GTEST_TARGET} KratosCoreTestUtilities "${TESTING_MPI_UTILITIES}" GTest::gtest_main GTest::gmock_main)
	set_target_properties("${KRATOS_ADD_GTEST_TARGET}Test" PROPERTIES COMPILE_DEFINITIONS "KRATOS_TEST_CORE=IMPORT,API")

	install(TARGETS ${KRATOS_ADD_GTEST_TARGET}Test DESTINATION test)

	if(DEFINED KRATOS_ADD_GTEST_WORKING_DIRECTORY)
		gtest_discover_tests(${KRATOS_ADD_GTEST_TARGET}Test DISCOVERY_MODE PRE_TEST WORKING_DIRECTORY "${KRATOS_ADD_GTEST_WORKING_DIRECTORY}")
	else()
		gtest_discover_tests(${KRATOS_ADD_GTEST_TARGET}Test DISCOVERY_MODE PRE_TEST )
	endif()
endmacro(kratos_add_gtests)
