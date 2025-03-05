if(EXISTS "/home/camarotti/Software/kratos/scripts/kratos/KratosCoreTest")
  if(NOT EXISTS "/home/camarotti/Software/kratos/scripts/kratos/KratosCoreTest[1]_tests.cmake" OR
     NOT "/home/camarotti/Software/kratos/scripts/kratos/KratosCoreTest[1]_tests.cmake" IS_NEWER_THAN "/home/camarotti/Software/kratos/scripts/kratos/KratosCoreTest" OR
     NOT "/home/camarotti/Software/kratos/scripts/kratos/KratosCoreTest[1]_tests.cmake" IS_NEWER_THAN "${CMAKE_CURRENT_LIST_FILE}")
    include("/usr/share/cmake-3.28/Modules/GoogleTestAddTests.cmake")
    gtest_discover_tests_impl(
      TEST_EXECUTABLE [==[/home/camarotti/Software/kratos/scripts/kratos/KratosCoreTest]==]
      TEST_EXECUTOR [==[]==]
      TEST_WORKING_DIR [==[/home/camarotti/Software/kratos/scripts/kratos]==]
      TEST_EXTRA_ARGS [==[]==]
      TEST_PROPERTIES [==[]==]
      TEST_PREFIX [==[]==]
      TEST_SUFFIX [==[]==]
      TEST_FILTER [==[]==]
      NO_PRETTY_TYPES [==[FALSE]==]
      NO_PRETTY_VALUES [==[FALSE]==]
      TEST_LIST [==[KratosCoreTest_TESTS]==]
      CTEST_FILE [==[/home/camarotti/Software/kratos/scripts/kratos/KratosCoreTest[1]_tests.cmake]==]
      TEST_DISCOVERY_TIMEOUT [==[5]==]
      TEST_XML_OUTPUT_DIR [==[]==]
    )
  endif()
  include("/home/camarotti/Software/kratos/scripts/kratos/KratosCoreTest[1]_tests.cmake")
else()
  add_test(KratosCoreTest_NOT_BUILT KratosCoreTest_NOT_BUILT)
endif()
