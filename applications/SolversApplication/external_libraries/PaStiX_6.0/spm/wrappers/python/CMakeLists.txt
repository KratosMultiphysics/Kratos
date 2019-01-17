###
#
#  @copyright 2017-2018 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
#                       Univ. Bordeaux. All rights reserved.
#
#  @version 6.0.0
#  @author Mathieu Faverge
#  @date 2018-05-14
#
###

# Configure enum.py
if (SPM_INT64)
  set(SPM_PYTHON_INTEGER c_int64)
else()
  set(SPM_PYTHON_INTEGER c_int)
endif()

configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/spm/enum.py.in"
  "${CMAKE_CURRENT_SOURCE_DIR}/spm/enum.py" @ONLY)

# Install python package
install(FILES
  ${CMAKE_CURRENT_SOURCE_DIR}/spm/__init__.py
  ${CMAKE_CURRENT_SOURCE_DIR}/spm/__spm__.py
  ${CMAKE_CURRENT_SOURCE_DIR}/spm/spm.py
  ${CMAKE_CURRENT_SOURCE_DIR}/spm/enum.py
  DESTINATION lib/python/spm )

# Install python examples
install(FILES
  ${CMAKE_CURRENT_SOURCE_DIR}/spm_driver.py
  ${CMAKE_CURRENT_SOURCE_DIR}/spm_scipy.py
  DESTINATION examples
  )

## CTest execution
find_package(PythonInterp QUIET)
if (PYTHONINTERP_FOUND)
  set( PYTHON_TESTS
    spm_driver spm_scipy )

  foreach(example ${PYTHON_TESTS} )
    add_test(python_${example} ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/${example}.py)
  endforeach()
endif()

