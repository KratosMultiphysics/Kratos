###
#
# @copyright (c) 2009-2014 The University of Tennessee and The University
#                          of Tennessee Research Foundation.
#                          All rights reserved.
# @copyright (c) 2012-2016 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria, Univ. Bordeaux. All rights reserved.
#
###
#
#  @file RulesPrecisions.cmake
#
#  @project MORSE
#  MORSE is a software package provided by:
#     Inria Bordeaux - Sud-Ouest,
#     Univ. of Tennessee,
#     King Abdullah Univesity of Science and Technology
#     Univ. of California Berkeley,
#     Univ. of Colorado Denver.
#
#  @version 0.9.0
#  @author Thomas Herault
#  @author George Bosilca
#  @author Aurelien Bouteiller
#  @author Mathieu Faverge
#  @author Florent Puvost
#  @date 13-07-2012
#
###
cmake_minimum_required(VERSION 2.8)
include(CMakeDependentOption)

set(RP_GENDEPENDENCIES ${MORSE_CMAKE_MODULE_PATH}/precision_generator/genDependencies.py)
set(RP_CODEGEN         ${MORSE_CMAKE_MODULE_PATH}/precision_generator/codegen.py)

# Dictionnary
# -----------
if( NOT DEFINED RP_${CMAKE_PROJECT_NAME}_DICTIONNARY )
  message( WARNING "RulesPrecisions included before RP_${CMAKE_PROJECT_NAME}_DICTIONNARY was defined (Default is used)" )
  set(RP_${CMAKE_PROJECT_NAME}_DICTIONNARY ${MORSE_CMAKE_MODULE_PATH}/precision_generator/subs.py
    CACHE INTERNAL "Dictionnary used for precision generation" )
else()
  set(RP_${CMAKE_PROJECT_NAME}_DICTIONNARY ${RP_${CMAKE_PROJECT_NAME}_DICTIONNARY}
    CACHE INTERNAL "Dictionnary used for precision generation" )
endif()

# Default Precisions
# ------------------
if( NOT DEFINED RP_${CMAKE_PROJECT_NAME}_PRECISIONS )
  message( WARNING "RulesPrecisions included before RP_${CMAKE_PROJECT_NAME}_PRECISIONS was defined (\"s;d;c;z\" is used)" )
  set(RP_${CMAKE_PROJECT_NAME}_PRECISIONS "s;d;c;z"
    CACHE INTERNAL "Set of available precisions for the project" )
else()
  set(RP_${CMAKE_PROJECT_NAME}_PRECISIONS ${RP_${CMAKE_PROJECT_NAME}_PRECISIONS}
    CACHE INTERNAL "Set of available precisions for the project" )
endif()

# Detect if compilation is done in or out of place
# ------------------------------------------------
string(COMPARE EQUAL "${CMAKE_SOURCE_DIR}" "${CMAKE_BINARY_DIR}" __inplace)
set( ${CMAKE_PROJECT_NAME}_COMPILE_INPLACE ${__inplace}
  CACHE INTERNAL "Defines if the project compilation is made inplace or not" )

# Detect default settings
# -----------------------
set( _prec_S  OFF )
set( _prec_D  OFF )
set( _prec_DS OFF )
set( _prec_C  OFF )
set( _prec_Z  OFF )
set( _prec_ZC OFF )
foreach(_prec ${RP_${CMAKE_PROJECT_NAME}_PRECISIONS})
  if ( ${_prec} STREQUAL zc )
    set( _prec_S  ON )
    set( _prec_D  ON )
    set( _prec_C  ON )
    set( _prec_Z  ON )
    set( _prec_ZC ON )
  elseif( ${_prec} STREQUAL z )
    set( _prec_D  ON )
    set( _prec_Z  ON )
  elseif( ${_prec} STREQUAL c )
    set( _prec_S  ON )
    set( _prec_C  ON )
  elseif( ${_prec} STREQUAL ds )
    set( _prec_S  ON )
    set( _prec_D  ON )
    set( _prec_DS ON )
  elseif( ${_prec} STREQUAL d )
    set( _prec_S  ON )
  elseif( ${_prec} STREQUAL s )
    set( _prec_S  ON )
  endif()
endforeach()

# set computed precisions
# -----------------------
if( ${_prec_S} )
  option(${PROJECT_NAME}_PREC_S
    "Build ${PROJECT_NAME} real single precision" ON)
endif()
if( ${_prec_D} )
  option(${PROJECT_NAME}_PREC_D
    "Build ${PROJECT_NAME} real double precision" ON)
endif()

if( ${_prec_C} )
  cmake_dependent_option(${PROJECT_NAME}_PREC_C
    "Build ${PROJECT_NAME} complex single precision" ON "${PROJECT_NAME}_PREC_S" OFF)
endif()
if( ${_prec_Z} )
  cmake_dependent_option(${PROJECT_NAME}_PREC_Z
    "Build ${PROJECT_NAME} complex double precision" ON "${PROJECT_NAME}_PREC_D" OFF)
endif()

if( ${_prec_DS} )
  cmake_dependent_option(${PROJECT_NAME}_PREC_DS
    "Build ${PROJECT_NAME} real mixed precision"    ON "${PROJECT_NAME}_PREC_S;${PROJECT_NAME}_PREC_D" OFF)
endif()
if( ${_prec_ZC} )
  cmake_dependent_option(${PROJECT_NAME}_PREC_ZC
    "Build ${PROJECT_NAME} complex mixed precision" ON "${PROJECT_NAME}_PREC_C;${PROJECT_NAME}_PREC_Z" OFF)
endif()

# Define precision supported by the project
# -----------------------------------------
set( ${CMAKE_PROJECT_NAME}_PRECISION "" )
if(${CMAKE_PROJECT_NAME}_PREC_S)
  list(APPEND ${CMAKE_PROJECT_NAME}_PRECISION "s")
endif()

if(${CMAKE_PROJECT_NAME}_PREC_D)
  list(APPEND ${CMAKE_PROJECT_NAME}_PRECISION "d")
endif()

if(${CMAKE_PROJECT_NAME}_PREC_DS)
  list(APPEND ${CMAKE_PROJECT_NAME}_PRECISION "ds")
endif()

if(${CMAKE_PROJECT_NAME}_PREC_C)
  list(APPEND ${CMAKE_PROJECT_NAME}_PRECISION "c")
endif()

if(${CMAKE_PROJECT_NAME}_PREC_Z)
  list(APPEND ${CMAKE_PROJECT_NAME}_PRECISION "z")
endif()

if(${CMAKE_PROJECT_NAME}_PREC_ZC)
  list(APPEND ${CMAKE_PROJECT_NAME}_PRECISION "zc")
endif()

#
# Generates a rule for every SOURCES file, to create the precisions in PRECISIONS. If TARGETDIR
# is not empty then all generated files will be prepended with the $TARGETDIR/.
# A new file is created, from a copy by default
# If the first precision is "/", all occurences of the basename in the file are remplaced by
# "pbasename" where p is the selected precision.
# the target receives a -DPRECISION_p in its cflags.
#
include(ParseArguments)

# Add a hint to help Cmake to find the correct python version:
# (see https://cmake.org/cmake/help/v3.0/module/FindPythonInterp.html)
set(Python_ADDITIONAL_VERSIONS 2)
find_package(PythonInterp REQUIRED)

MACRO(precisions_rules_py)
  PARSE_ARGUMENTS(PREC_RULE "TARGETDIR;PRECISIONS" "" ${ARGN})

  message(STATUS "Generate precision dependencies in ${CMAKE_CURRENT_SOURCE_DIR}")

  # The first is the output variable list
  CAR(OUTPUTLIST ${PREC_RULE_DEFAULT_ARGS})
  # Everything else should be source files.
  CDR(SOURCES ${PREC_RULE_DEFAULT_ARGS})

  if( NOT DEFINED ${CMAKE_PROJECT_NAME}_COMPILE_INPLACE )
    message( FATAL_ERROR "precisions_rules_init must be called before any call to precisions_rules_py" )
  endif()

  if( NOT DEFINED RP_${CMAKE_PROJECT_NAME}_DICTIONNARY )
    message( FATAL_ERROR "precisions_rules_init must be called before any call to precisions_rules_py" )
  endif()

  # By default the TARGETDIR is the current binary directory
  if( "${PREC_RULE_TARGETDIR}" STREQUAL "" )
    set(PREC_RULE_TARGETDIR "./")
    set(PRECISIONPP_prefix "./")
    set(PRECISIONPP_arg "-P")
  else( "${PREC_RULE_TARGETDIR}" STREQUAL "" )
    if(EXISTS ${CMAKE_CURRENT_BINARY_DIR}/${PREC_RULE_TARGETDIR})
    else(EXISTS ${CMAKE_CURRENT_BINARY_DIR}/${PREC_RULE_TARGETDIR})
      file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${PREC_RULE_TARGETDIR})
    endif(EXISTS ${CMAKE_CURRENT_BINARY_DIR}/${PREC_RULE_TARGETDIR})
    set(PRECISIONPP_arg "-P")
    set(PRECISIONPP_prefix "${PREC_RULE_TARGETDIR}")
  endif( "${PREC_RULE_TARGETDIR}" STREQUAL "" )

  set(options_list "")
  foreach(prec_rules_PREC ${PREC_RULE_PRECISIONS})
    set(options_list "${options_list} ${prec_rules_PREC}")
  endforeach()

  set(sources_list "")
  foreach(_src ${SOURCES})
    set(sources_list "${sources_list} ${_src}")
  endforeach()

  set(gencmd ${PYTHON_EXECUTABLE} ${RP_GENDEPENDENCIES} -f "${sources_list}" -p "${options_list}" -s "${CMAKE_CURRENT_SOURCE_DIR}" ${PRECISIONPP_arg} ${PRECISIONPP_prefix})
  EXECUTE_PROCESS(COMMAND ${gencmd} OUTPUT_VARIABLE dependencies_list)

  foreach(_dependency ${dependencies_list})

    string(STRIP "${_dependency}" _dependency)
    string(COMPARE NOTEQUAL "${_dependency}" "" not_empty)
    if( not_empty )

      string(REGEX REPLACE "^(.*),(.*),(.*)$" "\\1" _dependency_INPUT "${_dependency}")
      set(_dependency_PREC   "${CMAKE_MATCH_2}")
      set(_dependency_OUTPUT "${CMAKE_MATCH_3}")

      set(pythoncmd ${PYTHON_EXECUTABLE} ${RP_CODEGEN} -f ${CMAKE_CURRENT_SOURCE_DIR}/${_dependency_INPUT} -p ${_dependency_PREC} ${PRECISIONPP_arg} ${PRECISIONPP_prefix})

      string(STRIP "${_dependency_OUTPUT}" _dependency_OUTPUT)
      string(COMPARE NOTEQUAL "${_dependency_OUTPUT}" "" got_file)

      # Force the copy of the original files in the binary_dir
      # for VPATH compilation
      if( NOT ${CMAKE_PROJECT_NAME}_COMPILE_INPLACE )
	set(generate_out 1)
      else( NOT ${CMAKE_PROJECT_NAME}_COMPILE_INPLACE )
	string(COMPARE NOTEQUAL "${_dependency_OUTPUT}" "${_dependency_INPUT}" generate_out )
      endif()

      # We generate a dependency only if a file will be generated
      if( got_file )
        set( _compile_flags "-DPRECISION_${_dependency_PREC}" )
        set( _listtmp ${PREC_RULE_PRECISIONS})
        list(REMOVE_ITEM _listtmp ${_dependency_PREC})
        foreach( _prec  ${_listtmp})
          set( _compile_flags "${_compile_flags} -UPRECISION_${_prec}" )
        endforeach()
	if( generate_out )
	  # the custom command is executed in CMAKE_CURRENT_BINARY_DIR
	  ADD_CUSTOM_COMMAND(
	    OUTPUT ${_dependency_OUTPUT}
	    COMMAND ${CMAKE_COMMAND} -E remove -f ${_dependency_OUTPUT} && ${pythoncmd} && chmod a-w ${_dependency_OUTPUT}
	    DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/${_dependency_INPUT} ${RP_CODEGEN} ${RP_${CMAKE_PROJECT_NAME}_DICTIONNARY})
          set_source_files_properties(${_dependency_OUTPUT} PROPERTIES COMPILE_FLAGS ${_compile_flags} GENERATED 1 IS_IN_BINARY_DIR 1 )
	else( generate_out )
          set_source_files_properties(${_dependency_OUTPUT} PROPERTIES COMPILE_FLAGS ${_compile_flags} GENERATED 0 )
	endif( generate_out )

	list(APPEND ${OUTPUTLIST} ${_dependency_OUTPUT})
      endif( got_file )
    endif()
  endforeach()

  message(STATUS "Generate precision dependencies in ${CMAKE_CURRENT_SOURCE_DIR} - Done")

ENDMACRO(precisions_rules_py)

##
## @end file RulesPrecisions.cmake
##
