###############################################################################
#####
#####         Generation of a fortran header file
#####
###############################################################################

MACRO ( GENERATE_FORTRAN_HEADER name
    in_dir in_file include_dir out_dir out_file
    )
  # Wrap add_custom_command into add_custom target to remove dpendencies from
  # the custom command and thus allow parallel build.
  ADD_CUSTOM_COMMAND (
    OUTPUT ${out_dir}/${out_file}
    COMMAND genheader ${out_dir}/${out_file} ${in_dir}/${in_file} ${include_dir}
    ${PROJECT_SOURCE_DIR}/scripts/genfort.pl
    WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
    DEPENDS genheader ${in_dir}/${in_file}
    ${PROJECT_SOURCE_DIR}/scripts/genfort.pl
    COMMENT "Generating Fortran header for ${name}"
    )

  ADD_CUSTOM_TARGET (
    ${name}_fortran_header
    ALL
    DEPENDS ${out_dir}/${out_file}
    )

ENDMACRO ( )

###############################################################################
#####
#####         Copy an automatically generated header file to another place
#####
###############################################################################
MACRO ( COPY_FORTRAN_HEADER
    in_dir in_file out_dir out_file
    file_dependencies
    target_name
    )
  # Wrap add_custom_command into add_custom target to remove dpendencies from
  # the custom command and thus allow parallel build.
  ADD_CUSTOM_COMMAND (
    OUTPUT  ${out_dir}/${out_file}
    COMMAND ${CMAKE_COMMAND} -E copy  ${in_dir}/${in_file} ${out_dir}/${out_file}
    WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
    DEPENDS ${file_dependencies} ${in_dir}/${in_file}
    COMMENT "Copying ${in_dir}/${in_file} in ${out_dir}/${out_file}"
    )

  ADD_CUSTOM_TARGET ( ${target_name} ALL
    DEPENDS ${out_dir}/${out_file} )
  ADD_DEPENDENCIES ( ${target_name} ${file_dependencies} )

ENDMACRO ( )


###############################################################################
#####
#####         Copy an automatically generated header file to another place
#####         and create the associated target
#####
###############################################################################
MACRO ( COPY_FORTRAN_HEADER_AND_CREATE_TARGET
    binary_dir include_dir target_identifier )

  COPY_FORTRAN_HEADER (
    ${COMMON_BINARY_DIR} libmmgtypesf.h
    ${include_dir} libmmgtypesf.h
    mmg_fortran_header copy${target_identifier}_libmmgtypesf )

  COPY_FORTRAN_HEADER (
    ${binary_dir}
    libmmg${target_identifier}f.h ${include_dir}
    libmmg${target_identifier}f.h
    mmg${target_identifier}_fortran_header copy_libmmg${target_identifier}f
    )

  ADD_CUSTOM_TARGET(copy_${target_identifier}_headers ALL
    DEPENDS
    copy_libmmg${target_identifier}f copy${target_identifier}_libmmgtypesf
    ${include_dir}/libmmg${target_identifier}.h
    ${include_dir}/libmmgtypes.h )

ENDMACRO ( )

###############################################################################
#####
#####         Add a library to build and needed include dir, set its
#####         properties, add link dependencies and the install rule
#####
###############################################################################

MACRO ( ADD_AND_INSTALL_LIBRARY
    target_name target_type sources output_name )

  ADD_LIBRARY ( ${target_name} ${target_type} ${sources} )
  ADD_LIBRARY ( Mmg::${target_name} ALIAS ${target_name} )
  
  IF ( CMAKE_VERSION VERSION_LESS 2.8.12 )
    INCLUDE_DIRECTORIES ( ${target_name} PUBLIC
      ${COMMON_BINARY_DIR} ${COMMON_SOURCE_DIR} ${PROJECT_BINARY_DIR}/include )
  ELSE ( )
    target_include_directories( ${target_name} PUBLIC
      $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/include/>
      $<BUILD_INTERFACE:${COMMON_SOURCE_DIR}>
      $<BUILD_INTERFACE:${COMMON_BINARY_DIR}>
      $<BUILD_INTERFACE:${MMG3D_SOURCE_DIR}>
      $<BUILD_INTERFACE:${MMGS_SOURCE_DIR}>
      $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src/mmg/>
      $<BUILD_INTERFACE:${MMG2D_BINARY_DIR}>
      $<BUILD_INTERFACE:${MMG3D_BINARY_DIR}>
      $<BUILD_INTERFACE:${MMGS_BINARY_DIR}>
      $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}> )

  ENDIF ( )

  SET_TARGET_PROPERTIES ( ${target_name} PROPERTIES
    OUTPUT_NAME ${output_name}
    VERSION ${CMAKE_RELEASE_VERSION_MAJOR}.${CMAKE_RELEASE_VERSION_MINOR}.${CMAKE_RELEASE_VERSION_PATCH}
    SOVERSION ${CMAKE_RELEASE_VERSION_MAJOR} )

  SET_PROPERTY(TARGET ${target_name} PROPERTY C_STANDARD 99)

  TARGET_LINK_LIBRARIES ( ${target_name} ${LIBRARIES} )

  IF (NOT CMAKE_INSTALL_LIBDIR)
    SET(CMAKE_INSTALL_LIBDIR lib)
  ENDIF()

#  INSTALL ( TARGETS ${target_name}
#    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
#    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
#    COMPONENT lib
#    )


  SET ( MmgTargetsExported 1 )
  install(TARGETS ${target_name} EXPORT MmgTargets
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
    Component lib
    )
ENDMACRO ( )

###############################################################################
#####
#####         Add an executable to build and needed include dir, set its
#####         postfix, add link dependencies and the install rule
#####
###############################################################################

MACRO ( ADD_AND_INSTALL_EXECUTABLE
    exec_name lib_files main_file )

  IF ( NOT TARGET lib${exec_name}_a AND NOT TARGET lib${exec_name}_so )
    ADD_EXECUTABLE ( ${exec_name} ${lib_files} ${main_file} )
  ELSE ( )
    ADD_EXECUTABLE ( ${exec_name} ${main_file})

    SET_PROPERTY(TARGET ${exec_name} PROPERTY C_STANDARD 99)

    IF ( NOT TARGET lib${exec_name}_a )
      TARGET_LINK_LIBRARIES(${exec_name} lib${exec_name}_so)
    ELSE ( )
      TARGET_LINK_LIBRARIES(${exec_name} lib${exec_name}_a)
    ENDIF ( )

  ENDIF ( )

  IF ( WIN32 AND NOT MINGW AND USE_SCOTCH )
    my_add_link_flags ( ${exec_name} "/SAFESEH:NO")
  ENDIF ( )

 IF ( CMAKE_VERSION VERSION_LESS 2.8.12 )
   INCLUDE_DIRECTORIES ( ${exec_name} PUBLIC
     ${COMMON_BINARY_DIR} ${COMMON_SOURCE_DIR} ${PROJECT_BINARY_DIR}/include )
 ELSE ( )
   TARGET_INCLUDE_DIRECTORIES ( ${exec_name} PUBLIC
     ${COMMON_BINARY_DIR} ${COMMON_SOURCE_DIR} ${PROJECT_BINARY_DIR}/include )
 ENDIF ( )

  TARGET_LINK_LIBRARIES ( ${exec_name} ${LIBRARIES}  )

  INSTALL(TARGETS ${exec_name} RUNTIME DESTINATION bin COMPONENT appli)

  ADD_TARGET_POSTFIX(${exec_name})

ENDMACRO ( )


###############################################################################
#####
#####         Add a target postfix depending on the build type
#####
###############################################################################

MACRO ( ADD_TARGET_POSTFIX target_name )
  IF ( CMAKE_BUILD_TYPE MATCHES "Debug" )
    # in debug mode we name the executable mmgs_debug
    SET_TARGET_PROPERTIES(${target_name} PROPERTIES DEBUG_POSTFIX _debug)
  ELSEIF ( CMAKE_BUILD_TYPE MATCHES "Release" )
    # in Release mode we name the executable mmgs_O3
    SET_TARGET_PROPERTIES(${target_name} PROPERTIES RELEASE_POSTFIX _O3)
  ELSEIF ( CMAKE_BUILD_TYPE MATCHES "RelWithDebInfo" )
    # in RelWithDebInfo mode we name the executable mmgs_O3d
    SET_TARGET_PROPERTIES(${target_name} PROPERTIES RELWITHDEBINFO_POSTFIX _O3d)
  ELSEIF ( CMAKE_BUILD_TYPE MATCHES "MinSizeRel" )
    # in MinSizeRel mode we name the executable mmgs_O3
    SET_TARGET_PROPERTIES(${target_name} PROPERTIES MINSIZEREL_POSTFIX _Os)
  ENDIF ( )
ENDMACRO ( )

###############################################################################
#####
#####         Add Executable that must be tested by ci
#####
###############################################################################

MACRO ( ADD_EXEC_TO_CI_TESTS exec_name list_name )

  IF(${CMAKE_BUILD_TYPE} MATCHES "Debug")
    SET(EXECUT ${EXECUTABLE_OUTPUT_PATH}/${exec_name}_debug)
    SET(BUILDNAME ${BUILDNAME}_debug CACHE STRING "build name variable")
  ELSEIF(${CMAKE_BUILD_TYPE} MATCHES "Release")
    SET(EXECUT ${EXECUTABLE_OUTPUT_PATH}/${exec_name}_O3)
    SET(BUILDNAME ${BUILDNAME}_O3 CACHE STRING "build name variable")
  ELSEIF(${CMAKE_BUILD_TYPE} MATCHES "RelWithDebInfo")
    SET(EXECUT ${EXECUTABLE_OUTPUT_PATH}/${exec_name}_O3d)
    SET(BUILDNAME ${BUILDNAME}_O3d CACHE STRING "build name variable")
  ELSEIF(${CMAKE_BUILD_TYPE} MATCHES "MinSizeRel")
    SET(EXECUT ${EXECUTABLE_OUTPUT_PATH}/${exec_name}_Os)
    SET(BUILDNAME ${BUILDNAME}_Os CACHE STRING "build name variable")
  ELSE()
    SET(EXECUT ${EXECUTABLE_OUTPUT_PATH}/${exec_name})
    SET(BUILDNAME ${BUILDNAME} CACHE STRING "build name variable")
  ENDIF()

  SET ( ${list_name} ${EXECUT} )

ENDMACRO ( )


###############################################################################
#####
#####         Add a library test
#####
###############################################################################

MACRO ( ADD_LIBRARY_TEST target_name main_path target_dependency lib_name )
  ADD_EXECUTABLE ( ${target_name} ${main_path} )
  ADD_DEPENDENCIES( ${target_name} ${target_dependency} )

  IF ( CMAKE_VERSION VERSION_LESS 2.8.12 )
    INCLUDE_DIRECTORIES ( ${target_name} PUBLIC ${PROJECT_BINARY_DIR}/include )
  ELSE ( )
    TARGET_INCLUDE_DIRECTORIES ( ${target_name} PUBLIC ${PROJECT_BINARY_DIR}/include )
  ENDIF ( )

  IF ( WIN32 AND ((NOT MINGW) AND USE_SCOTCH) )
    MY_ADD_LINK_FLAGS ( ${target_name} "/SAFESEH:NO" )
  ENDIF ( )

  TARGET_LINK_LIBRARIES ( ${target_name}  ${lib_name} )
  INSTALL(TARGETS ${target_name} RUNTIME DESTINATION bin COMPONENT appli )

ENDMACRO ( )


###############################################################################
#####
#####         Add a test and run it again if RUN_AGAIN option is enabled
#####
###############################################################################

MACRO ( ADD_RUN_AGAIN_TESTS exec_name test_names args input_files )

  LIST(LENGTH test_names test_number)
  MATH(EXPR len2 "${test_number} - 1")

  FOREACH ( it RANGE ${len2} )
    LIST(GET test_names   ${it} test_name)
    LIST(GET input_files  ${it} input_file)
    LIST(GET args         ${it} _arg_)

    STRING(REPLACE " " ";" arg ${_arg_})

    ADD_TEST(NAME ${test_name}
      COMMAND ${exec_name} ${arg}
      ${input_file}
      -out ${CTEST_OUTPUT_DIR}/${test_name}-out.o.meshb )

    SET_TESTS_PROPERTIES ( ${test_name}
      PROPERTIES FIXTURES_SETUP ${test_name} )

    IF ( RUN_AGAIN )
      ADD_TEST(NAME ${test_name}_2
        COMMAND ${exec_name} ${arg} -hgrad -1
        ${CTEST_OUTPUT_DIR}/${test_name}-out.o.meshb
        -out ${CTEST_OUTPUT_DIR}/${test_name}_2-out.o.meshb
        )

      SET_TESTS_PROPERTIES ( ${test_name}_2
        PROPERTIES FIXTURES_REQUIRED ${test_name} )

    ENDIF ( )

  ENDFOREACH ( )


ENDMACRO ( )
