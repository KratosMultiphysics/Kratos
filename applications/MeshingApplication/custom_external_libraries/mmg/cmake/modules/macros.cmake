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
#####         Copy a header file to another place
#####
###############################################################################
MACRO ( COPY_HEADER
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
#####         Copy a header file to another place
#####         and create the associated target
#####
###############################################################################
MACRO ( COPY_1_HEADER_AND_CREATE_TARGET
    source_dir name include_dir target_identifier )

  ADD_CUSTOM_TARGET(mmg${target_identifier}_${name}_header ALL
    DEPENDS
    ${source_dir}/${name}.h )

  COPY_HEADER (
    ${source_dir} ${name}.h
    ${include_dir} ${name}.h
    mmg${target_identifier}_${name}_header copy${target_identifier}_${name} )


ENDMACRO()

###############################################################################
#####
#####         Copy an automatically generated header file to another place
#####         and create the associated target
#####
###############################################################################
MACRO ( COPY_HEADERS_AND_CREATE_TARGET
    source_dir binary_dir include_dir target_identifier )

  COPY_1_HEADER_AND_CREATE_TARGET(
    ${source_dir} libmmg${target_identifier}
    ${include_dir} ${target_identifier})

  COPY_1_HEADER_AND_CREATE_TARGET(
    ${source_dir} mmg${target_identifier}_export
    ${include_dir} ${target_identifier})

  SET ( tgt_list
    copy${target_identifier}_libmmg${target_identifier}
    copy${target_identifier}_mmg${target_identifier}_export )

  if (PERL_FOUND)
    COPY_HEADER (
      ${binary_dir} libmmg${target_identifier}f.h
      ${include_dir} libmmg${target_identifier}f.h
      mmg${target_identifier}_fortran_header copy_libmmg${target_identifier}f )

    LIST ( APPEND tgt_list copy_libmmg${target_identifier}f )
  endif (PERL_FOUND)

  IF ( MMG_INSTALL_PRIVATE_HEADERS )
    COPY_1_HEADER_AND_CREATE_TARGET(
      ${source_dir} libmmg${target_identifier}_private
      ${include_dir} ${target_identifier})
    COPY_1_HEADER_AND_CREATE_TARGET(
      ${source_dir} mmg${target_identifier}externs_private
      ${include_dir} ${target_identifier})

    LIST ( APPEND tgt_list copy${target_identifier}_libmmg${target_identifier}_private
      copy${target_identifier}_mmg${target_identifier}externs_private
      )
  ENDIF()

  ADD_CUSTOM_TARGET (copy_${target_identifier}_headers ALL
    DEPENDS ${tgt_list} ${tgt_opt_list} copy_mmgcommon_headers)

ENDMACRO ( )

###############################################################################
#####
#####         Add a library to build and needed include dir, set its
#####         properties, add link dependencies and the install rule
#####
###############################################################################

MACRO ( ADD_AND_INSTALL_LIBRARY
    target_name target_type target_dependencies sources output_name )

  ADD_LIBRARY ( ${target_name} ${target_type} ${sources} )
  ADD_LIBRARY ( Mmg::${target_name} ALIAS ${target_name} )

  IF ( ${CMAKE_C_COMPILER_ID} STREQUAL "Clang" AND DEFINED CMAKE_C_COMPILER_VERSION )
    IF ( ${CMAKE_C_COMPILER_VERSION} VERSION_GREATER 10 )
      target_compile_options(${target_name} PRIVATE "-fcommon")
    ENDIF()
  ENDIF()

  IF (NOT WIN32 OR MINGW)
    ADD_DEPENDENCIES(${target_name} GenerateGitHashMmg)
  ENDIF()
  ADD_DEPENDENCIES( ${target_name} ${target_dependencies})

  IF ( CMAKE_VERSION VERSION_LESS 2.8.12 )
    INCLUDE_DIRECTORIES ( BEFORE
      ${MMGCOMMON_BINARY_DIR} ${MMGCOMMON_SOURCE_DIR} ${PROJECT_BINARY_DIR}/include ${PROJECT_BINARY_DIR})
  ELSE ( )
    target_include_directories( ${target_name} BEFORE PUBLIC
      $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/include/>
      $<BUILD_INTERFACE:${MMGCOMMON_SOURCE_DIR}>
      $<BUILD_INTERFACE:${MMGCOMMON_BINARY_DIR}>
      $<BUILD_INTERFACE:${MMG3D_SOURCE_DIR}>
      $<BUILD_INTERFACE:${MMGS_SOURCE_DIR}>
      $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
      )

  ENDIF ( )

  if ( SCOTCH_FOUND AND NOT USE_SCOTCH MATCHES OFF )
    message(STATUS "[mmg] add include scotch directories ${SCOTCH_INCLUDE_DIRS}")
    IF ( CMAKE_VERSION VERSION_LESS 2.8.12 )
      INCLUDE_DIRECTORIES ( AFTER ${SCOTCH_INCLUDE_DIRS} )
    ELSE ( )
      target_include_directories( ${target_name} PUBLIC ${SCOTCH_INCLUDE_DIRS} )
    endif()
  endif( )

  SET_TARGET_PROPERTIES ( ${target_name} PROPERTIES
    OUTPUT_NAME ${output_name}
    VERSION ${CMAKE_RELEASE_VERSION_MAJOR}.${CMAKE_RELEASE_VERSION_MINOR}.${CMAKE_RELEASE_VERSION_PATCH}
    SOVERSION ${CMAKE_RELEASE_VERSION_MAJOR} )


  SET_PROPERTY(TARGET ${target_name} PROPERTY C_STANDARD 99)

  TARGET_LINK_LIBRARIES ( ${target_name} PRIVATE ${LIBRARIES} )

  IF (NOT CMAKE_INSTALL_LIBDIR)
    SET(CMAKE_INSTALL_LIBDIR lib)
  ENDIF()


  SET ( MmgTargetsExported 1 )
  install(TARGETS ${target_name} EXPORT MmgTargets
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
    )
ENDMACRO ( )

###############################################################################
#####
#####         Add an executable to build and needed include dir, set its
#####         postfix, add link dependencies and the install rule
#####
###############################################################################

MACRO ( ADD_AND_INSTALL_EXECUTABLE
    exec_name target_dependencies lib_files main_file )

  # if one of the Mmg lib is built, use it instead of compiling whole sources
  IF ( NOT TARGET lib${exec_name}_a AND NOT TARGET lib${exec_name}_so
      AND NOT TARGET libmmg_a AND NOT TARGET libmmg_so )
    ADD_EXECUTABLE ( ${exec_name} ${lib_files} ${main_file} )
  ELSE ( )
    ADD_EXECUTABLE ( ${exec_name} ${main_file})

    SET_PROPERTY(TARGET ${exec_name} PROPERTY C_STANDARD 99)

    # link libraries in order of simplicities and depending on built targets
    IF ( TARGET  lib${exec_name}_a )
      TARGET_LINK_LIBRARIES(${exec_name} PRIVATE lib${exec_name}_a)
    ELSEIF ( TARGET  libmmg_a )
      TARGET_LINK_LIBRARIES(${exec_name} PRIVATE libmmg_a)
    ELSEIF ( TARGET lib${exec_name}_so )
      TARGET_LINK_LIBRARIES(${exec_name} PRIVATE lib${exec_name}_so)
    ELSE ( )
      TARGET_LINK_LIBRARIES(${exec_name} PRIVATE libmmg_so)
    ENDIF ( )

  ENDIF ( )

  IF (NOT WIN32 OR MINGW)
    ADD_DEPENDENCIES(${exec_name} GenerateGitHashMmg)
  endif()
  ADD_DEPENDENCIES(${exec_name} ${target_dependencies})

  IF ( WIN32 AND NOT MINGW AND SCOTCH_FOUND AND NOT USE_SCOTCH MATCHES OFF )
    my_add_link_flags ( ${exec_name} "/SAFESEH:NO")
  ENDIF ( )

 IF ( CMAKE_VERSION VERSION_LESS 2.8.12 )
   INCLUDE_DIRECTORIES ( BEFORE
     ${MMGCOMMON_BINARY_DIR} ${MMGCOMMON_SOURCE_DIR} ${PROJECT_BINARY_DIR}/include ${PROJECT_BINARY_DIR} )
   if ( SCOTCH_FOUND AND NOT USE_SCOTCH MATCHES OFF )
     message(STATUS "[mmg] add include scotch directories ${SCOTCH_INCLUDE_DIRS}")
     INCLUDE_DIRECTORIES ( AFTER ${SCOTCH_INCLUDE_DIRS} )
   ENDIF()

 ELSE ( )
   TARGET_INCLUDE_DIRECTORIES ( ${exec_name} BEFORE PUBLIC
     ${MMGCOMMON_BINARY_DIR} ${MMGCOMMON_SOURCE_DIR} ${PROJECT_BINARY_DIR}/include ${PROJECT_BINARY_DIR} )
   if ( SCOTCH_FOUND AND NOT USE_SCOTCH MATCHES OFF )
     message(STATUS "[mmg] add include scotch directories ${SCOTCH_INCLUDE_DIRS}")
     target_include_directories( ${exec_name} BEFORE PUBLIC ${SCOTCH_INCLUDE_DIRS} )
   ENDIF()

 ENDIF ( )

  TARGET_LINK_LIBRARIES ( ${exec_name} PRIVATE ${LIBRARIES}  )

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
#####         Add a library test
#####
###############################################################################

MACRO ( ADD_LIBRARY_TEST target_name main_path target_dependency lib_name lib_type )
  ADD_EXECUTABLE ( ${target_name} ${main_path} )
  ADD_DEPENDENCIES( ${target_name} ${target_dependency} )

  IF ( CMAKE_VERSION VERSION_LESS 2.8.12 )
    INCLUDE_DIRECTORIES ( BEFORE ${PROJECT_BINARY_DIR}/include )
  ELSE ( )
    TARGET_INCLUDE_DIRECTORIES ( ${target_name} BEFORE PUBLIC ${PROJECT_BINARY_DIR}/include )
  ENDIF ( )

  IF ( WIN32 AND ((NOT MINGW) AND SCOTCH_FOUND AND NOT USE_SCOTCH MATCHES OFF) )
    MY_ADD_LINK_FLAGS ( ${target_name} "/SAFESEH:NO" )
  ENDIF ( )

  IF ( "${lib_type}" STREQUAL "SHARED" )
    ADD_DEFINITIONS(-D${lib_name}_IMPORTS)
  ENDIF()

  TARGET_LINK_LIBRARIES ( ${target_name}  PRIVATE ${lib_name} )
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
