ADD_EXECUTABLE(libmmg3d_example0_a_oldAPI
  ${CMAKE_SOURCE_DIR}/libexamples/mmg3d/adaptation_example0/example0_a/main.c)

ADD_EXECUTABLE(libmmg3d_example0_b_oldAPI
  ${CMAKE_SOURCE_DIR}/libexamples/mmg3d/adaptation_example0/example0_b/main.c)

ADD_EXECUTABLE(libmmg3d_example1_oldAPI
  ${CMAKE_SOURCE_DIR}/libexamples/mmg3d/adaptation_example1/main.c)

ADD_EXECUTABLE(libmmg3d_example2_oldAPI
  ${CMAKE_SOURCE_DIR}/libexamples/mmg3d/adaptation_example2/main.c)

ADD_EXECUTABLE(libmmg3d_example4_oldAPI
  ${CMAKE_SOURCE_DIR}/libexamples/mmg3d/LagrangianMotion_example0/main.c)

IF ( WIN32 AND ((NOT MINGW) AND USE_SCOTCH) )
  my_add_link_flags(libmmg3d_example0_a_oldAPI "/SAFESEH:NO")
  my_add_link_flags(libmmg3d_example0_b_oldAPI "/SAFESEH:NO")
  my_add_link_flags(libmmg3d_example1_oldAPI "/SAFESEH:NO")
  my_add_link_flags(libmmg3d_example2_oldAPI "/SAFESEH:NO")
  my_add_link_flags(libmmg3d_example4_oldAPI "/SAFESEH:NO")
ENDIF ( )

IF ( LIBMMG3D_STATIC )

  TARGET_LINK_LIBRARIES(libmmg3d_example0_a_oldAPI ${PROJECT_NAME}3d_a)
  TARGET_LINK_LIBRARIES(libmmg3d_example0_b_oldAPI ${PROJECT_NAME}3d_a)
  TARGET_LINK_LIBRARIES(libmmg3d_example1_oldAPI   ${PROJECT_NAME}3d_a)
  TARGET_LINK_LIBRARIES(libmmg3d_example2_oldAPI   ${PROJECT_NAME}3d_a)
  TARGET_LINK_LIBRARIES(libmmg3d_example4_oldAPI   ${PROJECT_NAME}3d_a)

ELSEIF ( LIBMMG3D_SHARED )
  TARGET_LINK_LIBRARIES(libmmg3d_example0_a_oldAPI ${PROJECT_NAME}3d_so)
  TARGET_LINK_LIBRARIES(libmmg3d_example0_b_oldAPI ${PROJECT_NAME}3d_so)
  TARGET_LINK_LIBRARIES(libmmg3d_example1_oldAPI   ${PROJECT_NAME}3d_so)
  TARGET_LINK_LIBRARIES(libmmg3d_example2_oldAPI   ${PROJECT_NAME}3d_so)
  TARGET_LINK_LIBRARIES(libmmg3d_example4_oldAPI   ${PROJECT_NAME}3d_so)

ELSE ()
  MESSAGE(WARNING "You must activate the compilation of the static or"
    " shared ${PROJECT_NAME} library to compile this tests." )
ENDIF ()

INSTALL(TARGETS libmmg3d_example0_a_oldAPI RUNTIME DESTINATION bin )
INSTALL(TARGETS libmmg3d_example0_b_oldAPI RUNTIME DESTINATION bin )
INSTALL(TARGETS libmmg3d_example1_oldAPI   RUNTIME DESTINATION bin )
INSTALL(TARGETS libmmg3d_example2_oldAPI   RUNTIME DESTINATION bin )
INSTALL(TARGETS libmmg3d_example4_oldAPI   RUNTIME DESTINATION bin )

###############################################################################
#####
#####         Fortran Tests
#####
###############################################################################
IF (CMAKE_Fortran_COMPILER)
  ENABLE_LANGUAGE (Fortran)

  ADD_EXECUTABLE(libmmg3d_fortran_a_oldAPI
    ${CMAKE_SOURCE_DIR}/libexamples/mmg3d/adaptation_example0_fortran/example0_a/main.F90)

  ADD_EXECUTABLE(libmmg3d_fortran_b_oldAPI
    ${CMAKE_SOURCE_DIR}/libexamples/mmg3d/adaptation_example0_fortran/example0_b/main.F90)

  IF ( WIN32 AND ((NOT MINGW) AND USE_SCOTCH) )
    my_add_link_flags(libmmg3d_fortran_a_oldAPI "/SAFESEH:NO")
    my_add_link_flags(libmmg3d_fortran_b_oldAPI "/SAFESEH:NO")
  ENDIF ( )

  IF ( LIBMMG3D_STATIC )

    TARGET_LINK_LIBRARIES(libmmg3d_fortran_a_oldAPI  ${PROJECT_NAME}3d_a)
    TARGET_LINK_LIBRARIES(libmmg3d_fortran_b_oldAPI  ${PROJECT_NAME}3d_a)

  ELSE ( )

    TARGET_LINK_LIBRARIES(libmmg3d_fortran_a_oldAPI  ${PROJECT_NAME}3d_so)
    TARGET_LINK_LIBRARIES(libmmg3d_fortran_b_oldAPI  ${PROJECT_NAME}3d_so)

  ENDIF ()

  INSTALL(TARGETS libmmg3d_fortran_b_oldAPI  RUNTIME DESTINATION bin )
  INSTALL(TARGETS libmmg3d_fortran_a_oldAPI  RUNTIME DESTINATION bin )

ENDIF()
