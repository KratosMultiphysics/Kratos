set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosExternalSolversApplication ****")

include_directories( ${CMAKE_SOURCE_DIR}/kratos )
include_directories( ${CMAKE_CURRENT_SOURCE_DIR}/custom_external_libraries/SuperLU_4.3 )
include_directories( ${CMAKE_CURRENT_SOURCE_DIR}/custom_external_libraries )

################### PYBIND11
include(pybind11Tools)

# ##add subdirectories for compilation of dependent libraries
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/custom_external_libraries/SuperLU_4.3)

if(${INCLUDE_PASTIX} MATCHES ON)
    add_definitions(-DINCLUDE_PASTIX)
    include_directories( ${PASTIX_INSTALL_DIR} )
    FIND_LIBRARY(PASTIX_LINK_LIB pastix  ${PASTIX_INSTALL_DIR})
    FIND_LIBRARY(SCOTCH_ERR_LIB scotcherr ${SCOTCH_INSTALL_DIR})
    FIND_LIBRARY(SCOTCH_LIB scotch  ${SCOTCH_INSTALL_DIR})
    message("pastix library found =" ${PASTIX_LINK_LIB})
    add_definitions( -DFORCE_NOMPI )
    add_library(kratos_pastix STATIC custom_external_libraries/solvePASTIX.c  )
    target_link_libraries(kratos_pastix  ${PASTIX_LINK_LIB} ${SCOTCH_LIB} ${SCOTCH_ERR_LIB} rt  ${BLAS_LIBRARIES} )

        install(FILES ${SCOTCH_LIB} ${SCOTCH_ERR_LIB}   ${BLAS_LIBRARIES} DESTINATION libs )
else(${INCLUDE_PASTIX} MATCHES ON)
    set(PASTIX_LINK_LIB "")
endif(${INCLUDE_PASTIX} MATCHES ON)

# if(${INCLUDE_FEAST} MATCHES ON)
#   if (NOT DEFINED FEAST_ROOT)
#     set(FEAST_ROOT ${CMAKE_CURRENT_SOURCE_DIR}/custom_external_libraries/FEAST/3.0 )
#   endif(NOT DEFINED FEAST_ROOT)
#   include_directories( ${FEAST_ROOT}/include )
#   link_directories( ${FEAST_ROOT}/lib/x64 )
#   add_definitions(-DINCLUDE_FEAST)
# endif(${INCLUDE_FEAST} MATCHES ON)
if(${INCLUDE_FEAST} MATCHES ON)
  include(CMakeAddFortranSubdirectory)
  cmake_add_fortran_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/custom_external_libraries/FEAST NO_EXTERNAL_INSTALL )
  if (NOT DEFINED FEAST_ROOT)
    set(FEAST_ROOT ${CMAKE_CURRENT_SOURCE_DIR}/custom_external_libraries/FEAST/3.0 )
  endif(NOT DEFINED FEAST_ROOT)
  include_directories( ${FEAST_ROOT}/include )
#   link_directories( ${FEAST_ROOT}/lib/x64 )
  add_definitions(-DINCLUDE_FEAST)
endif(${INCLUDE_FEAST} MATCHES ON)

## generate variables with the sources
set( KRATOS_EXTERNAL_SOLVERS_APPLICATION_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/externalsolvers_application.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_linear_solvers_to_python.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/externalsolvers_python_application.cpp
)
#
# ###############################################################
# ## define library Kratos which defines the basic python interface
# add_library(KratosExternalSolversApplication SHARED ${KRATOS_EXTERNAL_SOLVERS_APPLICATION_SOURCES})
pybind11_add_module(KratosExternalSolversApplication MODULE THIN_LTO ${KRATOS_EXTERNAL_SOLVERS_APPLICATION_SOURCES})


if(${INCLUDE_PASTIX} MATCHES ON)
        target_link_libraries(KratosExternalSolversApplication PRIVATE KratosCore superlu  kratos_pastix rt)
else(${INCLUDE_PASTIX} MATCHES ON)
        target_link_libraries(KratosExternalSolversApplication PRIVATE KratosCore superlu  )
endif(${INCLUDE_PASTIX} MATCHES ON)


if(${INCLUDE_FEAST} MATCHES ON)
  set(FEAST_LINK_LIBS feast ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} gfortran m)
  target_link_libraries(KratosExternalSolversApplication PRIVATE ${FEAST_LINK_LIBS})
endif(${INCLUDE_FEAST} MATCHES ON)

set_target_properties(KratosExternalSolversApplication PROPERTIES PREFIX "")

# changing the .dll suffix to .pyd
# if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
#     set_target_properties(KratosExternalSolversApplication PROPERTIES SUFFIX .pyd)
# endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
#
# # changing the .dylib suffix to .so
# if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
#     set_target_properties(KratosExternalSolversApplication PROPERTIES SUFFIX .so)
# endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

if(${INSTALL_PYTHON_FILES} MATCHES ON)
  get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
  install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/python_scripts DESTINATION applications/${CURRENT_DIR_NAME}  FILES_MATCHING PATTERN "*.py"  PATTERN ".svn" EXCLUDE)
endif(${INSTALL_PYTHON_FILES} MATCHES ON)

# Add Cotire
if(USE_COTIRE MATCHES ON)
    cotire(KratosExternalSolversApplication)
endif(USE_COTIRE MATCHES ON)

install(TARGETS KratosExternalSolversApplication DESTINATION libs )

# Add to the KratosMultiphisics Python module
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/ExternalSolversApplication.py" DESTINATION KratosMultiphysics )

#if(${INCLUDE_PASTIX} MATCHES ON)
#    install( TARGETS kratos_pastix DESTINATION libs)
#endif(${INCLUDE_PASTIX} MATCHES ON)
