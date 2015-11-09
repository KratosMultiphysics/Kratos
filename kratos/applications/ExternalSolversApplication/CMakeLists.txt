set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosExternalSolversApplication ****")

if(${INCLUDE_PASTIX} MATCHES ON)
	set(${EXCLUDE_ITSOL} ON)
	message("************************** ITSOL NOT COMPILED AS IT IS INCOMPATIBLE WITH PASTIX **************** ")
endif(${INCLUDE_PASTIX} MATCHES ON)




include_directories( ${CMAKE_SOURCE_DIR}/kratos )
include_directories( ${CMAKE_CURRENT_SOURCE_DIR}/custom_external_libraries/SuperLU_4.3 )
include_directories( ${CMAKE_CURRENT_SOURCE_DIR}/custom_external_libraries )

# ##add subdirectories for compilation of dependent libraries
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/custom_external_libraries/SuperLU_4.3)

if(${EXCLUDE_ITSOL} MATCHES ON)
    add_definitions(-DEXCLUDE_ITSOL)
else(${EXCLUDE_ITSOL} MATCHES ON)
    #for ITSOL
    add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/custom_external_libraries/ITSOL_2)
    include_directories( ${CMAKE_CURRENT_SOURCE_DIR}/custom_external_libraries/ITSOL_2 )
    include_directories( ${CMAKE_CURRENT_SOURCE_DIR}/custom_external_libraries/ITSOL_2/INC  )
endif(${EXCLUDE_ITSOL} MATCHES ON)

if(${EXCLUDE_AMGCL} MATCHES ON)
    add_definitions(-DEXCLUDE_AMGCL)
endif(${EXCLUDE_AMGCL} MATCHES ON)

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



## generate variables with the sources
set( KRATOS_EXTERNAL_SOLVERS_APPLICATION_SOURCES    
	${CMAKE_CURRENT_SOURCE_DIR}/externalsolvers_application.cpp 
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_linear_solvers_to_python.cpp 
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/externalsolvers_python_application.cpp 
)
# 
# ###############################################################
# ## define library Kratos which defines the basic python interface
add_library(KratosExternalSolversApplication SHARED ${KRATOS_EXTERNAL_SOLVERS_APPLICATION_SOURCES})

if(${EXCLUDE_ITSOL} MATCHES ON)
	if(${INCLUDE_PASTIX} MATCHES ON)
		target_link_libraries(KratosExternalSolversApplication KratosCore superlu kratos_pastix rt)
	else(${INCLUDE_PASTIX} MATCHES ON)
		target_link_libraries(KratosExternalSolversApplication KratosCore superlu )
	endif(${INCLUDE_PASTIX} MATCHES ON)
else(${EXCLUDE_ITSOL} MATCHES ON)
	if(${INCLUDE_PASTIX} MATCHES ON)
		target_link_libraries(KratosExternalSolversApplication KratosCore superlu itsol_solvers kratos_pastix rt)
	else(${INCLUDE_PASTIX} MATCHES ON)
		target_link_libraries(KratosExternalSolversApplication KratosCore superlu itsol_solvers )
	endif(${INCLUDE_PASTIX} MATCHES ON)
endif(${EXCLUDE_ITSOL} MATCHES ON)

set_target_properties(KratosExternalSolversApplication PROPERTIES PREFIX "")

# changing the .dll suffix to .pyd 
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
	set_target_properties(KratosExternalSolversApplication PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

# changing the .dylib suffix to .so 
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
	set_target_properties(KratosExternalSolversApplication PROPERTIES SUFFIX .so)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

if(${INSTALL_PYTHON_FILES} MATCHES ON)
  get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
  install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/python_scripts DESTINATION applications/${CURRENT_DIR_NAME}  FILES_MATCHING PATTERN "*.py"  PATTERN ".svn" EXCLUDE)
endif(${INSTALL_PYTHON_FILES} MATCHES ON)

install(TARGETS KratosExternalSolversApplication DESTINATION libs )

# Add to the KratosMultiphisics Python module
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/ExternalSolversApplication.py" DESTINATION KratosMultiphysics )

#if(${INCLUDE_PASTIX} MATCHES ON)
#	install( TARGETS kratos_pastix DESTINATION libs)
#endif(${INCLUDE_PASTIX} MATCHES ON)
