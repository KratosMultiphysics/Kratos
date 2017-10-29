set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosPfemApplication ****")

include_directories( ${CMAKE_SOURCE_DIR}/kratos )

## Set a default value for triangle and tetgen in case the user does not set it
## to avoid problems with the define linkage block. By default we will compile the lib
if(NOT DEFINED ${INCLUDE_TRIANGLE})
  set(INCLUDE_TRIANGLE ON)
  set(TRIANGLE_INSTALL_DIR ${CMAKE_SOURCE_DIR}/external_libraries/triangle/)
endif(NOT DEFINED ${INCLUDE_TRIANGLE})

if(NOT DEFINED ${INCLUDE_TETGEN})
  set(INCLUDE_TETGEN ON)
  #set(TETGEN_INSTALL_DIR ${CMAKE_SOURCE_DIR}/external_libraries/tetgen/tetgen1.5.0/)
  set(TETGEN_INSTALL_DIR ${CMAKE_CURRENT_SOURCE_DIR}/custom_external_libraries/tetgen/tetgen1.4.3/)
endif(NOT DEFINED ${INCLUDE_TETGEN})

if(${INCLUDE_TRIANGLE} MATCHES ON)
  add_definitions( -DINCLUDE_TRIANGLE )
  add_definitions( -DTRILIBRARY )
  add_definitions( -w )
  #perform platform specific defintions
  if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    add_definitions( -DCPU86 )
    add_definitions( -DNO_TIMER )
  else(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    add_definitions( -DLINUX )
  endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
  include_directories( ${TRIANGLE_INSTALL_DIR} )
  set( TRIANGLE_SOURCES ${TRIANGLE_INSTALL_DIR}/triangle.c )
  SET_SOURCE_FILES_PROPERTIES(${TRIANGLE_SOURCES} PROPERTIES COMPILE_FLAGS -O0)
  add_library(external_triangle STATIC ${TRIANGLE_INSTALL_DIR}/triangle.c)
  #install (FILES ${TRIANGLE_INSTALL_DIR}/triangle.h DESTINATION libs)
endif(${INCLUDE_TRIANGLE} MATCHES ON)

if(${INCLUDE_TETGEN} MATCHES ON)
  add_definitions( -DINCLUDE_TETGEN )
  add_definitions( -DTETLIBRARY )
  #add_definitions( -DTETGEN150 )
  add_definitions( -w )
  #perform platform specific defintions
  if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    add_definitions( -DCPU86 )
    add_definitions( -DNO_TIMER )
  else(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    add_definitions( -DLINUX )
  endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
  include_directories( ${TETGEN_INSTALL_DIR} )
  set( TETGEN_SOURCES ${TETGEN_INSTALL_DIR}/tetgen.cxx )
  set( UNOPTIMIZED_PREDICATES ${TETGEN_INSTALL_DIR}/predicates.cxx )
  SET_SOURCE_FILES_PROPERTIES(${UNOPTIMIZED_PREDICATES} PROPERTIES COMPILE_FLAGS -O0)
  SET_SOURCE_FILES_PROPERTIES(${TETGEN_SOURCES} PROPERTIES COMPILE_FLAGS -O2)
  add_library( external_tetgen STATIC ${TETGEN_SOURCES} ${UNOPTIMIZED_PREDICATES} )
  #install (FILES ${TETGEN_INSTALL_DIR}/tetgen.h DESTINATION libs)
endif(${INCLUDE_TETGEN} MATCHES ON)

#if(${INCLUDE_TETGEN} MATCHES ON)
#  add_definitions( -DINCLUDE_TETGEN )
#  add_definitions( -DTETLIBRARY )
#  add_definitions( -w )
#  #perform platform specific defintions
#  if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
#    add_definitions( -DCPU86 )
#    add_definitions( -DNO_TIMER )
#  else(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
#    add_definitions( -DLINUX )
#  endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
#  include_directories( ${TETGEN_INSTALL_DIR} )
#  link_directories( ${TETGEN_INSTALL_DIR} )
#endif(${INCLUDE_TETGEN} MATCHES ON)

##include directories for compilation of dependent libraries
#include_directories( ${CMAKE_CURRENT_SOURCE_DIR}/custom_external_libraries/triangle )
#include_directories( ${CMAKE_CURRENT_SOURCE_DIR}/custom_external_libraries/tetgen/tetgen1.4.3 )
##include_directories( ${CMAKE_CURRENT_SOURCE_DIR}/custom_external_libraries/msuite )
##include_directories( ${CMAKE_CURRENT_SOURCE_DIR}/custom_external_libraries/msuite/src/util )
##include_directories( ${CMAKE_CURRENT_SOURCE_DIR}/custom_external_libraries/msuite/src/malla )
##include_directories( ${CMAKE_CURRENT_SOURCE_DIR}/custom_external_libraries/tetMeshOpt )

##add subdirectories for compilation of dependent libraries
#add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/custom_external_libraries/triangle)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/custom_external_libraries/tetgen/tetgen1.4.3)
##add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/custom_external_libraries/msuite)
##add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/custom_external_libraries/tetMeshOpt)


## generate variables with the sources
set( KRATOS_PFEM_APPLICATION_SOURCES
	${CMAKE_CURRENT_SOURCE_DIR}/pfem_application.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/pfem_application_variables.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/pfem_python_application.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_processes_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_utilities_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_modelers_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_bounding_to_python.cpp

	##utilities
	${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/modeler_utilities.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/mesh_data_transfer_utilities.cpp

	##conditions
	${CMAKE_CURRENT_SOURCE_DIR}/custom_conditions/composite_condition.cpp

	##modelers
	${CMAKE_CURRENT_SOURCE_DIR}/custom_modelers/mesh_modeler.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_modelers/triangular_mesh_2D_modeler.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_modelers/tetrahedral_mesh_3D_modeler.cpp

)

###############################################################
## define library Kratos which defines the basic python interface
add_library(KratosPfemApplication SHARED ${KRATOS_PFEM_APPLICATION_SOURCES})

if(${MSVC})
  set(MATLIB "")
else(${MSVC})
  set(MATLIB "m")
endif(${MSVC})

## define the linkage of libraries
if( ${INCLUDE_TRIANGLE} MATCHES ON AND ${INCLUDE_TETGEN} MATCHES ON )
    set(LIBS KratosCore external_triangle external_tetgen ${MATLIB})
    #set(LIBS KratosCore external_triangle libextet.a ${MATLIB})
elseif(${INCLUDE_TRIANGLE} MATCHES ON)
    set(LIBS KratosCore external_triangle ${MATLIB})
elseif(${INCLUDE_TETGEN} MATCHES ON)
    set(LIBS KratosCore external_tetgen ${MATLIB})
    #set(LIBS KratosCore libextet.a ${MATLIB})
else(${INCLUDE_TRIANGLE} MATCHES ON AND ${INCLUDE_TETGEN} MATCHES ON)
    set(LIBS KratosCore)
endif(${INCLUDE_TRIANGLE} MATCHES ON AND ${INCLUDE_TETGEN} MATCHES ON)

target_link_libraries(KratosPfemApplication ${LIBS})
##target_link_libraries(KratosPfemApplication KratosCore triangle tetgen msuite)
##target_link_libraries(KratosPfemApplication KratosCore triangle tetgen msuite tetMeshOpt)

set_target_properties(KratosPfemApplication PROPERTIES PREFIX "")
set_target_properties(KratosPfemApplication PROPERTIES COMPILE_DEFINITIONS "PFEM_APPLICATION=EXPORT,API")

install(TARGETS KratosPfemApplication DESTINATION libs )


# get_property(inc_dirs DIRECTORY PROPERTY INCLUDE_DIRECTORIES)
# message("TestApplication subdir inc_dirs = ${inc_dirs}")

# changing the .dll suffix to .pyd
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
	set_target_properties(KratosPfemApplication PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

# changing the .dylib suffix to .so
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
	set_target_properties(KratosPfemApplication PROPERTIES SUFFIX .so)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

if(${INSTALL_PYTHON_FILES} MATCHES ON)
  get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
  install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/python_scripts DESTINATION applications/${CURRENT_DIR_NAME}  FILES_MATCHING PATTERN "*.py"  PATTERN ".svn" EXCLUDE)
endif(${INSTALL_PYTHON_FILES} MATCHES ON)

# get_property(inc_dirs DIRECTORY PROPERTY INCLUDE_DIRECTORIES)
# message("PfemApplication subdir inc_dirs = ${inc_dirs}")

# Add to the KratosMultiphisics Python module
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/PfemApplication.py" DESTINATION KratosMultiphysics )
