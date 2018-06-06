set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosPfemApplication ****")

include_directories( ${CMAKE_SOURCE_DIR}/kratos )

################### PYBIND11
include(pybind11Tools)

## Set a default value for triangle and tetgen in case the user does not set it
## to avoid problems with the define linkage block. By default we will compile the lib
if(NOT DEFINED ${INCLUDE_TRIANGLE})
  set(INCLUDE_TRIANGLE ON)
  set(TRIANGLE_INSTALL_DIR ${CMAKE_SOURCE_DIR}/external_libraries/triangle/)
endif(NOT DEFINED ${INCLUDE_TRIANGLE})

if(NOT DEFINED ${INCLUDE_TETGEN})
  set(INCLUDE_TETGEN ON)
  set(TETGEN_INSTALL_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external_modules/tetgen/)
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
  add_definitions( -DTETGEN150 )
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
##include_directories( ${CMAKE_CURRENT_SOURCE_DIR}/custom_external_libraries/msuite )
##include_directories( ${CMAKE_CURRENT_SOURCE_DIR}/custom_external_libraries/msuite/src/util )
##include_directories( ${CMAKE_CURRENT_SOURCE_DIR}/custom_external_libraries/msuite/src/malla )
##include_directories( ${CMAKE_CURRENT_SOURCE_DIR}/custom_external_libraries/tetMeshOpt )

##add subdirectories for compilation of dependent libraries
##add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/custom_external_libraries/msuite)
##add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/custom_external_libraries/tetMeshOpt)


## Generate variables with the sources
set( KRATOS_PFEM_APPLICATION_CORE
	${CMAKE_CURRENT_SOURCE_DIR}/pfem_application.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/pfem_application_variables.cpp

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

## Generate variables with the sources
set( KRATOS_PFEM_APPLICATION_PYTHON_INTERFACE
        # Custom python
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/pfem_python_application.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_processes_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_utilities_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_modelers_to_python.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_bounding_to_python.cpp
)

## matplot lib
if(${MSVC})
  set(MATLIB "")
else(${MSVC})
  set(MATLIB "m")
endif(${MSVC})

## define the linkage of libraries
if( ${INCLUDE_TRIANGLE} MATCHES ON AND ${INCLUDE_TETGEN} MATCHES ON )
    set(LIBS external_triangle external_tetgen ${MATLIB})
    #set(LIBS external_triangle libextet.a ${MATLIB})
elseif(${INCLUDE_TRIANGLE} MATCHES ON)
    set(LIBS external_triangle ${MATLIB})
elseif(${INCLUDE_TETGEN} MATCHES ON)
    set(LIBS external_tetgen ${MATLIB})
    #set(LIBS libextet.a ${MATLIB})
else(${INCLUDE_TRIANGLE} MATCHES ON AND ${INCLUDE_TETGEN} MATCHES ON)
    set(LIBS)
endif(${INCLUDE_TRIANGLE} MATCHES ON AND ${INCLUDE_TETGEN} MATCHES ON)

add_library(KratosPfemCore SHARED ${KRATOS_PFEM_APPLICATION_CORE})
target_link_libraries(KratosPfemCore PUBLIC KratosCore ${LIBS})
set_target_properties(KratosPfemCore PROPERTIES COMPILE_DEFINITIONS "PFEM_APPLICATION=EXPORT,API")

###############################################################
## define library Kratos which defines the basic python interface
pybind11_add_module(KratosPfemApplication MODULE THIN_LTO ${KRATOS_PFEM_APPLICATION_PYTHON_INTERFACE})
target_link_libraries(KratosPfemApplication PRIVATE KratosPfemCore)
set_target_properties(KratosPfemApplication PROPERTIES PREFIX "")

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

if(USE_COTIRE MATCHES ON)
    cotire(KratosPfemCore)
    cotire(KratosPfemApplication)
endif(USE_COTIRE MATCHES ON)

install(TARGETS KratosPfemCore DESTINATION libs )
install(TARGETS KratosPfemApplication DESTINATION libs )

# Add to the KratosMultiphisics Python module
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/PfemApplication.py" DESTINATION KratosMultiphysics )
