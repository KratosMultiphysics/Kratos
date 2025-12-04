set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosMedApplication ****")

################### PYBIND11
include(pybind11Tools)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}) # for FindMED
find_package(MED REQUIRED)

# MED uses HDF
# TODO check interoperability with the HDFapp, both in serial and parallel
if (NOT HDF5_FOUND)
    find_package(HDF5 REQUIRED COMPONENTS C)
endif()

include_directories(${MED_INCLUDE_DIR} ${HDF5_INCLUDE_DIRS})

include_directories( ${KRATOS_SOURCE_DIR}/kratos )

## Med Core sources
file(GLOB_RECURSE KRATOS_MED_APPLICATION_CORE
    ${CMAKE_CURRENT_SOURCE_DIR}/med_application.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_io/*.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/*.cpp
)

## Med testing sources
if(${KRATOS_BUILD_TESTING} MATCHES ON)
    file(GLOB_RECURSE KRATOS_MED_APPLICATION_TESTING_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/tests/*.cpp)
    kratos_add_gtests(TARGET KratosMedCore SOURCES "${KRATOS_MED_APPLICATION_TESTING_SOURCES}" USE_CUSTOM_MAIN)
endif(${KRATOS_BUILD_TESTING} MATCHES ON)

## Med python interface sources
file(GLOB_RECURSE KRATOS_MED_APPLICATION_PYTHON_INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/*.cpp)

add_library(KratosMedCore SHARED ${KRATOS_MED_APPLICATION_CORE})
target_link_libraries(KratosMedCore PUBLIC KratosCore ${MED_LIBRARIES} ${HDF5_C_LIBRARIES})
set_target_properties(KratosMedCore PROPERTIES COMPILE_DEFINITIONS "MED_APPLICATION=EXPORT,API")
target_include_directories(KratosMedCore SYSTEM PRIVATE hdf5::hdf5)


###############################################################
## define library Kratos which defines the basic python interface
pybind11_add_module(KratosMedApplication MODULE THIN_LTO ${KRATOS_MED_APPLICATION_PYTHON_INTERFACE})
target_link_libraries(KratosMedApplication PRIVATE KratosMedCore)
set_target_properties(KratosMedApplication PROPERTIES PREFIX "")

# changing the .dll suffix to .pyd (Windows)
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    set_target_properties(KratosMedApplication PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

# changing the .dylib suffix to .so (OS X)
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    set_target_properties(KratosMedApplication PROPERTIES SUFFIX .so)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

# Add to the KratosMultiphisics Python module
kratos_python_install(${INSTALL_PYTHON_USING_LINKS} ${CMAKE_CURRENT_SOURCE_DIR}/MedApplication.py KratosMultiphysics/MedApplication/__init__.py )
kratos_python_install(${INSTALL_PYTHON_USING_LINKS} ${CMAKE_CURRENT_SOURCE_DIR}/python_registry_lists.py KratosMultiphysics/MedApplication/python_registry_lists.py )

# Install python files
get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
kratos_python_install_directory(${INSTALL_PYTHON_USING_LINKS} ${CMAKE_CURRENT_SOURCE_DIR}/python_scripts KratosMultiphysics/${CURRENT_DIR_NAME} )

# Kratos Testing. Install everything except sources to ensure that reference and configuration files are copied.
if(${INSTALL_TESTING_FILES} MATCHES ON )
    get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
    install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/tests DESTINATION applications/${CURRENT_DIR_NAME}
        PATTERN "*.git" EXCLUDE
        PATTERN "*.c" EXCLUDE
        PATTERN "*.h" EXCLUDE
        PATTERN "*.cpp" EXCLUDE
        PATTERN "*.hpp" EXCLUDE
  )
endif(${INSTALL_TESTING_FILES} MATCHES ON)

# Fix Windows paths
file(TO_CMAKE_PATH "${HDF5_ROOT}/lib" HDF5_ROOT_DEP_LIB)
file(TO_CMAKE_PATH "${HDF5_ROOT}/bin" HDF5_ROOT_DEP_BIN)
file(TO_CMAKE_PATH "${MED_ROOT}/lib" MED_ROOT_DEP)

# Install targets
install(TARGETS KratosMedCore DESTINATION libs)
install(TARGETS KratosMedApplication DESTINATION libs)
if(WIN32)
    install(DIRECTORY "${HDF5_ROOT_DEP_LIB}/" DESTINATION libs FILES_MATCHING PATTERN "*.dll")
    install(DIRECTORY "${HDF5_ROOT_DEP_BIN}/" DESTINATION libs FILES_MATCHING PATTERN "*.dll")
    install(DIRECTORY "${MED_ROOT_DEP}/" DESTINATION libs FILES_MATCHING PATTERN "*.dll")
endif()
if(APPLE)
    install(DIRECTORY "${HDF5_ROOT_DEP_LIB}/" DESTINATION libs FILES_MATCHING PATTERN "*.dylib")
    install(DIRECTORY "${HDF5_ROOT_DEP_BIN}/" DESTINATION libs FILES_MATCHING PATTERN "*.dylib")
    install(DIRECTORY "${MED_ROOT_DEP}/" DESTINATION libs FILES_MATCHING PATTERN "*.dylib")
endif()
if(UNIX AND NOT APPLE)
    install(DIRECTORY "${HDF5_ROOT_DEP_LIB}/" DESTINATION libs FILES_MATCHING PATTERN "*.so*")
    install(DIRECTORY "${HDF5_ROOT_DEP_BIN}/" DESTINATION libs FILES_MATCHING PATTERN "*.so*")
    install(DIRECTORY "${MED_ROOT_DEP}/" DESTINATION libs FILES_MATCHING PATTERN "*.so*")
endif()

# Define custom targets
set(KRATOS_KERNEL "${KRATOS_KERNEL};KratosMedCore" PARENT_SCOPE)
set(KRATOS_PYTHON_INTERFACE "${KRATOS_PYTHON_INTERFACE};KratosMedApplication" PARENT_SCOPE)
