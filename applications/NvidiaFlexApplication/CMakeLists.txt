set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring KratosNvidiaFlexApplication ****")

################### PYBIND11
include(pybind11Tools)

include_directories( ${CMAKE_SOURCE_DIR}/kratos )
include_directories( ${CMAKE_SOURCE_DIR}/applications/DEMApplication )
include_directories( ${NVIDIA_FLEX_ROOT_DIR}/core )
include_directories( ${NVIDIA_FLEX_ROOT_DIR}/include )

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fpermissive  ")

## generate variables with the sources
set(KRATOS_NVIDIAFLEX_APPLICATION_CORE
    ${CMAKE_CURRENT_SOURCE_DIR}/NvidiaFlex_application.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_utilities/flex_wrapper.cpp
)

set(KRATOS_NVIDIAFLEX_APPLICATION_PYTHON_INTERFACE
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/NvidiaFlex_python_application.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/custom_python/add_custom_utilities_to_python.cpp
)

set(NVIDIA_FLEX_LIBRARY
    ${NVIDIA_FLEX_LIBRARY_DIR}/NvFlexReleaseCUDA_x64.a
    ${NVIDIA_FLEX_LIBRARY_DIR}/NvFlexExtReleaseCUDA_x64.a
    ${NVIDIA_FLEX_LIBRARY_DIR}/NvFlexDeviceRelease_x64.a
)

add_library(KratosNvidiaFlexCore SHARED ${KRATOS_NVIDIAFLEX_APPLICATION_CORE})
target_link_libraries(KratosNvidiaFlexCore PUBLIC KratosCore KratosDEMCore ${NVIDIA_FLEX_LIBRARY} ${CUDA_LIBRARIES})
set_target_properties(KratosNvidiaFlexCore PROPERTIES COMPILE_DEFINITIONS "NVIDIAFLEX_APPLICATION=EXPORT,API")

###############################################################
## define library Kratos which defines the basic python interface
pybind11_add_module(KratosNvidiaFlexApplication MODULE THIN_LTO ${KRATOS_NVIDIAFLEX_APPLICATION_PYTHON_INTERFACE})
target_link_libraries(KratosNvidiaFlexApplication PUBLIC KratosNvidiaFlexCore)
set_target_properties(KratosNvidiaFlexApplication PROPERTIES PREFIX "")


if(${ACTIVATE_DEBUG_MACRO} MATCHES ON) #MSI: Flag defined for debug Macro
    add_definitions(-DDEBUG_MACRO)
endif(${ACTIVATE_DEBUG_MACRO} MATCHES ON)

# changing the .dll suffix to .pyd
if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    set_target_properties(KratosNvidiaFlexApplication PROPERTIES SUFFIX .pyd)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

# changing the .dylib suffix to .so
if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
	set_target_properties(KratosNvidiaFlexApplication PROPERTIES SUFFIX .so)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

if(${INSTALL_TESTING_FILES} MATCHES ON)
  get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
  install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/tests DESTINATION applications/${CURRENT_DIR_NAME} FILES_MATCHING PATTERN "*.py" PATTERN  "*.json" PATTERN "*.mdpa" PATTERN ".svn" EXCLUDE)
endif(${INSTALL_TESTING_FILES} MATCHES ON)

if(${INSTALL_PYTHON_FILES} MATCHES ON)
  get_filename_component (CURRENT_DIR_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
  install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/python_scripts DESTINATION applications/${CURRENT_DIR_NAME}  FILES_MATCHING PATTERN "*.py"  PATTERN ".svn" EXCLUDE)
endif(${INSTALL_PYTHON_FILES} MATCHES ON)


# get_property(inc_dirs DIRECTORY PROPERTY INCLUDE_DIRECTORIES)
# message("KratosNvidiaFlexApplication subdir inc_dirs = ${inc_dirs}")
# Add Cotire
if(USE_COTIRE MATCHES ON)
  cotire(KratosNvidiaFlexCore)
  cotire(KratosNvidiaFlexApplication)
endif(USE_COTIRE MATCHES ON)

install(TARGETS KratosNvidiaFlexCore DESTINATION libs )
install(TARGETS KratosNvidiaFlexApplication DESTINATION libs )

# Add to the KratosMultiphysics Python module
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/NvidiaFlexApplication.py" DESTINATION KratosMultiphysics )
