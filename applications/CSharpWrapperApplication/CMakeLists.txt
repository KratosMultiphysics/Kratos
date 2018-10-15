set(CMAKE_INCLUDE_CURRENT_DIR ON)

message("**** configuring CSharpWrapperApplication ****")

include_directories( ${CMAKE_SOURCE_DIR}/kratos )
include_directories( ${CMAKE_SOURCE_DIR}/applications/StructuralMechanicsApplication )

set( KRATOS_WRAPPER_SOURCES
	${CMAKE_CURRENT_SOURCE_DIR}/external_interface.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/external_interface_bindings.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/mesh_converter.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/vector3.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/kratos_internals.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/id_translator.cpp
)

add_library(KratosWrapper SHARED ${KRATOS_WRAPPER_SOURCES})
target_link_libraries(KratosWrapper PUBLIC KratosCore KratosStructuralMechanicsCore)


install(TARGETS KratosWrapper DESTINATION libs)

