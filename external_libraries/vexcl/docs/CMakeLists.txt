cmake_minimum_required(VERSION 3.1)
project(docs)

find_package(Doxygen QUIET)

if (DOXYGEN_FOUND)
    configure_file(
	${CMAKE_CURRENT_SOURCE_DIR}/conf.py
	${CMAKE_CURRENT_BINARY_DIR}/conf.py
        COPYONLY
	)

    add_custom_target(doxygen
        COMMAND ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
	COMMENT "Build Doxygen-generated documentation"
	)

    add_custom_target(doc
        COMMAND sphinx-build -b html -c ${CMAKE_CURRENT_SOURCE_DIR} . ${CMAKE_CURRENT_BINARY_DIR}/html
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        DEPENDS doxygen
        COMMENT "Generate HTML documentation"
        )
endif()
