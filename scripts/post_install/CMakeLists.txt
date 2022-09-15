set(KRATOS_GENERATE_PYTHON_STUBS OFF CACHE BOOL "Generate hints for binary python packages")

if (${KRATOS_GENERATE_PYTHON_STUBS})
    install(CODE "message(STATUS \"Trying to install python stub files to enable python hinting support for cpp libraries in IDEs...\")")
    install(CODE "execute_process(COMMAND \"${PYTHON_EXECUTABLE}\" \"${CMAKE_CURRENT_SOURCE_DIR}/stub_generation.py\" \"${CMAKE_INSTALL_PREFIX}\" --quiet)")
endif()
