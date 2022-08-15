set(KRATOS_GENERATE_PYTHON_STUBS ON CACHE BOOL "Generate hints for binary python packages")

if (${CMAKE_SYSTEM_NAME} MATCHES "Linux") # todo: test the stub generation on other platforms too
    if (${KRATOS_GENERATE_PYTHON_STUBS})
        if (NOT DEFINED PYTHON_EXECUTABLE)
            set(PYTHON_EXECUTABLE "python")
        endif()

        install(CODE "message(STATUS \"Trying to install python stub files to enable python hinting support for cpp libraries in IDEs...\")")
        install(CODE "execute_process(COMMAND \"${PYTHON_EXECUTABLE}\" \"${CMAKE_CURRENT_SOURCE_DIR}/stub_generation.py\" \"${CMAKE_INSTALL_PREFIX}\" --quiet)")
    endif()
endif()