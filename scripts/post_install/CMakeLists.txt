set(KRATOS_GENERATE_PYTHON_STUBS OFF CACHE BOOL "Generate hints for binary python packages")

if (${KRATOS_GENERATE_PYTHON_STUBS})
    install(CODE "message(STATUS \"Trying to install python stub files to enable python hinting support for cpp libraries in IDEs...\")")
<<<<<<< HEAD
    install(CODE "execute_process(COMMAND \"${PYTHON_EXECUTABLE}\" \"${CMAKE_CURRENT_SOURCE_DIR}/stub_generation.py\" \"${CMAKE_INSTALL_PREFIX}\" --quiet)")
=======
    cmake_path(CONVERT ${PYTHON_EXECUTABLE} TO_CMAKE_PATH_LIST CONVERTED_PY_EXEC)
    install(CODE "execute_process(COMMAND \"${CONVERTED_PY_EXEC}\" \"${CMAKE_CURRENT_SOURCE_DIR}/stub_generation.py\" \"${CMAKE_INSTALL_PREFIX}\" --quiet)")
>>>>>>> 21c387f4469e81694616ffcfba50ef4788e0fb2a
endif()
