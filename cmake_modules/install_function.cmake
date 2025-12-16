# This functin installs a symbolic link
macro(install_symlink filepath sympath symfile)
    install(DIRECTORY DESTINATION ${sympath})
    install(CODE [[
        execute_process(
            COMMAND "${CMAKE_COMMAND}" -E create_symlink 
            "${filepath}" 
            "${sympath}/${symfile}"
        )
    ]])
endmacro(install_symlink)

# This function install a symlinc to an entire directory
macro(kratos_python_install_directory INSTALL_PYTHON_USING_LINKS origin_dir destination_dir)
    install(DIRECTORY DESTINATION ${CMAKE_INSTALL_PREFIX}/${destination_dir})
    file(GLOB_RECURSE SRC_FILES RELATIVE ${origin_dir} ${origin_dir}/*.py )
    if(${INSTALL_PYTHON_USING_LINKS} MATCHES ON )
        foreach(f ${SRC_FILES})
            get_filename_component(barename ${f} NAME)
            get_filename_component(relpath ${f} DIRECTORY)
            install_symlink(${origin_dir}/${f} ${CMAKE_INSTALL_PREFIX}/${destination_dir}/${relpath} ${barename} )
        endforeach(f ${SRC_FILES})
    else(${INSTALL_PYTHON_USING_LINKS} MATCHES ON )
        foreach(f ${SRC_FILES})
            get_filename_component(relpath ${f} DIRECTORY)
            install(FILES ${origin_dir}/${f} DESTINATION ${CMAKE_INSTALL_PREFIX}/${destination_dir}/${relpath} )
        endforeach(f ${SRC_FILES})
    endif(${INSTALL_PYTHON_USING_LINKS} MATCHES ON )
endmacro(kratos_python_install_directory)

# This function install a symlinc to an entire directory
macro(kratos_python_install_tests)
    set(options INSTALL_PYTHON_USING_LINKS)
    set(oneValueArgs DIRECTORY DESTINATION)
    set(multiValueArgs EXCLUSION_PATTERNS)

    cmake_parse_arguments(TEST "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

    # Install the destination directory (required for the file structure)
    install(DIRECTORY DESTINATION ${CMAKE_INSTALL_PREFIX}/${TEST_DESTINATION})

    # Recursively find all files in the origin directory
    file(GLOB_RECURSE SRC_FILES RELATIVE ${TEST_DIRECTORY} ${TEST_DIRECTORY}/* )

    # Filter the file list based on the exclusion patterns
    foreach(pattern IN LISTS TEST_EXCLUSION_PATTERNS)
        list(FILTER SRC_FILES EXCLUDE REGEX "${pattern}")
    endforeach()

    # Handle installation based on symlink preference
    if(${INSTALL_PYTHON_USING_LINKS} MATCHES ON )
        foreach(f ${SRC_FILES})
            get_filename_component(barename ${f} NAME)
            get_filename_component(relpath ${f} DIRECTORY)
            install_symlink(${TEST_DIRECTORY}/${f} ${CMAKE_INSTALL_PREFIX}/${TEST_DESTINATION}/tests/${relpath} ${barename} )
        endforeach(f ${SRC_FILES})
    else()
        foreach(f ${SRC_FILES})
            get_filename_component(relpath ${f} DIRECTORY)
            message("(D)=======> ${TEST_DIRECTORY}/${f} ==> ${CMAKE_INSTALL_PREFIX}/${TEST_DESTINATION}/${relpath}")
            install(FILES ${TEST_DIRECTORY}/${f} DESTINATION ${CMAKE_INSTALL_PREFIX}/${TEST_DESTINATION}/tests/${relpath} )
        endforeach(f ${SRC_FILES})
    endif()
endmacro(kratos_python_install_tests)

macro(kratos_python_install INSTALL_PYTHON_USING_LINKS origin_file destination_file)
    get_filename_component(destination_dir ${destination_file} DIRECTORY )
    get_filename_component(destination_filename ${destination_file} NAME )
    install(DIRECTORY DESTINATION ${CMAKE_INSTALL_PREFIX}/${destination_dir})
    if(${INSTALL_PYTHON_USING_LINKS} MATCHES ON )
        install_symlink(${origin_file} ${CMAKE_INSTALL_PREFIX} ${destination_file} )
    else(${INSTALL_PYTHON_USING_LINKS} MATCHES ON )
        install(FILES ${origin_file} DESTINATION ${CMAKE_INSTALL_PREFIX}/${destination_dir} RENAME ${destination_filename})
    endif(${INSTALL_PYTHON_USING_LINKS} MATCHES ON )
endmacro(kratos_python_install)
