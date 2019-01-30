# Tetgen automatic download
# Variables:  TETGEN_DIR

find_path(TETGEN_FOUND NAMES tetgen.h HINTS ${TETGEN_DIR})

if(TETGEN_FOUND)
  message(STATUS "Tetgen found in ${TETGEN_DIR}")
else(TETGEN_FOUND)
  message(STATUS "TETGEN_DIR not defined: ${TETGEN_DIR}")

  # Default version
  file(TO_NATIVE_PATH "${EXTERNAL_LIBRARIES_DIR}" tetgen_install_dir)

  # Default TETGEN_DIR
  file(TO_NATIVE_PATH "${tetgen_install_dir}/tetgen" DEFAULT_TETGEN_DIR)
  set(TETGEN_DIR ${DEFAULT_TETGEN_DIR})
  message(STATUS "set TETGEN_DIR: ${TETGEN_DIR}")

  set(tetgen_file_name "tetgen-1.5.0-delaunay")
  set(tetgen_packed_file "${tetgen_file_name}.zip")
  file(TO_NATIVE_PATH "${tetgen_install_dir}/${tetgen_packed_file}" tetgen_packed_dir)

  if(NOT EXISTS ${tetgen_packed_dir})
    message(STATUS "${tetgen_packed_file} not found in ${tetgen_install_dir}")
    message(STATUS "Downloading ${tetgen_packed_file} from https://github.com/PFEM/tetgen-1.5.0/archive/delaunay.zip to ${tetgen_install_dir} ...")
    file(DOWNLOAD https://github.com/PFEM/tetgen-1.5.0/archive/delaunay.zip ${tetgen_packed_dir} SHOW_PROGRESS EXPECTED_MD5 7f5b348b1edac952700d5635fa150faa)
    message(STATUS "${tetgen_packed_file} downloaded")
  endif(NOT EXISTS ${tetgen_packed_dir})

  find_path(TETGEN_FOUND NAMES tetgen.h HINTS ${TETGEN_DIR})

  if((NOT EXISTS ${TETGEN_DIR}) OR (NOT TETGEN_FOUND))
    message(STATUS "Unpacking tetgen-1.5.0-delaunay.zip in ${tetgen_install_dir} ...")
    execute_process(COMMAND ${CMAKE_COMMAND} -E tar xzf ${tetgen_packed_file}
      WORKING_DIRECTORY ${tetgen_install_dir})
    execute_process(COMMAND ${CMAKE_COMMAND} -E copy_directory ${tetgen_install_dir}/${tetgen_file_name} ${TETGEN_DIR}
      WORKING_DIRECTORY ${tetgen_install_dir})
    execute_process(COMMAND ${CMAKE_COMMAND} -E remove_directory ${tetgen_install_dir}/${tetgen_file_name}
      WORKING_DIRECTORY ${tetgen_install_dir})
  endif((NOT EXISTS ${TETGEN_DIR}) OR (NOT TETGEN_FOUND))
endif(TETGEN_FOUND)
