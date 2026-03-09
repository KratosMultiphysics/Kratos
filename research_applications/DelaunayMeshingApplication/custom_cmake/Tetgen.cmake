# Tetgen automatic download
# Variables:  TETGEN_DIR

find_path(TETGEN_FOUND NAMES tetgen.h HINTS ${TETGEN_DIR})

if(TETGEN_FOUND)
  message(STATUS "Tetgen found in ${TETGEN_DIR}")
else(TETGEN_FOUND)
  message(STATUS "TETGEN_DIR not defined: ${TETGEN_DIR}")

  # Default version
  set(tetgen_install_dir "${EXTERNAL_LIBRARIES_DIR}")

  # Default TETGEN_DIR
  set(TETGEN_DIR "${tetgen_install_dir}/tetgen")
  message(STATUS "set TETGEN_DIR: ${TETGEN_DIR}")

  set(tetgen_file_name "tetgen-1.5.0-kratos")
  set(tetgen_packed_file "${tetgen_file_name}.zip")
  set(tetgen_packed_dir "${tetgen_install_dir}/${tetgen_packed_file}")

  if(NOT EXISTS ${tetgen_packed_dir})
    message(STATUS "${tetgen_packed_file} not found in ${tetgen_install_dir}")
    message(STATUS "Downloading ${tetgen_packed_file} from https://github.com/PFEM/tetgen-1.5.0/archive/kratos.zip to ${tetgen_install_dir} ...")
    file(DOWNLOAD https://github.com/PFEM/tetgen-1.5.0/archive/kratos.zip ${tetgen_packed_dir} SHOW_PROGRESS EXPECTED_MD5 9038fbda461f087cd9d233dfa69e0816)
    message(STATUS "${tetgen_packed_file} downloaded")
  endif(NOT EXISTS ${tetgen_packed_dir})

  find_path(TETGEN_FOUND NAMES tetgen.h HINTS ${TETGEN_DIR})

  if((NOT EXISTS ${TETGEN_DIR}) OR (NOT TETGEN_FOUND))
    message(STATUS "Unpacking ${tetgen_packed_file} in ${tetgen_install_dir} ...")
    execute_process(COMMAND ${CMAKE_COMMAND} -E tar xzf ${tetgen_packed_file}
      WORKING_DIRECTORY ${tetgen_install_dir})
    execute_process(COMMAND ${CMAKE_COMMAND} -E copy_directory ${tetgen_install_dir}/${tetgen_file_name} ${TETGEN_DIR}
      WORKING_DIRECTORY ${tetgen_install_dir})
    execute_process(COMMAND ${CMAKE_COMMAND} -E remove_directory ${tetgen_install_dir}/${tetgen_file_name}
      WORKING_DIRECTORY ${tetgen_install_dir})
  endif((NOT EXISTS ${TETGEN_DIR}) OR (NOT TETGEN_FOUND))
endif(TETGEN_FOUND)
