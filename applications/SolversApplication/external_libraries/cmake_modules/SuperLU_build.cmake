message("**** compiling SuperLU ****")

set(superlu_cmd "cp")
set(superlu_linux "../cmake_modules/make.superlu")
set(superlu_inc "make.inc")
message(STATUS "superlu cmd: ${superlu_cmd} ${superlu_linux} ${superlu_inc}")
execute_process(COMMAND ${superlu_cmd} ${superlu_linux} ${superlu_inc}
  WORKING_DIRECTORY ${SUPERLU_DIR}
  RESULT_VARIABLE superlu_result
  OUTPUT_VARIABLE superlu_output)
#message(STATUS "superlu output[${superlu_result}]:\n${superlu_output}")

set(superlu_cmd "mkdir")
set(superlu_arg "-p")
set(superlu_lib "lib")
message(STATUS "superlu cmd: ${superlu_cmd} ${superlu_arg} ${superlu_lib}")
execute_process(COMMAND ${superlu_cmd} ${superlu_arg} ${superlu_lib}
  WORKING_DIRECTORY ${SUPERLU_DIR}
  RESULT_VARIABLE superlu_result
  OUTPUT_VARIABLE superlu_output)
#message(STATUS "superlu output[${superlu_result}]:\n${superlu_output}")

set(superlu_cmd "make")
set(superlu_omp "-j4")
message(STATUS "superlu cmd: ${superlu_cmd} ${superlu_omp}")
execute_process(COMMAND ${superlu_cmd} ${superlu_omp}
  WORKING_DIRECTORY ${SUPERLU_DIR}
  RESULT_VARIABLE superlu_result
  OUTPUT_VARIABLE superlu_output)
message(STATUS "superlu output[${superlu_result}]:\n${superlu_output}")


find_library(SUPERLU_LIBRARIES
  superlu_5.2.1
  "${SUPERLU_DIR}/lib"
  )
find_path(SUPERLU_INCLUDES
  slu_ddefs.h
  "${SUPERLU_DIR}/SRC"
  )

include_directories(${SUPERLU_INCLUDES})
message(STATUS "SuperLU lib found: ${SUPERLU_LIBRARIES}")
#add_library(${SUPERLU_LIBRARIES} STATIC IMPORTED)
#set_target_properties(${SUPERLU_LIBRARIES} PROPERTIES IMPORTED_LOCATION ${SUPERLU_DIR}/lib)

install(FILES ${SUPERLU_LIBRARIES} DESTINATION libs)
