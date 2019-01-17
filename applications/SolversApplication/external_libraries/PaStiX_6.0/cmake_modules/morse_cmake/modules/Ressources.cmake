###
#
# @copyright (c) 2009-2014 The University of Tennessee and The University
#                          of Tennessee Research Foundation.
#                          All rights reserved.
# @copyright (c) 2012-2014 Inria. All rights reserved.
# @copyright (c) 2012-2014 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria, Univ. Bordeaux. All rights reserved.
#
###
#
#  @file Ressources.cmake
#
#  @project MORSE
#  MORSE is a software package provided by:
#     Inria Bordeaux - Sud-Ouest,
#     Univ. of Tennessee,
#     King Abdullah Univesity of Science and Technology
#     Univ. of California Berkeley,
#     Univ. of Colorado Denver.
#
#  @version 0.9.0
#  @author Florent Pruvost
#  @date 13-05-2015
#
# Define variables for hardware ressources:
#  - NUMBER_OF_CPU.
#  - TODO: NUMBER_OF_CUDA.
###


if(NOT DEFINED PROCESSOR_COUNT)
  # Unknown:
  set(NUMBER_OF_CPU 0)

  # Linux:
  set(cpuinfo_file "/proc/cpuinfo")
  if(EXISTS "${cpuinfo_file}")
    file(STRINGS "${cpuinfo_file}" procs REGEX "^processor.: [0-9]+$")
    list(LENGTH procs NUMBER_OF_CPU)
  endif()

  # Mac:
  if(APPLE)
    find_program(cmd_sys_pro "system_profiler")
    if(cmd_sys_pro)
      execute_process(COMMAND ${cmd_sys_pro} SPHardwareDataType OUTPUT_VARIABLE info)
      string(REGEX REPLACE "^.*Total Number of Cores: ([0-9]+).*$" "\\1"
	NUMBER_OF_CPU "${info}")
    endif()
  endif()

  # Windows:
  if(WIN32)
    set(NUMBER_OF_CPU "$ENV{NUMBER_OF_PROCESSORS}")
  endif()
endif()
