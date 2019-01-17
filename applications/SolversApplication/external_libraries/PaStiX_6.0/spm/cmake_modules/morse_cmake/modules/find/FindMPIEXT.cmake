# - Find MPI EXTENDED for MORSE projects:
#
# This module allows to find MPI libraries by calling the official FindMPI module
# and adds an additional variable indicating the level of thread supported:
#  MPI_THREAD_SUPPORTED_LEVEL - MPI_THREAD_MULTIPLE, MPI_THREAD_SERIALIZED, MPI_THREAD_FUNNELED, or MPI_THREAD_SINGLE
#  MPIEXT_FOUND        - if a MPI has been found

#=============================================================================
# Copyright 2012-2013 Inria
# Copyright 2012-2013 Emmanuel Agullo
# Copyright 2012-2013 Mathieu Faverge
# Copyright 2012      Cedric Castagnede
# Copyright 2013      Florent Pruvost
#
# Distributed under the OSI-approved BSD License (the "License");
# see accompanying file MORSE-Copyright.txt for details.
#
# This software is distributed WITHOUT ANY WARRANTY; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the License for more information.
#=============================================================================
# (To distribute this file outside of Morse, substitute the full
#  License text for the above reference.)

# add a cache variable to let the user specify the BLAS vendor
if (NOT MPI_FOUND)
  if(MPIEXT_FIND_REQUIRED)
    find_package(MPI REQUIRED)
  else()
    find_package(MPI)
  endif()
endif ()
set(MPIEXT_FOUND ${MPI_FOUND})

if (MPI_FOUND)

  set(CMAKE_REQUIRED_INCLUDES_SAVE ${CMAKE_REQUIRED_INCLUDES})
  set(CMAKE_REQUIRED_LIBRARIES_SAVE ${CMAKE_REQUIRED_LIBRARIES})
  set(CMAKE_REQUIRED_FLAGS_SAVE ${CMAKE_REQUIRED_FLAGS})
  # Set flags for building test program
  set(CMAKE_REQUIRED_INCLUDES ${MPI_INCLUDE_PATH})
  set(CMAKE_REQUIRED_LIBRARIES ${MPI_LIBRARIES})
  set(CMAKE_REQUIRED_FLAGS ${MPI_COMPILE_FLAGS})

  include(CheckCSourceRuns)
  set(MPI_C_TEST_MULTIPLE_SOURCE "
#include <stdio.h>
#include <mpi.h>

int main(int argc, char **argv) {
  int provided;

  MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);

  if (MPI_THREAD_MULTIPLE > provided) {
    fprintf(stderr,\"MPI implementation is not supporting MPI_THREAD_MULTIPLE\");
    MPI_Finalize();
    return 1000+MPI_THREAD_MULTIPLE-provided;
  }
  MPI_Finalize();
  return 0;
}
")

  set(MPI_C_TEST_SERIALIZED_SOURCE "
#include <stdio.h>
#include <mpi.h>

int main(int argc, char **argv) {
  int provided;

  MPI_Init_thread(&argc, &argv, MPI_THREAD_SERIALIZED, &provided);

  if (MPI_THREAD_SERIALIZED > provided) {
    printf(\"MPI implementation is not supporting MPI_THREAD_SERIALIZED\");
    MPI_Finalize();
    return 1;
  }
  MPI_Finalize();
  return 0;
}
")

  set(MPI_C_TEST_FUNNELED_SOURCE "
#include <stdio.h>
#include <mpi.h>

int main(int argc, char **argv) {
  int provided;

  MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);

  if (MPI_THREAD_FUNNELED > provided) {
    printf(\"MPI implementation is not supporting MPI_THREAD_FUNNELED\");
    MPI_Finalize();
    return 1;
  }
  MPI_Finalize();
  return 0;
}
")

  check_c_source_runs("${MPI_C_TEST_MULTIPLE_SOURCE}" MULTIPLE_TEST_RUNS)
  if(NOT MULTIPLE_TEST_RUNS)
    check_c_source_runs("${MPI_C_TEST_SERIALIZED_SOURCE}" SERIALIZED_TEST_RUNS)
    if(NOT SERIALIZED_TEST_RUNS)
      check_c_source_runs("${MPI_C_TEST_FUNNELED_SOURCE}" FUNNELED_TEST_RUNS)
      if(NOT FUNNELED_TEST_RUNS)
        set(MPI_THREAD_SUPPORTED_LEVEL "MPI_THREAD_SINGLE")
      else(NOT FUNNELED_TEST_RUNS)
        set(MPI_THREAD_SUPPORTED_LEVEL "MPI_THREAD_FUNNELED")
      endif(NOT FUNNELED_TEST_RUNS)
    else(NOT SERIALIZED_TEST_RUNS)
      set(MPI_THREAD_SUPPORTED_LEVEL "MPI_THREAD_SERIALIZED")
    endif(NOT SERIALIZED_TEST_RUNS)
  else(NOT MULTIPLE_TEST_RUNS)
    set(MPI_THREAD_SUPPORTED_LEVEL "MPI_THREAD_MULTIPLE")
  endif(NOT MULTIPLE_TEST_RUNS)
  # Set flags for building test program
  set(CMAKE_REQUIRED_INCLUDES ${CMAKE_REQUIRED_INCLUDES_SAVE})
  set(CMAKE_REQUIRED_LIBRARIES ${CMAKE_REQUIRED_LIBRARIES_SAVE})
  set(CMAKE_REQUIRED_FLAGS ${CMAKE_REQUIRED_FLAGS_SAVE})

endif(MPI_FOUND)
