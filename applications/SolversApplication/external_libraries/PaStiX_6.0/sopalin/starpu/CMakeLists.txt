###
#
#  @copyright 2013-2018 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
#                       Univ. Bordeaux. All rights reserved.
#
#  @version 6.0.1
#  @author Mathieu Faverge
#  @date 2018-07-16
#
###
include(RulesPrecisions)
include(RulesJDF)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${CMAKE_CURRENT_BINARY_DIR})

# reset variables
set(generated_sources "")
set(generated_jdf     "")
set(generated_headers "")

set(EXTRA_SOURCES
  starpu.c
  starpu_sparse_matrix.c
  starpu_dense_matrix.c
  pastix_starpu_model.c
)

### Generate the headers in all precisions
set(HEADERS
  pastix_zstarpu.h
)

precisions_rules_py(generated_headers
  "${HEADERS}"
  PRECISIONS "s;d;c;z")

set(starpu_headers
  ${generated_headers}
  )

add_custom_target(starpu_headers_tgt
  DEPENDS ${starpu_headers} pastix_starpu.h )

### Generate the pastix_parsec wrappers for all required precisions
set(SOURCES
  codelet_zgemmsp.c
  codelet_ztrsmsp.c
  codelet_zscalo.c
  codelet_zgetrfsp.c
  codelet_zhetrfsp.c
  codelet_zpotrfsp.c
  codelet_zpxtrfsp.c
  codelet_zsytrfsp.c
  codelet_solve_zdiag.c
  codelet_solve_ztrsm.c
  codelet_solve_zgemm.c
  starpu_zdiag.c
  starpu_zgetrf.c
  starpu_zhetrf.c
  starpu_zpotrf.c
  starpu_zpxtrf.c
  starpu_zsytrf.c
  starpu_ztrsm.c
)

precisions_rules_py(generated_sources
  "${SOURCES}"
  PRECISIONS "s;d;c;z")

### Add generated documented files to the global property before generating JDFs
add_documented_files(
  DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  ${generated_headers}
  ${generated_sources}
  )

add_library(pastix_starpu
  ${generated_headers}
  ${generated_sources}
  ${EXTRA_SOURCES})

target_link_libraries(pastix_starpu
  pastix_kernels
  ${STARPU_LIBRARIES_DEP}
  )

add_dependencies(pastix_starpu
  sopalin_headers
  starpu_headers_tgt
)

### Generate the lib
if (MPI_C_FOUND)
  set_target_properties(pastix_starpu PROPERTIES COMPILE_FLAGS "${MPI_C_COMPILE_FLAGS}")
endif (MPI_C_FOUND)

install(TARGETS pastix_starpu
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib)

### Add non-generated documented files to the global property
add_documented_files(
  # Headers
  pastix_starpu.h
  # Source files
  starpu.c
  starpu_sparse_matrix.c
  starpu_dense_matrix.c
  )
