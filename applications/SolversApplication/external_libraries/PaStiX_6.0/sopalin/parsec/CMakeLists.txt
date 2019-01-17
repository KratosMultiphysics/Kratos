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
  parsec.c
  parsec_sparse_matrix.c
)

### Generate the headers in all precisions
set(HEADERS
  pastix_zparsec.h
)

precisions_rules_py(generated_headers
  "${HEADERS}"
  PRECISIONS "s;d;c;z")

set(parsec_headers
  ${generated_headers}
  )

add_custom_target(parsec_headers_tgt
  DEPENDS ${parsec_headers} pastix_parsec.h )

### Generate the pastix_parsec wrappers for all required precisions
set(SOURCES
  parsec_zgetrf.c
  parsec_zhetrf.c
  parsec_zpotrf.c
  parsec_zpxtrf.c
  parsec_zsytrf.c
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

### Generate .c files from .jdf for all required precisions
set(JDF
  zgetrf_sp1dplus.jdf
  zhetrf_sp1dplus.jdf
  zpotrf_sp1dplus.jdf
  zpxtrf_sp1dplus.jdf
  zsytrf_sp1dplus.jdf
  #
  zgetrf_sp2d.jdf
  zhetrf_sp2d.jdf
  zpotrf_sp2d.jdf
  zpxtrf_sp2d.jdf
  zsytrf_sp2d.jdf
)
precisions_rules_py(generated_jdf
  "${JDF}"
  PRECISIONS "s;d;c;z")

jdf_rules(generated_sources
	 "${generated_jdf}"
	)

add_library(pastix_parsec
  ${generated_headers}
  ${generated_sources}
  ${EXTRA_SOURCES})

if ( ${PARSEC_parsec_distribution_LIBRARY} )
  target_link_libraries(pastix_parsec
    pastix_kernels
    ${PARSEC_LIBRARIES_DEP}
    parsec_distribution
    parsec_distribution_matrix
    ${PARSEC_LIBRARIES_DEP}
    parsec-base
    )
else()
  target_link_libraries(pastix_parsec
    pastix_kernels
    ${PARSEC_LIBRARIES_DEP}
    )
endif()

add_dependencies(pastix_parsec
  sopalin_headers
  parsec_headers_tgt
)

### Generate the lib
if (MPI_C_FOUND)
  set_target_properties(pastix_parsec PROPERTIES COMPILE_FLAGS "${MPI_C_COMPILE_FLAGS}")
endif (MPI_C_FOUND)

install(TARGETS pastix_parsec
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib)

### Add non-generated documented files to the global property
add_documented_files(
  # Headers
  pastix_parsec.h
  # Source files
  parsec.c
  parsec_sparse_matrix.c
  )
