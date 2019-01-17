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
include(AddSourceFiles)

### reset variables
set(generated_sources "")
set(generated_headers "")

include_directories("${CMAKE_CURRENT_BINARY_DIR}")

### Generate the headers in all precisions
set(HEADERS
  z_refine_functions.h
)

precisions_rules_py(generated_headers
  "${HEADERS}"
  PRECISIONS "p;s;d;c;z")

set(refinement_headers
  ${generated_headers}
  )

### Generate the sources in all precisions
# set(SOURCES
#   z_refine_functions.c
#   z_refine_gmres.c
#   z_refine_grad.c
#   z_refine_pivot.c
#   z_refine_bicgstab.c
# )

# precisions_rules_py(generated_sources
#   "${SOURCES}"
#   PRECISIONS "p;s;d;c;z")

add_custom_target(refinement_headers_tgt
  DEPENDS ${refinement_headers} ${generated_sources})

# set(refinement_sources
#   ${generated_sources}
# )

### Add documented files to the global property
add_documented_files(
  DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  ${generated_headers}
  ${generated_sources}
  )
