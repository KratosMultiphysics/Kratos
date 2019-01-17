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

# reset variables
set(generated_sources "")
set(generated_headers "")

include_directories("${CMAKE_CURRENT_BINARY_DIR}")

### Add trace generation
if(PASTIX_WITH_EZTRACE)
  add_library(eztrace-convert-kernels SHARED
        eztrace_module/eztrace_convert_kernels.c
  )
endif(PASTIX_WITH_EZTRACE)

### Generate the headers in all precisions
set(HEADERS
  pastix_zcores.h
  pastix_zcuda.h
  pastix_zlrcores.h
  z_nan_check.h
)

precisions_rules_py(generated_headers
  "${HEADERS}"
  PRECISIONS "s;d;c;z")

set(kernels_headers
  ${generated_headers}
  )

add_custom_target(kernels_headers_tgt
  DEPENDS ${kernels_headers} )

### Generate the sources in all precisions
set(SOURCES
  # CPU extra kernel
  core_zgemdm.c
  core_zgetro.c
  core_zgeadd.c
  core_zplrnt.c
  core_ztradd.c
  core_zgemmsp.c
  core_ztrsmsp.c
  core_zscalo.c
  # Kernels
  core_zsytrfsp.c
  core_zhetrfsp.c
  core_zpotrfsp.c
  core_zpxtrfsp.c
  core_zgetrfsp.c
  # Low-rank kernels
  core_zlrdbg.c
  core_zlr2xx.c
  core_zxx2lr.c
  core_zxx2fr.c
  core_zlrmm.c
  core_zgelrops_rrqr.c
  core_zgelrops_svd.c
  core_zgelrops.c
  # cblk operations
  cpucblk_zinit.c
  cpucblk_zcompress.c
  cpucblk_zdiff.c
  cpucblk_zschur.c
  # Solve operation
  solve_ztrsmsp.c
)

precisions_rules_py(generated_sources
  "${SOURCES}"
  PRECISIONS "s;d;c;z")

set(kernels_sources
  ${generated_sources}
  kernels_trace.c
  )

if(PASTIX_WITH_CUDA)
   set( CUDA_SOURCES
     # GPU kernel
     gpu_zgemmsp.c
    )
   set(CUDA_generated_files "")
   precisions_rules_py(CUDA_generated_files
        "${CUDA_SOURCES}"
        PRECISIONS "s;d;c;z")

   set(kernels_sources
       ${kernels_sources}
       ${CUDA_generated_files})

   add_subdirectory(gpus)

endif(PASTIX_WITH_CUDA)

add_library(pastix_kernels
  ${kernels_sources}
  )

add_dependencies(pastix_kernels
  kernels_headers_tgt
)

if(PASTIX_WITH_EZTRACE)
   target_link_libraries(pastix_kernels
     ${EZTRACE_LIBRARIES}
     litl
   )
endif(PASTIX_WITH_EZTRACE)

target_link_libraries( pastix_kernels
  ${LAPACKE_LIBRARIES_DEP}
  ${CBLAS_LIBRARIES_DEP}
  ${EXTRA_LIBS}
)

if(PASTIX_WITH_CUDA)
    # if(PASTIX_CUDA_FERMI)
    #   target_link_libraries( pastix_kernels
    #                          pastix_cucores_sm20 )
    # else()
    #   target_link_libraries( pastix_kernels
    #                          pastix_cucores_sm35 )
    # endif()
    target_link_libraries( pastix_kernels
      pastix_kernels_cuda )
endif()

### Install library
install(TARGETS pastix_kernels
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib)

### Add documented files to the global property
add_documented_files(
  DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  ${generated_headers}
  ${generated_sources}
  )

# add_documented_files(
#   # Headers
#   # Source files
#   )
