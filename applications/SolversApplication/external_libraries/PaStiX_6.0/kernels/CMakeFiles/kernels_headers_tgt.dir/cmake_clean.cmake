file(REMOVE_RECURSE
  "CMakeFiles/kernels_headers_tgt"
  "pastix_ccores.h"
  "pastix_scores.h"
  "pastix_dcores.h"
  "pastix_ccuda.h"
  "pastix_scuda.h"
  "pastix_dcuda.h"
  "pastix_clrcores.h"
  "pastix_slrcores.h"
  "pastix_dlrcores.h"
  "c_nan_check.h"
  "s_nan_check.h"
  "d_nan_check.h"
)

# Per-language clean rules from dependency scanning.
foreach(lang )
  include(CMakeFiles/kernels_headers_tgt.dir/cmake_clean_${lang}.cmake OPTIONAL)
endforeach()
