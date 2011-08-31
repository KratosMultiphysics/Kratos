# Find the PARDISO that comes with MKL includes and libraries
#

# MKLSOLVER_INCLUDE_DIR - where to find autopack.h
# MKLSOLVER_LIBRARIES   - List of fully qualified libraries to link against.
# MKLSOLVER_FOUND       - Do not attempt to use if "no" or undefined.

FIND_PATH(MKLSOLVER_INCLUDE_DIR mkl_blas.h
  ${MKLSOLVER_INCLUDE_DIR}
)

# set(LIB_MKL_LIST mkl_intel_ilp64 mkl_intel_thread mkl_core mkl_solver_ilp64 mkl_mc mkl_mc3 mkl_lapack pthread mkl_p4n iomp5)

set(MY_LIB_SEARCH_DIR ${MKLSOLVER_LIB_DIR})

FIND_LIBRARY(AUX1 mkl_intel_ilp64  
  ${MY_LIB_SEARCH_DIR}
)

if(CMAKE_COMPILER_IS_GNUCXX )
    FIND_LIBRARY(AUX2 mkl_gnu_thread  
      ${MY_LIB_SEARCH_DIR}
    )
else(CMAKE_COMPILER_IS_GNUCXX )
    FIND_LIBRARY(AUX2 mkl_intel_thread  
      ${MY_LIB_SEARCH_DIR}
    )
endif(CMAKE_COMPILER_IS_GNUCXX )

FIND_LIBRARY(AUX3 mkl_core  
  ${MY_LIB_SEARCH_DIR}
)

FIND_LIBRARY(AUX4 mkl_solver_ilp64  
  ${MY_LIB_SEARCH_DIR}
)

FIND_LIBRARY(AUX5 mkl_mc  
  ${MY_LIB_SEARCH_DIR}
)

FIND_LIBRARY(AUX6 mkl_mc3  
  ${MY_LIB_SEARCH_DIR}
)

FIND_LIBRARY(AUX7 mkl_lapack  
  ${MY_LIB_SEARCH_DIR}
)

FIND_LIBRARY(AUX8 pthread  
  ${MY_LIB_SEARCH_DIR}
)

FIND_LIBRARY(AUX9 mkl_p4n  
  ${MY_LIB_SEARCH_DIR}
)







##ATTENTION! the order of libraries is FUNDAMENTAL!! take care!!!
IF(MKLSOLVER_INCLUDE_DIR)
  IF(AUX1)
    SET( MKLSOLVER_LIBRARIES ${AUX9} ${AUX8} ${AUX6} ${AUX5}  ${AUX7} ${AUX4} ${AUX3} ${AUX2}  ${AUX1}  )
#     SET( MKLSOLVER_LIBRARIES ${AUX1} ${AUX2} ${AUX3} ${AUX4}  ${AUX7} ${AUX5} ${AUX6} ${AUX8}  ${AUX9}  )
    SET( MKLSOLVER_FOUND "YES" )
    message("************ mkl solver libraries found  ${MKLSOLVER_LIBRARIES}") 
  ELSE(AUX1)
    message("finding MKLSOLVER library failed, please try to set the var MKLSOLVER_ROOT_DIR")
  ENDIF(AUX1)
ELSE(MKLSOLVER_INCLUDE_DIR)
  message("finding MKLSOLVER includes failed, please try to set the var MKLSOLVER_ROOT_DIR")
ENDIF(MKLSOLVER_INCLUDE_DIR)


