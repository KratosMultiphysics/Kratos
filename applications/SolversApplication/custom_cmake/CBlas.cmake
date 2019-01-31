include_directories( ${SUPERLU_DIR}/CBLAS )

set(CBLAS_HEADERS
  ${SUPERLU_DIR}/CBLAS/slu_Cnames.h
  ${SUPERLU_DIR}/CBLAS/f2c.h
  )

set(CBLAS_SOURCES "")

if(enable_single)
  list(APPEND CBLAS_SOURCES
    ${SUPERLU_DIR}/CBLAS/isamax.c
    ${SUPERLU_DIR}/CBLAS/sasum.c
    ${SUPERLU_DIR}/CBLAS/saxpy.c
    ${SUPERLU_DIR}/CBLAS/scopy.c
    ${SUPERLU_DIR}/CBLAS/sdot.c
    ${SUPERLU_DIR}/CBLAS/snrm2.c
    ${SUPERLU_DIR}/CBLAS/srot.c
    ${SUPERLU_DIR}/CBLAS/sscal.c
    ${SUPERLU_DIR}/CBLAS/sswap.c
    ${SUPERLU_DIR}/CBLAS/sgemv.c
    ${SUPERLU_DIR}/CBLAS/ssymv.c
    ${SUPERLU_DIR}/CBLAS/strsv.c
    ${SUPERLU_DIR}/CBLAS/sger.c
    ${SUPERLU_DIR}/CBLAS/ssyr2.c
    )
endif(enable_single)

if(enable_double)
  list(APPEND CBLAS_SOURCES
    ${SUPERLU_DIR}/CBLAS/idamax.c
    ${SUPERLU_DIR}/CBLAS/dasum.c
    ${SUPERLU_DIR}/CBLAS/daxpy.c
    ${SUPERLU_DIR}/CBLAS/dcopy.c
    ${SUPERLU_DIR}/CBLAS/ddot.c
    ${SUPERLU_DIR}/CBLAS/dnrm2.c
    ${SUPERLU_DIR}/CBLAS/drot.c
    ${SUPERLU_DIR}/CBLAS/dscal.c
    ${SUPERLU_DIR}/CBLAS/dswap.c
    ${SUPERLU_DIR}/CBLAS/dgemv.c
    ${SUPERLU_DIR}/CBLAS/dsymv.c
    ${SUPERLU_DIR}/CBLAS/dtrsv.c
    ${SUPERLU_DIR}/CBLAS/dger.c
    ${SUPERLU_DIR}/CBLAS/dsyr2.c
    )
endif(enable_double)

if(enable_complex)
  list(APPEND CBLAS_SOURCES
    ${SUPERLU_DIR}/CBLAS/icamax.c
    ${SUPERLU_DIR}/CBLAS/scasum.c
    ${SUPERLU_DIR}/CBLAS/caxpy.c
    ${SUPERLU_DIR}/CBLAS/ccopy.c
    ${SUPERLU_DIR}/CBLAS/scnrm2.c
    ${SUPERLU_DIR}/CBLAS/cscal.c
    ${SUPERLU_DIR}/CBLAS/cswap.c
    ${SUPERLU_DIR}/CBLAS/cdotc.c
    ${SUPERLU_DIR}/CBLAS/cgemv.c
    ${SUPERLU_DIR}/CBLAS/chemv.c
    ${SUPERLU_DIR}/CBLAS/ctrsv.c
    ${SUPERLU_DIR}/CBLAS/cgerc.c
    ${SUPERLU_DIR}/CBLAS/cher2.c
    )
endif(enable_complex)

if(enable_complex16)
  list(APPEND CBLAS_SOURCES
    ${SUPERLU_DIR}/CBLAS/izamax.c
    ${SUPERLU_DIR}/CBLAS/dzasum.c
    ${SUPERLU_DIR}/CBLAS/zaxpy.c
    ${SUPERLU_DIR}/CBLAS/zcopy.c
    ${SUPERLU_DIR}/CBLAS/dznrm2.c
    ${SUPERLU_DIR}/CBLAS/zscal.c
    ${SUPERLU_DIR}/CBLAS/dcabs1.c
    ${SUPERLU_DIR}/CBLAS/zswap.c
    ${SUPERLU_DIR}/CBLAS/zdotc.c
    ${SUPERLU_DIR}/CBLAS/zgemv.c
    ${SUPERLU_DIR}/CBLAS/zhemv.c
    ${SUPERLU_DIR}/CBLAS/ztrsv.c
    ${SUPERLU_DIR}/CBLAS/zgerc.c
    ${SUPERLU_DIR}/CBLAS/zher2.c
    )
endif(enable_complex16)

add_library(external_libblas STATIC ${CBLAS_SOURCES} ${CBLAS_HEADERS})
