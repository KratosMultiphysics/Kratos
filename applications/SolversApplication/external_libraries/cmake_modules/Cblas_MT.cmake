include_directories( ${SUPERLU_MT_DIR}/CBLAS )

set(CBLAS_HEADERS
  ${SUPERLU_MT_DIR}/CBLAS/slu_mt_Cnames.h
  ${SUPERLU_MT_DIR}/CBLAS/superlu_f2c.h
  ${SUPERLU_MT_DIR}/CBLAS/f2c.h
  )

set(CBLAS_SOURCES "")

if(enable_single)
  list(APPEND CBLAS_SOURCES
    ${SUPERLU_MT_DIR}/CBLAS/isamax.c
    ${SUPERLU_MT_DIR}/CBLAS/sasum.c
    ${SUPERLU_MT_DIR}/CBLAS/saxpy.c
    ${SUPERLU_MT_DIR}/CBLAS/scopy.c
    ${SUPERLU_MT_DIR}/CBLAS/sdot.c
    ${SUPERLU_MT_DIR}/CBLAS/snrm2.c
    ${SUPERLU_MT_DIR}/CBLAS/srot.c
    ${SUPERLU_MT_DIR}/CBLAS/sscal.c
    ${SUPERLU_MT_DIR}/CBLAS/sgemv.c
    ${SUPERLU_MT_DIR}/CBLAS/ssymv.c
    ${SUPERLU_MT_DIR}/CBLAS/strsv.c
    ${SUPERLU_MT_DIR}/CBLAS/sger.c
    ${SUPERLU_MT_DIR}/CBLAS/ssyr2.c
    )
endif(enable_single)

if(enable_double)
  list(APPEND CBLAS_SOURCES
    ${SUPERLU_MT_DIR}/CBLAS/idamax.c
    ${SUPERLU_MT_DIR}/CBLAS/dasum.c
    ${SUPERLU_MT_DIR}/CBLAS/daxpy.c
    ${SUPERLU_MT_DIR}/CBLAS/dcopy.c
    ${SUPERLU_MT_DIR}/CBLAS/ddot.c
    ${SUPERLU_MT_DIR}/CBLAS/dnrm2.c
    ${SUPERLU_MT_DIR}/CBLAS/drot.c
    ${SUPERLU_MT_DIR}/CBLAS/dscal.c
    ${SUPERLU_MT_DIR}/CBLAS/dgemv.c
    ${SUPERLU_MT_DIR}/CBLAS/dsymv.c
    ${SUPERLU_MT_DIR}/CBLAS/dtrsv.c
    ${SUPERLU_MT_DIR}/CBLAS/dger.c
    ${SUPERLU_MT_DIR}/CBLAS/dsyr2.c
    )
endif(enable_double)

if(enable_complex)
  list(APPEND CBLAS_SOURCES
    ${SUPERLU_MT_DIR}/CBLAS/icamax.c
    ${SUPERLU_MT_DIR}/CBLAS/scasum.c
    ${SUPERLU_MT_DIR}/CBLAS/caxpy.c
    ${SUPERLU_MT_DIR}/CBLAS/ccopy.c
    ${SUPERLU_MT_DIR}/CBLAS/scnrm2.c
    ${SUPERLU_MT_DIR}/CBLAS/cscal.c
    ${SUPERLU_MT_DIR}/CBLAS/cgemv.c
    ${SUPERLU_MT_DIR}/CBLAS/chemv.c
    ${SUPERLU_MT_DIR}/CBLAS/ctrsv.c
    ${SUPERLU_MT_DIR}/CBLAS/cgerc.c
    ${SUPERLU_MT_DIR}/CBLAS/cher2.c
    )
endif(enable_complex)

if(enable_complex16)
  list(APPEND CBLAS_SOURCES
    ${SUPERLU_MT_DIR}/CBLAS/izamax.c
    ${SUPERLU_MT_DIR}/CBLAS/dzasum.c
    ${SUPERLU_MT_DIR}/CBLAS/zaxpy.c
    ${SUPERLU_MT_DIR}/CBLAS/zcopy.c
    ${SUPERLU_MT_DIR}/CBLAS/dznrm2.c
    ${SUPERLU_MT_DIR}/CBLAS/zscal.c
    ${SUPERLU_MT_DIR}/CBLAS/dcabs1.c
    ${SUPERLU_MT_DIR}/CBLAS/zgemv.c
    ${SUPERLU_MT_DIR}/CBLAS/zhemv.c
    ${SUPERLU_MT_DIR}/CBLAS/ztrsv.c
    ${SUPERLU_MT_DIR}/CBLAS/zgerc.c
    ${SUPERLU_MT_DIR}/CBLAS/zher2.c
    )
endif(enable_complex16)

add_library(libblas_mt STATIC ${CBLAS_SOURCES} ${CBLAS_HEADERS})
