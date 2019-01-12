include_directories( ${SUPER_LU_MT_DIR}/CBLAS )

set(CBLAS_HEADERS
  ${SUPER_LU_MT_DIR}/CBLAS/slu_mt_Cnames.h
  ${SUPER_LU_MT_DIR}/CBLAS/superlu_f2c.h
  ${SUPER_LU_MT_DIR}/CBLAS/f2c.h
  )

set(CBLAS_SOURCES "")

if(enable_single)
  list(APPEND CBLAS_SOURCES
    ${SUPER_LU_MT_DIR}/CBLAS/isamax.c
    ${SUPER_LU_MT_DIR}/CBLAS/sasum.c
    ${SUPER_LU_MT_DIR}/CBLAS/saxpy.c
    ${SUPER_LU_MT_DIR}/CBLAS/scopy.c
    ${SUPER_LU_MT_DIR}/CBLAS/sdot.c
    ${SUPER_LU_MT_DIR}/CBLAS/snrm2.c
    ${SUPER_LU_MT_DIR}/CBLAS/srot.c
    ${SUPER_LU_MT_DIR}/CBLAS/sscal.c
    ${SUPER_LU_MT_DIR}/CBLAS/sgemv.c
    ${SUPER_LU_MT_DIR}/CBLAS/ssymv.c
    ${SUPER_LU_MT_DIR}/CBLAS/strsv.c
    ${SUPER_LU_MT_DIR}/CBLAS/sger.c
    ${SUPER_LU_MT_DIR}/CBLAS/ssyr2.c
    )
endif(enable_single)

if(enable_double)
  list(APPEND CBLAS_SOURCES
    ${SUPER_LU_MT_DIR}/CBLAS/idamax.c
    ${SUPER_LU_MT_DIR}/CBLAS/dasum.c
    ${SUPER_LU_MT_DIR}/CBLAS/daxpy.c
    ${SUPER_LU_MT_DIR}/CBLAS/dcopy.c
    ${SUPER_LU_MT_DIR}/CBLAS/ddot.c
    ${SUPER_LU_MT_DIR}/CBLAS/dnrm2.c
    ${SUPER_LU_MT_DIR}/CBLAS/drot.c
    ${SUPER_LU_MT_DIR}/CBLAS/dscal.c
    ${SUPER_LU_MT_DIR}/CBLAS/dgemv.c
    ${SUPER_LU_MT_DIR}/CBLAS/dsymv.c
    ${SUPER_LU_MT_DIR}/CBLAS/dtrsv.c
    ${SUPER_LU_MT_DIR}/CBLAS/dger.c
    ${SUPER_LU_MT_DIR}/CBLAS/dsyr2.c
    )
endif(enable_double)

if(enable_complex)
  list(APPEND CBLAS_SOURCES
    ${SUPER_LU_MT_DIR}/CBLAS/icamax.c
    ${SUPER_LU_MT_DIR}/CBLAS/scasum.c
    ${SUPER_LU_MT_DIR}/CBLAS/caxpy.c
    ${SUPER_LU_MT_DIR}/CBLAS/ccopy.c
    ${SUPER_LU_MT_DIR}/CBLAS/scnrm2.c
    ${SUPER_LU_MT_DIR}/CBLAS/cscal.c
    ${SUPER_LU_MT_DIR}/CBLAS/cgemv.c
    ${SUPER_LU_MT_DIR}/CBLAS/chemv.c
    ${SUPER_LU_MT_DIR}/CBLAS/ctrsv.c
    ${SUPER_LU_MT_DIR}/CBLAS/cgerc.c
    ${SUPER_LU_MT_DIR}/CBLAS/cher2.c
    )
endif(enable_complex)

if(enable_complex16)
  list(APPEND CBLAS_SOURCES
    ${SUPER_LU_MT_DIR}/CBLAS/izamax.c
    ${SUPER_LU_MT_DIR}/CBLAS/dzasum.c
    ${SUPER_LU_MT_DIR}/CBLAS/zaxpy.c
    ${SUPER_LU_MT_DIR}/CBLAS/zcopy.c
    ${SUPER_LU_MT_DIR}/CBLAS/dznrm2.c
    ${SUPER_LU_MT_DIR}/CBLAS/zscal.c
    ${SUPER_LU_MT_DIR}/CBLAS/dcabs1.c
    ${SUPER_LU_MT_DIR}/CBLAS/zgemv.c
    ${SUPER_LU_MT_DIR}/CBLAS/zhemv.c
    ${SUPER_LU_MT_DIR}/CBLAS/ztrsv.c
    ${SUPER_LU_MT_DIR}/CBLAS/zgerc.c
    ${SUPER_LU_MT_DIR}/CBLAS/zher2.c
    )
endif(enable_complex16)

add_library(libblas_mt STATIC ${CBLAS_SOURCES} ${CBLAS_HEADERS})
