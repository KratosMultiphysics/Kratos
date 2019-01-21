project (CBLAS C)

include_directories( ${CMAKE_CURRENT_SOURCE_DIR} )

set(CBLAS_HEADERS
  ${CMAKE_CURRENT_SOURCE_DIR}/slu_mt_Cnames.h
  ${CMAKE_CURRENT_SOURCE_DIR}/superlu_f2c.h
  ${CMAKE_CURRENT_SOURCE_DIR}/f2c.h
  )

set(CBLAS_SOURCES "")

if(enable_single)
  list(APPEND CBLAS_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/isamax.c
    ${CMAKE_CURRENT_SOURCE_DIR}/sasum.c
    ${CMAKE_CURRENT_SOURCE_DIR}/saxpy.c
    ${CMAKE_CURRENT_SOURCE_DIR}/scopy.c
    ${CMAKE_CURRENT_SOURCE_DIR}/sdot.c
    ${CMAKE_CURRENT_SOURCE_DIR}/snrm2.c
    ${CMAKE_CURRENT_SOURCE_DIR}/srot.c
    ${CMAKE_CURRENT_SOURCE_DIR}/sscal.c
    ${CMAKE_CURRENT_SOURCE_DIR}/sgemv.c
    ${CMAKE_CURRENT_SOURCE_DIR}/ssymv.c
    ${CMAKE_CURRENT_SOURCE_DIR}/strsv.c
    ${CMAKE_CURRENT_SOURCE_DIR}/sger.c
    ${CMAKE_CURRENT_SOURCE_DIR}/ssyr2.c
    )
endif(enable_single)

if(enable_double)
  list(APPEND CBLAS_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/idamax.c
    ${CMAKE_CURRENT_SOURCE_DIR}/dasum.c
    ${CMAKE_CURRENT_SOURCE_DIR}/daxpy.c
    ${CMAKE_CURRENT_SOURCE_DIR}/dcopy.c
    ${CMAKE_CURRENT_SOURCE_DIR}/ddot.c
    ${CMAKE_CURRENT_SOURCE_DIR}/dnrm2.c
    ${CMAKE_CURRENT_SOURCE_DIR}/drot.c
    ${CMAKE_CURRENT_SOURCE_DIR}/dscal.c
    ${CMAKE_CURRENT_SOURCE_DIR}/dgemv.c
    ${CMAKE_CURRENT_SOURCE_DIR}/dsymv.c
    ${CMAKE_CURRENT_SOURCE_DIR}/dtrsv.c
    ${CMAKE_CURRENT_SOURCE_DIR}/dger.c
    ${CMAKE_CURRENT_SOURCE_DIR}/dsyr2.c
    )
endif(enable_double)

if(enable_complex)
  list(APPEND CBLAS_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/icamax.c
    ${CMAKE_CURRENT_SOURCE_DIR}/scasum.c
    ${CMAKE_CURRENT_SOURCE_DIR}/caxpy.c
    ${CMAKE_CURRENT_SOURCE_DIR}/ccopy.c
    ${CMAKE_CURRENT_SOURCE_DIR}/scnrm2.c
    ${CMAKE_CURRENT_SOURCE_DIR}/cscal.c
    ${CMAKE_CURRENT_SOURCE_DIR}/cgemv.c
    ${CMAKE_CURRENT_SOURCE_DIR}/chemv.c
    ${CMAKE_CURRENT_SOURCE_DIR}/ctrsv.c
    ${CMAKE_CURRENT_SOURCE_DIR}/cgerc.c
    ${CMAKE_CURRENT_SOURCE_DIR}/cher2.c
    )
endif(enable_complex)

if(enable_complex16)
  list(APPEND CBLAS_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/izamax.c
    ${CMAKE_CURRENT_SOURCE_DIR}/dzasum.c
    ${CMAKE_CURRENT_SOURCE_DIR}/zaxpy.c
    ${CMAKE_CURRENT_SOURCE_DIR}/zcopy.c
    ${CMAKE_CURRENT_SOURCE_DIR}/dznrm2.c
    ${CMAKE_CURRENT_SOURCE_DIR}/zscal.c
    ${CMAKE_CURRENT_SOURCE_DIR}/dcabs1.c
    ${CMAKE_CURRENT_SOURCE_DIR}/zgemv.c
    ${CMAKE_CURRENT_SOURCE_DIR}/zhemv.c
    ${CMAKE_CURRENT_SOURCE_DIR}/ztrsv.c
    ${CMAKE_CURRENT_SOURCE_DIR}/zgerc.c
    ${CMAKE_CURRENT_SOURCE_DIR}/zher2.c
    )
endif(enable_complex16)

#add_definitions( -D_LONGINT )
add_definitions( -DNoChange )
add_definitions( -D__OPENMP )
add_definitions( -DPRNTlevel=0 )
add_definitions( -DDEBUGlevel=0 )
add_definitions( -w )

add_library(external_libblas_mt STATIC ${CBLAS_SOURCES} ${CBLAS_HEADERS})
