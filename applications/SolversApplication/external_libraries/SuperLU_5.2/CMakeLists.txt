project (SUPERLU C)

include_directories( ${CMAKE_CURRENT_SOURCE_DIR}/SRC )

# setup options
option(enable_single    "Enable single precision library" ON)
option(enable_double    "Enable double precision library" ON)
option(enable_complex   "Enable complex precision library" ON)
option(enable_complex16 "Enable complex16 precision library" ON)

set(SUPERLU_HEADERS
  ${CMAKE_CURRENT_SOURCE_DIR}/SRC/supermatrix.h
  ${CMAKE_CURRENT_SOURCE_DIR}/SRC/slu_Cnames.h
  ${CMAKE_CURRENT_SOURCE_DIR}/SRC/slu_dcomplex.h
  ${CMAKE_CURRENT_SOURCE_DIR}/SRC/slu_scomplex.h
  ${CMAKE_CURRENT_SOURCE_DIR}/SRC/slu_util.h
  ${CMAKE_CURRENT_SOURCE_DIR}/SRC/superlu_enum_consts.h
  )

set(SUPERLU_SOURCES
  ${CMAKE_CURRENT_SOURCE_DIR}/SRC/superlu_timer.c
  ${CMAKE_CURRENT_SOURCE_DIR}/SRC/util.c
  ${CMAKE_CURRENT_SOURCE_DIR}/SRC/memory.c
  ${CMAKE_CURRENT_SOURCE_DIR}/SRC/get_perm_c.c
  ${CMAKE_CURRENT_SOURCE_DIR}/SRC/mmd.c
  ${CMAKE_CURRENT_SOURCE_DIR}/SRC/sp_coletree.c
  ${CMAKE_CURRENT_SOURCE_DIR}/SRC/sp_preorder.c
  ${CMAKE_CURRENT_SOURCE_DIR}/SRC/sp_ienv.c
  ${CMAKE_CURRENT_SOURCE_DIR}/SRC/relax_snode.c
  ${CMAKE_CURRENT_SOURCE_DIR}/SRC/heap_relax_snode.c
  ${CMAKE_CURRENT_SOURCE_DIR}/SRC/colamd.c
  ${CMAKE_CURRENT_SOURCE_DIR}/SRC/ilu_relax_snode.c
  ${CMAKE_CURRENT_SOURCE_DIR}/SRC/ilu_heap_relax_snode.c
  ${CMAKE_CURRENT_SOURCE_DIR}/SRC/mark_relax.c
  ${CMAKE_CURRENT_SOURCE_DIR}/SRC/mc64ad.c
  ${CMAKE_CURRENT_SOURCE_DIR}/SRC/qselect.c
  ${CMAKE_CURRENT_SOURCE_DIR}/SRC/input_error.c
  )
set_source_files_properties(superlu_timer.c PROPERTIES COMPILE_FLAGS -O0)

if(enable_single)
  list(APPEND SUPERLU_HEADERS
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/slu_sdefs.h
    )

  list(APPEND SUPERLU_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/slacon2.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/smach.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/sgssv.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/sgssvx.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/ssp_blas2.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/ssp_blas3.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/sgscon.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/slangs.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/sgsequ.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/slaqgs.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/spivotgrowth.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/sgsrfs.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/sgstrf.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/sgstrs.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/scopy_to_ucol.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/ssnode_dfs.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/ssnode_bmod.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/spanel_dfs.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/spanel_bmod.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/sreadhb.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/sreadrb.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/sreadtriple.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/scolumn_dfs.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/scolumn_bmod.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/spivotL.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/spruneL.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/smemory.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/sutil.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/smyblas2.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/sgsisx.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/sgsitrf.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/sldperm.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/ilu_sdrop_row.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/ilu_ssnode_dfs.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/ilu_scolumn_dfs.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/ilu_spanel_dfs.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/ilu_scopy_to_ucol.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/ilu_spivotL.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/sdiagonal.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/sreadMM.c
    )
  set_source_files_properties(smach.c PROPERTIES COMPILE_FLAGS -O0)
endif(enable_single)

if(enable_double)
  list(APPEND SUPERLU_HEADERS
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/slu_ddefs.h
  )

  list(APPEND SUPERLU_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/dlacon2.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/dmach.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/dgssv.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/dgssvx.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/dsp_blas2.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/dsp_blas3.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/dgscon.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/dlangs.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/dgsequ.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/dlaqgs.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/dpivotgrowth.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/dgsrfs.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/dgstrf.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/dgstrs.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/dcopy_to_ucol.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/dsnode_dfs.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/dsnode_bmod.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/dpanel_dfs.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/dpanel_bmod.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/dreadhb.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/dreadrb.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/dreadtriple.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/dcolumn_dfs.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/dcolumn_bmod.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/dpivotL.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/dpruneL.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/dmemory.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/dutil.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/dmyblas2.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/dgsisx.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/dgsitrf.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/dldperm.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/ilu_ddrop_row.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/ilu_dsnode_dfs.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/ilu_dcolumn_dfs.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/ilu_dpanel_dfs.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/ilu_dcopy_to_ucol.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/ilu_dpivotL.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/ddiagonal.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/dreadMM.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/dGetDiagU.c
  )
  set_source_files_properties(dmach.c PROPERTIES COMPILE_FLAGS -O0)
endif(enable_double)

if(enable_complex)
  list(APPEND SUPERLU_HEADERS
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/slu_cdefs.h
  )

  list(APPEND SUPERLU_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/clacon2.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/scsum1.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/icmax1.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/scomplex.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/cgssv.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/cgssvx.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/csp_blas2.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/csp_blas3.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/cgscon.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/clangs.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/cgsequ.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/claqgs.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/cpivotgrowth.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/cgsrfs.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/cgstrf.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/cgstrs.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/ccopy_to_ucol.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/csnode_dfs.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/csnode_bmod.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/cpanel_dfs.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/cpanel_bmod.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/creadhb.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/creadrb.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/creadtriple.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/ccolumn_dfs.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/ccolumn_bmod.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/cpivotL.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/cpruneL.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/cmemory.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/cutil.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/cmyblas2.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/cgsisx.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/cgsitrf.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/cldperm.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/ilu_cdrop_row.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/ilu_csnode_dfs.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/ilu_ccolumn_dfs.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/ilu_cpanel_dfs.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/ilu_ccopy_to_ucol.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/ilu_cpivotL.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/cdiagonal.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/creadMM.c
  )
endif(enable_complex)

if(enable_complex16)
  list(APPEND SUPERLU_HEADERS
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/slu_zdefs.h
  )

  list(APPEND SUPERLU_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/zlacon2.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/dzsum1.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/izmax1.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/dcomplex.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/zgssv.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/zgssvx.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/zsp_blas2.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/zsp_blas3.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/zgscon.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/zlangs.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/zgsequ.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/zlaqgs.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/zpivotgrowth.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/zgsrfs.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/zgstrf.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/zgstrs.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/zcopy_to_ucol.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/zsnode_dfs.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/zsnode_bmod.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/zpanel_dfs.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/zpanel_bmod.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/zreadhb.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/zreadrb.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/zreadtriple.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/zcolumn_dfs.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/zcolumn_bmod.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/zpivotL.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/zpruneL.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/zmemory.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/zutil.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/zmyblas2.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/zgsisx.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/zgsitrf.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/zldperm.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/ilu_zdrop_row.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/ilu_zsnode_dfs.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/ilu_zcolumn_dfs.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/ilu_zpanel_dfs.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/ilu_zcopy_to_ucol.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/ilu_zpivotL.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/zdiagonal.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/zreadMM.c
  )
endif(enable_complex16)

#add_definitions( -D_LONGINT )
add_definitions( -DAdd_ )
add_definitions( -fPIC )
add_definitions( -w )

if(NOT BLAS_FOUND)
  find_package(BLAS)
endif(NOT BLAS_FOUND)

if(BLAS_FOUND)
  set(CBLAS_LIBRARIES ${BLAS_LIBRARIES})
  add_definitions( -DUSE_VENDOR_BLAS )
else(BLAS_FOUND)
  add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/CBLAS")
  set(CBLAS_LIBRARIES external_libblas)
endif(BLAS_FOUND)

message(STATUS "SUPERLU cblas: ${CBLAS_LIBRARIES}")

add_library(external_superlu STATIC ${SUPERLU_SOURCES} ${SUPERLU_HEADERS})

message(STATUS "BLAS_LIBRARIES : ${CBLAS_LIBRARIES}" )
target_link_libraries(external_superlu ${CBLAS_LIBRARIES} )
