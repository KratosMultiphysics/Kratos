project (SUPERLU C)
include_directories( ${SUPERLU_DIR}/SRC )

# setup options
option(enable_single    "Enable single precision library" ON)
option(enable_double    "Enable double precision library" ON)
option(enable_complex   "Enable complex precision library" ON)
option(enable_complex16 "Enable complex16 precision library" ON)

set(SUPERLU_HEADERS
  ${SUPERLU_DIR}/SRC/supermatrix.h
  ${SUPERLU_DIR}/SRC/slu_Cnames.h
  ${SUPERLU_DIR}/SRC/slu_dcomplex.h
  ${SUPERLU_DIR}/SRC/slu_scomplex.h
  ${SUPERLU_DIR}/SRC/slu_util.h
  ${SUPERLU_DIR}/SRC/superlu_enum_consts.h
  )

set(SUPERLU_SOURCES
  ${SUPERLU_DIR}/SRC/superlu_timer.c
  ${SUPERLU_DIR}/SRC/util.c
  ${SUPERLU_DIR}/SRC/memory.c
  ${SUPERLU_DIR}/SRC/get_perm_c.c
  ${SUPERLU_DIR}/SRC/mmd.c
  ${SUPERLU_DIR}/SRC/sp_coletree.c
  ${SUPERLU_DIR}/SRC/sp_preorder.c
  ${SUPERLU_DIR}/SRC/sp_ienv.c
  ${SUPERLU_DIR}/SRC/relax_snode.c
  ${SUPERLU_DIR}/SRC/heap_relax_snode.c
  ${SUPERLU_DIR}/SRC/colamd.c
  ${SUPERLU_DIR}/SRC/ilu_relax_snode.c
  ${SUPERLU_DIR}/SRC/ilu_heap_relax_snode.c
  ${SUPERLU_DIR}/SRC/mark_relax.c
  ${SUPERLU_DIR}/SRC/mc64ad.c
  ${SUPERLU_DIR}/SRC/qselect.c
  ${SUPERLU_DIR}/SRC/input_error.c
  ${SUPERLU_DIR}/EXAMPLE/fgmr.c
  )
set_source_files_properties(superlu_timer.c PROPERTIES COMPILE_FLAGS -O0)

if(enable_single)
  list(APPEND SUPERLU_HEADERS
    ${SUPERLU_DIR}/SRC/slu_sdefs.h
    )

  list(APPEND SUPERLU_SOURCES
    ${SUPERLU_DIR}/SRC/slacon2.c
    ${SUPERLU_DIR}/SRC/smach.c
    ${SUPERLU_DIR}/SRC/sgssv.c
    ${SUPERLU_DIR}/SRC/sgssvx.c
    ${SUPERLU_DIR}/SRC/ssp_blas2.c
    ${SUPERLU_DIR}/SRC/ssp_blas3.c
    ${SUPERLU_DIR}/SRC/sgscon.c
    ${SUPERLU_DIR}/SRC/slangs.c
    ${SUPERLU_DIR}/SRC/sgsequ.c
    ${SUPERLU_DIR}/SRC/slaqgs.c
    ${SUPERLU_DIR}/SRC/spivotgrowth.c
    ${SUPERLU_DIR}/SRC/sgsrfs.c
    ${SUPERLU_DIR}/SRC/sgstrf.c
    ${SUPERLU_DIR}/SRC/sgstrs.c
    ${SUPERLU_DIR}/SRC/scopy_to_ucol.c
    ${SUPERLU_DIR}/SRC/ssnode_dfs.c
    ${SUPERLU_DIR}/SRC/ssnode_bmod.c
    ${SUPERLU_DIR}/SRC/spanel_dfs.c
    ${SUPERLU_DIR}/SRC/spanel_bmod.c
    ${SUPERLU_DIR}/SRC/sreadhb.c
    ${SUPERLU_DIR}/SRC/sreadrb.c
    ${SUPERLU_DIR}/SRC/sreadtriple.c
    ${SUPERLU_DIR}/SRC/scolumn_dfs.c
    ${SUPERLU_DIR}/SRC/scolumn_bmod.c
    ${SUPERLU_DIR}/SRC/spivotL.c
    ${SUPERLU_DIR}/SRC/spruneL.c
    ${SUPERLU_DIR}/SRC/smemory.c
    ${SUPERLU_DIR}/SRC/sutil.c
    ${SUPERLU_DIR}/SRC/smyblas2.c
    ${SUPERLU_DIR}/SRC/sgsisx.c
    ${SUPERLU_DIR}/SRC/sgsitrf.c
    ${SUPERLU_DIR}/SRC/sldperm.c
    ${SUPERLU_DIR}/SRC/ilu_sdrop_row.c
    ${SUPERLU_DIR}/SRC/ilu_ssnode_dfs.c
    ${SUPERLU_DIR}/SRC/ilu_scolumn_dfs.c
    ${SUPERLU_DIR}/SRC/ilu_spanel_dfs.c
    ${SUPERLU_DIR}/SRC/ilu_scopy_to_ucol.c
    ${SUPERLU_DIR}/SRC/ilu_spivotL.c
    ${SUPERLU_DIR}/SRC/sdiagonal.c
    ${SUPERLU_DIR}/SRC/sreadMM.c
    ${SUPERLU_DIR}/EXAMPLE/sfgmr.c
    )
  set_source_files_properties(smach.c PROPERTIES COMPILE_FLAGS -O0)
endif(enable_single)

if(enable_double)
  list(APPEND SUPERLU_HEADERS
    ${SUPERLU_DIR}/SRC/slu_ddefs.h
  )

  list(APPEND SUPERLU_SOURCES
    ${SUPERLU_DIR}/SRC/dlacon2.c
    ${SUPERLU_DIR}/SRC/dmach.c
    ${SUPERLU_DIR}/SRC/dgssv.c
    ${SUPERLU_DIR}/SRC/dgssvx.c
    ${SUPERLU_DIR}/SRC/dsp_blas2.c
    ${SUPERLU_DIR}/SRC/dsp_blas3.c
    ${SUPERLU_DIR}/SRC/dgscon.c
    ${SUPERLU_DIR}/SRC/dlangs.c
    ${SUPERLU_DIR}/SRC/dgsequ.c
    ${SUPERLU_DIR}/SRC/dlaqgs.c
    ${SUPERLU_DIR}/SRC/dpivotgrowth.c
    ${SUPERLU_DIR}/SRC/dgsrfs.c
    ${SUPERLU_DIR}/SRC/dgstrf.c
    ${SUPERLU_DIR}/SRC/dgstrs.c
    ${SUPERLU_DIR}/SRC/dcopy_to_ucol.c
    ${SUPERLU_DIR}/SRC/dsnode_dfs.c
    ${SUPERLU_DIR}/SRC/dsnode_bmod.c
    ${SUPERLU_DIR}/SRC/dpanel_dfs.c
    ${SUPERLU_DIR}/SRC/dpanel_bmod.c
    ${SUPERLU_DIR}/SRC/dreadhb.c
    ${SUPERLU_DIR}/SRC/dreadrb.c
    ${SUPERLU_DIR}/SRC/dreadtriple.c
    ${SUPERLU_DIR}/SRC/dcolumn_dfs.c
    ${SUPERLU_DIR}/SRC/dcolumn_bmod.c
    ${SUPERLU_DIR}/SRC/dpivotL.c
    ${SUPERLU_DIR}/SRC/dpruneL.c
    ${SUPERLU_DIR}/SRC/dmemory.c
    ${SUPERLU_DIR}/SRC/dutil.c
    ${SUPERLU_DIR}/SRC/dmyblas2.c
    ${SUPERLU_DIR}/SRC/dgsisx.c
    ${SUPERLU_DIR}/SRC/dgsitrf.c
    ${SUPERLU_DIR}/SRC/dldperm.c
    ${SUPERLU_DIR}/SRC/ilu_ddrop_row.c
    ${SUPERLU_DIR}/SRC/ilu_dsnode_dfs.c
    ${SUPERLU_DIR}/SRC/ilu_dcolumn_dfs.c
    ${SUPERLU_DIR}/SRC/ilu_dpanel_dfs.c
    ${SUPERLU_DIR}/SRC/ilu_dcopy_to_ucol.c
    ${SUPERLU_DIR}/SRC/ilu_dpivotL.c
    ${SUPERLU_DIR}/SRC/ddiagonal.c
    ${SUPERLU_DIR}/SRC/dreadMM.c
    ${SUPERLU_DIR}/SRC/dGetDiagU.c
    ${SUPERLU_DIR}/EXAMPLE/dfgmr.c
  )
  set_source_files_properties(dmach.c PROPERTIES COMPILE_FLAGS -O0)
endif(enable_double)

if(enable_complex)
  list(APPEND SUPERLU_HEADERS
    ${SUPERLU_DIR}/SRC/slu_cdefs.h
  )

  list(APPEND SUPERLU_SOURCES
    ${SUPERLU_DIR}/SRC/clacon2.c
    ${SUPERLU_DIR}/SRC/scsum1.c
    ${SUPERLU_DIR}/SRC/icmax1.c
    ${SUPERLU_DIR}/SRC/scomplex.c
    ${SUPERLU_DIR}/SRC/cgssv.c
    ${SUPERLU_DIR}/SRC/cgssvx.c
    ${SUPERLU_DIR}/SRC/csp_blas2.c
    ${SUPERLU_DIR}/SRC/csp_blas3.c
    ${SUPERLU_DIR}/SRC/cgscon.c
    ${SUPERLU_DIR}/SRC/clangs.c
    ${SUPERLU_DIR}/SRC/cgsequ.c
    ${SUPERLU_DIR}/SRC/claqgs.c
    ${SUPERLU_DIR}/SRC/cpivotgrowth.c
    ${SUPERLU_DIR}/SRC/cgsrfs.c
    ${SUPERLU_DIR}/SRC/cgstrf.c
    ${SUPERLU_DIR}/SRC/cgstrs.c
    ${SUPERLU_DIR}/SRC/ccopy_to_ucol.c
    ${SUPERLU_DIR}/SRC/csnode_dfs.c
    ${SUPERLU_DIR}/SRC/csnode_bmod.c
    ${SUPERLU_DIR}/SRC/cpanel_dfs.c
    ${SUPERLU_DIR}/SRC/cpanel_bmod.c
    ${SUPERLU_DIR}/SRC/creadhb.c
    ${SUPERLU_DIR}/SRC/creadrb.c
    ${SUPERLU_DIR}/SRC/creadtriple.c
    ${SUPERLU_DIR}/SRC/ccolumn_dfs.c
    ${SUPERLU_DIR}/SRC/ccolumn_bmod.c
    ${SUPERLU_DIR}/SRC/cpivotL.c
    ${SUPERLU_DIR}/SRC/cpruneL.c
    ${SUPERLU_DIR}/SRC/cmemory.c
    ${SUPERLU_DIR}/SRC/cutil.c
    ${SUPERLU_DIR}/SRC/cmyblas2.c
    ${SUPERLU_DIR}/SRC/cgsisx.c
    ${SUPERLU_DIR}/SRC/cgsitrf.c
    ${SUPERLU_DIR}/SRC/cldperm.c
    ${SUPERLU_DIR}/SRC/ilu_cdrop_row.c
    ${SUPERLU_DIR}/SRC/ilu_csnode_dfs.c
    ${SUPERLU_DIR}/SRC/ilu_ccolumn_dfs.c
    ${SUPERLU_DIR}/SRC/ilu_cpanel_dfs.c
    ${SUPERLU_DIR}/SRC/ilu_ccopy_to_ucol.c
    ${SUPERLU_DIR}/SRC/ilu_cpivotL.c
    ${SUPERLU_DIR}/SRC/cdiagonal.c
    ${SUPERLU_DIR}/SRC/creadMM.c
    ${SUPERLU_DIR}/EXAMPLE/cfgmr.c
  )
endif(enable_complex)

if(enable_complex16)
  list(APPEND SUPERLU_HEADERS
    ${SUPERLU_DIR}/SRC/slu_zdefs.h
  )

  list(APPEND SUPERLU_SOURCES
    ${SUPERLU_DIR}/SRC/zlacon2.c
    ${SUPERLU_DIR}/SRC/dzsum1.c
    ${SUPERLU_DIR}/SRC/izmax1.c
    ${SUPERLU_DIR}/SRC/dcomplex.c
    ${SUPERLU_DIR}/SRC/zgssv.c
    ${SUPERLU_DIR}/SRC/zgssvx.c
    ${SUPERLU_DIR}/SRC/zsp_blas2.c
    ${SUPERLU_DIR}/SRC/zsp_blas3.c
    ${SUPERLU_DIR}/SRC/zgscon.c
    ${SUPERLU_DIR}/SRC/zlangs.c
    ${SUPERLU_DIR}/SRC/zgsequ.c
    ${SUPERLU_DIR}/SRC/zlaqgs.c
    ${SUPERLU_DIR}/SRC/zpivotgrowth.c
    ${SUPERLU_DIR}/SRC/zgsrfs.c
    ${SUPERLU_DIR}/SRC/zgstrf.c
    ${SUPERLU_DIR}/SRC/zgstrs.c
    ${SUPERLU_DIR}/SRC/zcopy_to_ucol.c
    ${SUPERLU_DIR}/SRC/zsnode_dfs.c
    ${SUPERLU_DIR}/SRC/zsnode_bmod.c
    ${SUPERLU_DIR}/SRC/zpanel_dfs.c
    ${SUPERLU_DIR}/SRC/zpanel_bmod.c
    ${SUPERLU_DIR}/SRC/zreadhb.c
    ${SUPERLU_DIR}/SRC/zreadrb.c
    ${SUPERLU_DIR}/SRC/zreadtriple.c
    ${SUPERLU_DIR}/SRC/zcolumn_dfs.c
    ${SUPERLU_DIR}/SRC/zcolumn_bmod.c
    ${SUPERLU_DIR}/SRC/zpivotL.c
    ${SUPERLU_DIR}/SRC/zpruneL.c
    ${SUPERLU_DIR}/SRC/zmemory.c
    ${SUPERLU_DIR}/SRC/zutil.c
    ${SUPERLU_DIR}/SRC/zmyblas2.c
    ${SUPERLU_DIR}/SRC/zgsisx.c
    ${SUPERLU_DIR}/SRC/zgsitrf.c
    ${SUPERLU_DIR}/SRC/zldperm.c
    ${SUPERLU_DIR}/SRC/ilu_zdrop_row.c
    ${SUPERLU_DIR}/SRC/ilu_zsnode_dfs.c
    ${SUPERLU_DIR}/SRC/ilu_zcolumn_dfs.c
    ${SUPERLU_DIR}/SRC/ilu_zpanel_dfs.c
    ${SUPERLU_DIR}/SRC/ilu_zcopy_to_ucol.c
    ${SUPERLU_DIR}/SRC/ilu_zpivotL.c
    ${SUPERLU_DIR}/SRC/zdiagonal.c
    ${SUPERLU_DIR}/SRC/zreadMM.c
    ${SUPERLU_DIR}/EXAMPLE/zfgmr.c
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
  INCLUDE("${CMAKE_CURRENT_SOURCE_DIR}/external_libraries/CMakeFiles/CBlas.cmake")
  set(CBLAS_LIBRARIES external_libblas)
endif(BLAS_FOUND)

message(STATUS "SUPERLU cblas: ${CBLAS_LIBRARIES}")

add_library(external_superlu STATIC ${SUPERLU_SOURCES} ${SUPERLU_HEADERS})

message(STATUS "BLAS_LIBRARIES : ${CBLAS_LIBRARIES}" )
target_link_libraries(external_superlu ${CBLAS_LIBRARIES} )
set(SUPERLU_LIBRARIES external_superlu)
