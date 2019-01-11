include_directories( ${SUPER_LU_DIR}/SRC )

# setup options
option(enable_single    "Enable single precision library" ON)
option(enable_double    "Enable double precision library" ON)
option(enable_complex   "Enable complex precision library" ON)
option(enable_complex16 "Enable complex16 precision library" ON)

set(SUPER_LU_HEADERS
  ${SUPER_LU_DIR}/SRC/supermatrix.h
  ${SUPER_LU_DIR}/SRC/slu_Cnames.h
  ${SUPER_LU_DIR}/SRC/slu_dcomplex.h
  ${SUPER_LU_DIR}/SRC/slu_scomplex.h
  ${SUPER_LU_DIR}/SRC/slu_util.h
  ${SUPER_LU_DIR}/SRC/superlu_enum_consts.h
  )

set(SUPER_LU_SOURCES
  ${SUPER_LU_DIR}/SRC/superlu_timer.c
  ${SUPER_LU_DIR}/SRC/util.c
  ${SUPER_LU_DIR}/SRC/memory.c
  ${SUPER_LU_DIR}/SRC/get_perm_c.c
  ${SUPER_LU_DIR}/SRC/mmd.c
  ${SUPER_LU_DIR}/SRC/sp_coletree.c
  ${SUPER_LU_DIR}/SRC/sp_preorder.c
  ${SUPER_LU_DIR}/SRC/sp_ienv.c
  ${SUPER_LU_DIR}/SRC/relax_snode.c
  ${SUPER_LU_DIR}/SRC/heap_relax_snode.c
  ${SUPER_LU_DIR}/SRC/colamd.c
  ${SUPER_LU_DIR}/SRC/ilu_relax_snode.c
  ${SUPER_LU_DIR}/SRC/ilu_heap_relax_snode.c
  ${SUPER_LU_DIR}/SRC/mark_relax.c
  ${SUPER_LU_DIR}/SRC/mc64ad.c
  ${SUPER_LU_DIR}/SRC/qselect.c
  ${SUPER_LU_DIR}/SRC/input_error.c
  )
set_source_files_properties(superlu_timer.c PROPERTIES COMPILE_FLAGS -O0)

if(enable_single)
  list(APPEND SUPER_LU_HEADERS
    ${SUPER_LU_DIR}/SRC/slu_sdefs.h
    )

  list(APPEND SUPER_LU_SOURCES
    ${SUPER_LU_DIR}/SRC/slacon2.c
    ${SUPER_LU_DIR}/SRC/smach.c
    ${SUPER_LU_DIR}/SRC/sgssv.c
    ${SUPER_LU_DIR}/SRC/sgssvx.c
    ${SUPER_LU_DIR}/SRC/ssp_blas2.c
    ${SUPER_LU_DIR}/SRC/ssp_blas3.c
    ${SUPER_LU_DIR}/SRC/sgscon.c
    ${SUPER_LU_DIR}/SRC/slangs.c
    ${SUPER_LU_DIR}/SRC/sgsequ.c
    ${SUPER_LU_DIR}/SRC/slaqgs.c
    ${SUPER_LU_DIR}/SRC/spivotgrowth.c
    ${SUPER_LU_DIR}/SRC/sgsrfs.c
    ${SUPER_LU_DIR}/SRC/sgstrf.c
    ${SUPER_LU_DIR}/SRC/sgstrs.c
    ${SUPER_LU_DIR}/SRC/scopy_to_ucol.c
    ${SUPER_LU_DIR}/SRC/ssnode_dfs.c
    ${SUPER_LU_DIR}/SRC/ssnode_bmod.c
    ${SUPER_LU_DIR}/SRC/spanel_dfs.c
    ${SUPER_LU_DIR}/SRC/spanel_bmod.c
    ${SUPER_LU_DIR}/SRC/sreadhb.c
    ${SUPER_LU_DIR}/SRC/sreadrb.c
    ${SUPER_LU_DIR}/SRC/sreadtriple.c
    ${SUPER_LU_DIR}/SRC/scolumn_dfs.c
    ${SUPER_LU_DIR}/SRC/scolumn_bmod.c
    ${SUPER_LU_DIR}/SRC/spivotL.c
    ${SUPER_LU_DIR}/SRC/spruneL.c
    ${SUPER_LU_DIR}/SRC/smemory.c
    ${SUPER_LU_DIR}/SRC/sutil.c
    ${SUPER_LU_DIR}/SRC/smyblas2.c
    ${SUPER_LU_DIR}/SRC/sgsisx.c
    ${SUPER_LU_DIR}/SRC/sgsitrf.c
    ${SUPER_LU_DIR}/SRC/sldperm.c
    ${SUPER_LU_DIR}/SRC/ilu_sdrop_row.c
    ${SUPER_LU_DIR}/SRC/ilu_ssnode_dfs.c
    ${SUPER_LU_DIR}/SRC/ilu_scolumn_dfs.c
    ${SUPER_LU_DIR}/SRC/ilu_spanel_dfs.c
    ${SUPER_LU_DIR}/SRC/ilu_scopy_to_ucol.c
    ${SUPER_LU_DIR}/SRC/ilu_spivotL.c
    ${SUPER_LU_DIR}/SRC/sdiagonal.c
    )
  set_source_files_properties(smach.c PROPERTIES COMPILE_FLAGS -O0)
endif(enable_single)

if(enable_double)
  list(APPEND SUPER_LU_HEADERS
    ${SUPER_LU_DIR}/SRC/slu_ddefs.h
  )

  list(APPEND SUPER_LU_SOURCES
    ${SUPER_LU_DIR}/SRC/dlacon2.c
    ${SUPER_LU_DIR}/SRC/dmach.c
    ${SUPER_LU_DIR}/SRC/dgssv.c
    ${SUPER_LU_DIR}/SRC/dgssvx.c
    ${SUPER_LU_DIR}/SRC/dsp_blas2.c
    ${SUPER_LU_DIR}/SRC/dsp_blas3.c
    ${SUPER_LU_DIR}/SRC/dgscon.c
    ${SUPER_LU_DIR}/SRC/dlangs.c
    ${SUPER_LU_DIR}/SRC/dgsequ.c
    ${SUPER_LU_DIR}/SRC/dlaqgs.c
    ${SUPER_LU_DIR}/SRC/dpivotgrowth.c
    ${SUPER_LU_DIR}/SRC/dgsrfs.c
    ${SUPER_LU_DIR}/SRC/dgstrf.c
    ${SUPER_LU_DIR}/SRC/dgstrs.c
    ${SUPER_LU_DIR}/SRC/dcopy_to_ucol.c
    ${SUPER_LU_DIR}/SRC/dsnode_dfs.c
    ${SUPER_LU_DIR}/SRC/dsnode_bmod.c
    ${SUPER_LU_DIR}/SRC/dpanel_dfs.c
    ${SUPER_LU_DIR}/SRC/dpanel_bmod.c
    ${SUPER_LU_DIR}/SRC/dreadhb.c
    ${SUPER_LU_DIR}/SRC/dreadrb.c
    ${SUPER_LU_DIR}/SRC/dreadtriple.c
    ${SUPER_LU_DIR}/SRC/dcolumn_dfs.c
    ${SUPER_LU_DIR}/SRC/dcolumn_bmod.c
    ${SUPER_LU_DIR}/SRC/dpivotL.c
    ${SUPER_LU_DIR}/SRC/dpruneL.c
    ${SUPER_LU_DIR}/SRC/dmemory.c
    ${SUPER_LU_DIR}/SRC/dutil.c
    ${SUPER_LU_DIR}/SRC/dmyblas2.c
    ${SUPER_LU_DIR}/SRC/dgsisx.c
    ${SUPER_LU_DIR}/SRC/dgsitrf.c
    ${SUPER_LU_DIR}/SRC/dldperm.c
    ${SUPER_LU_DIR}/SRC/ilu_ddrop_row.c
    ${SUPER_LU_DIR}/SRC/ilu_dsnode_dfs.c
    ${SUPER_LU_DIR}/SRC/ilu_dcolumn_dfs.c
    ${SUPER_LU_DIR}/SRC/ilu_dpanel_dfs.c
    ${SUPER_LU_DIR}/SRC/ilu_dcopy_to_ucol.c
    ${SUPER_LU_DIR}/SRC/ilu_dpivotL.c
    ${SUPER_LU_DIR}/SRC/ddiagonal.c
  )
  set_source_files_properties(dmach.c PROPERTIES COMPILE_FLAGS -O0)
endif(enable_double)

if(enable_complex)
  list(APPEND SUPER_LU_HEADERS
    ${SUPER_LU_DIR}/SRC/slu_cdefs.h
  )

  list(APPEND SUPER_LU_SOURCES
    ${SUPER_LU_DIR}/SRC/clacon2.c
    ${SUPER_LU_DIR}/SRC/scsum1.c
    ${SUPER_LU_DIR}/SRC/icmax1.c
    ${SUPER_LU_DIR}/SRC/scomplex.c
    ${SUPER_LU_DIR}/SRC/cgssv.c
    ${SUPER_LU_DIR}/SRC/cgssvx.c
    ${SUPER_LU_DIR}/SRC/csp_blas2.c
    ${SUPER_LU_DIR}/SRC/csp_blas3.c
    ${SUPER_LU_DIR}/SRC/cgscon.c
    ${SUPER_LU_DIR}/SRC/clangs.c
    ${SUPER_LU_DIR}/SRC/cgsequ.c
    ${SUPER_LU_DIR}/SRC/claqgs.c
    ${SUPER_LU_DIR}/SRC/cpivotgrowth.c
    ${SUPER_LU_DIR}/SRC/cgsrfs.c
    ${SUPER_LU_DIR}/SRC/cgstrf.c
    ${SUPER_LU_DIR}/SRC/cgstrs.c
    ${SUPER_LU_DIR}/SRC/ccopy_to_ucol.c
    ${SUPER_LU_DIR}/SRC/csnode_dfs.c
    ${SUPER_LU_DIR}/SRC/csnode_bmod.c
    ${SUPER_LU_DIR}/SRC/cpanel_dfs.c
    ${SUPER_LU_DIR}/SRC/cpanel_bmod.c
    ${SUPER_LU_DIR}/SRC/creadhb.c
    ${SUPER_LU_DIR}/SRC/creadrb.c
    ${SUPER_LU_DIR}/SRC/creadtriple.c
    ${SUPER_LU_DIR}/SRC/ccolumn_dfs.c
    ${SUPER_LU_DIR}/SRC/ccolumn_bmod.c
    ${SUPER_LU_DIR}/SRC/cpivotL.c
    ${SUPER_LU_DIR}/SRC/cpruneL.c
    ${SUPER_LU_DIR}/SRC/cmemory.c
    ${SUPER_LU_DIR}/SRC/cutil.c
    ${SUPER_LU_DIR}/SRC/cmyblas2.c
    ${SUPER_LU_DIR}/SRC/cgsisx.c
    ${SUPER_LU_DIR}/SRC/cgsitrf.c
    ${SUPER_LU_DIR}/SRC/cldperm.c
    ${SUPER_LU_DIR}/SRC/ilu_cdrop_row.c
    ${SUPER_LU_DIR}/SRC/ilu_csnode_dfs.c
    ${SUPER_LU_DIR}/SRC/ilu_ccolumn_dfs.c
    ${SUPER_LU_DIR}/SRC/ilu_cpanel_dfs.c
    ${SUPER_LU_DIR}/SRC/ilu_ccopy_to_ucol.c
    ${SUPER_LU_DIR}/SRC/ilu_cpivotL.c
    ${SUPER_LU_DIR}/SRC/cdiagonal.c
  )
endif(enable_complex)

if(enable_complex16)
  list(APPEND SUPER_LU_HEADERS
    ${SUPER_LU_DIR}/SRC/slu_zdefs.h
  )

  list(APPEND SUPER_LU_SOURCES
    ${SUPER_LU_DIR}/SRC/zlacon2.c
    ${SUPER_LU_DIR}/SRC/dzsum1.c
    ${SUPER_LU_DIR}/SRC/izmax1.c
    ${SUPER_LU_DIR}/SRC/dcomplex.c
    ${SUPER_LU_DIR}/SRC/zgssv.c
    ${SUPER_LU_DIR}/SRC/zgssvx.c
    ${SUPER_LU_DIR}/SRC/zsp_blas2.c
    ${SUPER_LU_DIR}/SRC/zsp_blas3.c
    ${SUPER_LU_DIR}/SRC/zgscon.c
    ${SUPER_LU_DIR}/SRC/zlangs.c
    ${SUPER_LU_DIR}/SRC/zgsequ.c
    ${SUPER_LU_DIR}/SRC/zlaqgs.c
    ${SUPER_LU_DIR}/SRC/zpivotgrowth.c
    ${SUPER_LU_DIR}/SRC/zgsrfs.c
    ${SUPER_LU_DIR}/SRC/zgstrf.c
    ${SUPER_LU_DIR}/SRC/zgstrs.c
    ${SUPER_LU_DIR}/SRC/zcopy_to_ucol.c
    ${SUPER_LU_DIR}/SRC/zsnode_dfs.c
    ${SUPER_LU_DIR}/SRC/zsnode_bmod.c
    ${SUPER_LU_DIR}/SRC/zpanel_dfs.c
    ${SUPER_LU_DIR}/SRC/zpanel_bmod.c
    ${SUPER_LU_DIR}/SRC/zreadhb.c
    ${SUPER_LU_DIR}/SRC/zreadrb.c
    ${SUPER_LU_DIR}/SRC/zreadtriple.c
    ${SUPER_LU_DIR}/SRC/zcolumn_dfs.c
    ${SUPER_LU_DIR}/SRC/zcolumn_bmod.c
    ${SUPER_LU_DIR}/SRC/zpivotL.c
    ${SUPER_LU_DIR}/SRC/zpruneL.c
    ${SUPER_LU_DIR}/SRC/zmemory.c
    ${SUPER_LU_DIR}/SRC/zutil.c
    ${SUPER_LU_DIR}/SRC/zmyblas2.c
    ${SUPER_LU_DIR}/SRC/zgsisx.c
    ${SUPER_LU_DIR}/SRC/zgsitrf.c
    ${SUPER_LU_DIR}/SRC/zldperm.c
    ${SUPER_LU_DIR}/SRC/ilu_zdrop_row.c
    ${SUPER_LU_DIR}/SRC/ilu_zsnode_dfs.c
    ${SUPER_LU_DIR}/SRC/ilu_zcolumn_dfs.c
    ${SUPER_LU_DIR}/SRC/ilu_zpanel_dfs.c
    ${SUPER_LU_DIR}/SRC/ilu_zcopy_to_ucol.c
    ${SUPER_LU_DIR}/SRC/ilu_zpivotL.c
    ${SUPER_LU_DIR}/SRC/zdiagonal.c
  )
endif(enable_complex16)

add_definitions( -DUSE_VENDOR_BLAS )
#add_definitions( -D_LONGINT )
add_definitions( -DAdd_ )
add_definitions( -w )

add_library(super_lu STATIC ${SUPER_LU_SOURCES} ${SUPER_LU_HEADERS})
target_link_libraries(super_lu ${BLAS_LIBRARIES} )
