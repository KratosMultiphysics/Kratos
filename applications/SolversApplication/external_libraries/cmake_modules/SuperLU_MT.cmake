include_directories( ${SUPERLU_MT_DIR}/SRC )

# setup options
option(enable_single    "Enable single precision library" ON)
option(enable_double    "Enable double precision library" ON)
option(enable_complex   "Enable complex precision library" ON)
option(enable_complex16 "Enable complex16 precision library" ON)

add_definitions( -DNoChange )
add_definitions( -D__OPENMP )
add_definitions( -DPRNTlevel=0 )
add_definitions( -DDEBUGlevel=0 )
add_definitions( -w )

message(STATUS "****compiling cblas*****")
INCLUDE("${CMAKE_CURRENT_SOURCE_DIR}/external_libraries/cmake_modules/Cblas_MT.cmake")

set(SUPERLU_HEADERS
  ${SUPERLU_MT_DIR}/SRC/supermatrix.h
  ${SUPERLU_MT_DIR}/SRC/slu_dcomplex.h
  ${SUPERLU_MT_DIR}/SRC/slu_scomplex.h
  ${SUPERLU_MT_DIR}/SRC/slu_mt_Cnames.h
  ${SUPERLU_MT_DIR}/SRC/slu_mt_util.h
  ${SUPERLU_MT_DIR}/SRC/slu_mt_machines.h
  )

set(SUPERLU_SOURCES
  ${SUPERLU_MT_DIR}/SRC/superlu_timer.c
  ${SUPERLU_MT_DIR}/SRC/dclock.c
  ${SUPERLU_MT_DIR}/SRC/sp_ienv.c
  ${SUPERLU_MT_DIR}/SRC/lsame.c
  ${SUPERLU_MT_DIR}/SRC/xerbla.c
  ${SUPERLU_MT_DIR}/SRC/util.c
  ${SUPERLU_MT_DIR}/SRC/pmemory.c
  ${SUPERLU_MT_DIR}/SRC/qrnzcnt.c
  ${SUPERLU_MT_DIR}/SRC/cholnzcnt.c
  ${SUPERLU_MT_DIR}/SRC/await.c
  ${SUPERLU_MT_DIR}/SRC/get_perm_c.c
  ${SUPERLU_MT_DIR}/SRC/mmd.c
  ${SUPERLU_MT_DIR}/SRC/colamd.c
  ${SUPERLU_MT_DIR}/SRC/sp_coletree.c
  ${SUPERLU_MT_DIR}/SRC/pxgstrf_scheduler.c
  ${SUPERLU_MT_DIR}/SRC/sp_colorder.c
  ${SUPERLU_MT_DIR}/SRC/pxgstrf_mark_busy_descends.c
  ${SUPERLU_MT_DIR}/SRC/pxgstrf_pruneL.c
  ${SUPERLU_MT_DIR}/SRC/pxgstrf_super_bnd_dfs.c
  ${SUPERLU_MT_DIR}/SRC/pxgstrf_relax_snode.c
  ${SUPERLU_MT_DIR}/SRC/heap_relax_snode.c
  ${SUPERLU_MT_DIR}/SRC/pxgstrf_synch.c
  ${SUPERLU_MT_DIR}/SRC/pxgstrf_finalize.c
  )
set_source_files_properties(superlu_timer.c PROPERTIES COMPILE_FLAGS -O0)
set_source_files_properties(dclock.c PROPERTIES COMPILE_FLAGS -O0)
set_source_files_properties(await.c PROPERTIES COMPILE_FLAGS -O0)

if(enable_single)
  list(APPEND SUPERLU_HEADERS
    ${SUPERLU_MT_DIR}/SRC/slu_mt_sdefs.h
    )

  list(APPEND SUPERLU_SOURCES
    ${SUPERLU_MT_DIR}/SRC/sreadhb.c
    ${SUPERLU_MT_DIR}/SRC/sreadrb.c
    ${SUPERLU_MT_DIR}/SRC/smatgen.c
    ${SUPERLU_MT_DIR}/SRC/psgssv.c
    ${SUPERLU_MT_DIR}/SRC/psgssvx.c
    ${SUPERLU_MT_DIR}/SRC/sgstrs.c
    ${SUPERLU_MT_DIR}/SRC/sgsrfs.c
    ${SUPERLU_MT_DIR}/SRC/sgscon.c
    ${SUPERLU_MT_DIR}/SRC/slacon.c
    ${SUPERLU_MT_DIR}/SRC/slangs.c
    ${SUPERLU_MT_DIR}/SRC/sgsequ.c
    ${SUPERLU_MT_DIR}/SRC/slaqgs.c
    ${SUPERLU_MT_DIR}/SRC/spivotgrowth.c
    ${SUPERLU_MT_DIR}/SRC/psmemory.c
    ${SUPERLU_MT_DIR}/SRC/psutil.c
    ${SUPERLU_MT_DIR}/SRC/smyblas2.c
    ${SUPERLU_MT_DIR}/SRC/psgstrf.c
    ${SUPERLU_MT_DIR}/SRC/psgstrf_init.c
    ${SUPERLU_MT_DIR}/SRC/psgstrf_thread.c
    ${SUPERLU_MT_DIR}/SRC/psgstrf_thread_init.c
    ${SUPERLU_MT_DIR}/SRC/psgstrf_thread_finalize.c
    ${SUPERLU_MT_DIR}/SRC/psgstrf_factor_snode.c
    ${SUPERLU_MT_DIR}/SRC/psgstrf_snode_dfs.c
    ${SUPERLU_MT_DIR}/SRC/psgstrf_snode_bmod.c
    ${SUPERLU_MT_DIR}/SRC/psgstrf_panel_dfs.c
    ${SUPERLU_MT_DIR}/SRC/psgstrf_panel_bmod.c
    ${SUPERLU_MT_DIR}/SRC/psgstrf_copy_to_ucol.c
    ${SUPERLU_MT_DIR}/SRC/psgstrf_pivotL.c
    ${SUPERLU_MT_DIR}/SRC/psgstrf_column_dfs.c
    ${SUPERLU_MT_DIR}/SRC/psgstrf_column_bmod.c
    ${SUPERLU_MT_DIR}/SRC/psgstrf_bmod1D.c
    ${SUPERLU_MT_DIR}/SRC/psgstrf_bmod2D.c
    ${SUPERLU_MT_DIR}/SRC/psgstrf_bmod1D_mv2.c
    ${SUPERLU_MT_DIR}/SRC/psgstrf_bmod2D_mv2.c
    ${SUPERLU_MT_DIR}/SRC/ssp_blas2.c
    ${SUPERLU_MT_DIR}/SRC/ssp_blas3.c
    )
endif(enable_single)

if(enable_double)
  list(APPEND SUPERLU_HEADERS
    ${SUPERLU_MT_DIR}/SRC/slu_mt_ddefs.h
    )

  list(APPEND SUPERLU_SOURCES
    ${SUPERLU_MT_DIR}/SRC/dreadhb.c
    ${SUPERLU_MT_DIR}/SRC/dreadrb.c
    ${SUPERLU_MT_DIR}/SRC/dmatgen.c
    ${SUPERLU_MT_DIR}/SRC/pdgssv.c
    ${SUPERLU_MT_DIR}/SRC/pdgssvx.c
    ${SUPERLU_MT_DIR}/SRC/dgstrs.c
    ${SUPERLU_MT_DIR}/SRC/dgsrfs.c
    ${SUPERLU_MT_DIR}/SRC/dgscon.c
    ${SUPERLU_MT_DIR}/SRC/dlacon.c
    ${SUPERLU_MT_DIR}/SRC/dlangs.c
    ${SUPERLU_MT_DIR}/SRC/dgsequ.c
    ${SUPERLU_MT_DIR}/SRC/dlaqgs.c
    ${SUPERLU_MT_DIR}/SRC/dpivotgrowth.c
    ${SUPERLU_MT_DIR}/SRC/pdmemory.c
    ${SUPERLU_MT_DIR}/SRC/pdutil.c
    ${SUPERLU_MT_DIR}/SRC/dmyblas2.c
    ${SUPERLU_MT_DIR}/SRC/pdgstrf.c
    ${SUPERLU_MT_DIR}/SRC/pdgstrf_init.c
    ${SUPERLU_MT_DIR}/SRC/pdgstrf_thread.c
    ${SUPERLU_MT_DIR}/SRC/pdgstrf_thread_init.c
    ${SUPERLU_MT_DIR}/SRC/pdgstrf_thread_finalize.c
    ${SUPERLU_MT_DIR}/SRC/pdgstrf_factor_snode.c
    ${SUPERLU_MT_DIR}/SRC/pdgstrf_snode_dfs.c
    ${SUPERLU_MT_DIR}/SRC/pdgstrf_snode_bmod.c
    ${SUPERLU_MT_DIR}/SRC/pdgstrf_panel_dfs.c
    ${SUPERLU_MT_DIR}/SRC/pdgstrf_panel_bmod.c
    ${SUPERLU_MT_DIR}/SRC/pdgstrf_copy_to_ucol.c
    ${SUPERLU_MT_DIR}/SRC/pdgstrf_pivotL.c
    ${SUPERLU_MT_DIR}/SRC/pdgstrf_column_dfs.c
    ${SUPERLU_MT_DIR}/SRC/pdgstrf_column_bmod.c
    ${SUPERLU_MT_DIR}/SRC/pdgstrf_bmod1D.c
    ${SUPERLU_MT_DIR}/SRC/pdgstrf_bmod2D.c
    ${SUPERLU_MT_DIR}/SRC/pdgstrf_bmod1D_mv2.c
    ${SUPERLU_MT_DIR}/SRC/pdgstrf_bmod2D_mv2.c
    ${SUPERLU_MT_DIR}/SRC/dsp_blas2.c
    ${SUPERLU_MT_DIR}/SRC/dsp_blas3.c
    )
  set_source_files_properties(dmach.c PROPERTIES COMPILE_FLAGS -O0)
endif(enable_double)

if(enable_complex)
  list(APPEND SUPERLU_HEADERS
    ${SUPERLU_MT_DIR}/SRC/slu_mt_cdefs.h
    )

  list(APPEND SUPERLU_SOURCES
    ${SUPERLU_MT_DIR}/SRC/scomplex.c
    ${SUPERLU_MT_DIR}/SRC/creadhb.c
    ${SUPERLU_MT_DIR}/SRC/creadrb.c
    ${SUPERLU_MT_DIR}/SRC/cmatgen.c
    ${SUPERLU_MT_DIR}/SRC/scsum1.c
    ${SUPERLU_MT_DIR}/SRC/icmax1.c
    ${SUPERLU_MT_DIR}/SRC/pcgssv.c
    ${SUPERLU_MT_DIR}/SRC/pcgssvx.c
    ${SUPERLU_MT_DIR}/SRC/cgstrs.c
    ${SUPERLU_MT_DIR}/SRC/cgsrfs.c
    ${SUPERLU_MT_DIR}/SRC/cgscon.c
    ${SUPERLU_MT_DIR}/SRC/clacon.c
    ${SUPERLU_MT_DIR}/SRC/clangs.c
    ${SUPERLU_MT_DIR}/SRC/cgsequ.c
    ${SUPERLU_MT_DIR}/SRC/claqgs.c
    ${SUPERLU_MT_DIR}/SRC/cpivotgrowth.c
    ${SUPERLU_MT_DIR}/SRC/pcmemory.c
    ${SUPERLU_MT_DIR}/SRC/pcutil.c
    ${SUPERLU_MT_DIR}/SRC/cmyblas2.c
    ${SUPERLU_MT_DIR}/SRC/pcgstrf.c
    ${SUPERLU_MT_DIR}/SRC/pcgstrf_init.c
    ${SUPERLU_MT_DIR}/SRC/pcgstrf_thread.c
    ${SUPERLU_MT_DIR}/SRC/pcgstrf_thread_init.c
    ${SUPERLU_MT_DIR}/SRC/pcgstrf_thread_finalize.c
    ${SUPERLU_MT_DIR}/SRC/pcgstrf_factor_snode.c
    ${SUPERLU_MT_DIR}/SRC/pcgstrf_snode_dfs.c
    ${SUPERLU_MT_DIR}/SRC/pcgstrf_snode_bmod.c
    ${SUPERLU_MT_DIR}/SRC/pcgstrf_panel_dfs.c
    ${SUPERLU_MT_DIR}/SRC/pcgstrf_panel_bmod.c
    ${SUPERLU_MT_DIR}/SRC/pcgstrf_copy_to_ucol.c
    ${SUPERLU_MT_DIR}/SRC/pcgstrf_pivotL.c
    ${SUPERLU_MT_DIR}/SRC/pcgstrf_column_dfs.c
    ${SUPERLU_MT_DIR}/SRC/pcgstrf_column_bmod.c
    ${SUPERLU_MT_DIR}/SRC/pcgstrf_bmod1D.c
    ${SUPERLU_MT_DIR}/SRC/pcgstrf_bmod2D.c
    ${SUPERLU_MT_DIR}/SRC/pcgstrf_bmod1D_mv2.c
    ${SUPERLU_MT_DIR}/SRC/pcgstrf_bmod2D_mv2.c
    ${SUPERLU_MT_DIR}/SRC/csp_blas2.c
    ${SUPERLU_MT_DIR}/SRC/csp_blas3.c
    )
endif(enable_complex)

if(enable_complex16)
  list(APPEND SUPERLU_HEADERS
    ${SUPERLU_MT_DIR}/SRC/slu_mt_zdefs.h
    )

  list(APPEND SUPERLU_SOURCES
    ${SUPERLU_MT_DIR}/SRC/dcomplex.c
    ${SUPERLU_MT_DIR}/SRC/zreadhb.c
    ${SUPERLU_MT_DIR}/SRC/zreadrb.c
    ${SUPERLU_MT_DIR}/SRC/zmatgen.c
    ${SUPERLU_MT_DIR}/SRC/dzsum1.c
    ${SUPERLU_MT_DIR}/SRC/izmax1.c
    ${SUPERLU_MT_DIR}/SRC/pzgssv.c
    ${SUPERLU_MT_DIR}/SRC/pzgssvx.c
    ${SUPERLU_MT_DIR}/SRC/zgstrs.c
    ${SUPERLU_MT_DIR}/SRC/zgsrfs.c
    ${SUPERLU_MT_DIR}/SRC/zgscon.c
    ${SUPERLU_MT_DIR}/SRC/zlacon.c
    ${SUPERLU_MT_DIR}/SRC/zlangs.c
    ${SUPERLU_MT_DIR}/SRC/zgsequ.c
    ${SUPERLU_MT_DIR}/SRC/zlaqgs.c
    ${SUPERLU_MT_DIR}/SRC/zpivotgrowth.c
    ${SUPERLU_MT_DIR}/SRC/pzmemory.c
    ${SUPERLU_MT_DIR}/SRC/pzutil.c
    ${SUPERLU_MT_DIR}/SRC/zmyblas2.c
    ${SUPERLU_MT_DIR}/SRC/pzgstrf.c
    ${SUPERLU_MT_DIR}/SRC/pzgstrf_init.c
    ${SUPERLU_MT_DIR}/SRC/pzgstrf_thread.c
    ${SUPERLU_MT_DIR}/SRC/pzgstrf_thread_init.c
    ${SUPERLU_MT_DIR}/SRC/pzgstrf_thread_finalize.c
    ${SUPERLU_MT_DIR}/SRC/pzgstrf_factor_snode.c
    ${SUPERLU_MT_DIR}/SRC/pzgstrf_snode_dfs.c
    ${SUPERLU_MT_DIR}/SRC/pzgstrf_snode_bmod.c
    ${SUPERLU_MT_DIR}/SRC/pzgstrf_panel_dfs.c
    ${SUPERLU_MT_DIR}/SRC/pzgstrf_panel_bmod.c
    ${SUPERLU_MT_DIR}/SRC/pzgstrf_copy_to_ucol.c
    ${SUPERLU_MT_DIR}/SRC/pzgstrf_pivotL.c
    ${SUPERLU_MT_DIR}/SRC/pzgstrf_column_dfs.c
    ${SUPERLU_MT_DIR}/SRC/pzgstrf_column_bmod.c
    ${SUPERLU_MT_DIR}/SRC/pzgstrf_bmod1D.c
    ${SUPERLU_MT_DIR}/SRC/pzgstrf_bmod2D.c
    ${SUPERLU_MT_DIR}/SRC/pzgstrf_bmod1D_mv2.c
    ${SUPERLU_MT_DIR}/SRC/pzgstrf_bmod2D_mv2.c
    ${SUPERLU_MT_DIR}/SRC/zsp_blas2.c
    ${SUPERLU_MT_DIR}/SRC/zsp_blas3.c
    )
endif(enable_complex16)

if(enable_single OR enable_complex)
  list(APPEND SUPERLU_SOURCES
    ${SUPERLU_MT_DIR}/SRC/slamch.c
    )
  set_source_files_properties(smach.c PROPERTIES COMPILE_FLAGS -O0)
endif(enable_single OR enable_complex)

if(enable_double OR enable_complex16)
  list(APPEND SUPERLU_SOURCES
    ${SUPERLU_MT_DIR}/SRC/dlamch.c
    )
  set_source_files_properties(dmach.c PROPERTIES COMPILE_FLAGS -O0)
endif(enable_double OR enable_complex16)

message(STATUS "****compiling super_lu_mt*****")

#add_definitions( -DUSE_VENDOR_BLAS )
#add_definitions( -D_LONGINT )

add_library(external_superlu_mt STATIC ${SUPERLU_SOURCES} ${SUPERLU_HEADERS})
#target_link_libraries(external_superlu_mt ${BLAS_LIBRARIES})
target_link_libraries(external_superlu_mt external_libblas_mt)

set(SUPERLU_LIBRARIES external_superlu_mt)
