include_directories( ${SUPER_LU_MT_DIR}/SRC )

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

set(SUPER_LU_HEADERS
  ${SUPER_LU_MT_DIR}/SRC/supermatrix.h
  ${SUPER_LU_MT_DIR}/SRC/slu_dcomplex.h
  ${SUPER_LU_MT_DIR}/SRC/slu_scomplex.h
  ${SUPER_LU_MT_DIR}/SRC/slu_mt_Cnames.h
  ${SUPER_LU_MT_DIR}/SRC/slu_mt_util.h
  ${SUPER_LU_MT_DIR}/SRC/slu_mt_machines.h
  )

set(SUPER_LU_SOURCES
  ${SUPER_LU_MT_DIR}/SRC/superlu_timer.c
  ${SUPER_LU_MT_DIR}/SRC/dclock.c
  ${SUPER_LU_MT_DIR}/SRC/sp_ienv.c
  ${SUPER_LU_MT_DIR}/SRC/lsame.c
  ${SUPER_LU_MT_DIR}/SRC/xerbla.c
  ${SUPER_LU_MT_DIR}/SRC/util.c
  ${SUPER_LU_MT_DIR}/SRC/pmemory.c
  ${SUPER_LU_MT_DIR}/SRC/qrnzcnt.c
  ${SUPER_LU_MT_DIR}/SRC/cholnzcnt.c
  ${SUPER_LU_MT_DIR}/SRC/await.c
  ${SUPER_LU_MT_DIR}/SRC/get_perm_c.c
  ${SUPER_LU_MT_DIR}/SRC/mmd.c
  ${SUPER_LU_MT_DIR}/SRC/colamd.c
  ${SUPER_LU_MT_DIR}/SRC/sp_coletree.c
  ${SUPER_LU_MT_DIR}/SRC/pxgstrf_scheduler.c
  ${SUPER_LU_MT_DIR}/SRC/sp_colorder.c
  ${SUPER_LU_MT_DIR}/SRC/pxgstrf_mark_busy_descends.c
  ${SUPER_LU_MT_DIR}/SRC/pxgstrf_pruneL.c
  ${SUPER_LU_MT_DIR}/SRC/pxgstrf_super_bnd_dfs.c
  ${SUPER_LU_MT_DIR}/SRC/pxgstrf_relax_snode.c
  ${SUPER_LU_MT_DIR}/SRC/heap_relax_snode.c
  ${SUPER_LU_MT_DIR}/SRC/pxgstrf_synch.c
  ${SUPER_LU_MT_DIR}/SRC/pxgstrf_finalize.c
  )
set_source_files_properties(superlu_timer.c PROPERTIES COMPILE_FLAGS -O0)
set_source_files_properties(dclock.c PROPERTIES COMPILE_FLAGS -O0)
set_source_files_properties(await.c PROPERTIES COMPILE_FLAGS -O0)

if(enable_single)
  list(APPEND SUPER_LU_HEADERS
    ${SUPER_LU_MT_DIR}/SRC/slu_mt_sdefs.h
    )

  list(APPEND SUPER_LU_SOURCES
    ${SUPER_LU_MT_DIR}/SRC/sreadhb.c
    ${SUPER_LU_MT_DIR}/SRC/sreadrb.c
    ${SUPER_LU_MT_DIR}/SRC/smatgen.c
    ${SUPER_LU_MT_DIR}/SRC/psgssv.c
    ${SUPER_LU_MT_DIR}/SRC/psgssvx.c
    ${SUPER_LU_MT_DIR}/SRC/sgstrs.c
    ${SUPER_LU_MT_DIR}/SRC/sgsrfs.c
    ${SUPER_LU_MT_DIR}/SRC/sgscon.c
    ${SUPER_LU_MT_DIR}/SRC/slacon.c
    ${SUPER_LU_MT_DIR}/SRC/slangs.c
    ${SUPER_LU_MT_DIR}/SRC/sgsequ.c
    ${SUPER_LU_MT_DIR}/SRC/slaqgs.c
    ${SUPER_LU_MT_DIR}/SRC/spivotgrowth.c
    ${SUPER_LU_MT_DIR}/SRC/psmemory.c
    ${SUPER_LU_MT_DIR}/SRC/psutil.c
    ${SUPER_LU_MT_DIR}/SRC/smyblas2.c
    ${SUPER_LU_MT_DIR}/SRC/psgstrf.c
    ${SUPER_LU_MT_DIR}/SRC/psgstrf_init.c
    ${SUPER_LU_MT_DIR}/SRC/psgstrf_thread.c
    ${SUPER_LU_MT_DIR}/SRC/psgstrf_thread_init.c
    ${SUPER_LU_MT_DIR}/SRC/psgstrf_thread_finalize.c
    ${SUPER_LU_MT_DIR}/SRC/psgstrf_factor_snode.c
    ${SUPER_LU_MT_DIR}/SRC/psgstrf_snode_dfs.c
    ${SUPER_LU_MT_DIR}/SRC/psgstrf_snode_bmod.c
    ${SUPER_LU_MT_DIR}/SRC/psgstrf_panel_dfs.c
    ${SUPER_LU_MT_DIR}/SRC/psgstrf_panel_bmod.c
    ${SUPER_LU_MT_DIR}/SRC/psgstrf_copy_to_ucol.c
    ${SUPER_LU_MT_DIR}/SRC/psgstrf_pivotL.c
    ${SUPER_LU_MT_DIR}/SRC/psgstrf_column_dfs.c
    ${SUPER_LU_MT_DIR}/SRC/psgstrf_column_bmod.c
    ${SUPER_LU_MT_DIR}/SRC/psgstrf_bmod1D.c
    ${SUPER_LU_MT_DIR}/SRC/psgstrf_bmod2D.c
    ${SUPER_LU_MT_DIR}/SRC/psgstrf_bmod1D_mv2.c
    ${SUPER_LU_MT_DIR}/SRC/psgstrf_bmod2D_mv2.c
    ${SUPER_LU_MT_DIR}/SRC/ssp_blas2.c
    ${SUPER_LU_MT_DIR}/SRC/ssp_blas3.c
    )
endif(enable_single)

if(enable_double)
  list(APPEND SUPER_LU_HEADERS
    ${SUPER_LU_MT_DIR}/SRC/slu_mt_ddefs.h
    )

  list(APPEND SUPER_LU_SOURCES
    ${SUPER_LU_MT_DIR}/SRC/dreadhb.c
    ${SUPER_LU_MT_DIR}/SRC/dreadrb.c
    ${SUPER_LU_MT_DIR}/SRC/dmatgen.c
    ${SUPER_LU_MT_DIR}/SRC/pdgssv.c
    ${SUPER_LU_MT_DIR}/SRC/pdgssvx.c
    ${SUPER_LU_MT_DIR}/SRC/dgstrs.c
    ${SUPER_LU_MT_DIR}/SRC/dgsrfs.c
    ${SUPER_LU_MT_DIR}/SRC/dgscon.c
    ${SUPER_LU_MT_DIR}/SRC/dlacon.c
    ${SUPER_LU_MT_DIR}/SRC/dlangs.c
    ${SUPER_LU_MT_DIR}/SRC/dgsequ.c
    ${SUPER_LU_MT_DIR}/SRC/dlaqgs.c
    ${SUPER_LU_MT_DIR}/SRC/dpivotgrowth.c
    ${SUPER_LU_MT_DIR}/SRC/pdmemory.c
    ${SUPER_LU_MT_DIR}/SRC/pdutil.c
    ${SUPER_LU_MT_DIR}/SRC/dmyblas2.c
    ${SUPER_LU_MT_DIR}/SRC/pdgstrf.c
    ${SUPER_LU_MT_DIR}/SRC/pdgstrf_init.c
    ${SUPER_LU_MT_DIR}/SRC/pdgstrf_thread.c
    ${SUPER_LU_MT_DIR}/SRC/pdgstrf_thread_init.c
    ${SUPER_LU_MT_DIR}/SRC/pdgstrf_thread_finalize.c
    ${SUPER_LU_MT_DIR}/SRC/pdgstrf_factor_snode.c
    ${SUPER_LU_MT_DIR}/SRC/pdgstrf_snode_dfs.c
    ${SUPER_LU_MT_DIR}/SRC/pdgstrf_snode_bmod.c
    ${SUPER_LU_MT_DIR}/SRC/pdgstrf_panel_dfs.c
    ${SUPER_LU_MT_DIR}/SRC/pdgstrf_panel_bmod.c
    ${SUPER_LU_MT_DIR}/SRC/pdgstrf_copy_to_ucol.c
    ${SUPER_LU_MT_DIR}/SRC/pdgstrf_pivotL.c
    ${SUPER_LU_MT_DIR}/SRC/pdgstrf_column_dfs.c
    ${SUPER_LU_MT_DIR}/SRC/pdgstrf_column_bmod.c
    ${SUPER_LU_MT_DIR}/SRC/pdgstrf_bmod1D.c
    ${SUPER_LU_MT_DIR}/SRC/pdgstrf_bmod2D.c
    ${SUPER_LU_MT_DIR}/SRC/pdgstrf_bmod1D_mv2.c
    ${SUPER_LU_MT_DIR}/SRC/pdgstrf_bmod2D_mv2.c
    ${SUPER_LU_MT_DIR}/SRC/dsp_blas2.c
    ${SUPER_LU_MT_DIR}/SRC/dsp_blas3.c
    )
  set_source_files_properties(dmach.c PROPERTIES COMPILE_FLAGS -O0)
endif(enable_double)

if(enable_complex)
  list(APPEND SUPER_LU_HEADERS
    ${SUPER_LU_MT_DIR}/SRC/slu_mt_cdefs.h
    )

  list(APPEND SUPER_LU_SOURCES
    ${SUPER_LU_MT_DIR}/SRC/scomplex.c
    ${SUPER_LU_MT_DIR}/SRC/creadhb.c
    ${SUPER_LU_MT_DIR}/SRC/creadrb.c
    ${SUPER_LU_MT_DIR}/SRC/cmatgen.c
    ${SUPER_LU_MT_DIR}/SRC/scsum1.c
    ${SUPER_LU_MT_DIR}/SRC/icmax1.c
    ${SUPER_LU_MT_DIR}/SRC/pcgssv.c
    ${SUPER_LU_MT_DIR}/SRC/pcgssvx.c
    ${SUPER_LU_MT_DIR}/SRC/cgstrs.c
    ${SUPER_LU_MT_DIR}/SRC/cgsrfs.c
    ${SUPER_LU_MT_DIR}/SRC/cgscon.c
    ${SUPER_LU_MT_DIR}/SRC/clacon.c
    ${SUPER_LU_MT_DIR}/SRC/clangs.c
    ${SUPER_LU_MT_DIR}/SRC/cgsequ.c
    ${SUPER_LU_MT_DIR}/SRC/claqgs.c
    ${SUPER_LU_MT_DIR}/SRC/cpivotgrowth.c
    ${SUPER_LU_MT_DIR}/SRC/pcmemory.c
    ${SUPER_LU_MT_DIR}/SRC/pcutil.c
    ${SUPER_LU_MT_DIR}/SRC/cmyblas2.c
    ${SUPER_LU_MT_DIR}/SRC/pcgstrf.c
    ${SUPER_LU_MT_DIR}/SRC/pcgstrf_init.c
    ${SUPER_LU_MT_DIR}/SRC/pcgstrf_thread.c
    ${SUPER_LU_MT_DIR}/SRC/pcgstrf_thread_init.c
    ${SUPER_LU_MT_DIR}/SRC/pcgstrf_thread_finalize.c
    ${SUPER_LU_MT_DIR}/SRC/pcgstrf_factor_snode.c
    ${SUPER_LU_MT_DIR}/SRC/pcgstrf_snode_dfs.c
    ${SUPER_LU_MT_DIR}/SRC/pcgstrf_snode_bmod.c
    ${SUPER_LU_MT_DIR}/SRC/pcgstrf_panel_dfs.c
    ${SUPER_LU_MT_DIR}/SRC/pcgstrf_panel_bmod.c
    ${SUPER_LU_MT_DIR}/SRC/pcgstrf_copy_to_ucol.c
    ${SUPER_LU_MT_DIR}/SRC/pcgstrf_pivotL.c
    ${SUPER_LU_MT_DIR}/SRC/pcgstrf_column_dfs.c
    ${SUPER_LU_MT_DIR}/SRC/pcgstrf_column_bmod.c
    ${SUPER_LU_MT_DIR}/SRC/pcgstrf_bmod1D.c
    ${SUPER_LU_MT_DIR}/SRC/pcgstrf_bmod2D.c
    ${SUPER_LU_MT_DIR}/SRC/pcgstrf_bmod1D_mv2.c
    ${SUPER_LU_MT_DIR}/SRC/pcgstrf_bmod2D_mv2.c
    ${SUPER_LU_MT_DIR}/SRC/csp_blas2.c
    ${SUPER_LU_MT_DIR}/SRC/csp_blas3.c
    )
endif(enable_complex)

if(enable_complex16)
  list(APPEND SUPER_LU_HEADERS
    ${SUPER_LU_MT_DIR}/SRC/slu_mt_zdefs.h
    )

  list(APPEND SUPER_LU_SOURCES
    ${SUPER_LU_MT_DIR}/SRC/dcomplex.c
    ${SUPER_LU_MT_DIR}/SRC/zreadhb.c
    ${SUPER_LU_MT_DIR}/SRC/zreadrb.c
    ${SUPER_LU_MT_DIR}/SRC/zmatgen.c
    ${SUPER_LU_MT_DIR}/SRC/dzsum1.c
    ${SUPER_LU_MT_DIR}/SRC/izmax1.c
    ${SUPER_LU_MT_DIR}/SRC/pzgssv.c
    ${SUPER_LU_MT_DIR}/SRC/pzgssvx.c
    ${SUPER_LU_MT_DIR}/SRC/zgstrs.c
    ${SUPER_LU_MT_DIR}/SRC/zgsrfs.c
    ${SUPER_LU_MT_DIR}/SRC/zgscon.c
    ${SUPER_LU_MT_DIR}/SRC/zlacon.c
    ${SUPER_LU_MT_DIR}/SRC/zlangs.c
    ${SUPER_LU_MT_DIR}/SRC/zgsequ.c
    ${SUPER_LU_MT_DIR}/SRC/zlaqgs.c
    ${SUPER_LU_MT_DIR}/SRC/zpivotgrowth.c
    ${SUPER_LU_MT_DIR}/SRC/pzmemory.c
    ${SUPER_LU_MT_DIR}/SRC/pzutil.c
    ${SUPER_LU_MT_DIR}/SRC/zmyblas2.c
    ${SUPER_LU_MT_DIR}/SRC/pzgstrf.c
    ${SUPER_LU_MT_DIR}/SRC/pzgstrf_init.c
    ${SUPER_LU_MT_DIR}/SRC/pzgstrf_thread.c
    ${SUPER_LU_MT_DIR}/SRC/pzgstrf_thread_init.c
    ${SUPER_LU_MT_DIR}/SRC/pzgstrf_thread_finalize.c
    ${SUPER_LU_MT_DIR}/SRC/pzgstrf_factor_snode.c
    ${SUPER_LU_MT_DIR}/SRC/pzgstrf_snode_dfs.c
    ${SUPER_LU_MT_DIR}/SRC/pzgstrf_snode_bmod.c
    ${SUPER_LU_MT_DIR}/SRC/pzgstrf_panel_dfs.c
    ${SUPER_LU_MT_DIR}/SRC/pzgstrf_panel_bmod.c
    ${SUPER_LU_MT_DIR}/SRC/pzgstrf_copy_to_ucol.c
    ${SUPER_LU_MT_DIR}/SRC/pzgstrf_pivotL.c
    ${SUPER_LU_MT_DIR}/SRC/pzgstrf_column_dfs.c
    ${SUPER_LU_MT_DIR}/SRC/pzgstrf_column_bmod.c
    ${SUPER_LU_MT_DIR}/SRC/pzgstrf_bmod1D.c
    ${SUPER_LU_MT_DIR}/SRC/pzgstrf_bmod2D.c
    ${SUPER_LU_MT_DIR}/SRC/pzgstrf_bmod1D_mv2.c
    ${SUPER_LU_MT_DIR}/SRC/pzgstrf_bmod2D_mv2.c
    ${SUPER_LU_MT_DIR}/SRC/zsp_blas2.c
    ${SUPER_LU_MT_DIR}/SRC/zsp_blas3.c
    )
endif(enable_complex16)

if(enable_single OR enable_complex)
  list(APPEND SUPER_LU_SOURCES
    ${SUPER_LU_MT_DIR}/SRC/slamch.c
    )
  set_source_files_properties(smach.c PROPERTIES COMPILE_FLAGS -O0)
endif(enable_single OR enable_complex)

if(enable_double OR enable_complex16)
  list(APPEND SUPER_LU_SOURCES
    ${SUPER_LU_MT_DIR}/SRC/dlamch.c
    )
  set_source_files_properties(dmach.c PROPERTIES COMPILE_FLAGS -O0)
endif(enable_double OR enable_complex16)

message(STATUS "****compiling super_lu_mt*****")

#add_definitions( -DUSE_VENDOR_BLAS )
#add_definitions( -D_LONGINT )

add_library(super_lu_mt STATIC ${SUPER_LU_SOURCES} ${SUPER_LU_HEADERS})
#target_link_libraries(super_lu_mt ${BLAS_LIBRARIES})
target_link_libraries(super_lu_mt libblas_mt)
