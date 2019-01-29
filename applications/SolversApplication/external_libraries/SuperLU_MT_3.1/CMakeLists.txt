project (SUPERLU_MT C)

include_directories( ${CMAKE_CURRENT_SOURCE_DIR}/SRC )

# setup options
option(enable_single    "Enable single precision library" ON)
option(enable_double    "Enable double precision library" ON)
option(enable_complex   "Enable complex precision library" ON)
option(enable_complex16 "Enable complex16 precision library" ON)

set(SUPERLU_HEADERS
  ${CMAKE_CURRENT_SOURCE_DIR}/SRC/supermatrix.h
  ${CMAKE_CURRENT_SOURCE_DIR}/SRC/slu_dcomplex.h
  ${CMAKE_CURRENT_SOURCE_DIR}/SRC/slu_scomplex.h
  ${CMAKE_CURRENT_SOURCE_DIR}/SRC/slu_mt_Cnames.h
  ${CMAKE_CURRENT_SOURCE_DIR}/SRC/slu_mt_util.h
  ${CMAKE_CURRENT_SOURCE_DIR}/SRC/slu_mt_machines.h
  )

set(SUPERLU_SOURCES
  ${CMAKE_CURRENT_SOURCE_DIR}/SRC/superlu_timer.c
  ${CMAKE_CURRENT_SOURCE_DIR}/SRC/dclock.c
  ${CMAKE_CURRENT_SOURCE_DIR}/SRC/sp_ienv.c
  ${CMAKE_CURRENT_SOURCE_DIR}/SRC/lsame.c
  ${CMAKE_CURRENT_SOURCE_DIR}/SRC/xerbla.c
  ${CMAKE_CURRENT_SOURCE_DIR}/SRC/util.c
  ${CMAKE_CURRENT_SOURCE_DIR}/SRC/pmemory.c
  ${CMAKE_CURRENT_SOURCE_DIR}/SRC/qrnzcnt.c
  ${CMAKE_CURRENT_SOURCE_DIR}/SRC/cholnzcnt.c
  ${CMAKE_CURRENT_SOURCE_DIR}/SRC/await.c
  ${CMAKE_CURRENT_SOURCE_DIR}/SRC/get_perm_c.c
  ${CMAKE_CURRENT_SOURCE_DIR}/SRC/mmd.c
  ${CMAKE_CURRENT_SOURCE_DIR}/SRC/colamd.c
  ${CMAKE_CURRENT_SOURCE_DIR}/SRC/sp_coletree.c
  ${CMAKE_CURRENT_SOURCE_DIR}/SRC/pxgstrf_scheduler.c
  ${CMAKE_CURRENT_SOURCE_DIR}/SRC/sp_colorder.c
  ${CMAKE_CURRENT_SOURCE_DIR}/SRC/pxgstrf_mark_busy_descends.c
  ${CMAKE_CURRENT_SOURCE_DIR}/SRC/pxgstrf_pruneL.c
  ${CMAKE_CURRENT_SOURCE_DIR}/SRC/pxgstrf_super_bnd_dfs.c
  ${CMAKE_CURRENT_SOURCE_DIR}/SRC/pxgstrf_relax_snode.c
  ${CMAKE_CURRENT_SOURCE_DIR}/SRC/heap_relax_snode.c
  ${CMAKE_CURRENT_SOURCE_DIR}/SRC/pxgstrf_synch.c
  ${CMAKE_CURRENT_SOURCE_DIR}/SRC/pxgstrf_finalize.c
  )
set_source_files_properties(superlu_timer.c PROPERTIES COMPILE_FLAGS -O0)
set_source_files_properties(dclock.c PROPERTIES COMPILE_FLAGS -O0)
set_source_files_properties(await.c PROPERTIES COMPILE_FLAGS -O0)

if(enable_single)
  list(APPEND SUPERLU_HEADERS
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/slu_mt_sdefs.h
    )

  list(APPEND SUPERLU_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/sreadhb.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/sreadrb.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/smatgen.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/psgssv.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/psgssvx.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/sgstrs.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/sgsrfs.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/sgscon.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/slacon.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/slangs.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/sgsequ.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/slaqgs.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/spivotgrowth.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/psmemory.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/psutil.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/smyblas2.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/psgstrf.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/psgstrf_init.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/psgstrf_thread.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/psgstrf_thread_init.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/psgstrf_thread_finalize.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/psgstrf_factor_snode.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/psgstrf_snode_dfs.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/psgstrf_snode_bmod.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/psgstrf_panel_dfs.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/psgstrf_panel_bmod.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/psgstrf_copy_to_ucol.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/psgstrf_pivotL.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/psgstrf_column_dfs.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/psgstrf_column_bmod.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/psgstrf_bmod1D.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/psgstrf_bmod2D.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/psgstrf_bmod1D_mv2.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/psgstrf_bmod2D_mv2.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/ssp_blas2.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/ssp_blas3.c
    )
endif(enable_single)

if(enable_double)
  list(APPEND SUPERLU_HEADERS
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/slu_mt_ddefs.h
    )

  list(APPEND SUPERLU_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/dreadhb.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/dreadrb.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/dmatgen.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/pdgssv.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/pdgssvx.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/dgstrs.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/dgsrfs.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/dgscon.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/dlacon.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/dlangs.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/dgsequ.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/dlaqgs.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/dpivotgrowth.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/pdmemory.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/pdutil.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/dmyblas2.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/pdgstrf.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/pdgstrf_init.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/pdgstrf_thread.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/pdgstrf_thread_init.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/pdgstrf_thread_finalize.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/pdgstrf_factor_snode.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/pdgstrf_snode_dfs.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/pdgstrf_snode_bmod.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/pdgstrf_panel_dfs.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/pdgstrf_panel_bmod.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/pdgstrf_copy_to_ucol.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/pdgstrf_pivotL.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/pdgstrf_column_dfs.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/pdgstrf_column_bmod.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/pdgstrf_bmod1D.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/pdgstrf_bmod2D.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/pdgstrf_bmod1D_mv2.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/pdgstrf_bmod2D_mv2.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/dsp_blas2.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/dsp_blas3.c
    )
  set_source_files_properties(dmach.c PROPERTIES COMPILE_FLAGS -O0)
endif(enable_double)

if(enable_complex)
  list(APPEND SUPERLU_HEADERS
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/slu_mt_cdefs.h
    )

  list(APPEND SUPERLU_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/scomplex.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/creadhb.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/creadrb.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/cmatgen.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/scsum1.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/icmax1.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/pcgssv.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/pcgssvx.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/cgstrs.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/cgsrfs.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/cgscon.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/clacon.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/clangs.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/cgsequ.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/claqgs.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/cpivotgrowth.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/pcmemory.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/pcutil.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/cmyblas2.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/pcgstrf.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/pcgstrf_init.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/pcgstrf_thread.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/pcgstrf_thread_init.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/pcgstrf_thread_finalize.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/pcgstrf_factor_snode.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/pcgstrf_snode_dfs.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/pcgstrf_snode_bmod.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/pcgstrf_panel_dfs.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/pcgstrf_panel_bmod.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/pcgstrf_copy_to_ucol.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/pcgstrf_pivotL.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/pcgstrf_column_dfs.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/pcgstrf_column_bmod.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/pcgstrf_bmod1D.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/pcgstrf_bmod2D.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/pcgstrf_bmod1D_mv2.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/pcgstrf_bmod2D_mv2.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/csp_blas2.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/csp_blas3.c
    )
endif(enable_complex)

if(enable_complex16)
  list(APPEND SUPERLU_HEADERS
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/slu_mt_zdefs.h
    )

  list(APPEND SUPERLU_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/dcomplex.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/zreadhb.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/zreadrb.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/zmatgen.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/dzsum1.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/izmax1.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/pzgssv.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/pzgssvx.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/zgstrs.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/zgsrfs.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/zgscon.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/zlacon.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/zlangs.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/zgsequ.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/zlaqgs.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/zpivotgrowth.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/pzmemory.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/pzutil.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/zmyblas2.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/pzgstrf.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/pzgstrf_init.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/pzgstrf_thread.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/pzgstrf_thread_init.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/pzgstrf_thread_finalize.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/pzgstrf_factor_snode.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/pzgstrf_snode_dfs.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/pzgstrf_snode_bmod.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/pzgstrf_panel_dfs.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/pzgstrf_panel_bmod.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/pzgstrf_copy_to_ucol.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/pzgstrf_pivotL.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/pzgstrf_column_dfs.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/pzgstrf_column_bmod.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/pzgstrf_bmod1D.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/pzgstrf_bmod2D.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/pzgstrf_bmod1D_mv2.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/pzgstrf_bmod2D_mv2.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/zsp_blas2.c
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/zsp_blas3.c
    )
endif(enable_complex16)

if(enable_single OR enable_complex)
  list(APPEND SUPERLU_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/slamch.c
    )
  set_source_files_properties(smach.c PROPERTIES COMPILE_FLAGS -O0)
endif(enable_single OR enable_complex)

if(enable_double OR enable_complex16)
  list(APPEND SUPERLU_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/SRC/dlamch.c
    )
  set_source_files_properties(dmach.c PROPERTIES COMPILE_FLAGS -O0)
endif(enable_double OR enable_complex16)

#add_definitions( -D_LONGINT )
add_definitions( -DNoChange )
add_definitions( -D__OPENMP )
add_definitions( -DPRNTlevel=0 )
add_definitions( -DDEBUGlevel=0 )
add_definitions( -w )

#if(NOT BLAS_FOUND)
#  find_package(BLAS)
#endif(NOT BLAS_FOUND)

#if(BLAS_FOUND)
#  set(CBLAS_LIBRARIES ${BLAS_LIBRARIES})
#  add_definitions( -DUSE_VENDOR_BLAS )
#else(BLAS_FOUND)
#  add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/CBLAS")
#  set(CBLAS_LIBRARIES external_libblas_mt)
#endif(BLAS_FOUND)

add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/CBLAS")
set(CBLAS_LIBRARIES external_libblas_mt)

message(STATUS "SUPER_LU_MT cblas: ${CBLAS_LIBRARIES}")

add_library(external_superlu_mt STATIC ${SUPERLU_SOURCES} ${SUPERLU_HEADERS})
target_link_libraries(external_superlu_mt ${CBLAS_LIBRARIES})
