set(CMAKE_INCLUDE_CURRENT_DIR ON)

include_directories( ${CMAKE_CURRENT_SOURCE_DIR}/SRC )

## generate variables with the sources
set( LAAUX     
 	${CMAKE_CURRENT_SOURCE_DIR}/SRC/lsame.c
#         ${CMAKE_CURRENT_SOURCE_DIR}/SRC/xerbla.c
     )
##shall we add xerbla.c here?

set( DLASRC     
 	${CMAKE_CURRENT_SOURCE_DIR}/SRC/dlacon.c
     )

set( DZLAUX     
 	${CMAKE_CURRENT_SOURCE_DIR}/SRC/dlamch.c
     )
SET_SOURCE_FILES_PROPERTIES(${DZLAUX} PROPERTIES COMPILE_FLAGS -O0)

set( ALLAUX   
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
     )
SET_SOURCE_FILES_PROPERTIES(superlu_timer.c PROPERTIES COMPILE_FLAGS -O0)

set( DLUSRC  
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
 	${CMAKE_CURRENT_SOURCE_DIR}/EXAMPLE/dfgmr.c
    )

##removing optimizations
SET_SOURCE_FILES_PROPERTIES(${TRIANGLE_SOURCES} PROPERTIES COMPILE_FLAGS -O0)

add_definitions( -DUSE_VENDOR_BLAS )
add_definitions( -DAdd_ )
add_definitions( -D_LONGINT )
add_definitions( -w )

if(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    add_definitions( -DNO_TIMER )
endif(${CMAKE_SYSTEM_NAME} MATCHES "Windows")


###############################################################
add_library(superlu STATIC ${DLUSRC} ${ALLAUX} ${LAAUX} ${DLASRC} ${DZLAUX} )
target_link_libraries(superlu ${BLAS_LIBRARIES} )
install(TARGETS superlu DESTINATION libs )