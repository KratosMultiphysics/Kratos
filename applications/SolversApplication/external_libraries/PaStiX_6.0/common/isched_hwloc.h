/**
 *
 * @file isched_hwloc.h
 *
 * @copyright 2008-2018 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 * @copyright 2010-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 *
 * PaStiX thread binding routines with hwloc header.
 *
 * @version 6.0.1
 * @author Mathieu Faverge
 * @date 2018-07-16
 *
 **/
#ifndef _isched_hwloc_h_
#define _isched_hwloc_h_

/**
 * This file is available only if hwloc is available
 */
#if !defined(HAVE_HWLOC)
#error "This file should not be included if HwLoc is not available"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <hwloc.h>

BEGIN_C_DECLS

typedef struct {
    int lvl;
    int processor_id;
    int master_id;
    int id1;
    int id2;
    int set;
} hwloc_info;

/**
 * Find the master for the processor_id at n level
 */
extern int isched_hwloc_master_id( int level, int processor_id );

/**
 * Find the number of core for master_id at n level
 */
extern unsigned int isched_hwloc_nb_cores( int level, int master_id );

/**
 * Find the number of level from the computer architecture
 */
extern int isched_hwloc_nb_levels( void );

/**
 * Find the cache size for master at n level
 */
extern size_t isched_hwloc_cache_size( unsigned int level, int master_id );

/**
 * Find the distance between id1 and id2
 */
extern int isched_hwloc_distance( int id1, int id2 );

/**
 * load the HWLOC topology.
 */
extern int isched_hwloc_init(void);

/**
 * unload the HWLOC topology.
 */
extern int isched_hwloc_fini(void);

/**
 * Find the number of core of the architecture.
 */
extern int isched_hwloc_nb_real_cores();

/**
 * Bind the current thread on the core of index cpu_index.
 */
int isched_hwloc_bind_on_core_index(int cpu_index);

/**
 * Unbind the current thread.
 */
int isched_hwloc_unbind();

/**
 * Return the logical socket index for a core index (hwloc numbering).
 */
int isched_hwloc_socket_id(int core_id);

/**
 * Return the logical NUMA node index for a core index (hwloc numbering).
 */
int isched_hwloc_numa_id(int core_id);

/**
 * Return the depth of the first core hardware ancestor: NUMA node or socket.
 */
int isched_hwloc_core_first_hrwd_ancestor_depth();

/**
 * Return the number of hwloc objects at the "level" depth.
 */
int isched_hwloc_get_nb_objects(int level);

/**
 * Return the number of hwloc objects at the "level" depth.
 */
int isched_hwloc_get_nb_objects(int level);

/**
 * Find the number of core under the object number index at the topology depth
 * level.
 */
unsigned int isched_hwloc_nb_cores_per_obj( hwloc_obj_type_t level, int index );

/**
 * Return the number of thread on the machine.
 */
unsigned int isched_world_size();

/**
 * Bind the current thread according the mask of index mask_index.
 */
int isched_hwloc_bind_on_mask_index(hwloc_cpuset_t mask_index);

/**
 * Allow serial thread binding per core to use the SMT/HT capabilities of the
 * processor
 */
int isched_hwloc_allow_ht(int htnb);
int isched_hwloc_get_ht(void);

END_C_DECLS

#endif /* _isched_hwloc_h_ */
