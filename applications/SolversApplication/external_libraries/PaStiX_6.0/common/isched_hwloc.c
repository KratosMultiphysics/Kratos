/**
 *
 * @file isched_hwloc.c
 *
 * @copyright 2008-2018 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 * @copyright 2010-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 *
 * PaStiX thread binding routines.
 *
 * @version 6.0.1
 * @author Mathieu Faverge
 * @date 2018-07-16
 *
 */
#include "common.h"
#include "isched_hwloc.h"

#if defined(HAVE_HWLOC)

static hwloc_topology_t topology;
static int first_init = 0;
static volatile pastix_atomic_lock_t topo_lock = PASTIX_ATOMIC_UNLOCKED;

#if defined(HAVE_HWLOC_PARENT_MEMBER)
#define HWLOC_GET_PARENT(OBJ)  (OBJ)->parent
#else
#define HWLOC_GET_PARENT(OBJ)  (OBJ)->father
#endif  /* defined(HAVE_HWLOC_PARENT_MEMBER) */

int isched_hwloc_init(void)
{
    pastix_atomic_lock( &topo_lock );
    if ( first_init == 0 ) {
        hwloc_topology_init(&topology);
        hwloc_topology_load(topology);
    }
    first_init++;
    pastix_atomic_unlock( &topo_lock );
    return 0;
}

int isched_hwloc_destroy(void)
{
    pastix_atomic_lock( &topo_lock );
    first_init--;
    if ( first_init == 0 ) {
        hwloc_topology_destroy(topology);
    }
    pastix_atomic_unlock( &topo_lock );
    return 0;
}

unsigned int isched_hwloc_nb_cores_per_obj( hwloc_obj_type_t type, int index )
{
    hwloc_obj_t obj = hwloc_get_obj_by_type(topology, type, index);
    assert( obj != NULL );
    return hwloc_get_nbobjs_inside_cpuset_by_type(topology, obj->cpuset, HWLOC_OBJ_CORE);
}

int isched_hwloc_world_size()
{
    return isched_hwloc_nb_cores_per_obj( HWLOC_OBJ_MACHINE, 0 );
}

int isched_hwloc_bind_on_core_index(int cpu_index)
{
    hwloc_obj_t      core;     /* Hwloc object    */
    hwloc_cpuset_t   cpuset;   /* Hwloc cpuset    */

    /* Get the core of index cpu_index */
    core = hwloc_get_obj_by_type(topology, HWLOC_OBJ_CORE, cpu_index);
    if (!core) {
        fprintf(stderr,
                "isched_hwloc_bind_on_core_index: unable to get the core of index %i (nb physical cores = %i )\n",
                cpu_index, isched_hwloc_world_size());
        return -1;
    }

    /* Get a copy of its cpuset that we may modify.  */
#if !defined(HAVE_HWLOC_BITMAP)
    cpuset = hwloc_cpuset_dup(core->cpuset);
    hwloc_cpuset_singlify(cpuset);
#else
    cpuset = hwloc_bitmap_dup(core->cpuset);
    hwloc_bitmap_singlify(cpuset);
#endif

    /* And try to bind ourself there.  */
    if (hwloc_set_cpubind(topology, cpuset, HWLOC_CPUBIND_THREAD)) {
        char *str = NULL;
#if !defined(HAVE_HWLOC_BITMAP)
        hwloc_cpuset_asprintf(&str, core->cpuset);
#else
        hwloc_bitmap_asprintf(&str, core->cpuset);
#endif
        fprintf(stderr, "isched_hwloc: couldn't bind to cpuset %s\n", str);
        free(str);

        /* Free our cpuset copy */
#if !defined(HAVE_HWLOC_BITMAP)
        hwloc_cpuset_free(cpuset);
#else
        hwloc_bitmap_free(cpuset);
#endif
        return -1;
    }

    /* Get the number at Proc level*/
    cpu_index = core->os_index;

    /* Free our cpuset copy */
#if !defined(HAVE_HWLOC_BITMAP)
    hwloc_cpuset_free(cpuset);
#else
    hwloc_bitmap_free(cpuset);
#endif
    return cpu_index;
}

int isched_hwloc_unbind()
{
#if defined(HAVE_HWLOC_BITMAP)
    hwloc_obj_t      obj;      /* Hwloc object    */
    assert( first_init > 0 );

    /* Get last one.  */
    obj = hwloc_get_obj_by_type(topology, HWLOC_OBJ_MACHINE, 0);
    if (!obj) {
        fprintf(stderr, "isched_hwloc_unbind: Could not get object\n");
        return PASTIX_ERR_UNKNOWN;
    }

    if (hwloc_set_cpubind(topology, obj->cpuset, HWLOC_CPUBIND_THREAD)) {
        char *str = NULL;
        hwloc_bitmap_asprintf(&str, obj->cpuset);
        fprintf(stderr, "isched_hwloc_unbind: Couldn't unbind with cpuset %s\n", str);
        free(str);
        return -1;
    }
#endif
    return PASTIX_SUCCESS;
}

#endif /* defined(HAVE_HWLOC) */
