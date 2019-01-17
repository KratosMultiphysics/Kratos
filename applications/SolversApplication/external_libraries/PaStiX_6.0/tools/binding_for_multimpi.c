/**
 *
 * Small test to see what's the bindtab array should look like in case of
 * multiple MPI processes per node.
 * This test can be compiled with:
 *    mpicc binding_for_multimpi.c -o binding_for_multimpi -Wall -lpthread -lhwloc
 *
 */
#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include <unistd.h>
#include <string.h>
#include <hwloc.h>

/**
 * Extract of PaStiX internal scheduler for simpler compilation of this small test
 */
static hwloc_topology_t topology;

int isched_topo_init(void)
{
    hwloc_topology_init(&topology);
    hwloc_topology_load(topology);
    return 0;
}

int isched_topo_world_size()
{
    hwloc_obj_t obj = hwloc_get_obj_by_type( topology, HWLOC_OBJ_MACHINE, 0 );
    return hwloc_get_nbobjs_inside_cpuset_by_type(topology, obj->cpuset, HWLOC_OBJ_CORE);
}
/**
 * End of the extract
 */

#define BUF_MAX 256

int main( int argc, char *argv[] )
{
    MPI_Comm intra_comm;
    int     i, len;
    char    procname[BUF_MAX];
    int     rc, key;
    int64_t color;
    int world_size, world_rank;
    int intra_size, intra_rank;
    int nbthread, intra_nbthread;

    MPI_Init( &argc, &argv );
    MPI_Comm_size( MPI_COMM_WORLD, &world_size );
    MPI_Comm_rank( MPI_COMM_WORLD, &world_rank );
    key = world_rank;


    /**
     * Section that need to be added to the code to generate the bindtab arrays
     */
    {
        /**
         * Get hostname to generate a hash that will be the color of each node
         * MPI_Get_processor_name is not used as it can returned different
         * strings for processes of a same physical node.
         */
        rc = gethostname(procname, BUF_MAX-1);
        procname[BUF_MAX-1] = '\0';
        len = strlen( procname );

        /**
         * Compute hash of the procname
         */
        color = 0;
        for (i = 0; i < len; i++) {
            color = color*256*sizeof(char) + procname[i];
        }

        /**
         * Create intra-node communicator
         */
        MPI_Comm_split( MPI_COMM_WORLD, color, key, &intra_comm );
        MPI_Comm_size( intra_comm, &intra_size );
        MPI_Comm_rank( intra_comm, &intra_rank );

        /**
         * Initialize hwloc topology of pastix
         * Can be done before pastixInit, and must be done in that case to get
         * the number of cores available
         */
        isched_topo_init();
        nbthread = isched_topo_world_size();

        intra_nbthread = nbthread / intra_size;
        /* Make sure it's at least 1 */
        intra_nbthread = (intra_nbthread < 1) ? 1 : intra_nbthread;
    }

    if ( world_rank == 0 ) {
        printf( " Number of MPI processes:     %d\n"
                " Number of nodes:             %d\n"
                " Number of cores per node:    %d\n"
                " Number of cores per process: %d\n",
                world_size, world_size / intra_size,
                nbthread, intra_nbthread );
    }

    /**
     * Print the bindtab arrays per node
     */
    {
        char corelists[2*BUF_MAX];
        char *c = corelists;
        rc = sprintf( c, "[%2d - %s] :", world_rank, procname );
        c += rc;
        for(i=0; i<intra_nbthread; i++) {
            rc = sprintf( c, "%2d ", intra_nbthread * intra_rank + i );
            c += rc;
        }
        printf( "%s\n", corelists );
    }

    MPI_Finalize();

    return EXIT_SUCCESS;
}
