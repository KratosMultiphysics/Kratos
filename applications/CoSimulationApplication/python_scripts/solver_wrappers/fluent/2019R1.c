#include "udf.h"
#include <math.h>


/* dynamic memory allocation */
#define DECLARE_MEMORY(name,type) type *name = NULL

#define RELEASE_MEMORY(name)										\
if (NNULLP(name)) {													\
	free(name); 													\
	name = NULL;													\
}

#define ASSIGN_MEMORY(name,size,type) 								\
if (size) {															\
	if (NNULLP(name)) { 											\
		name = (type *)realloc(name, size*sizeof(type));			\
	} else {														\
		name = (type *)malloc(size*sizeof(type));					\
	} 																\
	if (NULLP(name)) {												\
		Error("\nUDF-error: Memory assignment failed for name."); 	\
		exit(1);													\
	}																\
}

#define pi 3.1415926535
#define e 2.7182818284


/* Make UDF compatible with 2D and 3D cases, use ND_ND etc. Start with only 2D. */

int n_threads;
DECLARE_MEMORY(thread_ids, int);

#if !RP_HOST
	DECLARE_MEMORY(n_nodes, int); /* most nodes appear 2, 3, or 4 times! */
    DECLARE_MEMORY(n_faces, int);
#endif /* !RP_HOST */


  /*----------------*/
 /* get_thread_ids */
/*----------------*/

DEFINE_ON_DEMAND(get_thread_ids) {
    /* read in thread thread ids, should be called early on;
    expand this explanation */

#if !RP_NODE
    char tmp;
    int k;
    FILE *file;
    file = fopen("bcs.txt", "r");
	fscanf(file, "%i", &n_threads);
#endif /* !RP_NODE */

	host_to_node_int_1(n_threads);
	ASSIGN_MEMORY(thread_ids, n_threads, int);

#if !RP_NODE
	for (k = 0; k < n_threads; k++) {
		fscanf(file, "%s %i", &tmp, &thread_ids[k]);
	}
	fclose(file);
#endif /* !RP_NODE */

	host_to_node_int(thread_ids, n_threads);

#if !RP_HOST
    ASSIGN_MEMORY(n_nodes, n_threads, int);
    ASSIGN_MEMORY(n_faces, n_threads, int);
#endif /* !RP_HOST */

    /* test UDF
    int j;
	Message0("\nUDF read_thread_ids, result:");
	for (j=0;j<n_threads;j++) {
	    Message0(" %i", thread_ids[j]);
	} */

}


  /*-------------*/
 /* store_nodes */
/*-------------*/

DEFINE_ON_DEMAND(store_nodes) {

    /*** annotate this function after it's finished */

    /* printf("\nPOS 1, myid = %i", myid); fflush(stdout); */

    /*** define a bunch of things, such as domain, thread, ints, reals, arrays... */
    /*** check for warnings, some variables are not used on every node/host/... */
#if !RP_HOST
	Domain *domain;
	Thread *face_thread;
	face_t face;
	Node *node;
	int node_number;
#endif /* !RP_HOST */

#if PARALLEL
	int compute_node;
	int n_tmp;
#endif /* PARALLEL */

    int i, d;
    int thread;
    DECLARE_MEMORY(node_ids, int);
	real *node_coords[ND_ND] = {NULL};

    for (thread=0; thread<n_threads; thread++) {

#if !RP_HOST
        domain = Get_Domain(1);
		face_thread = Lookup_Thread(domain, thread_ids[thread]);

		n_nodes[thread] = 0;
		begin_f_loop(face, face_thread) {
			n_nodes[thread] += F_NNODES(face, face_thread);
		} end_f_loop(face, face_thread)

        ASSIGN_MEMORY(node_ids, n_nodes[thread], int);
        for (d = 0; d < ND_ND; d++) {
            ASSIGN_MEMORY(node_coords[d], n_nodes[thread], real);
        }

        i = 0;
        begin_f_loop(face, face_thread) {
            f_node_loop(face, face_thread, node_number) {
                if (i >= n_nodes[thread]) {
                    Error("\nIndex %i larger than array size %i.", i, n_nodes[thread]);
                }
                else {
                    node = F_NODE(face, face_thread, node_number);
                    node_ids[i] = NODE_DM_ID(node);
                    for (d = 0; d < ND_ND; d++) {
                        node_coords[d][i] = NODE_COORD(node)[d];
                    }
                }
                i++;
            }
        } end_f_loop(face, face_thread);

#endif /* !RP_HOST */

#if RP_NODE
        compute_node = (I_AM_NODE_ZERO_P) ? node_host : node_zero;
        /*
        - if compute_node = node_zero: node-zero sends data to the host
        - if compute_node is not node_zero: the node sends data to node_zero
        */

        PRF_CSEND_INT(compute_node, &n_nodes[thread], 1, myid); /* send pointer to n_faces */
        PRF_CSEND_INT(compute_node, node_ids, n_nodes[thread], myid);

        /*for (d = 0; d < ND_ND; d++) {
            PRF_CSEND_REAL(compute_node, node_coords[d], n_nodes[thread], myid);
        } */

        RELEASE_MEMORY(node_ids);

        /*
        Now node_zero has to receive the data from the other nodes, and
        redirect it to the host. Therefore, node_zero loops over all other
        nodes. First, the size of the p_face is stored in n_faces, so that
        enough space can be allocated. Then the received data is stored in
        p_face, and immediately send onward to the host. At the end, p_face
        is freed.
        */
        if(I_AM_NODE_ZERO_P){
            compute_node_loop_not_zero(compute_node) {
                PRF_CRECV_INT(compute_node, &n_tmp, 1, compute_node);

                ASSIGN_MEMORY(node_ids, n_tmp, int);
                PRF_CRECV_INT(compute_node, node_ids, n_tmp, compute_node);

                /*printf("\nFor thread %i, node %i receives from node %i the value n = %i", thread, myid, compute_node, n_tmp); fflush(stdout);*/

                PRF_CSEND_INT(node_host, &n_tmp, 1, compute_node);
                PRF_CSEND_INT(node_host, node_ids, n_tmp, compute_node);

                RELEASE_MEMORY(node_ids);
            }
        }


#endif /* RP_NODE */


#if RP_HOST	/* Receive data on host before printing */
	compute_node_loop(compute_node) {
		PRF_CRECV_INT(node_zero, &n_tmp, 1, compute_node);

		ASSIGN_MEMORY(node_ids, n_tmp, int);
        PRF_CRECV_INT(node_zero, node_ids, n_tmp, compute_node);

		printf("\nFor thread %i, host receives from node %i the value n = %i", thread, compute_node, n_tmp); fflush(stdout);
    }
#endif /* RP_HOST */



    } /* close loop over threads */


    /*** then send data to node 0 and host, and all the receiving... */

    /*** define a filename (per thread, timestep), open it, write data to it */

    printf("\nSuccessfully finished UDF, myid = %i", myid); fflush(stdout);

}

/*
find UDF properties:
grep -r NODE_COORD /apps/SL6.3/ANSYS/2019R1/ansys_inc/v193/fluent/fluent19.3.0/src/
*/


  /*-------------*/
 /* store_faces */
/*-------------*/

DEFINE_ON_DEMAND(store_faces) {
    /* write faces data (ID, nodes),
    loop over all threads */


    /*** find a way to defina a unique ID for face, probably based on node unique ids */
}


  /*------------*/
 /* move_nodes */
/*------------*/

DEFINE_GRID_MOTION(move_nodes, domain, thread, time, dtime) {
    /* moves the nodes on the given thread */

}

