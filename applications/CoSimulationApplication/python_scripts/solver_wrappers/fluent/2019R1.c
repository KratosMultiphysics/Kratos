#include "udf.h"
#include <math.h>


/* dynamic memory allocation */
#define DECLARE_MEMORY(name, type) type *name = NULL

#define DECLARE_MEMORY_ND(name, type) type *name[ND_ND] = {NULL}

#define RELEASE_MEMORY(name)										\
if (NNULLP(name)) {													\
	free(name); 													\
	name = NULL;													\
}

#define RELEASE_MEMORY_ND(name)										\
for (_d = 0; _d < ND_ND; _d++) {                                    \
    RELEASE_MEMORY(name[_d]);                                       \
}

#define ASSIGN_MEMORY(name, size, type) 							\
if (size) {															\
	if (NNULLP(name)) { 											\
		name = (type *)realloc(name, size * sizeof(type));			\
	} else {														\
		name = (type *)malloc(size * sizeof(type));					\
	} 																\
	if (NULLP(name)) {												\
		Error("\nUDF-error: Memory assignment failed for name."); 	\
		exit(1);													\
	}																\
}

#define ASSIGN_MEMORY_ND(name, size, type) 							\
for (_d = 0; _d < ND_ND; _d++) {                                    \
    ASSIGN_MEMORY(name[_d], size, type);                             \
}

/* sending and receiving */
#define PRF_CSEND_REAL_ND(to, name, n, tag)                         \
for (_d = 0; _d < ND_ND; _d++) {                                    \
    PRF_CSEND_REAL(to, name[_d], n, tag);                           \
}

#define PRF_CRECV_REAL_ND(from, name, n, tag)                       \
for (_d = 0; _d < ND_ND; _d++) {                                    \
    PRF_CRECV_REAL(from, name[_d], n, tag);                         \
}


/*** perhaps add fnction DECLARE_MEMORY for 2D? or specifically ND_ND? */

#define pi 3.1415926535
#define e 2.7182818284


/* Make UDF compatible with 2D and 3D cases, use ND_ND etc. Start with only 2D. */

int _d; /* don't use in UDFs! */

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
    /* read in thread thread ids, should be called early on; */
    /* expand this explanation */

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
    /*** look which checks Joris built into his UDF */

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

#if !RP_NODE
	char file_name[256];
	FILE *file = NULL;
#endif /* !RP_NODE */

#if PARALLEL
	int compute_node;
#endif /* PARALLEL */

    int i, d, n_tmp;
    int thread;
    DECLARE_MEMORY(node_ids, int);
    DECLARE_MEMORY_ND(node_coords, real);

	/*real *node_coords[ND_ND] = {NULL};*/

    for (thread=0; thread<n_threads; thread++) {

#if !RP_NODE
        sprintf(file_name, "nodes_thread%i.dat", thread_ids[thread]); /*** temp name */

        if (NULLP(file = fopen(file_name, "w"))) {
			Error("\nUDF-error: Unable to open %s for writing\n", file_name);
			exit(1);
		}

#if RP_2D
		fprintf(file, "node-id\tx-coordinate\ty-coordinate\n");
#else /* RP_2D */
		fprintf(file, "node-id\tx-coordinate\ty-coordinate\tz-coordinate\n");
#endif /* RP_2D */

#endif /* !RP_NODE */

#if !RP_HOST
        domain = Get_Domain(1);
		face_thread = Lookup_Thread(domain, thread_ids[thread]);

		n_nodes[thread] = 0;
		begin_f_loop(face, face_thread) {
			n_nodes[thread] += F_NNODES(face, face_thread);
		} end_f_loop(face, face_thread)

        ASSIGN_MEMORY(node_ids, n_nodes[thread], int);
        ASSIGN_MEMORY_ND(node_coords, n_nodes[thread], real);

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
        PRF_CSEND_REAL_ND(compute_node, node_coords, n_nodes[thread], myid);

        RELEASE_MEMORY(node_ids);
        RELEASE_MEMORY_ND(node_coords);

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
                ASSIGN_MEMORY_ND(node_coords, n_tmp, real);

                PRF_CRECV_INT(compute_node, node_ids, n_tmp, compute_node);
                PRF_CRECV_REAL_ND(compute_node, node_coords, n_tmp, compute_node);

                PRF_CSEND_INT(node_host, &n_tmp, 1, compute_node);
                PRF_CSEND_INT(node_host, node_ids, n_tmp, compute_node);
                PRF_CSEND_REAL_ND(node_host, node_coords, n_tmp, compute_node);

                RELEASE_MEMORY(node_ids);
                RELEASE_MEMORY_ND(node_coords);
            }
        }
#endif /* RP_NODE */


#if RP_HOST	/* Receive data on host before printing */
        compute_node_loop(compute_node) {
            PRF_CRECV_INT(node_zero, &n_tmp, 1, compute_node);

            ASSIGN_MEMORY(node_ids, n_tmp, int);
            ASSIGN_MEMORY_ND(node_coords, n_tmp, real);

            PRF_CRECV_INT(node_zero, node_ids, n_tmp, compute_node);
            PRF_CRECV_REAL_ND(node_zero, node_coords, n_tmp, compute_node);
#endif /* RP_HOST */

#if !PARALLEL
            n_tmp = n_nodes[thread];
#endif /* !PARALLEL */

#if !RP_NODE
            for (i = 0; i < n_tmp; i++) {
                fprintf(file, "%i", node_ids[i]);

                for (d = 0; d < ND_ND; d++) {
                    fprintf(file, "\t%27.17e", node_coords[d][i]);
                }
                fprintf(file, "\n");

            }

            RELEASE_MEMORY(node_ids);
            RELEASE_MEMORY_ND(node_coords);
#endif /* !RP_NODE */

#if RP_HOST
        } /* close compute_node_loop */
#endif /* RP_HOST */

#if !RP_NODE
        fclose(file);
#endif /* !RP_NODE */

    } /* close loop over threads */

    printf("\nNode %i: Finished UDF store_nodes", myid); fflush(stdout);

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

