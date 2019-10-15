#include "udf.h"
#include <math.h>


/* dynamic memory allocation */
#define DECLARE_MEMORY(name, type) type *name = NULL

#define DECLARE_MEMORY_N(name, type, dim) type *name[dim] = {NULL}

#define RELEASE_MEMORY(name)										\
if (NNULLP(name)) {													\
	free(name); 													\
	name = NULL;													\
}

#define RELEASE_MEMORY_N(name, dim)							        \
for (_d = 0; _d < dim; _d++) {                                      \
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

#define ASSIGN_MEMORY_N(name, size, type, dim)						\
for (_d = 0; _d < dim; _d++) {                                      \
    ASSIGN_MEMORY(name[_d], size, type);                            \
}

/* sending and receiving arrays in parallel */
#define PRF_CSEND_INT_N(to, name, n, tag, dim)                      \
for (_d = 0; _d < dim; _d++) {                                      \
    PRF_CSEND_INT(to, name[_d], n, tag);                            \
}

#define PRF_CSEND_REAL_N(to, name, n, tag, dim)                     \
for (_d = 0; _d < dim; _d++) {                                      \
    PRF_CSEND_REAL(to, name[_d], n, tag);                           \
}

#define PRF_CRECV_INT_N(from, name, n, tag, dim)                    \
for (_d = 0; _d < dim; _d++) {                                      \
    PRF_CRECV_INT(from, name[_d], n, tag);                          \
}

#define PRF_CRECV_REAL_N(from, name, n, tag, dim)                   \
for (_d = 0; _d < dim; _d++) {                                      \
    PRF_CRECV_REAL(from, name[_d], n, tag);                         \
}

#define mnpf |max_nodes_per_face|

int _d; /* don't use in UDFs! */
int n_threads;
DECLARE_MEMORY(thread_ids, int);
int iteration = 0;
int timestep = 0;


  /*----------------*/
 /* get_thread_ids */
/*----------------*/

DEFINE_ON_DEMAND(get_thread_ids) {
    /* read in thread thread ids, should be called early on; */
    /*** expand this explanation */

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
}


  /*----------------------*/
 /* store_coordinates_id */
/*----------------------*/

DEFINE_ON_DEMAND(store_coordinates_id) {

    /*** annotate this function after it's finished */

    /*** look which checks Joris built into his UDF, add those later */

    /*** principle faces?? */

    int thread, n_nodes, n_faces, i_n, i_f, d;
    DECLARE_MEMORY_N(node_coords, real, ND_ND);
    DECLARE_MEMORY(node_ids, int);
    DECLARE_MEMORY_N(face_coords, real, ND_ND);
    DECLARE_MEMORY_N(face_ids, int, mnpf);

#if !RP_HOST
	Domain *domain;
	Thread *face_thread;
	face_t face;
	Node *node;
	int node_number, j;
	real centroid[ND_ND];
#endif /* !RP_HOST */

#if !RP_NODE
	char file_nodes_name[256];
	char file_faces_name[256];
	FILE *file_nodes = NULL;
	FILE *file_faces = NULL;
#endif /* !RP_NODE */

#if PARALLEL
	int compute_node;
#endif /* PARALLEL */

    for (thread=0; thread<n_threads; thread++) {

#if !RP_NODE
        sprintf(file_nodes_name, "nodes_thread%i.dat", thread_ids[thread]); /*** temp name */
        sprintf(file_faces_name, "faces_thread%i.dat", thread_ids[thread]); /*** temp name */

        if (NULLP(file_nodes = fopen(file_nodes_name, "w"))) {
			Error("\nUDF-error: Unable to open %s for writing\n", file_nodes_name);
			exit(1);
		}
        if (NULLP(file_faces = fopen(file_faces_name, "w"))) {
			Error("\nUDF-error: Unable to open %s for writing\n", file_faces_name);
			exit(1);
		}

#if RP_2D
		fprintf(file_nodes, "%27s %27s %10s\n", "x-coordinate", "y-coordinate", "unique-id");
        fprintf(file_faces, "%27s %27s  %10s\n", "x-coordinate", "y-coordinate", "unique-ids");
#else /* RP_2D */
		fprintf(file_nodes, "%27s %27s %27s %10s\n", "x-coordinate", "y-coordinate", "z-coordinate", "unique-id");
		fprintf(file_faces, "%27s %27s %27s  %10s\n", "x-coordinate", "y-coordinate", "z-coordinate", "unique-ids");
#endif /* RP_2D */
#endif /* !RP_NODE */

#if !RP_HOST
        domain = Get_Domain(1);
		face_thread = Lookup_Thread(domain, thread_ids[thread]);

		n_nodes = 0;
		begin_f_loop(face, face_thread) {
			n_nodes += F_NNODES(face, face_thread);
		} end_f_loop(face, face_thread)
		n_faces = THREAD_N_ELEMENTS_INT(face_thread);

        ASSIGN_MEMORY_N(node_coords, n_nodes, real, ND_ND);
        ASSIGN_MEMORY(node_ids, n_nodes, int);
        ASSIGN_MEMORY_N(face_coords, n_faces, real, ND_ND);
        ASSIGN_MEMORY_N(face_ids, n_faces, int, mnpf);

        i_n = 0;
        i_f = 0;
        begin_f_loop(face, face_thread) {
            if (i_f >= n_faces) {Error("\nIndex %i >= array size %i.", i_f, n_faces);}

            F_CENTROID(centroid, face, face_thread);
            for (d = 0; d < ND_ND; d++) {
                face_coords[d][i_f] = centroid[d];
            }

            for (j = 0; j < mnpf; j++) {
                face_ids[j][i_f] = -1;
            }

            j = 0;
            f_node_loop(face, face_thread, node_number) {
                node = F_NODE(face, face_thread, node_number);

                if (j >= mnpf) {Error("\nIndex %i >= array size %i.", j, mnpf);}
                face_ids[j][i_f] = NODE_DM_ID(node);
                j++;

                if (i_n >= n_nodes) {Error("\nIndex %i >= array size %i.", i_n, n_nodes);}
                node_ids[i_n] = NODE_DM_ID(node);
                for (d = 0; d < ND_ND; d++) {
                    node_coords[d][i_n] = NODE_COORD(node)[d];
                }
                i_n++;
            }
            i_f++;
        } end_f_loop(face, face_thread);
#endif /* !RP_HOST */

#if RP_NODE
        compute_node = (I_AM_NODE_ZERO_P) ? node_host : node_zero;

        PRF_CSEND_INT(compute_node, &n_nodes, 1, myid);
        PRF_CSEND_INT(compute_node, &n_faces, 1, myid);

        PRF_CSEND_REAL_N(compute_node, node_coords, n_nodes, myid, ND_ND);
        PRF_CSEND_INT(compute_node, node_ids, n_nodes, myid);
        PRF_CSEND_REAL_N(compute_node, face_coords, n_faces, myid, ND_ND);
        PRF_CSEND_INT_N(compute_node, face_ids, n_faces, myid, mnpf);

        RELEASE_MEMORY_N(node_coords, ND_ND);
        RELEASE_MEMORY(node_ids);
        RELEASE_MEMORY_N(face_coords, ND_ND);
        RELEASE_MEMORY_N(face_ids, mnpf);

        if(I_AM_NODE_ZERO_P){
            compute_node_loop_not_zero(compute_node) {
                PRF_CRECV_INT(compute_node, &n_nodes, 1, compute_node);
                PRF_CRECV_INT(compute_node, &n_faces, 1, compute_node);

                ASSIGN_MEMORY_N(node_coords, n_nodes, real, ND_ND);
                ASSIGN_MEMORY(node_ids, n_nodes, int);
                ASSIGN_MEMORY_N(face_coords, n_faces, real, ND_ND);
                ASSIGN_MEMORY_N(face_ids, n_faces, int, mnpf);

                PRF_CRECV_REAL_N(compute_node, node_coords, n_nodes, compute_node, ND_ND);
                PRF_CRECV_INT(compute_node, node_ids, n_nodes, compute_node);
                PRF_CRECV_REAL_N(compute_node, face_coords, n_faces, compute_node, ND_ND);
                PRF_CRECV_INT_N(compute_node, face_ids, n_faces, compute_node, mnpf);

                PRF_CSEND_INT(node_host, &n_nodes, 1, compute_node);
                PRF_CSEND_INT(node_host, &n_faces, 1, compute_node);

                PRF_CSEND_REAL_N(node_host, node_coords, n_nodes, compute_node, ND_ND);
                PRF_CSEND_INT(node_host, node_ids, n_nodes, compute_node);
                PRF_CSEND_REAL_N(node_host, face_coords, n_faces, compute_node, ND_ND);
                PRF_CSEND_INT_N(node_host, face_ids, n_faces, compute_node, mnpf);

                RELEASE_MEMORY_N(node_coords, ND_ND);
                RELEASE_MEMORY(node_ids);
                RELEASE_MEMORY_N(face_coords, ND_ND);
                RELEASE_MEMORY_N(face_ids, mnpf);
            }
        }
#endif /* RP_NODE */

#if RP_HOST
        compute_node_loop(compute_node) {
            PRF_CRECV_INT(node_zero, &n_nodes, 1, compute_node);
            PRF_CRECV_INT(node_zero, &n_faces, 1, compute_node);

            ASSIGN_MEMORY_N(node_coords, n_nodes, real, ND_ND);
            ASSIGN_MEMORY(node_ids, n_nodes, int);
            ASSIGN_MEMORY_N(face_coords, n_faces, real, ND_ND);
            ASSIGN_MEMORY_N(face_ids, n_faces, int, mnpf);

            PRF_CRECV_REAL_N(node_zero, node_coords, n_nodes, compute_node, ND_ND);
            PRF_CRECV_INT(node_zero, node_ids, n_nodes, compute_node);
            PRF_CRECV_REAL_N(node_zero, face_coords, n_faces, compute_node, ND_ND);
            PRF_CRECV_INT_N(node_zero, face_ids, n_faces, compute_node, mnpf);
#endif /* RP_HOST */

#if !RP_NODE
            for (i_n = 0; i_n < n_nodes; i_n++) {
                for (d = 0; d < ND_ND; d++) {
                    fprintf(file_nodes, "%27.17e ", node_coords[d][i_n]);
                }
                fprintf(file_nodes, "%10d\n", node_ids[i_n]);
            }

            for (i_f = 0; i_f < n_faces; i_f++) {
                for (d = 0; d < ND_ND; d++) {
                    fprintf(file_faces, "%27.17e ", face_coords[d][i_f]);
                }
                for (d = 0; d < mnpf; d++) {
                    fprintf(file_faces, " %10d", face_ids[d][i_f]);
                }
                fprintf(file_faces, "\n");
            }

            RELEASE_MEMORY_N(node_coords, ND_ND);
            RELEASE_MEMORY(node_ids);
            RELEASE_MEMORY_N(face_coords, ND_ND);
            RELEASE_MEMORY_N(face_ids, mnpf);
#endif /* !RP_NODE */

#if RP_HOST
        } /* close compute_node_loop */
#endif /* RP_HOST */

#if !RP_NODE
        fclose(file_nodes);
        fclose(file_faces);
#endif /* !RP_NODE */

    } /* close loop over threads */

    printf("\n\nNode %i: Finished UDF store_coordinates_id", myid); fflush(stdout);
}


  /*-------------------------*/
 /* store_pressure_traction */
/*-------------------------*/


DEFINE_ON_DEMAND(store_pressure_traction) {

    /*** add inviscid option? */

    int thread, n, i, d;
    DECLARE_MEMORY_N(array, real, ND_ND + 1);
    DECLARE_MEMORY_N(ids, int, mnpf);

#if !RP_HOST
	Domain *domain;
	Thread *face_thread;
	face_t face;
	Node *node;
	int node_number, j;
	real traction[ND_ND], area[ND_ND];
#endif /* !RP_HOST */

#if !RP_NODE
	char file_name[256];
	FILE *file = NULL;
	iteration = RP_Get_Integer("udf/iteration");
	timestep = RP_Get_Integer("udf/timestep");
#endif /* !RP_NODE */

    host_to_node_int_1(iteration);
	host_to_node_int_1(timestep);

#if PARALLEL
	int compute_node;
#endif /* PARALLEL */

    for (thread=0; thread<n_threads; thread++) {

#if !RP_NODE
        sprintf(file_name, "pressure_traction_timestep%i_thread%i.dat",
                timestep, thread_ids[thread]);

        if (NULLP(file = fopen(file_name, "w"))) {
			Error("\nUDF-error: Unable to open %s for writing\n", file_name);
			exit(1);
		}

#if RP_2D
		fprintf(file, "%27s %27s %27s  %10s\n",
		    "x-shear", "y-shear", "pressure", "unique-ids");
#else /* RP_2D */
        fprintf(file, "%27s %27s %27s %27s  %10s\n",
            "x-shear", "y-shear", "z-shear", "pressure", "unique-ids");

#endif /* RP_2D */
#endif /* !RP_NODE */

#if !RP_HOST
        domain = Get_Domain(1);
		face_thread = Lookup_Thread(domain, thread_ids[thread]);

		n = THREAD_N_ELEMENTS_INT(face_thread);

        ASSIGN_MEMORY_N(array, n, real, ND_ND + 1);
        ASSIGN_MEMORY_N(ids, n, int, mnpf);

        i = 0;
        begin_f_loop(face, face_thread) {
            if (i >= n) {Error("\nIndex %i >= array size %i.", i, n);}

            F_AREA(area, face, face_thread);
		    NV_VS(traction, =, F_STORAGE_R_N3V(face, face_thread, SV_WALL_SHEAR), *, -1.0 / NV_MAG(area));
            for (d = 0; d < ND_ND; d++) {
                array[d][i] = traction[d];
            }
            array[ND_ND][i] = F_P(face, face_thread);

            for (j = 0; j < mnpf; j++) {
                ids[j][i] = -1;
            }

            j = 0;
            f_node_loop(face, face_thread, node_number) {
                if (j >= mnpf) {Error("\nIndex %i >= array size %i.", j, mnpf);}
                node = F_NODE(face, face_thread, node_number);
                ids[j][i] = NODE_DM_ID(node);
                j++;
            }
            i++;
        } end_f_loop(face, face_thread);
#endif /* !RP_HOST */

#if RP_NODE
        compute_node = (I_AM_NODE_ZERO_P) ? node_host : node_zero;

        PRF_CSEND_INT(compute_node, &n, 1, myid);

        PRF_CSEND_REAL_N(compute_node, array, n, myid, ND_ND + 1);
        PRF_CSEND_INT_N(compute_node, ids, n, myid, mnpf);

        RELEASE_MEMORY_N(array, ND_ND + 1);
        RELEASE_MEMORY_N(ids, mnpf);

        if(I_AM_NODE_ZERO_P){
            compute_node_loop_not_zero(compute_node) {
                PRF_CRECV_INT(compute_node, &n, 1, compute_node);

                ASSIGN_MEMORY_N(array, n, real, ND_ND + 1);
                ASSIGN_MEMORY_N(ids, n, int, mnpf);

                PRF_CRECV_REAL_N(compute_node, array, n, compute_node, ND_ND + 1);
                PRF_CRECV_INT_N(compute_node, ids, n, compute_node, mnpf);

                PRF_CSEND_INT(node_host, &n, 1, compute_node);

                PRF_CSEND_REAL_N(node_host, array, n, compute_node, ND_ND + 1);
                PRF_CSEND_INT_N(node_host, ids, n, compute_node, mnpf);

                RELEASE_MEMORY_N(array, ND_ND + 1);
                RELEASE_MEMORY_N(ids, mnpf);
            }
        }
#endif /* RP_NODE */

#if RP_HOST
        compute_node_loop(compute_node) {
            PRF_CRECV_INT(node_zero, &n, 1, compute_node);

            ASSIGN_MEMORY_N(array, n, real, ND_ND + 1);
            ASSIGN_MEMORY_N(ids, n, int, mnpf);

            PRF_CRECV_REAL_N(node_zero, array, n, compute_node, ND_ND + 1);
            PRF_CRECV_INT_N(node_zero, ids, n, compute_node, mnpf);
#endif /* RP_HOST */

#if !RP_NODE
            for (i = 0; i < n; i++) {
                for (d = 0; d < ND_ND + 1; d++) {
                    fprintf(file, "%27.17e ", array[d][i]);
                }
                for (d = 0; d < mnpf; d++) {
                    fprintf(file, " %10d", ids[d][i]);
                }
                fprintf(file, "\n");
            }

            RELEASE_MEMORY_N(array, ND_ND + 1);
            RELEASE_MEMORY_N(ids, mnpf);
#endif /* !RP_NODE */

#if RP_HOST
        } /* close compute_node_loop */
#endif /* RP_HOST */

#if !RP_NODE
        fclose(file);
#endif /* !RP_NODE */

    } /* close loop over threads */

    if (myid == 0) {printf("\nFinished UDF store_pressure_traction.\n"); fflush(stdout);}
    /*printf("\n\nNode %i: Finished UDF store_pressure_traction", myid); fflush(stdout);*/
}


  /*------------*/
 /* move_nodes */
/*------------*/

DEFINE_GRID_MOTION(move_nodes, domain, dynamic_thread, time, dtime) {

    /*** add checks and stuff from Joris */

    char file_name[256];
    Thread *face_thread = DT_THREAD(dynamic_thread);
    int thread_id = THREAD_ID(face_thread);

#if !RP_HOST
	face_t face;
	Node *node;
    int i, d, n, node_number;
	DECLARE_MEMORY_N(coords, real, ND_ND);
    DECLARE_MEMORY(ids, int);
    FILE *file = NULL;
#endif /* !RP_HOST */

#if !RP_NODE
	iteration = RP_Get_Integer("udf/iteration");
	timestep = RP_Get_Integer("udf/timestep");
#endif /* !RP_NODE */

	host_to_node_int_1(iteration);
	host_to_node_int_1(timestep);



#if !RP_NODE
        sprintf(file_name, "nodes_update_timestep%i_thread%i.dat",
                timestep, thread_id);
#else
        sprintf(file_name, "nodes_update_timestep%i_thread%i.dat",
                timestep, thread_id);
        host_to_node_sync_file("/tmp");
#endif /* !RP_NODE */

#if RP_HOST
        host_to_node_sync_file(file_name);
#endif /* RP_HOST */

#if !RP_HOST
        if (NULLP(file = fopen(file_name, "r"))) {
            Error("\nUDF-error: Unable to open %s for reading\n", file_name);
            exit(1);
        }

        fscanf(file, "%i", &n);

        ASSIGN_MEMORY_N(coords, n, real, ND_ND);
        ASSIGN_MEMORY(ids, n, int);

        for (i=0; i < n; i++) {
            for (d = 0; d < ND_ND; d++) {
                fscanf(file, "%lf", &coords[d][i]);
            }
            fscanf(file, "%i", &ids[i]);
        }

        fclose(file);

        begin_f_loop(face, face_thread) {
            f_node_loop(face, face_thread, node_number) {
                node = F_NODE(face, face_thread, node_number);
                if NODE_POS_NEED_UPDATE(node) {
                    for (i=0; i < n; i++) {
                        if (NODE_DM_ID(node) == ids[i]) {
                            for (d = 0; d < ND_ND; d++) {
                                NODE_COORD(node)[d] = coords[d][i];
                            }
                            NODE_POS_UPDATED(node);
                            break;
                        }
                    }
                }
            }
        } end_f_loop(face, face_thread);
#endif /* !RP_HOST */
        if (myid == 0) {printf("\nFinished UDF move_nodes.\n"); fflush(stdout);}

}

