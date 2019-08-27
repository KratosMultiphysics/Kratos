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


/*
This is the UDF.
What should I put in the UDF?


Make UDF compatible with 2D and 3D cases, use ND_ND etc. Start with only 2D.



Create UDF to which reads in surface ids: read bcs.txt, store the IDs in a global, dynamic array of ints.


*/

int n_threads;
DECLARE_MEMORY(thread_ids, int);


DEFINE_ON_DEMAND(get_thread_ids){
    /* read in surface thread ids, should be called early on;
    expand this explanation */

#if !RP_NODE
    char tmp;
    int k;
    FILE *file;
    file = fopen("bcs.txt","r");
	fscanf(file,"%i",&n_threads);
#endif /* !RP_NODE */

	host_to_node_int_1(n_threads);
	ASSIGN_MEMORY(thread_ids,n_threads,int);

#if !RP_NODE
	for (k=0;k<n_threads;k++) {
		fscanf(file,"%s %i",&tmp,&thread_ids[k]);
	}
	fclose(file);
#endif /* !RP_NODE */

	host_to_node_int(thread_ids,n_threads);

    /* test UDF
    int j;
	Message0("\nUDF read_thread_ids, result:");
	for (j=0;j<n_threads;j++) {
	    Message0(" %i", thread_ids[j]);
	} */

}


DEFINE_ON_DEMAND(store_faces){
    /* write faces data (ID, nodes), loop over all threads */
}

DEFINE_ON_DEMAND(store_nodes){
    /* write nodes data (ID, coords), loop over all threads */
}

DEFINE_GRID_MOTION(move_nodes, domain, thread, time, dtime) {
    /* moves the nodes on the given thread */
}

