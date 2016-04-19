//file for testing the use of octree mesher for Jorge Lopez (CIMAT), in collaboration with Abel Coll and Pooyan Dadvand (CIMNE)

#include "octree_driver.h"
#include "octree_binary_cell.h"
#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"
#include <iostream>
#include <fstream>
int main( int argc, char **argv) {
	int rank,size;
	MPI_Init(&argc,&argv);	
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	MPI_Comm_size(MPI_COMM_WORLD,&size);
	OctreeDriver octree_driver;
	Octree* octree = octree_driver.GetOctree();
	octree->RefineWithUniformSizeNormalized(0.4);
	Node nod1[1],nod2[1],nod3[1];
	nod1->MakeTest(0.1,0.1,0.1);
	nod2->MakeTest(0.9,0.9,0.9);
	nod3->MakeTest(0.1,0.9,0.5);
	Triangle *tri=new Triangle(nod1,nod2,nod3);
	const int n_levels=3;
	octree_driver.SubdivideCellsIntersected(n_levels,tri,octree);
	octree_driver.BalanceOctree();
	octree_driver.PrintMeshGiD();
	delete tri;
	MPI_Finalize();
	return 0;
}
