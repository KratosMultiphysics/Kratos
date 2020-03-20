#include "Epetra_FECrsGraph.h"
#include "Epetra_FECrsMatrix.h"
#include "Epetra_Import.h"
#include "Epetra_IntSerialDenseVector.h"
#include "Epetra_IntVector.h"
#include "Epetra_Map.h"
#include "Epetra_MpiComm.h"
#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_SerialDenseVector.h"
#include "Epetra_Vector.h"
#include "EpetraExt_MatrixMatrix.h"
#include "Epetra_FEVector.h"

#include "EpetraExt_CrsMatrixIn.h"
#include <EpetraExt_VectorIn.h>
#include <EpetraExt_RowMatrixOut.h>
#include <EpetraExt_MultiVectorOut.h>

#include <iostream>
#include <vector>
#include <sstream>
#include <fstream>

#define PRINT() std::cout<<"Rank : "<<mpi_rank

int main(int argc, char *argv[])
{
    // Initializing MPI and making epetra comm
    (void)MPI_Init(&argc, &argv);
    Epetra_MpiComm comm(MPI_COMM_WORLD);

    // Get current rank and a total num process
    const int mpi_rank = comm.MyPID();
    const int mpi_size = comm.NumProc();
    const int local_elements_per_rank = 5;
    const int num_over_lapping_elems = 0;
    const int total_num_local_elems = local_elements_per_rank+num_over_lapping_elems;
    std::vector<int> my_local_elements;
    for(int i=0; i<local_elements_per_rank; ++i)
    {
        my_local_elements.emplace_back( (mpi_rank*local_elements_per_rank)+i );
    }
    //defining a map as needed
    const Epetra_Map local_map(-1, my_local_elements.size(), my_local_elements.data(), 0, comm);

    // Making the EpetraFE vector 10x1
    Epetra_FEVector epetra_vector(local_map);

    std::vector<double> values(total_num_local_elems, 0.0);
    for(int i=0; i<total_num_local_elems; ++i)
    {
        values[i] = (mpi_rank+1)*total_num_local_elems+i;
    }
    std::vector<int> ids;
    for(int i=0; i<total_num_local_elems; ++i)
    {
        ids.emplace_back( (mpi_rank*local_elements_per_rank)+i );
    }
    comm.Barrier();

    // Assembling the vector
    int ierr = epetra_vector.SumIntoGlobalValues(ids.size(), ids.data(), values.data());
    if(ierr != 0) PRINT()<<" Epetra failure found : "<< ierr<<std::endl;
    ierr = epetra_vector.GlobalAssemble();
    if(ierr != 0) PRINT()<<" Epetra failure found : "<< ierr<<std::endl;

    // Making dummy epetra matrix (square) 10x10
    Epetra_FECrsGraph a_graph(Copy, local_map, 10);
    // Fill the graph
    ierr = a_graph.InsertGlobalIndices(ids.size(), ids.data(), ids.size(), ids.data());
    ierr = a_graph.GlobalAssemble();
    // use the graph to make the matrix
    Epetra_FECrsMatrix mat_a(Copy,a_graph);
    if(mpi_rank == 0){
        std::vector<double> diag_values{2, 5.50, 2.111};
        std::vector<int> col_indices{1,8,5};
        mat_a.ReplaceGlobalValues(5, col_indices.size(), diag_values.data(), col_indices.data());
    }
    // mat_a.PutScalar(1.0);
    mat_a.GlobalAssemble();

    {
        EpetraExt::RowMatrixToMatrixMarketFile("amat.mm", mat_a);
    }

    // Making dummy epetra matrix (rectangle) 10x2
    std::vector<int> row_elems = my_local_elements;
    std::vector<int> col_elems{0,1};
    const Epetra_Map row_map(-1, row_elems.size(), row_elems.data(), 0, comm);
    const Epetra_Map col_map(-1, col_elems.size(), col_elems.data(), 0, comm);

    Epetra_FECrsGraph b_graph(Copy, row_map, 10);
    ierr = b_graph.InsertGlobalIndices(ids.size(), ids.data(), col_elems.size(), col_elems.data());
    ierr = b_graph.GlobalAssemble();
    std::cout<<"b_graph "<<b_graph<<std::endl;

    Epetra_FECrsMatrix mat_b(Copy,b_graph);
    mat_b.FillComplete();
    if(mpi_rank == 1){
        std::vector<double> diag_valuesb{2, 0.50};
        std::vector<int> col_indicesb{0,1};
        mat_b.ReplaceGlobalValues(5, 2, diag_valuesb.data(), col_indicesb.data());
    }
    mat_b.GlobalAssemble();

    {
        EpetraExt::RowMatrixToMatrixMarketFile("bmat.mm", mat_b);
    }

    // Now multiplying c = a*b
    Epetra_FECrsMatrix mat_c(Copy,local_map, 1);
    // c <- A*B
    EpetraExt::MatrixMatrix::Multiply(mat_a, false, mat_b, false, mat_c, true);
    mat_c.FillComplete();
    mat_c.GlobalAssemble();

    {
        EpetraExt::RowMatrixToMatrixMarketFile("cmat.mm", mat_c);
    }

    Epetra_FEVector vec_x(local_map);
    Epetra_FEVector vec_y(local_map);
    if(mpi_rank == 0){
        std::vector<double> vals = {1.0,1.0,1.0,1.0};
        std::vector<int> ids = {1,4,5,8};
        vec_x.ReplaceGlobalValues(ids.size(), ids.data(), vals.data());
    }
    vec_x.GlobalAssemble();

    EpetraExt::MultiVectorToMatrixMarketFile("xvec.mm", vec_x);

    mat_a.Multiply(false, vec_x, vec_y);

    EpetraExt::MultiVectorToMatrixMarketFile("yvec.mm", vec_y);

    vec_y.GlobalAssemble();

    // MPI_Finalize after you are done using MPI.
    (void) MPI_Finalize ();
    return 0;
}