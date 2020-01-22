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

#include <iostream>
#include <vector>

#define PRINT() std::cout<<"Rank : "<<mpi_rank

int main(int argc, char *argv[])
{
    // Initializing MPI and making epetra comm
    (void)MPI_Init(&argc, &argv);
    Epetra_MpiComm comm(MPI_COMM_WORLD);

    // Get current rank and a total num process
    const int mpi_rank = comm.MyPID();
    const int mpi_size = comm.NumProc();

    /*
    * Creating the intial EpetraFE vector
    */
    const int local_elements_per_rank = 10;
    const int num_over_lapping_elems = 3;
    const int total_num_local_elems = local_elements_per_rank+num_over_lapping_elems;
    std::vector<int> my_local_elements;
    for(int i=0; i<local_elements_per_rank; ++i)
    {
        my_local_elements.emplace_back( (mpi_rank*local_elements_per_rank)+i );
    }

    //defining a map as needed
    const Epetra_Map local_map(-1, my_local_elements.size(), my_local_elements.data(), 0, comm);

    // Making the EpetraFE vector
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
    std::cout << epetra_vector<<std::endl;

    PRINT()<<"NumMyElements() "<<epetra_vector.Map().NumMyElements()<<std::endl;
    PRINT()<<"MyLength() "<<epetra_vector.MyLength()<<std::endl;




    // Now make the map which is to be used as source in import
    std::vector<int> target_map_ids;
    // std::vector<int> target_map_ids(my_local_elements.size(),0);
    // epetra_vector.Map().MyGlobalElements(&target_map_ids[0]);
    if(mpi_rank == 0){
        target_map_ids.push_back(11);
        target_map_ids.push_back(12);
        target_map_ids.push_back(2);
        const Epetra_Map target_map(-1, target_map_ids.size(), target_map_ids.data(), 0, comm);
        Epetra_Import importer(target_map, epetra_vector.Map());

        // Making the EpetraFE vector
        Epetra_FEVector import_epetra_vector(target_map);

        //importing in the new temp vector the values
        ierr = import_epetra_vector.Import(epetra_vector, importer, Insert);
        if(ierr != 0) PRINT()<<" Epetra failure found : "<< ierr<<std::endl;
        std::cout << import_epetra_vector<<std::endl;
    }
    if(mpi_rank == 1){
        // target_map_ids.push_back(4);
        target_map_ids.push_back(5);
        target_map_ids.push_back(6);
        const Epetra_Map target_map(-1, target_map_ids.size(), target_map_ids.data(), 0, comm);
        Epetra_Import importer(target_map, epetra_vector.Map());

        // Making the EpetraFE vector
        Epetra_FEVector import_epetra_vector(target_map);

        //importing in the new temp vector the values
        ierr = import_epetra_vector.Import(epetra_vector, importer, Insert);
        if(ierr != 0) PRINT()<<" Epetra failure found : "<< ierr<<std::endl;
        std::cout << import_epetra_vector<<std::endl;
    }

    std::cout << "Current MPI rank      : "<<mpi_rank << std::endl;
    comm.Barrier();
    if(mpi_rank == 0){
        std::cout << "Test Trilinos successfully complete !" << std::endl;
        std::cout << "Total MPI ranks       : "<<mpi_size << std::endl;
    }

    // MPI_Finalize after you are done using MPI.
    (void) MPI_Finalize ();
    return 0;
}