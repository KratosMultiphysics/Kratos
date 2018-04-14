//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher, Jordi Cotela
//
// See Master-Thesis P.Bucher
// "Development and Implementation of a Parallel
//  Framework for Non-Matching Grid Mapping"

// System includes

// External includes

// Project includes
#include "matrix_based_mapping_operation_utility_mpi.h"

namespace Kratos
{
    /***********************************************************************************/
    /* PUBLIC Methods */
    /***********************************************************************************/
    MatrixBasedMappingOperationUtilityMPI::MatrixBasedMappingOperationUtilityMPI(
        ModelPartPointerType pInterfaceModelPart)
        : MappingOperationUtility(pInterfaceModelPart)
    {
        KRATOS_ERROR_IF_NOT(SparseSpaceType::IsDistributed())
            << "Wrong SparseSpace used!" << std::endl;

        //creating a work array


        Epetra_MpiComm epetra_comm(MPI_COMM_WORLD);


        int* temp = new int[1000]; //

        Epetra_Map my_map(-1, 1000, temp, 0, epetra_comm);


        Epetra_FECrsGraph Agraph(Copy, my_map, 50);

        //////


        int ierr = Agraph.GlobalAssemble();

        // TSystemMatrixPointerType pNewA = TSystemMatrixPointerType(new TSystemMatrixType(Copy,Agraph) );

        SystemMatrixPointerType pNewA = Kratos::make_shared<SystemMatrixType>(Copy,Agraph);

        KRATOS_WATCH("After the Trilinos Stuff")


    }

    /***********************************************************************************/
    /* PROTECTED Methods */
    /***********************************************************************************/


    /***********************************************************************************/
    /* PRIVATE Methods */
    /***********************************************************************************/


}  // namespace Kratos.
