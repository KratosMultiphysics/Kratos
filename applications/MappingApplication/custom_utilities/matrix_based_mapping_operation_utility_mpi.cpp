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


        // big TODO what if the rank doesn't have local nodes ... ?

        const int num_local_nodes = 100; //pInterfaceModelPart->GetCommunicator().LocalMesh().NumberOfNodes();

        const Epetra_MpiComm epetra_comm(MPI_COMM_WORLD);

        std::vector<int> global_elements(num_local_nodes);

        // fill the
        // #pragma omp parallel for
        // for local_nodes:
        //     global_elements[i] = node.GetValue(INTERFACE_EQUATION_ID); // Is it that easy?

        // Epetra_Map (long long NumGlobalElements, int NumMyElements, const long long *MyGlobalElements, int IndexBase, const Epetra_Comm &Comm)

        const int num_global_elements = -1; // computed by Epetra_Map
        const int index_base = 0; // for C/C++

        Epetra_Map epetra_map(num_global_elements,
                              num_local_nodes,
                              global_elements.data(), // taken as const
                              index_base,
                              epetra_comm);


        const int num_indices_per_row = 25; // TODO this is to be tested

        Epetra_FECrsGraph epetra_graph(Epetra_DataAccess::Copy,
                                       epetra_map,
                                       num_indices_per_row);

        // //////

        // now loop the mapperconditions and fill the graph with the equation ids
        // for mapper_conditions:
        //     // InsertGlobalIndices (int numRows, const int *rows, int numCols, const int *cols)
        //     ierr = epetra_graph.InsertGlobalIndices(...)


        int ierr = epetra_graph.GlobalAssemble(); // TODO check if it should call "FillComplete"
        KRATOS_ERROR_IF( ierr != 0 ) << "Epetra failure in Epetra_FECrsGraph.GlobalAssemble. "
            << "Error code: " << ierr << std::endl;

        // // TSystemMatrixPointerType pNewA = TSystemMatrixPointerType(new TSystemMatrixType(Copy,Agraph) );
        // // https://trilinos.org/docs/dev/packages/epetra/doc/html/Epetra__DataAccess_8h.html#ad1a985e79f94ad63030815a0d7d90928
        SystemMatrixPointerType p_new_mapping_matrix = Kratos::make_shared<SystemMatrixType>(Epetra_DataAccess::Copy, epetra_graph);
        mpMappingMatrix.swap(p_new_mapping_matrix);

        SystemVectorPointerType p_new_vector_destination = Kratos::make_shared<SystemVectorType>(epetra_map);
        mpVectorDestination.swap(p_new_vector_destination);

        SystemVectorPointerType p_new_vector_origin = Kratos::make_shared<SystemVectorType>(epetra_map/*_origin*/);
        mpVectorOrigin.swap(p_new_vector_origin);

        // SparseSpaceType::Mult(*mpMappingMatrix, *mpVectorOrigin, *mpVectorDestination);


        KRATOS_WATCH("After the Trilinos Stuff")


    }

    /***********************************************************************************/
    /* PROTECTED Methods */
    /***********************************************************************************/


    /***********************************************************************************/
    /* PRIVATE Methods */
    /***********************************************************************************/


}  // namespace Kratos.
