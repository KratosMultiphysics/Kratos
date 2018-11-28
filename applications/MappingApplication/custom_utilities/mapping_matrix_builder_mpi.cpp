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
#include<set>

// External includes

// Project includes
#include "mapping_matrix_builder.h"
#include "custom_utilities/mapper_typedefs.h"
#include "custom_utilities/mapper_utilities.h"
#include "mapping_application_variables.h"

namespace Kratos
{
typedef typename MapperDefinitions::MPISparseSpaceType SparseSpaceType;
typedef typename MapperDefinitions::DenseSpaceType DenseSpaceType;

typedef MappingMatrixBuilder<SparseSpaceType, DenseSpaceType> BuilderType;

typedef typename MapperLocalSystem::MatrixType MatrixType;
typedef typename MapperLocalSystem::EquationIdVectorType EquationIdVectorType;

void ConstructRowColIdSets(BuilderType::MapperLocalSystemPointerVector& rMapperLocalSystems,
                           std::set<int>& rRowEquationIds,
                           std::set<int>& rColEquationIds)
{
    EquationIdVectorType origin_ids;
    EquationIdVectorType destination_ids;

    for (auto& rp_local_sys : rMapperLocalSystems) {
        rp_local_sys->EquationIdVectors(origin_ids, destination_ids);

        rRowEquationIds.insert(destination_ids.begin(), destination_ids.end());
        rColEquationIds.insert(origin_ids.begin(), origin_ids.end());
    }
}

void ConstructMatrixStructure(Epetra_FECrsGraph& rGraph,
                              BuilderType::MapperLocalSystemPointerVector& rMapperLocalSystems)
{
    EquationIdVectorType origin_ids;
    EquationIdVectorType destination_ids;
    int ierr;

    for (auto& rp_local_sys : rMapperLocalSystems) {
        rp_local_sys->EquationIdVectors(origin_ids, destination_ids);

        if (origin_ids.size() > 0) {
            ierr = rGraph.InsertGlobalIndices(
                destination_ids.size(), destination_ids.data(),
                origin_ids.size(),      origin_ids.data());

            // TODO maybe change this to a debug error only
            KRATOS_ERROR_IF( ierr < 0 ) << "Epetra failure in Epetra_FECrsGraph.InsertGlobalIndices. "
                << "Error code: " << ierr << std::endl;
        }
    }
}

void BuildMatrix(BuilderType::TMappingMatrixUniquePointerType& rpMdo,
                       BuilderType::MapperLocalSystemPointerVector& rMapperLocalSystems)
{
    MatrixType local_mapping_matrix;
    EquationIdVectorType origin_ids;
    EquationIdVectorType destination_ids;
    int ierr;

    for (auto& rp_local_sys : rMapperLocalSystems)
    {
        rp_local_sys->CalculateLocalSystem(local_mapping_matrix, origin_ids, destination_ids);

        KRATOS_DEBUG_ERROR_IF(local_mapping_matrix.size1() != destination_ids.size())
            << "DestinationID vector size mismatch" << std::endl;
        KRATOS_DEBUG_ERROR_IF(local_mapping_matrix.size2() != origin_ids.size())
            << "OriginID vector size mismatch" << std::endl;

        if (local_mapping_matrix.size1() > 0) {
            ierr = rpMdo->SumIntoGlobalValues(
                destination_ids.size(), destination_ids.data(),
                origin_ids.size(), origin_ids.data(),
                local_mapping_matrix.data().begin(), // TODO I think this changes with AMatrix
                Epetra_FECrsMatrix::ROW_MAJOR ); // same for Ublas ad AMatrix

            // TODO maybe change this to a debug error only
            KRATOS_ERROR_IF( ierr != 0 ) << "Epetra failure in Epetra_FECrsMatrix.SumIntoGlobalValues. "
                << "Error code: " << ierr << std::endl;
        }

        // The local-systems are always cleared since they would be recomputed
        // to fill a new MappingMatrix
        rp_local_sys->Clear();
    }
}

/***********************************************************************************/
/* PUBLIC Methods */
/***********************************************************************************/
template<>
BuilderType::MappingMatrixBuilder()
{
    KRATOS_ERROR_IF_NOT(SparseSpaceType::IsDistributed())
        << "Using a non-distributed Space!" << std::endl;
}

template<>
void BuilderType::BuildMappingMatrix(
        const InterfaceVectorContainerPointerType& rpVectorContainerOrigin,
        const InterfaceVectorContainerPointerType& rpVectorContainerDestination,
        MapperLocalSystemPointerVector& rMapperLocalSystems)
{
    // ***** Creating vectors with information abt which IDs are local *****
    const auto& r_local_mesh_origin = rpVectorContainerOrigin->GetModelPart().GetCommunicator().LocalMesh();
    const auto& r_local_mesh_destination = rpVectorContainerDestination->GetModelPart().GetCommunicator().LocalMesh();

    const int num_local_nodes_orig = r_local_mesh_origin.NumberOfNodes();
    const int num_local_nodes_dest = r_local_mesh_destination.NumberOfNodes();

    std::vector<int> global_equation_ids_origin(num_local_nodes_orig);
    std::vector<int> global_equation_ids_destination(num_local_nodes_dest);

    const auto nodes_begin_orig = r_local_mesh_origin.NodesBegin();
    #pragma omp parallel for
    for (int i=0; i<num_local_nodes_orig; ++i) {
        global_equation_ids_origin[i] = ( nodes_begin_orig + i )->GetValue(INTERFACE_EQUATION_ID);
    }

    const auto nodes_begin_dest = r_local_mesh_destination.NodesBegin();
    #pragma omp parallel for
    for (int i=0; i<num_local_nodes_dest; ++i) {
        global_equation_ids_destination[i] = ( nodes_begin_dest + i )->GetValue(INTERFACE_EQUATION_ID);
    }

    // Construct vectors containing all the equation ids of rows and columns this processor contributes to
    std::set<int> row_equation_ids_set;
    std::set<int> col_equation_ids_set;
    ConstructRowColIdSets(rMapperLocalSystems, row_equation_ids_set, col_equation_ids_set);
    std::vector<int> row_equation_ids(row_equation_ids_set.begin(), row_equation_ids_set.end());
    std::vector<int> col_equation_ids(col_equation_ids_set.begin(), col_equation_ids_set.end());

    // ***** Creating the maps for the MappingMatrix and the SystemVectors *****
    const Epetra_MpiComm epetra_comm(MPI_COMM_WORLD);

    // Epetra_Map (long long NumGlobalElements, int NumMyElements, const long long *MyGlobalElements, int IndexBase, const Epetra_Comm &Comm)

    const int num_global_elements = -1; // this means its gonna be computed by Epetra_Map // TODO I think I know this...
    const int index_base = 0; // for C/C++

    Epetra_Map epetra_col_map(num_global_elements,
                              col_equation_ids.size(),
                              col_equation_ids.data(), // taken as const
                              index_base,
                              epetra_comm);

    Epetra_Map epetra_row_map(num_global_elements,
                              row_equation_ids.size(),
                              row_equation_ids.data(), // taken as const
                              index_base,
                              epetra_comm);

    Epetra_Map epetra_domain_map(num_global_elements,
                                 global_equation_ids_origin.size(),
                                 global_equation_ids_origin.data(), // taken as const
                                 index_base,
                                 epetra_comm);

    Epetra_Map epetra_range_map(num_global_elements,
                                global_equation_ids_destination.size(),
                                global_equation_ids_destination.data(), // taken as const
                                index_base,
                                epetra_comm);

    // ***** Creating the graph for the MappingMatrix *****
    // explanation in here: https://trilinos.org/docs/dev/packages/epetra/doc/html/classEpetra__CrsGraph.html
    const int num_indices_per_row = 5; // TODO this is to be tested => set to zero maybe ... => also applicable for the serial version

    // Performance optimization see https://trilinos.org/docs/dev/packages/epetra/doc/html/classEpetra__CrsGraph.html
    Epetra_FECrsGraph epetra_graph(Epetra_DataAccess::Copy,
                                   epetra_row_map,
                                   epetra_col_map,
                                   num_indices_per_row);

    ConstructMatrixStructure(epetra_graph, rMapperLocalSystems);

    // range- and domain-map have to be passed since the matrix is rectangular
    int ierr = epetra_graph.GlobalAssemble(epetra_domain_map, epetra_range_map); // TODO check if it should call "FillComplete"

    KRATOS_ERROR_IF( ierr != 0 ) << "Epetra failure in Epetra_FECrsGraph.GlobalAssemble. "
        << "Error code: " << ierr << std::endl;

    epetra_graph.OptimizeStorage(); // TODO is an extra-call needed?

    // ***** Creating the MappingMatrix *****
    TMappingMatrixUniquePointerType p_Mdo =
        Kratos::make_unique<TMappingMatrixType>(Epetra_DataAccess::Copy, epetra_graph);

    BuildMatrix(p_Mdo, rMapperLocalSystems);

    // range- and domain-map have to be passed since the matrix is rectangular
    ierr = p_Mdo->GlobalAssemble(epetra_domain_map, epetra_range_map);

    KRATOS_ERROR_IF( ierr != 0 ) << "Epetra failure in Epetra_FECrsMatrix.GlobalAssemble. "
        << "Error code: " << ierr << std::endl;

    if (mEchoLevel > 2) {
        SparseSpaceType::WriteMatrixMarketMatrix("TrilinosMappingMatrix.mm", *p_Mdo, false);
    }

    mpMappingMatrix.swap(p_Mdo);

    // ***** Creating the SystemVectors *****
    TSystemVectorUniquePointerType p_new_vector_destination =
        Kratos::make_unique<TSystemVectorType>(epetra_range_map);
    TSystemVectorUniquePointerType p_new_vector_origin =
        Kratos::make_unique<TSystemVectorType>(epetra_domain_map);
    rpVectorContainerDestination->pGetVector().swap(p_new_vector_destination);
    rpVectorContainerOrigin->pGetVector().swap(p_new_vector_origin);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Class template instantiation
template class MappingMatrixBuilder< SparseSpaceType, DenseSpaceType >;

}  // namespace Kratos.
