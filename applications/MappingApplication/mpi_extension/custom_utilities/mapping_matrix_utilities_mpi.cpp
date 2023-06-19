//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher, Jordi Cotela
//
// See Master-Thesis P.Bucher
// "Development and Implementation of a Parallel
//  Framework for Non-Matching Grid Mapping"

// System includes
#include <set>

// External includes
#include "Epetra_FECrsGraph.h"

// Project includes
#include "utilities/parallel_utilities.h"
#include "custom_utilities/mapping_matrix_utilities.h"
#include "mapper_mpi_define.h"
#include "custom_utilities/mapper_utilities.h"
#include "mapping_application_variables.h"

namespace Kratos {

namespace {

typedef typename MPIMapperDefinitions::SparseSpaceType MappingSparseSpaceType;
typedef typename MPIMapperDefinitions::DenseSpaceType  DenseSpaceType;

typedef MappingMatrixUtilities<MappingSparseSpaceType, DenseSpaceType> MappingMatrixUtilitiesType;

typedef typename MapperLocalSystem::MatrixType MatrixType;
typedef typename MapperLocalSystem::EquationIdVectorType EquationIdVectorType;

void ConstructRowColIdSets(std::vector<Kratos::unique_ptr<MapperLocalSystem>>& rMapperLocalSystems,
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
                              std::vector<Kratos::unique_ptr<MapperLocalSystem>>& rMapperLocalSystems)
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

void BuildMatrix(Kratos::unique_ptr<typename MappingSparseSpaceType::MatrixType>& rpMdo,
                 std::vector<Kratos::unique_ptr<MapperLocalSystem>>& rMapperLocalSystems)
{
    MatrixType local_mapping_matrix;
    EquationIdVectorType origin_ids;
    EquationIdVectorType destination_ids;
    int ierr;

    for (auto& rp_local_sys : rMapperLocalSystems) {
        rp_local_sys->CalculateLocalSystem(local_mapping_matrix, origin_ids, destination_ids);

        KRATOS_DEBUG_ERROR_IF(local_mapping_matrix.size1() != destination_ids.size()) << "MPI-MappingMatrixAssembly: DestinationID vector size mismatch: LocalMappingMatrix-Size1: " << local_mapping_matrix.size1() << " | DestinationIDs-size: " << destination_ids.size() << std::endl;
        KRATOS_DEBUG_ERROR_IF(local_mapping_matrix.size2() != origin_ids.size())<< "MPI-MappingMatrixAssembly: OriginID vector size mismatch: LocalMappingMatrix-Size2: " << local_mapping_matrix.size2() << " | OriginIDs-size: " << origin_ids.size() << std::endl;

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

} // anonymous namespace

template<>
void MappingMatrixUtilitiesType::InitializeSystemVector(
    Kratos::unique_ptr<typename MappingSparseSpaceType::VectorType>& rpVector,
    const std::size_t VectorSize)
{
    KRATOS_ERROR << "this function was not yet implemented in Trilinos!" << std::endl;
}

template<>
void MappingMatrixUtilitiesType::BuildMappingMatrix(
    Kratos::unique_ptr<typename MappingSparseSpaceType::MatrixType>& rpMappingMatrix,
    Kratos::unique_ptr<typename MappingSparseSpaceType::VectorType>& rpInterfaceVectorOrigin,
    Kratos::unique_ptr<typename MappingSparseSpaceType::VectorType>& rpInterfaceVectorDestination,
    const ModelPart& rModelPartOrigin,
    const ModelPart& rModelPartDestination,
    std::vector<Kratos::unique_ptr<MapperLocalSystem>>& rMapperLocalSystems,
    const int EchoLevel)
{
    KRATOS_TRY

    static_assert(MappingSparseSpaceType::IsDistributed(), "Using a non-distributed Space!");

    // ***** Creating vectors with information abt which IDs are local *****
    const auto& r_local_mesh_origin = rModelPartOrigin.GetCommunicator().LocalMesh();
    const auto& r_local_mesh_destination = rModelPartDestination.GetCommunicator().LocalMesh();

    const int num_local_nodes_orig = r_local_mesh_origin.NumberOfNodes();
    const int num_local_nodes_dest = r_local_mesh_destination.NumberOfNodes();

    std::vector<int> global_equation_ids_origin(num_local_nodes_orig);
    std::vector<int> global_equation_ids_destination(num_local_nodes_dest);

    const auto nodes_begin_orig = r_local_mesh_origin.NodesBegin();
    IndexPartition<std::size_t>(num_local_nodes_orig).for_each([&global_equation_ids_origin, &nodes_begin_orig](const std::size_t Index){
        global_equation_ids_origin[Index] = (nodes_begin_orig+Index)->GetValue(INTERFACE_EQUATION_ID);
    });

    const auto nodes_begin_dest = r_local_mesh_destination.NodesBegin();
    IndexPartition<std::size_t>(num_local_nodes_dest).for_each([&global_equation_ids_destination, &nodes_begin_dest](const std::size_t Index){
        global_equation_ids_destination[Index] = (nodes_begin_dest+Index)->GetValue(INTERFACE_EQUATION_ID);
    });

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
    Kratos::unique_ptr<typename MappingSparseSpaceType::MatrixType> p_Mdo =
        Kratos::make_unique<typename MappingSparseSpaceType::MatrixType>(Epetra_DataAccess::Copy, epetra_graph);

    BuildMatrix(p_Mdo, rMapperLocalSystems);

    // range- and domain-map have to be passed since the matrix is rectangular
    ierr = p_Mdo->GlobalAssemble(epetra_domain_map, epetra_range_map);

    KRATOS_ERROR_IF( ierr != 0 ) << "Epetra failure in Epetra_FECrsMatrix.GlobalAssemble. "
        << "Error code: " << ierr << std::endl;

    if (EchoLevel > 2) {
        const std::string file_name = "TrilinosMappingMatrix_O_" + rModelPartOrigin.Name() + "__D_" + rModelPartDestination.Name() +".mm";
        MappingSparseSpaceType::WriteMatrixMarketMatrix(file_name.c_str(), *p_Mdo, false);
    }

    rpMappingMatrix.swap(p_Mdo);

    // ***** Creating the SystemVectors *****
    Kratos::unique_ptr<typename MappingSparseSpaceType::VectorType> p_new_vector_destination =
        Kratos::make_unique<typename MappingSparseSpaceType::VectorType>(epetra_range_map);
    Kratos::unique_ptr<typename MappingSparseSpaceType::VectorType> p_new_vector_origin =
        Kratos::make_unique<typename MappingSparseSpaceType::VectorType>(epetra_domain_map);
    rpInterfaceVectorDestination.swap(p_new_vector_destination);
    rpInterfaceVectorOrigin.swap(p_new_vector_origin);

    KRATOS_CATCH("")
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Class template instantiation
template class MappingMatrixUtilities< MappingSparseSpaceType, DenseSpaceType >;

}  // namespace Kratos.
