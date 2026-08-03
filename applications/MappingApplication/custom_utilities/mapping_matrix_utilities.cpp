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

// External includes

// Project includes
#include "mapping_matrix_utilities.h"
#include "mappers/mapper_define.h"
#include "custom_utilities/mapper_utilities.h"
#include "utilities/parallel_utilities.h"
#include "containers/sparse_contiguous_row_graph.h"

namespace Kratos {

namespace {

using MappingSparseSpaceType = typename MapperDefinitions::SparseSpaceType;
using DenseSpaceType = typename MapperDefinitions::DenseSpaceType;

using MappingMatrixUtilitiesType = MappingMatrixUtilities<MappingSparseSpaceType, DenseSpaceType>;

using MatrixType = typename MapperLocalSystem::MatrixType;
using EquationIdVectorType = typename MapperLocalSystem::EquationIdVectorType;

using IndexType = std::size_t;
using SizeType = std::size_t;

/***********************************************************************************/
/* Functions for internal use in this file */
/***********************************************************************************/

void ConstructMatrixStructure(Kratos::unique_ptr<typename MappingSparseSpaceType::MatrixType>& rpMdo,
                              std::vector<Kratos::unique_ptr<MapperLocalSystem>>& rMapperLocalSystems,
                              const SizeType NumNodesOrigin,
                              const SizeType NumNodesDestination)
{
    // Single shared graph with per-row locks — one unordered_set per destination
    // row, locked independently. Memory is O(NumNodesDestination + nnz) regardless
    // of thread count or number of local systems.
    SparseContiguousRowGraph<IndexType> graph(NumNodesDestination);

    struct TLSData
    {
        EquationIdVectorType OriginIds;
        EquationIdVectorType DestinationIds;
    };

    block_for_each(rMapperLocalSystems.begin(), rMapperLocalSystems.end(), TLSData(),
        [&graph](auto& rpLocalSys, TLSData& rTLS) {
            rTLS.OriginIds.clear();
            rTLS.DestinationIds.clear();

            rpLocalSys->EquationIdVectors(rTLS.OriginIds, rTLS.DestinationIds);

            for (const auto id_dest : rTLS.DestinationIds) {
                graph.AddEntries(id_dest, rTLS.OriginIds);
            }
        }
    );

    // Export sorted CSR arrays (column indices are sorted per row by ExportCSRArrays).
    // ExportCSRArrays allocates with new[] and transfers ownership to the caller.
    IndexType* p_row_data = nullptr;
    IndexType row_data_size = 0;
    IndexType* p_col_data = nullptr;
    IndexType col_data_size = 0;
    graph.ExportCSRArrays(p_row_data, row_data_size, p_col_data, col_data_size);

    const SizeType num_non_zero_entries = col_data_size;

    auto p_Mdo = Kratos::make_unique<typename MappingSparseSpaceType::MatrixType>(
        NumNodesDestination,
        NumNodesOrigin,
        num_non_zero_entries);

    IndexType* p_matrix_row_indices = p_Mdo->index1_data().begin();
    IndexType* p_matrix_col_indices = p_Mdo->index2_data().begin();
    double*    p_matrix_values       = p_Mdo->value_data().begin();

    IndexPartition<IndexType>(NumNodesDestination + 1).for_each([&](IndexType i) {
        p_matrix_row_indices[i] = p_row_data[i];
    });
    IndexPartition<IndexType>(num_non_zero_entries).for_each([&](IndexType i) {
        p_matrix_col_indices[i] = p_col_data[i];
        p_matrix_values[i] = 0.0;
    });

    delete[] p_row_data;
    delete[] p_col_data;

    p_Mdo->set_filled(NumNodesDestination + 1, num_non_zero_entries);

    rpMdo.swap(p_Mdo);
}

void BuildMatrix(Kratos::unique_ptr<typename MappingSparseSpaceType::MatrixType>& rpMdo,
                 std::vector<Kratos::unique_ptr<MapperLocalSystem>>& rMapperLocalSystems)
{
    struct TLS {
        MatrixType local_mapping_matrix;
        EquationIdVectorType origin_ids;
        EquationIdVectorType destination_ids;
    };

    block_for_each(rMapperLocalSystems, TLS(), [&rpMdo] (auto& r_local_sys, TLS& rTls){

        r_local_sys->CalculateLocalSystem(rTls.local_mapping_matrix, rTls.origin_ids, rTls.destination_ids);

        KRATOS_DEBUG_ERROR_IF(rTls.local_mapping_matrix.size1() != rTls.destination_ids.size()) << "MappingMatrixAssembly: DestinationID vector size mismatch: LocalMappingMatrix-Size1: " << rTls.local_mapping_matrix.size1() << " | DestinationIDs-size: " << rTls.destination_ids.size() << std::endl;
        KRATOS_DEBUG_ERROR_IF(rTls.local_mapping_matrix.size2() != rTls.origin_ids.size()) << "MappingMatrixAssembly: OriginID vector size mismatch: LocalMappingMatrix-Size2: " << rTls.local_mapping_matrix.size2() << " | OriginIDs-size: " << rTls.origin_ids.size() << std::endl;

        for (IndexType i = 0; i < rTls.destination_ids.size(); ++i) {
            for (IndexType j = 0; j < rTls.origin_ids.size(); ++j) {
                AtomicAdd((*rpMdo)(rTls.destination_ids[i], rTls.origin_ids[j]).ref(), rTls.local_mapping_matrix(i,j));
            }
        }

        r_local_sys->Clear();
    });

} 

} // anonymous namespace

template<>
void MappingMatrixUtilitiesType::CheckRowSum(
    const typename MappingSparseSpaceType::MatrixType& rM,
    const std::string& rBaseFileName,
    const bool ThrowError,
    const double Tolerance)
{
    typename MappingSparseSpaceType::VectorType unit_vector(MappingSparseSpaceType::Size2(rM));
    MappingSparseSpaceType::Set(unit_vector, 1.0);

    typename MappingSparseSpaceType::VectorType row_sums_vector(MappingSparseSpaceType::Size1(rM));

    MappingSparseSpaceType::Mult(rM, unit_vector, row_sums_vector);

    bool write_mm_file = false;
    for (std::size_t i = 0; i < MappingSparseSpaceType::Size(row_sums_vector); ++i) {
        if (std::abs(row_sums_vector[i] - 1.0) > Tolerance) {
            KRATOS_WARNING("MappingMatrixAssembly") << "The row sum in row " << i << " is unequal 1.0: " << row_sums_vector[i] << std::endl;
            write_mm_file = true;
        }
    }

    if (write_mm_file) {
        MappingSparseSpaceType::WriteMatrixMarketVector(("RowSumVector_" + rBaseFileName).c_str(), row_sums_vector);
        KRATOS_ERROR_IF(ThrowError) << "Mapping matrix does not sum to unity. Please check file " << rBaseFileName << " in your project directory for row sums\n";
    }
}

template<>
void MappingMatrixUtilitiesType::InitializeSystemVector(
    Kratos::unique_ptr<typename MappingSparseSpaceType::VectorType>& rpVector,
    const std::size_t VectorSize)
{
    // The vectors dont have graphs, that why we don't always have to reinitialize them
    if (rpVector == nullptr || rpVector->size() != VectorSize) { //if the pointer is not initialized initialize it to an empty vector
        Kratos::unique_ptr<typename MappingSparseSpaceType::VectorType> p_new_vector = Kratos::make_unique<typename MappingSparseSpaceType::VectorType>(VectorSize);
        rpVector.swap(p_new_vector);

        // TODO do I also have to set to zero the contents?
    }
    else {
        MappingSparseSpaceType::SetToZero(*rpVector);
    }
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

    static_assert(!MappingSparseSpaceType::IsDistributed(), "Using a distributed Space!");

    const SizeType num_nodes_origin = rModelPartOrigin.NumberOfNodes();
    const SizeType num_nodes_destination = rModelPartDestination.NumberOfNodes();

    // Initialize the Matrix
    // This has to be done always since the Graph has changed if the Interface is updated!
    ConstructMatrixStructure(rpMappingMatrix, rMapperLocalSystems,
                             num_nodes_origin, num_nodes_destination);

    BuildMatrix(rpMappingMatrix, rMapperLocalSystems);

    MappingMatrixUtilitiesType::InitializeSystemVector(rpInterfaceVectorOrigin, num_nodes_origin);
    MappingMatrixUtilitiesType::InitializeSystemVector(rpInterfaceVectorDestination, num_nodes_destination);

    KRATOS_CATCH("")
}

template<>
void MappingMatrixUtilitiesType::BuildMappingMatrixRBFMapper(
    Kratos::unique_ptr<typename MappingSparseSpaceType::MatrixType>& rpMappingMatrix,
    Kratos::unique_ptr<typename MappingSparseSpaceType::VectorType>& rpInterfaceVectorOrigin,
    Kratos::unique_ptr<typename MappingSparseSpaceType::VectorType>& rpInterfaceVectorDestination,
    const ModelPart& rModelPartOrigin,
    const ModelPart& rModelPartDestination,
    std::vector<Kratos::unique_ptr<MapperLocalSystem>>& rMapperLocalSystems,
    const IndexType NumberOfPolynomialTerms,
    const bool BuildOriginInterpolationMatrix,
    const bool OriginIsIga,
    const int EchoLevel)
{
    KRATOS_TRY

    static_assert(!MappingSparseSpaceType::IsDistributed(), "Using a distributed Space!");

    const SizeType num_nodes_origin = rModelPartOrigin.NumberOfNodes();
    const SizeType num_conditions_origin = rModelPartOrigin.NumberOfConditions();
    const SizeType num_nodes_destination = rModelPartDestination.NumberOfNodes();
    const SizeType num_conditions_destination = rModelPartDestination.NumberOfConditions();

    IndexType origin_size;
    IndexType destination_size;

    if (!OriginIsIga){
        origin_size = num_nodes_origin;
        destination_size = num_nodes_destination;
    } else if (OriginIsIga && BuildOriginInterpolationMatrix) {
        origin_size = num_conditions_origin;
        destination_size = num_conditions_destination;
    } else if (OriginIsIga && !BuildOriginInterpolationMatrix) {
        origin_size = num_conditions_origin;
        destination_size = num_nodes_destination;
    }
    
    // Initialize the Matrix
    // This has to be done always since the Graph has changed if the Interface is updated!
    if (BuildOriginInterpolationMatrix){
        ConstructMatrixStructure(rpMappingMatrix, rMapperLocalSystems,
                                origin_size + NumberOfPolynomialTerms, destination_size + NumberOfPolynomialTerms);
    } else {
        ConstructMatrixStructure(rpMappingMatrix, rMapperLocalSystems,
                                origin_size + NumberOfPolynomialTerms, destination_size);
    }

    BuildMatrix(rpMappingMatrix, rMapperLocalSystems);

    MappingMatrixUtilitiesType::InitializeSystemVector(rpInterfaceVectorOrigin, num_nodes_origin);
    MappingMatrixUtilitiesType::InitializeSystemVector(rpInterfaceVectorDestination, num_nodes_destination);

    KRATOS_CATCH("")
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Class template instantiation
template class MappingMatrixUtilities< MappingSparseSpaceType, DenseSpaceType >;

}  // namespace Kratos.
