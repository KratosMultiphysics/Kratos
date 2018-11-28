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
#include <unordered_set>

// External includes

// Project includes
#include "mapping_matrix_builder.h"
#include "custom_utilities/mapper_typedefs.h"
#include "custom_utilities/mapper_utilities.h"

namespace Kratos
{
typedef typename MapperDefinitions::SparseSpaceType SparseSpaceType;
typedef typename MapperDefinitions::DenseSpaceType DenseSpaceType;

typedef MappingMatrixBuilder<SparseSpaceType, DenseSpaceType> BuilderType;

typedef typename MapperLocalSystem::MatrixType MatrixType;
typedef typename MapperLocalSystem::EquationIdVectorType EquationIdVectorType;

typedef std::size_t IndexType;
typedef std::size_t SizeType;

/***********************************************************************************/
/* Functions for internal use in this file */
/***********************************************************************************/
void InitializeSystemVector(BuilderType::TSystemVectorUniquePointerType& rpVector,
                            const SizeType VectorSize)
{
    // The vectors dont have graphs, that why we don't always have to reinitialize them
    if (rpVector == nullptr || rpVector->size() != VectorSize) { //if the pointer is not initialized initialize it to an empty vector
        BuilderType::TSystemVectorUniquePointerType p_new_vector = Kratos::make_unique<BuilderType::TSystemVectorType>(VectorSize);
        rpVector.swap(p_new_vector);

        // TODO do I also have to set to zero the contents?
    }
    else {
        SparseSpaceType::SetToZero(*rpVector);
    }
}

void ConstructMatrixStructure(BuilderType::TMappingMatrixUniquePointerType& rpMdo,
                              BuilderType::MapperLocalSystemPointerVector& rMapperLocalSystems,
                              const SizeType NumNodesOrigin,
                              const SizeType NumNodesDestination)
{
    // one set for each row storing the corresponding col-IDs
    std::vector<std::unordered_set<IndexType> > indices(NumNodesDestination);

    // preallocate memory for the column indices
    for (IndexType i=0; i<NumNodesDestination; ++i) {
        // TODO I guess this can be optimized...
        // this highly depends on the used mapper => same goes for the Graph in Trilinos
        indices[i].reserve(3);
    }

    EquationIdVectorType origin_ids;
    EquationIdVectorType destination_ids;

    // Looping the local-systems to get the entries for the matrix
    // TODO omp
    for (/*const*/auto& r_local_sys : rMapperLocalSystems) { // TODO I think this can be const bcs it is the ptr
        r_local_sys->EquationIdVectors(origin_ids, destination_ids);
        for (const auto dest_idx : destination_ids) {
            indices[dest_idx].insert(origin_ids.begin(), origin_ids.end());
        }
    }

    // computing the number of non-zero entries
    SizeType num_non_zero_entries = 0;
    for (const auto& r_row_indices : indices) { // looping the indices per row
        num_non_zero_entries += r_row_indices.size(); // adding the number of col-indices
    }

    auto p_Mdo = Kratos::make_unique<BuilderType::TMappingMatrixType>(
        NumNodesDestination,
        NumNodesOrigin,
        num_non_zero_entries);

    double* p_matrix_values = p_Mdo->value_data().begin();
    IndexType* p_matrix_row_indices = p_Mdo->index1_data().begin();
    IndexType* p_matrix_col_indices = p_Mdo->index2_data().begin();

    // filling the index1 vector - do NOT make parallel!
    p_matrix_row_indices[0] = 0;
    for (IndexType i=0; i<NumNodesDestination; ++i) {
        p_matrix_row_indices[i+1] = p_matrix_row_indices[i] + indices[i].size();
    }

    for (IndexType i=0; i<NumNodesDestination; ++i) {
        const IndexType row_begin = p_matrix_row_indices[i];
        const IndexType row_end = p_matrix_row_indices[i+1];
        IndexType j = row_begin;
        for (const auto index : indices[i]) {
            p_matrix_col_indices[j] = index;
            p_matrix_values[j] = 0.0;
            ++j;
        }

        indices[i].clear(); //deallocating the memory // TODO necessary?

        std::sort(&p_matrix_col_indices[row_begin], &p_matrix_col_indices[row_end]);
    }

    p_Mdo->set_filled(indices.size()+1, num_non_zero_entries);

    rpMdo.swap(p_Mdo);
}

void BuildMatrix(BuilderType::TMappingMatrixUniquePointerType& rpMdo,
                       BuilderType::MapperLocalSystemPointerVector& rMapperLocalSystems)
{
    MatrixType local_mapping_matrix;
    EquationIdVectorType origin_ids;
    EquationIdVectorType destination_ids;

    for (auto& r_local_sys : rMapperLocalSystems) { // TODO omp

        r_local_sys->CalculateLocalSystem(local_mapping_matrix, origin_ids, destination_ids);

        KRATOS_DEBUG_ERROR_IF(local_mapping_matrix.size1() != destination_ids.size())
            << "DestinationID vector size mismatch" << std::endl;
        KRATOS_DEBUG_ERROR_IF(local_mapping_matrix.size2() != origin_ids.size())
            << "OriginID vector size mismatch" << std::endl;

        for (IndexType i=0; i<destination_ids.size(); ++i) {
            for (IndexType j=0; j<origin_ids.size(); ++j) {
                // #pragma omp atomic
                (*rpMdo)(destination_ids[i], origin_ids[j]) += local_mapping_matrix(i,j);
            }
        }

        r_local_sys->Clear();
    }
}

/***********************************************************************************/
/* PUBLIC Methods */
/***********************************************************************************/
template<>
BuilderType::MappingMatrixBuilder()
{
    KRATOS_ERROR_IF(SparseSpaceType::IsDistributed())
        << "Using a distributed Space!" << std::endl;
}

template<>
void BuilderType::BuildMappingMatrix(
    const InterfaceVectorContainerPointerType& rpVectorContainerOrigin,
    const InterfaceVectorContainerPointerType& rpVectorContainerDestination,
    MapperLocalSystemPointerVector& rMapperLocalSystems)
{
    KRATOS_TRY

    const SizeType num_nodes_origin = rpVectorContainerOrigin->GetModelPart().NumberOfNodes();
    const SizeType num_nodes_destination = rpVectorContainerDestination->GetModelPart().NumberOfNodes();

    // Initialize the Matrix
    // This has to be done always since the Graph has changed if the Interface is updated!
    ConstructMatrixStructure(mpMappingMatrix, rMapperLocalSystems,
                             num_nodes_origin, num_nodes_destination);

    BuildMatrix(mpMappingMatrix, rMapperLocalSystems);

    if (mEchoLevel > 2) {
        SparseSpaceType::WriteMatrixMarketMatrix("MappingMatrix.mm", *mpMappingMatrix, false);
    }

    InitializeSystemVector(rpVectorContainerOrigin->pGetVector(), num_nodes_origin);
    InitializeSystemVector(rpVectorContainerDestination->pGetVector(), num_nodes_destination);

    KRATOS_CATCH("")
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Class template instantiation
template class MappingMatrixBuilder< SparseSpaceType, DenseSpaceType >;


}  // namespace Kratos.
