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
#include "matrix_based_mapping_operation_utility.h"
#include "custom_utilities/mapper_typedefs.h"
#include "custom_utilities/mapper_utilities.h"

namespace Kratos
{
typedef typename MapperDefinitions::SparseSpaceType SparseSpaceType;
typedef typename MapperDefinitions::DenseSpaceType DenseSpaceType;

typedef MatrixBasedMappingOperationUtility<SparseSpaceType, DenseSpaceType> UtilityType;

typedef typename MapperLocalSystem::MatrixType MatrixType;
typedef typename MapperLocalSystem::EquationIdVectorType EquationIdVectorType;

typedef std::size_t IndexType;
typedef std::size_t SizeType;

/***********************************************************************************/
/* Functions for internal use in this file */
/***********************************************************************************/
void InitializeSystemVector(UtilityType::TSystemVectorUniquePointerType& rpVector,
                            const SizeType VectorSize)
{
    // The vectors dont have graphs, that why we don't always have to reinitialize them
    if (rpVector == nullptr || rpVector->size() != VectorSize) { //if the pointer is not initialized initialize it to an empty vector
        UtilityType::TSystemVectorUniquePointerType p_new_vector = Kratos::make_unique<UtilityType::TSystemVectorType>(VectorSize);
        rpVector.swap(p_new_vector);

        // TODO do I also have to set to zero the contents?
    }
    else {
        SparseSpaceType::SetToZero(*rpVector);
    }
}

void ConstructMatrixStructure(UtilityType::TSystemMatrixUniquePointerType& rpMdo,
                              UtilityType::MapperLocalSystemPointerVector& rMapperLocalSystems,
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

    auto p_Mdo = Kratos::make_unique<UtilityType::TSystemMatrixType>(
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

void BuildMappingMatrix(UtilityType::TSystemMatrixUniquePointerType& rpMdo,
                       UtilityType::MapperLocalSystemPointerVector& rMapperLocalSystems)
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
UtilityType::MatrixBasedMappingOperationUtility(Parameters Settings)
    : MappingOperationUtility<SparseSpaceType, DenseSpaceType>(Settings)
{
    KRATOS_ERROR_IF(SparseSpaceType::IsDistributed())
        << "Using a distributed Space!" << std::endl;
}

template<>
void UtilityType::BuildMappingSystem(
    TSystemMatrixUniquePointerType& rpMdo,
    TSystemVectorUniquePointerType& rpQo,
    TSystemVectorUniquePointerType& rpQd,
    ModelPart& rModelPartOrigin,
    ModelPart& rModelPartDestination,
    MapperLocalSystemPointerVector& rMapperLocalSystems) const
{
    KRATOS_TRY

    const SizeType num_nodes_origin = rModelPartOrigin.NumberOfNodes();
    const SizeType num_nodes_destination = rModelPartDestination.NumberOfNodes();

    // Initialize the Matrix
    // This has to be done always since the Graph has changed if the Interface is updated!
    ConstructMatrixStructure(rpMdo, rMapperLocalSystems,
                             num_nodes_origin, num_nodes_destination);

    BuildMappingMatrix(rpMdo, rMapperLocalSystems);

    if (GetEchoLevel() > 2) {
        SparseSpaceType::WriteMatrixMarketMatrix("MappingMatrix.mm", *rpMdo, false);
    }

    InitializeSystemVector(rpQo, num_nodes_origin);
    InitializeSystemVector(rpQd, num_nodes_destination);

    KRATOS_CATCH("")
}

template<>
void UtilityType::InitializeMappingStep(
    TSystemMatrixType& rMdo,
    TSystemVectorType& rQo,
    TSystemVectorType& rQd,
    ModelPart& rModelPartOrigin,
    ModelPart& rModelPartDestination,
    const DoubleVariableType& rOriginVariable,
    const DoubleVariableType& rDestinationVariable,
    const Kratos::Flags MappingOptions,
    const bool UseTranspose) const
{
    TInitializeMappingStep(rQo, rQd,
                           rModelPartOrigin, rModelPartDestination,
                           rOriginVariable, rDestinationVariable,
                           MappingOptions, UseTranspose);
    if (GetEchoLevel() > 2) {
        SparseSpaceType::WriteMatrixMarketVector("InitializeMappingStep_Qo.mm", rQo);
        SparseSpaceType::WriteMatrixMarketVector("InitializeMappingStep_Qd.mm", rQd);
    }
}

template<>
void UtilityType::InitializeMappingStep(
    TSystemMatrixType& rMdo,
    TSystemVectorType& rQo,
    TSystemVectorType& rQd,
    ModelPart& rModelPartOrigin,
    ModelPart& rModelPartDestination,
    const ComponentVariableType& rOriginVariable,
    const ComponentVariableType& rDestinationVariable,
    const Kratos::Flags MappingOptions,
    const bool UseTranspose) const
{
    TInitializeMappingStep(rQo, rQd,
                           rModelPartOrigin, rModelPartDestination,
                           rOriginVariable, rDestinationVariable,
                           MappingOptions, UseTranspose);
    if (GetEchoLevel() > 2) {
        SparseSpaceType::WriteMatrixMarketVector("InitializeMappingStep_Qo.mm", rQo);
        SparseSpaceType::WriteMatrixMarketVector("InitializeMappingStep_Qd.mm", rQd);
    }
}

template<>
void UtilityType::FinalizeMappingStep(
    TSystemMatrixType& rMdo,
    TSystemVectorType& rQo,
    TSystemVectorType& rQd,
    ModelPart& rModelPartOrigin,
    ModelPart& rModelPartDestination,
    const DoubleVariableType& rOriginVariable,
    const DoubleVariableType& rDestinationVariable,
    const Kratos::Flags MappingOptions,
    const bool UseTranspose) const
{
    TFinalizeMappingStep(rQo, rQd,
                         rModelPartOrigin, rModelPartDestination,
                         rOriginVariable, rDestinationVariable,
                         MappingOptions, UseTranspose);
    if (GetEchoLevel() > 2) {
        SparseSpaceType::WriteMatrixMarketVector("FinalizeMappingStep_Qo.mm", rQo);
        SparseSpaceType::WriteMatrixMarketVector("FinalizeMappingStep_Qd.mm", rQd);
    }
}

template<>
void UtilityType::FinalizeMappingStep(
    TSystemMatrixType& rMdo,
    TSystemVectorType& rQo,
    TSystemVectorType& rQd,
    ModelPart& rModelPartOrigin,
    ModelPart& rModelPartDestination,
    const ComponentVariableType& rOriginVariable,
    const ComponentVariableType& rDestinationVariable,
    const Kratos::Flags MappingOptions,
    const bool UseTranspose) const
{
    TFinalizeMappingStep(rQo, rQd,
                         rModelPartOrigin, rModelPartDestination,
                         rOriginVariable, rDestinationVariable,
                         MappingOptions, UseTranspose);
    if (GetEchoLevel() > 2) {
        SparseSpaceType::WriteMatrixMarketVector("FinalizeMappingStep_Qo.mm", rQo);
        SparseSpaceType::WriteMatrixMarketVector("FinalizeMappingStep_Qd.mm", rQd);
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Class template instantiation
template class MatrixBasedMappingOperationUtility< SparseSpaceType, DenseSpaceType >;


}  // namespace Kratos.
