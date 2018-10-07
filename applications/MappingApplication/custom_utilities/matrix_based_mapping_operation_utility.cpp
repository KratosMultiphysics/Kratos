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
#include "matrix_based_mapping_operation_utility.h"
#include "custom_utilities/mapper_typedefs.h"
#include "custom_utilities/mapper_utilities.h"

namespace Kratos
{
using SparseSpaceType = MapperDefinitions::SparseSpaceType;
using DenseSpaceType = MapperDefinitions::DenseSpaceType;

using UtilityType = MatrixBasedMappingOperationUtility<SparseSpaceType, DenseSpaceType>;

using EquationIdVectorType = typename MapperLocalSystem::EquationIdVectorType;
typedef typename MapperLocalSystem::MatrixType MatrixType;

using SizeType = std::size_t;
using IndexType = std::size_t;


/***********************************************************************************/
/* Functions for internal use in this file */
/***********************************************************************************/
void InitializeVector(UtilityType::TSystemVectorUniquePointerType& rpVector,
                        const SizeType VectorSize)
{
    // The vectors dont have graphs, that why we don't always have to reinitialize them
    if (rpVector == nullptr || rpVector->size() != VectorSize) //if the pointer is not initialized initialize it to an empty vector
    {
        UtilityType::TSystemVectorUniquePointerType p_new_vector = Kratos::make_unique<UtilityType::TSystemVectorType>(VectorSize);
        rpVector.swap(p_new_vector);

        // TODO do I also have to set to zero the contents?
    }
    else
    {
        SparseSpaceType::SetToZero(*rpVector);
    }
}

void ConstructMatrixStructure(UtilityType::MapperLocalSystemPointerVector& rMapperLocalSystems,
                              UtilityType::TSystemMatrixType& rMdo)
{
    // A = boost::numeric::ublas::compressed_matrix<double>(indices.size(), indices.size(), nnz);
    EquationIdVectorType origin_ids;
    EquationIdVectorType destination_ids;

    // TODO omp
    for (/*const*/auto& r_local_sys : rMapperLocalSystems) // TODO I think this can be const bcs it is the ptr
    {
        r_local_sys->EquationIdVectors(origin_ids, destination_ids);

    }
}

template< class TVarType >
void FillSystemVector(UtilityType::TSystemVectorType& rVector,
                        ModelPart& rModelPart,
                        const TVarType& rVariable,
                        const Kratos::Flags& rMappingOptions)
{
    // Here we construct a function pointer to not have the if all the time inside the loop
    const auto fill_fct = MapperUtilities::GetFillFunction<TVarType>(rMappingOptions);

    const auto nodes_begin = rModelPart.NodesBegin();

    #pragma omp parallel for
    for (int i = 0; i<static_cast<int>(rModelPart.NumberOfNodes()); i++)
        fill_fct(*(nodes_begin + i), rVariable, rVector[i]);
}

template< class TVarType >
void Update(UtilityType::TSystemVectorType& rVector,
            ModelPart& rModelPart,
            const TVarType& rVariable,
            const Kratos::Flags& rMappingOptions)
{
    const double factor = rMappingOptions.Is(MapperFlags::SWAP_SIGN) ? -1.0 : 1.0;

    // Here we construct a function pointer to not have the if all the time inside the loop
    const auto update_fct = std::bind(MapperUtilities::GetUpdateFunction<TVarType>(rMappingOptions),
                                        std::placeholders::_1,
                                        std::placeholders::_2,
                                        std::placeholders::_3,
                                        factor);

    const auto nodes_begin = rModelPart.NodesBegin();

    #pragma omp parallel for
    for (int i = 0; i<static_cast<int>(rModelPart.NumberOfNodes()); i++)
        update_fct(*(nodes_begin + i), rVariable, rVector[i]);
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
void UtilityType::ResizeAndInitializeVectors(
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
    const SizeType num_non_zeros = 100; // TODO this should be computed

    // ConstructMatrixStructure(rpMdo, rMapperLocalSystems);

    TSystemMatrixUniquePointerType p_Mdo = Kratos::make_unique<TSystemMatrixType>(
        num_nodes_destination, num_nodes_origin, num_non_zeros);
    rpMdo.swap(p_Mdo);

    // TODO do I also have to set to zero the contents?
    // SparseSpaceType::SetToZero(*rpMdo);

    InitializeVector(rpQo, num_nodes_origin);
    InitializeVector(rpQd, num_nodes_destination);

    KRATOS_CATCH("")
}

// The "Build" function
template<>
void UtilityType::BuildMappingMatrix(
    const MapperLocalSystemPointerVector& rMapperLocalSystems,
    TSystemMatrixType& rMdo) const
{
    MatrixType local_mapping_matrix;

    EquationIdVectorType origin_ids;
    EquationIdVectorType destination_ids;

    for (auto& r_local_sys : rMapperLocalSystems) // TODO omp
    {
        r_local_sys->CalculateLocalSystem(local_mapping_matrix, origin_ids, destination_ids);
        KRATOS_DEBUG_ERROR_IF(local_mapping_matrix.size1() != destination_ids.size())
            << "DestinationID vector size mismatch" << std::endl;
        KRATOS_DEBUG_ERROR_IF(local_mapping_matrix.size2() != origin_ids.size())
            << "OriginID vector size mismatch" << std::endl;

        // Insert the mapping weights from the local_systems into the mapping matrix
        for (IndexType i=0; i<destination_ids.size(); ++i) {
            for (IndexType j=0; j<origin_ids.size(); ++j) {
                rMdo(destination_ids[i], origin_ids[j]) += local_mapping_matrix(i,j);
            }
        }

        r_local_sys->Clear();
    }

    std::cout << "BuildMappingMatrix, non-mpi: " << "Leaving" << std::endl;

    if (GetEchoLevel() > 2)
        SparseSpaceType::WriteMatrixMarketMatrix("MappingMatrix", rMdo, false);
}

template<class TVarType>
void TInitializeMappingStep(UtilityType::TSystemMatrixType& rMdo,
    UtilityType::TSystemVectorType& rQo,
    UtilityType::TSystemVectorType& rQd,
    ModelPart& rModelPartOrigin,
    ModelPart& rModelPartDestination,
    const TVarType& rOriginVariable,
    const TVarType& rDestinationVariable,
    const Kratos::Flags MappingOptions,
    const bool UseTranspose)
{
    if (UseTranspose)
        FillSystemVector(rQd, rModelPartDestination, rDestinationVariable, MappingOptions);
    else
        FillSystemVector(rQo, rModelPartOrigin, rOriginVariable, MappingOptions);
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
    TInitializeMappingStep(rMdo, rQo, rQd,
                            rModelPartOrigin, rModelPartDestination,
                            rOriginVariable, rDestinationVariable,
                            MappingOptions, UseTranspose);
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
    TInitializeMappingStep(rMdo, rQo, rQd,
                            rModelPartOrigin, rModelPartDestination,
                            rOriginVariable, rDestinationVariable,
                            MappingOptions, UseTranspose);
}

void ExecuteMapping(UtilityType::TSystemMatrixType& rMdo,
                    UtilityType::TSystemVectorType& rQo,
                    UtilityType::TSystemVectorType& rQd,
                    const bool UseTranspose)
{
    if (UseTranspose)
        SparseSpaceType::TransposeMult(rMdo, rQd, rQo); // rQo = rMdo^T * rQo
    else
        SparseSpaceType::Mult(rMdo, rQo, rQd); // rQd = rMdo * rQo
}

template<>
void UtilityType::ExecuteMappingStep(
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
    ExecuteMapping(rMdo, rQo, rQd, UseTranspose);
}

template<>
void UtilityType::ExecuteMappingStep(
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
    ExecuteMapping(rMdo, rQo, rQd, UseTranspose);
}

template<class TVarType>
void TFinalizeMappingStep(UtilityType::TSystemMatrixType& rMdo,
    UtilityType::TSystemVectorType& rQo,
    UtilityType::TSystemVectorType& rQd,
    ModelPart& rModelPartOrigin,
    ModelPart& rModelPartDestination,
    const TVarType& rOriginVariable,
    const TVarType& rDestinationVariable,
    const Kratos::Flags MappingOptions,
    const bool UseTranspose)
{
    if (UseTranspose)
        Update(rQo, rModelPartOrigin, rOriginVariable, MappingOptions);
    else
        Update(rQd, rModelPartDestination, rDestinationVariable, MappingOptions);
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
    TFinalizeMappingStep(rMdo, rQo, rQd,
                            rModelPartOrigin, rModelPartDestination,
                            rOriginVariable, rDestinationVariable,
                            MappingOptions, UseTranspose);
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
    TFinalizeMappingStep(rMdo, rQo, rQd,
                            rModelPartOrigin, rModelPartDestination,
                            rOriginVariable, rDestinationVariable,
                            MappingOptions, UseTranspose);
}

/***********************************************************************************/
/* PROTECTED Methods */
/***********************************************************************************/


/***********************************************************************************/
/* PRIVATE Methods */
/***********************************************************************************/

///////////////////////////////////////////////////////////////////////////////////////////////////
// Class template instantiation
template class MatrixBasedMappingOperationUtility< SparseSpaceType, DenseSpaceType >;


}  // namespace Kratos.
