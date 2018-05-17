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

namespace Kratos
{
    using SparseSpaceType = MapperDefinitions::SparseSpaceType;
    using DenseSpaceType = MapperDefinitions::DenseSpaceType;

    using UtilityType = MatrixBasedMappingOperationUtility<SparseSpaceType, DenseSpaceType>;

    using EquationIdVectorType = typename MapperLocalSystem::EquationIdVectorType;

    using SizeType = std::size_t;
    using IndexType = std::size_t;



    void InitializeVector(UtilityType::TSystemVectorUniquePointerType& rpVector,
                         const SizeType VectorSize)
    {
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

        for (/*const*/auto& r_local_sys : rMapperLocalSystems) // TODO I think this can be const bcs it is the ptr
        {
            r_local_sys->EquationIdVectors(origin_ids, destination_ids);

        }
    }

    /***********************************************************************************/
    /* PUBLIC Methods */
    /***********************************************************************************/
    template<>
    UtilityType::MatrixBasedMappingOperationUtility(Parameters Settings)
        : MappingOperationUtility<SparseSpaceType, DenseSpaceType>(Settings)
    {
        KRATOS_WATCH("Non-MPI-Ctor")
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

        if (rpMdo == nullptr || rpMdo->size1() != num_nodes_origin || rpMdo->size2() != num_nodes_destination) //if the pointer is not initialized initialize it to an empty matrix
        {
            const SizeType num_non_zeros = 100; // TODO this should be computed

            // ConstructMatrixStructure(rpMdo, rMapperLocalSystems);

            TSystemMatrixUniquePointerType p_Mdo = Kratos::make_unique<TSystemMatrixType>(
                num_nodes_origin,num_nodes_destination, num_non_zeros);
            rpMdo.swap(p_Mdo);

            // TODO do I also have to set to zero the contents?
        }
        else
        {
            TSparseSpace::SetToZero(*rpMdo);
        }

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

    }

    // The "Solve" function
    template<>
    void UtilityType::ExecuteMapping(
        const TSystemMatrixType& rMdo,
        TSystemVectorType& rQo,
        TSystemVectorType& rQd,
        ModelPart& rModelPartOrigin,
        ModelPart& rModelPartDestination,
        const Variable<double>& rOriginVariable,
        const Variable<double>& rDestinationVariable,
        const Kratos::Flags MappingOptions,
        const bool UseTranspose) const
    {

    }

    // The "Solve" function
    template<>
    void UtilityType::ExecuteMapping(
        const TSystemMatrixType& rMdo,
        TSystemVectorType& rQo,
        TSystemVectorType& rQd,
        ModelPart& rModelPartOrigin,
        ModelPart& rModelPartDestination,
        const Variable<array_1d<double, 3>>& rOriginVariable,
        const Variable<array_1d<double, 3>>& rDestinationVariable,
        const Kratos::Flags MappingOptions,
        const bool UseTranspose) const
    {

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
