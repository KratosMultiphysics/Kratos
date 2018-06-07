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
    using MappingWeightsVector = typename MapperLocalSystem::MappingWeightsVector;

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

    template< class TVarType>
    void FillSystemVector(UtilityType::TSystemVectorType& rVector,
                          ModelPart& rModelPart,
                          const TVarType& rVariable,
                          const Kratos::Flags& rMappingOptions)
    {
        #pragma omp parallel for
        for (int i = 0; i<static_cast<int>(rModelPart.NumberOfNodes()); i++)
        {
            auto it = rModelPart.NodesBegin() + i;
            rVector[i] = it->FastGetSolutionStepValue(rVariable);
        }
    }

    template< class TVarType>
    void Update(UtilityType::TSystemVectorType& rVector,
                ModelPart& rModelPart,
                const TVarType& rVariable,
                const Kratos::Flags& rMappingOptions)
    {
        // TODO create a function ptr to not have the if all the time
        #pragma omp parallel for
        for (int i = 0; i<static_cast<int>(rModelPart.NumberOfNodes()); i++)
        {
            auto it = rModelPart.NodesBegin() + i;

            if (rMappingOptions.Is(MapperFlags::ADD_VALUES))
                it->FastGetSolutionStepValue(rVariable) += rVector[i] /** Factor*/;
            else
                it->FastGetSolutionStepValue(rVariable) = rVector[i] /** Factor*/;
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
        MappingWeightsVector mapping_weights;

        EquationIdVectorType origin_ids;
        EquationIdVectorType destination_ids;

        KRATOS_INFO("BuildMappingMatrix, non-mpi") << "Entering" << std::endl;

        for (auto& r_local_sys : rMapperLocalSystems) // TODO omp
        {
            r_local_sys->CalculateLocalSystem(mapping_weights, origin_ids, destination_ids);
            KRATOS_DEBUG_ERROR_IF(mapping_weights.size() != origin_ids.size()) << "OriginID vector size mismatch" << std::endl;
            KRATOS_DEBUG_ERROR_IF(mapping_weights.size() != destination_ids.size()) << "DestinationID vector size mismatch" << std::endl;

            // Insert the mapping weights from the local_systems into the mapping matrix
            for (IndexType i=0; i<mapping_weights.size(); ++i)
                rMdo(origin_ids[i], destination_ids[i]) += mapping_weights[i];

            r_local_sys->Clear();
        }

        if (GetEchoLevel() > 2)
            SparseSpaceType::WriteMatrixMarketMatrix("MappingMatrix", rMdo, false);
    }

    // The "Solve" function
    template<>
    void UtilityType::ExecuteMapping(
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
        if (UseTranspose)
        {
            FillSystemVector(rQd, rModelPartDestination, rDestinationVariable, MappingOptions);
            SparseSpaceType::TransposeMult(rMdo, rQd, rQo); // rQo = rMdo^T * rQo
            Update(rQo, rModelPartOrigin, rOriginVariable, MappingOptions);
        }
        else
        {
            KRATOS_INFO("BuildMappingMatrix, non-mpi") << "Doing the mapping steps" << std::endl;
            FillSystemVector(rQo, rModelPartOrigin, rOriginVariable, MappingOptions);
            SparseSpaceType::Mult(rMdo, rQo, rQd); // rQd = rMdo * rQo
            Update(rQd, rModelPartDestination, rDestinationVariable, MappingOptions);
        }
    }

    // The "Solve" function
    template<>
    void UtilityType::ExecuteMapping(
        TSystemMatrixType& rMdo,
        TSystemVectorType& rQo,
        TSystemVectorType& rQd,
        ModelPart& rModelPartOrigin,
        ModelPart& rModelPartDestination,
        const Array3VariableType& rOriginVariable,
        const Array3VariableType& rDestinationVariable,
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
