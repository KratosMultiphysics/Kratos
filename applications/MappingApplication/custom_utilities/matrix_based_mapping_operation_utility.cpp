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

    using UtilityType = MatrixBasedMappingOperationUtility<MapperDefinitions::SparseSpaceType,
        MapperDefinitions::DenseSpaceType>;
    /***********************************************************************************/
    /* PUBLIC Methods */
    /***********************************************************************************/
    template<>
    UtilityType::MatrixBasedMappingOperationUtility(Parameters Settings)
        : MappingOperationUtility<MapperDefinitions::SparseSpaceType,
          MapperDefinitions::DenseSpaceType>(Settings)
    {
        KRATOS_WATCH("Non-MPI-Ctor")
    }

    template<>
    void UtilityType::ResizeAndInitializeVectors(
        TSystemMatrixTypeUniquePointerType& rpMdo,
        TSystemVectorTypeUniquePointerType& rpQo,
        TSystemVectorTypeUniquePointerType& rpQd,
        ModelPart& rModelPartOrigin,
        ModelPart& rModelPartDestination) const
    {
        KRATOS_TRY

        if (rpMdo == nullptr) //if the pointer is not initialized initialize it to an empty matrix
        {
            TSystemMatrixTypeUniquePointerType p_Mdo = Kratos::make_unique<TSystemMatrixType>(0,0);
            rpMdo.swap(p_Mdo);
        }
        if (rpQo == nullptr) //if the pointer is not initialized initialize it to an empty vector
        {
            TSystemVectorTypeUniquePointerType p_Do = Kratos::make_unique<TSystemVectorType>(0);
            rpQo.swap(p_Do);
        }
        if (rpQd == nullptr) //if the pointer is not initialized initialize it to an empty vector
        {
            TSystemVectorTypeUniquePointerType p_Dd = Kratos::make_unique<TSystemVectorType>(0);
            rpQd.swap(p_Dd);
        }

        TSystemMatrixType& r_Mdo = *rpMdo;
        TSystemVectorType& r_Qo = *rpQo;
        TSystemVectorType& r_Qd = *rpQd;

        const std::size_t num_nodes_origin = rModelPartOrigin.NumberOfNodes();
        const std::size_t num_nodes_destination = rModelPartDestination.NumberOfNodes();

        // TODO CHECK THIS!!!
        if (r_Mdo.size1() != num_nodes_origin || r_Mdo.size2() != num_nodes_destination)
            r_Mdo.resize(num_nodes_origin, num_nodes_destination, false);

        if (r_Qo.size() != num_nodes_origin)
            r_Qo.resize(num_nodes_origin, false);

        if (r_Qd.size() != num_nodes_destination)
            r_Qd.resize(num_nodes_destination, false);

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
template class MatrixBasedMappingOperationUtility< MapperDefinitions::SparseSpaceType, MapperDefinitions::DenseSpaceType >;


}  // namespace Kratos.
