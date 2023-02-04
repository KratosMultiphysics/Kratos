//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 license: HDF5Application/license.txt
//
//  Main author:     Suneth Warnakulasuriya
//

#pragma once

// System includes
#include <string>
#include <functional>
#include <vector>

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "spatial_containers/spatial_containers.h"

// Application includes
#include "custom_utilities/mappers/container_variable_data_mapper.h"
#include "custom_utilities/mappers/entity_point.h"
#include "custom_utilities/mappers/mapping_filter_functions.h"

namespace Kratos {

///@name Kratos Classes
///@{

template<class TContainerType>
class KRATOS_API(OPTIMIZATION_APPLICATION) VertexMorphingContainerVariableDataMapper : public ContainerVariableDataMapper<TContainerType>
{
public:
    ///@name Type definitions
    ///@{

    using BaseType = ContainerVariableDataMapper<TContainerType>;

    using ContainerVariableDataHolderType = ContainerVariableDataHolderBase<TContainerType>;

    using TEntityType = typename TContainerType::data_type;

    using EntityPointType = EntityPoint<TEntityType>;

    using EntityPointTypePointer = typename EntityPointType::Pointer;

    using EntityPointVector = std::vector<EntityPointTypePointer>;

    using EntityPointIterator = typename std::vector<EntityPointTypePointer>::iterator;

    using DoubleVectorIterator = std::vector<double>::iterator;

    using BucketType = Bucket<3, EntityPointType, EntityPointVector, EntityPointTypePointer, EntityPointIterator, DoubleVectorIterator>;

    using KDTree = Tree<KDTreePartition<BucketType>>;

    using SparseSpaceType = UblasSpace<double, CompressedMatrix, Vector>;

    using SparseMatrixType = SparseSpaceType::MatrixType;

    /// Pointer definition of ContainerMapper
    KRATOS_CLASS_POINTER_DEFINITION(VertexMorphingContainerVariableDataMapper);

    ///@}
    ///@name LifeCycle
    ///@{

    VertexMorphingContainerVariableDataMapper(
        const ModelPart& rOriginModelPart,
        const ModelPart& rDestinationModelPart,
        Parameters Params);

    ~VertexMorphingContainerVariableDataMapper() override = default;

    ///@}
    ///@name Public operations

    void Update() override;

    void Map(
        const ContainerVariableDataHolderType& rOriginDataContainer,
        ContainerVariableDataHolderType& rDestinationDataContainer) const override;

    void InverseMap(
        ContainerVariableDataHolderType& rOriginDataContainer,
        const ContainerVariableDataHolderType& rDestinationDataContainer) const override;

    std::string Info() const override;

    ///@}
private:
    ///@name Private member variables
    ///@{

    const ModelPart& mrOriginModelPart;

    const ModelPart& mrDestinationModelPart;

    IndexType mMaxEntitiesInFilterRadius;

    std::string mFilterFunctionType;

    double mFilterRadius;

    bool mIsConsistentMapping;

    EntityPointVector mOriginEntityPointsVector;

    EntityPointVector mDestinationEntityPointsVector;

    MappingFilterFunction::Pointer mpMappingFilterFunction;

    SparseMatrixType mMappingMatrix;

    IndexType mBucketSize = 100;

    typename KDTree::Pointer mpSearchTree;

    ///@}
    ///@name Private operations
    ///@{

    void InitializeEntityPoints(
        EntityPointVector& rEntityPointsVector,
        const ModelPart& rModelPart) const;

    void CreateSearchTreeWithAllNodesInOriginModelPart();

    void ComputeWeightForAllNeighbors(
        double& rSumOfWeights,
        std::vector<double>& rListOfWeights,
        const IndexType NumberOfNeighbours,
        const EntityPointType& rEntityPoint,
        const EntityPointVector& rNeighbourEntityPoints) const;

    void FillMappingMatrixWithWeights(
        const EntityPointType& rEntityPoint,
        const EntityPointVector& rNeighbourEntityPoints,
        const std::vector<double>& rListOfWeights,
        const double WeightSum,
        const IndexType NumberOfNeighbours);

    void ComputeMappingMatrix(const EntityPointVector& rDestinationEntityPointsVector);

    ///@}
};

///@}
///@name Input and output
///@{

/// output stream function
template<class TContainerType>
inline std::ostream& operator<<(
    std::ostream& rOStream,
    const VertexMorphingContainerVariableDataMapper<TContainerType>& rThis)
{
    return rOStream << rThis.Info();
}

///@}
} // namespace Kratos