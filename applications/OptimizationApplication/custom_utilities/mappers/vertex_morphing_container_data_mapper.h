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
#include <variant>
#include <functional>

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "spatial_containers/spatial_containers.h"

// Application includes
#include "custom_utilities/mappers/container_data_mapper.h"
#include "custom_utilities/mappers/entity_point.h"

namespace Kratos {

///@name Kratos Classes
///@{

class KRATOS_API(OPTIMIZATION_APPLICATION) VertexMorphingContainerDataMapper : public ContainerDataMapper
{
private:
    ///@name Private classes
    ///@{

    template<class TContainerType>
    class ContainerMapper
    {
    public:
        ///@name Type definitions
        ///@{

        using TEntityType = typename TContainerType::data_type;

        using EntityPointType = EntityPoint<TEntityType>;

        using EntityPointTypePointer = typename EntityPointType::Pointer;

        using EntityPointVector = std::vector<EntityPointTypePointer>;

        using EntityPointIterator = typename std::vector<EntityPointTypePointer>::iterator;

        using DoubleVectorIterator = std::vector<double>::iterator;

        using  BucketType = Bucket<3, EntityPointType, EntityPointVector, EntityPointTypePointer, EntityPointIterator, DoubleVectorIterator>;

        using KDTree = Tree<KDTreePartition<BucketType>>;

        /// Pointer definition of ContainerMapper
        KRATOS_CLASS_POINTER_DEFINITION(ContainerMapper);

        ///@}
        ///@name LifeCycle
        ///@{

        ContainerMapper(
            TContainerType& rContainer,
            const std::function<double(const array_1d<double, 3>&, const array_1d<double, 3>&)>& rFilterFunction,
            const double FilterRadius,
            const IndexType MaxEntitiesInFilterRadius,
            const IndexType BucketSize);

        ///@}
    private:
        ///@name Private member variables
        ///@{

        TContainerType& mrContainer;

        const std::function<double(const array_1d<double, 3>&, const array_1d<double, 3>&)>& mrFilterFunction;

        typename KDTree::UniquePointer mpKDTree = nullptr;

        EntityPointVector mEntityPointsVector;

        const double mFilterRadius;

        const IndexType mMaxEntitiesInFilterRadius;

        ///@}
        ///@name Private operations
        ///@{

        void InitializeEntityPoints();

        void ComputeMappingMatrix();

        void ComputeWeightForAllNeighbors(
            double& rSumOfWeights,
            std::vector<double>& rListOfWeights,
            const EntityPointType& rEntityPoint,
            const EntityPointVector& rNeighbourEntityPoints);

        ///@}
    };

    ///@}

public:
    ///@name Type definitions
    ///@{

    using IndexType = std::size_t;

    /// Pointer definition of ContainerData
    KRATOS_CLASS_POINTER_DEFINITION(VertexMorphingContainerDataMapper);

    ///@}
    ///@name Life cycle
    ///@{

    VertexMorphingContainerDataMapper(
        ModelPart& rOriginModelPart,
        ModelPart& rDestinationModelPart,
        const ContainerData::ContainerDataType& rContainerDataType);

    ~VertexMorphingContainerDataMapper() override = default;

    ///@}
    ///@name Public operations

    void Update() override;

    void Map(
        ContainerData& rOutputDataContainer,
        const ContainerData& rInputDataContainer) override;

    void InverseMap(
        ContainerData& rOutputDataContainer,
        const ContainerData& rInputDataContainer) override;

    ///@}
private:
    ///@name Private member variables
    ///@{

    ModelPart& mrOriginModelPart;
    ModelPart& mrDestinationModelPart;
    const ContainerData::ContainerDataType& mrContainerDataType;
    const IndexType mBucketSize = 100;

    std::variant<
        ContainerMapper<ModelPart::NodesContainerType>::UniquePointer,
        ContainerMapper<ModelPart::ConditionsContainerType>::UniquePointer,
        ContainerMapper<ModelPart::ElementsContainerType>::UniquePointer> mpContainerMapper;

    ///@}
    ///@name Private operations
    ///@{

    ///@}
};

///@}
} // namespace Kratos