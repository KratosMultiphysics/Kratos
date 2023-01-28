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

// System includes

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"

// Application includes
#include "custom_utilities/mappers/container_data_mapper.h"
#include "utilities/parallel_utilities.h"

// Include base h
#include "vertex_morphing_container_data_mapper.h"

namespace Kratos {

template<class TContainerType>
void VertexMorphingContainerDataMapper::ContainerMapper<TContainerType>::InitializeEntityPoints()
{
    if (mEntityPointsVector.size() != mrContainer.size()) {
        mEntityPointsVector.resize(mrContainer.size());
    }

    IndexPartition<IndexType>(mEntityPointsVector.size()).for_each([&](const IndexType Index) {
        mEntityPointsVector[Index] = Kratos::make_shared<EntityPointType>(*(mrContainer.begin() + Index), Index);
    });
}

template<class TContainerType>
void VertexMorphingContainerDataMapper::ContainerMapper<TContainerType>::ComputeWeightForAllNeighbors(
    double& rSumOfWeights,
    std::vector<double>& rListOfWeights,
    const EntityPointType& rEntityPoint,
    const EntityPointVector& rNeighbourEntityPoints)
{
    IndexType index = 0;
    for(const auto& r_neighbour_entity : rNeighbourEntityPoints)
    {
        const double weight = mrFilterFunction(rEntityPoint.Coordinates(), r_neighbour_entity->Coordinates());

        rListOfWeights[index++] = weight;
        rSumOfWeights += weight;
    }
}

template<class TContainerType>
void VertexMorphingContainerDataMapper::ContainerMapper<TContainerType>::ComputeMappingMatrix()
{
    for(const auto& r_entity : mEntityPointsVector)
    {
        EntityPointVector neighbor_entity_points(mMaxEntitiesInFilterRadius);
        std::vector<double> resulting_squared_distances(mMaxEntitiesInFilterRadius);
        unsigned int number_of_neighbors = mpKDTree->SearchInRadius(*r_entity,
                                                                    mFilterRadius,
                                                                    neighbor_entity_points.begin(),
                                                                    resulting_squared_distances.begin(),
                                                                    mMaxEntitiesInFilterRadius);



        std::vector<double> list_of_weights( number_of_neighbors, 0.0 );
        double sum_of_weights = 0.0;

        if(number_of_neighbors >= mMaxEntitiesInFilterRadius) {
            KRATOS_WARNING("OptApp::VertexMorphingContainerDataMapper")
                << "For entity " << r_entity->GetEntity().Id() << " and specified filter radius, maximum number of neighbor nodes (="
                << mMaxEntitiesInFilterRadius << " entities) reached!" << std::endl;
        }

        ComputeWeightForAllNeighbors(sum_of_weights, list_of_weights, *r_entity, neighbor_entity_points);
        // FillMappingMatrixWithWeights( node_i, neighbor_nodes, number_of_neighbors, list_of_weights, sum_of_weights );
    }
}

template<class TContainerType>
VertexMorphingContainerDataMapper::ContainerMapper<TContainerType>::ContainerMapper(
    TContainerType& rContainer,
    const std::function<double(const array_1d<double, 3>&, const array_1d<double, 3>&)>& rFilterFunction,
    const double FilterRadius,
    const IndexType MaxEntitiesInFilterRadius,
    const IndexType BucketSize)
    : mrContainer(rContainer),
      mrFilterFunction(rFilterFunction),
      mFilterRadius(FilterRadius),
      mMaxEntitiesInFilterRadius(MaxEntitiesInFilterRadius)
{
    InitializeEntityPoints();
    mpKDTree = Kratos::make_unique<KDTree>(mEntityPointsVector.begin(), mEntityPointsVector.end(), BucketSize);

    ComputeMappingMatrix();
}

VertexMorphingContainerDataMapper::VertexMorphingContainerDataMapper(
    ModelPart& rOriginModelPart,
    ModelPart& rDestinationModelPart,
    const ContainerData::ContainerDataType& rContainerDataType)
    : ContainerDataMapper(),
      mrOriginModelPart(rOriginModelPart),
      mrDestinationModelPart(rDestinationModelPart),
      mrContainerDataType(rContainerDataType)
{
    switch (mrContainerDataType)
    {
        case ContainerData::ContainerDataType::NodalHistorical:
        case ContainerData::ContainerDataType::NodalNonHistorical:
            mpContainerMapper = Kratos::make_unique<ContainerMapper<ModelPart::NodesContainerType>>(rOriginModelPart.Nodes(), [](const array_1d<double, 3>&, const array_1d<double, 3>&) {return 0.0;}, 10, 100, mBucketSize);
            break;
        case ContainerData::ContainerDataType::ConditionNonHistorical:
        case ContainerData::ContainerDataType::ConditionProperties:
            mpContainerMapper = Kratos::make_unique<ContainerMapper<ModelPart::ConditionsContainerType>>(rOriginModelPart.Conditions(), [](const array_1d<double, 3>&, const array_1d<double, 3>&) {return 0.0;}, 10, 100, mBucketSize);
            break;
        case ContainerData::ContainerDataType::ElementNonHistorical:
        case ContainerData::ContainerDataType::ElementProperties:
            mpContainerMapper = Kratos::make_unique<ContainerMapper<ModelPart::ElementsContainerType>>(rOriginModelPart.Elements(), [](const array_1d<double, 3>&, const array_1d<double, 3>&) {return 0.0;}, 10, 100, mBucketSize);
            break;
    }
}


void VertexMorphingContainerDataMapper::Update()
{

}

void VertexMorphingContainerDataMapper::Map(
    ContainerData& rOutputDataContainer,
    const ContainerData& rInputDataContainer)
{

}

void VertexMorphingContainerDataMapper::InverseMap(
    ContainerData& rOutputDataContainer,
    const ContainerData& rInputDataContainer)
{
}

} // namespace Kratos