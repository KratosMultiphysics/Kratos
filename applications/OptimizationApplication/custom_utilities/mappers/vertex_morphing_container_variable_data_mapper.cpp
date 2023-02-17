//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   license: OptimizationApplication/license.txt
//
//  Main authors:    Baumgaertner Daniel
//                   Aditya Ghantasala
//                   Suneth Warnakulasuriya
//

// System includes
#include <sstream>

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "utilities/parallel_utilities.h"

// Application includes
#include "custom_utilities/container_variable_data_holder_utils.h"

// Include base h
#include "vertex_morphing_container_variable_data_mapper.h"

namespace Kratos {

namespace VertexMorphingContainerVariableDataMapperHelpers
{
template<class TContainerType>
const TContainerType& GetContainer(const ModelPart& rModelPart);

template <>
const ModelPart::NodesContainerType& GetContainer<ModelPart::NodesContainerType>(const ModelPart& rModelPart)
{
    return rModelPart.GetCommunicator().LocalMesh().Nodes();
}

template <>
const ModelPart::ConditionsContainerType& GetContainer<ModelPart::ConditionsContainerType>(const ModelPart& rModelPart)
{
    return rModelPart.GetCommunicator().LocalMesh().Conditions();
}

template <>
const ModelPart::ElementsContainerType& GetContainer<ModelPart::ElementsContainerType>(const ModelPart& rModelPart)
{
    return rModelPart.GetCommunicator().LocalMesh().Elements();
}
}

template<class TContainerType>
VertexMorphingContainerVariableDataMapper<TContainerType>::VertexMorphingContainerVariableDataMapper(
    const ModelPart& rOriginModelPart,
    const ModelPart& rDestinationModelPart,
    Parameters Params)
    : mrOriginModelPart(rOriginModelPart),
      mrDestinationModelPart(rDestinationModelPart)
{
    KRATOS_TRY

    Parameters default_parameters = Parameters(R"(
    {
        "filter_function_type"         : "linear",
        "filter_radius"                : 1.0,
        "max_entities_in_filter_radius": 10000
    })" );

    Params.ValidateAndAssignDefaults(default_parameters);

    mFilterFunctionType = Params["filter_function_type"].GetString();
    mpMappingFilterFunction = Kratos::make_shared<MappingFilterFunction>(mFilterFunctionType);
    mFilterRadius = Params["filter_radius"].GetDouble();
    mMaxEntitiesInFilterRadius = Params["max_entities_in_filter_radius"].GetInt();
    mIsConsistentMapping = (mrOriginModelPart.FullName() == mrDestinationModelPart.FullName());

    KRATOS_CATCH("");
}

template<class TContainerType>
void VertexMorphingContainerVariableDataMapper<TContainerType>::InitializeEntityPoints(
    EntityPointVector& rEntityPointsVector,
    const ModelPart& rModelPart) const
{
    auto& r_container = VertexMorphingContainerVariableDataMapperHelpers::GetContainer<TContainerType>(rModelPart);

    if (rEntityPointsVector.size() != r_container.size()) {
        rEntityPointsVector.resize(r_container.size());
    }

    IndexPartition<IndexType>(rEntityPointsVector.size()).for_each([&](const IndexType Index) {
        rEntityPointsVector[Index] = Kratos::make_shared<EntityPointType>(*(r_container.begin() + Index), Index);
    });
}

template<class TContainerType>
void VertexMorphingContainerVariableDataMapper<TContainerType>::CreateSearchTreeWithAllNodesInOriginModelPart()
{
    mpSearchTree = Kratos::shared_ptr<KDTree>(new KDTree(mOriginEntityPointsVector.begin(), mOriginEntityPointsVector.end(), mBucketSize));
}

template<class TContainerType>
void VertexMorphingContainerVariableDataMapper<TContainerType>::ComputeWeightForAllNeighbors(
    double& rWieghtsSum,
    std::vector<double>& rListOfWeights,
    const IndexType NumberOfNeighbours,
    const EntityPointType& rEntityPoint,
    const EntityPointVector& rNeighbourEntityPoints) const
{
    for (IndexType neighbour_index = 0; neighbour_index < NumberOfNeighbours; ++neighbour_index) {
        const auto& r_neighbour_entity = *(*(rNeighbourEntityPoints.begin() + neighbour_index));
        const double weight = mpMappingFilterFunction->ComputeWeight(rEntityPoint.Coordinates(), r_neighbour_entity.Coordinates(), mFilterRadius);

        rListOfWeights[neighbour_index] = weight;
        rWieghtsSum += weight;
    }
}

template<class TContainerType>
void VertexMorphingContainerVariableDataMapper<TContainerType>::FillMappingMatrixWithWeights(
    const EntityPointType& rEntityPoint,
    const EntityPointVector& rNeighbourEntityPoints,
    const std::vector<double>& rListOfWeights,
    const double WeightSum,
    const IndexType NumberOfNeighbours)
{
    const IndexType row_id = rEntityPoint.Id();
    for (IndexType neighbour_index = 0; neighbour_index < NumberOfNeighbours; ++neighbour_index) {
        const auto& r_neighbour_entity = *(rNeighbourEntityPoints.begin() + neighbour_index);
        const IndexType column_id = r_neighbour_entity->Id();

        double weight = rListOfWeights[neighbour_index] / WeightSum;
        mMappingMatrix.insert_element(row_id, column_id, weight);
    }
}

template<class TContainerType>
void VertexMorphingContainerVariableDataMapper<TContainerType>::ComputeMappingMatrix(
    const EntityPointVector& rDestinationEntityPointsVector)
{
    for(const auto& r_entity : rDestinationEntityPointsVector) {
        EntityPointVector neighbor_entity_points(mMaxEntitiesInFilterRadius);
        std::vector<double> resulting_squared_distances(mMaxEntitiesInFilterRadius);
        unsigned int number_of_neighbors = mpSearchTree->SearchInRadius(*r_entity,
                                                                        mFilterRadius,
                                                                        neighbor_entity_points.begin(),
                                                                        resulting_squared_distances.begin(),
                                                                        mMaxEntitiesInFilterRadius);



        std::vector<double> list_of_weights( number_of_neighbors, 0.0 );
        double sum_of_weights = 0.0;

        if(number_of_neighbors >= mMaxEntitiesInFilterRadius) {
            KRATOS_WARNING("OptApp::VertexMorphingContainerDataMapper")
                << "Reached maximum number of neighbor nodes for radius (="
                << mMaxEntitiesInFilterRadius << "!" << std::endl;
        }

        ComputeWeightForAllNeighbors(sum_of_weights, list_of_weights, number_of_neighbors, *r_entity, neighbor_entity_points);
        FillMappingMatrixWithWeights(*r_entity, neighbor_entity_points, list_of_weights, sum_of_weights, number_of_neighbors);
    }
}

template<class TContainerType>
void VertexMorphingContainerVariableDataMapper<TContainerType>::Update()
{
    KRATOS_TRY

    if (!mIsConsistentMapping) {
        // if we are mapping from/to different model parts
        InitializeEntityPoints(mOriginEntityPointsVector, mrOriginModelPart);
        InitializeEntityPoints(mDestinationEntityPointsVector, mrDestinationModelPart);
        mMappingMatrix.resize(mDestinationEntityPointsVector.size(), mOriginEntityPointsVector.size(), false);

        // search tree is always computed on the origin model part
        CreateSearchTreeWithAllNodesInOriginModelPart();
        ComputeMappingMatrix(mDestinationEntityPointsVector);
    } else {
        // if we are mapping from/to the same model part
        InitializeEntityPoints(mOriginEntityPointsVector, mrOriginModelPart);
        mMappingMatrix.resize(mOriginEntityPointsVector.size(), mOriginEntityPointsVector.size(), false);

        CreateSearchTreeWithAllNodesInOriginModelPart();
        ComputeMappingMatrix(mOriginEntityPointsVector);
    }

    KRATOS_CATCH("");
}

template<class TContainerType>
void VertexMorphingContainerVariableDataMapper<TContainerType>::Map(
    const ContainerVariableDataHolderType& rOriginDataContainer,
    ContainerVariableDataHolderType& rDestinationDataContainer) const
{
    KRATOS_TRY

    KRATOS_ERROR_IF_NOT(rOriginDataContainer.GetModelPart() == mrOriginModelPart)
        << "Origin model part ["
        << mrOriginModelPart.FullName() << "] of the mapper and model part of the origin data container are mismatching. "
        << "Followings are the provided containers: "
        << "\n\tDestination container = " << rDestinationDataContainer
        << "\n\tOrigin container      = " << rOriginDataContainer << "\n";

    KRATOS_ERROR_IF_NOT(rDestinationDataContainer.GetModelPart() == mrDestinationModelPart)
        << "Destination model part ["
        << mrDestinationModelPart.FullName() << "] of the mapper and model part of the destination data container are mismatching. "
        << "Followings are the provided containers: "
        << "\n\tDestination container = " << rDestinationDataContainer
        << "\n\tOrigin container      = " << rDestinationDataContainer << "\n";

    ContainerVariableDataHolderUtils::ProductWithEntityMatrix(rDestinationDataContainer, mMappingMatrix, rOriginDataContainer);

    KRATOS_CATCH("");
}

template<class TContainerType>
void VertexMorphingContainerVariableDataMapper<TContainerType>::InverseMap(
    ContainerVariableDataHolderType& rOriginDataContainer,
    const ContainerVariableDataHolderType& rDestinationDataContainer) const
{
    KRATOS_TRY

    KRATOS_ERROR_IF_NOT(rOriginDataContainer.GetModelPart() == mrOriginModelPart)
        << "Origin model part ["
        << mrOriginModelPart.FullName() << "] of the mapper and model part of the origin data container are mismatching. "
        << "Followings are the provided containers: "
        << "\n\tDestination container = " << rDestinationDataContainer
        << "\n\tOrigin container      = " << rOriginDataContainer << "\n";

    KRATOS_ERROR_IF_NOT(rDestinationDataContainer.GetModelPart() == mrDestinationModelPart)
        << "Destination model part ["
        << mrDestinationModelPart.FullName() << "] of the mapper and model part of the destination data container are mismatching. "
        << "Followings are the provided containers: "
        << "\n\tDestination container = " << rDestinationDataContainer
        << "\n\tOrigin container      = " << rDestinationDataContainer << "\n";

    if(mIsConsistentMapping) {
        ContainerVariableDataHolderUtils::ProductWithEntityMatrix(rOriginDataContainer, mMappingMatrix, rDestinationDataContainer);
    }
    else {
        SparseMatrixType transpose;
        ContainerVariableDataHolderUtils::Transpose(transpose, mMappingMatrix);
        ContainerVariableDataHolderUtils::ProductWithEntityMatrix(rOriginDataContainer, transpose, rDestinationDataContainer);
    }

    KRATOS_CATCH("");
}

template<class TContainerType>
std::string VertexMorphingContainerVariableDataMapper<TContainerType>::Info() const
{
    std::stringstream msg;

    msg << "VertexMorphing";

    if constexpr(std::is_same_v<TContainerType, ModelPart::NodesContainerType>) {
        msg << "Nodal";
    } else if constexpr(std::is_same_v<TContainerType, ModelPart::ConditionsContainerType>) {
        msg << "Condition";
    } else if constexpr(std::is_same_v<TContainerType, ModelPart::ElementsContainerType>) {
        msg << "Element";
    }

    msg << "ContainerVariableDataMapper [ ";
    msg << "Origin model part = " << mrOriginModelPart.FullName() << ", ";
    msg << "Destination model part = " << mrDestinationModelPart.FullName() << ", ";
    msg << "Is consistent mapping = " << (mIsConsistentMapping ? "yes" : "no") << ", ";
    msg << "Filter function type = " << mFilterFunctionType << ", ";
    msg << "Filter radius = " << mFilterRadius << ", ";
    msg << "Max num. neighbours = " << mMaxEntitiesInFilterRadius << " ]" << std::endl;

    return msg.str();
}

// template instantiations
template class VertexMorphingContainerVariableDataMapper<ModelPart::NodesContainerType>;
template class VertexMorphingContainerVariableDataMapper<ModelPart::ElementsContainerType>;
template class VertexMorphingContainerVariableDataMapper<ModelPart::ConditionsContainerType>;

} // namespace Kratos