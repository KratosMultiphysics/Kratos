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
#include <type_traits>
#include <tuple>

// Project includes
#include "expression/literal_flat_expression.h"
#include "utilities/builtin_timer.h"
#include "utilities/parallel_utilities.h"
#include "utilities/atomic_utilities.h"

// Application includes

// Include base h
#include "explicit_vertex_morphing_filter.h"

namespace Kratos {

namespace ExplicitVertexMorphingFilterHelperUtilities
{

template<class TContainerType>
const TContainerType& GetContainer(const ModelPart& rModelPart);

template<>
const ModelPart::NodesContainerType& GetContainer(const ModelPart& rModelPart)
{
    return rModelPart.Nodes();
}

template<>
const ModelPart::ConditionsContainerType& GetContainer(const ModelPart& rModelPart)
{
    return rModelPart.Conditions();
}

template<>
const ModelPart::ElementsContainerType& GetContainer(const ModelPart& rModelPart)
{
    return rModelPart.Elements();
}

template <class TEntityType>
void ComputeWeightForAllNeighbors(
    double& rSumOfWeights,
    std::vector<double>& rListOfWeights,
    const FilterFunction& rFilterFunction,
    const double VertexMorphingRadius,
    const EntityPoint<TEntityType>& rDesignPoint,
    const std::vector<typename EntityPoint<TEntityType>::Pointer>& rNeighbourNodes,
    const IndexType NumberOfNeighbours)
{
    for (IndexType i = 0; i < NumberOfNeighbours; ++i) {
        const double weight = rFilterFunction.ComputeWeight(rDesignPoint.Coordinates(), rNeighbourNodes[i]->Coordinates(), VertexMorphingRadius);
        rListOfWeights[i] = weight;
        rSumOfWeights += weight;
    }
}

}; // namespace ExplicitVertexMorphingFilterHelperUtilities

template<class TContainerType>
ExplicitVertexMorphingFilter<TContainerType>::ExplicitVertexMorphingFilter(
    const ModelPart& rModelPart,
    const std::string& rKernelFunctionType,
    const IndexType MaxNumberOfNeighbours)
    : mrModelPart(rModelPart),
      mMaxNumberOfNeighbors(MaxNumberOfNeighbours)
{
    mpKernelFunction = Kratos::make_unique<FilterFunction>(rKernelFunctionType);
    Update();
}

template<class TContainerType>
void ExplicitVertexMorphingFilter<TContainerType>::SetFilterRadius(const ContainerExpression<TContainerType>& rContainerExpression)
{
    KRATOS_ERROR_IF_NOT(rContainerExpression.GetItemComponentCount() == 1)
        << "Only scalar values are allowed for the filter radius container expression. "
        << "Provided contaienr expression = " << rContainerExpression << ".\n";

    KRATOS_ERROR_IF_NOT(&rContainerExpression.GetModelPart() == &mrModelPart)
        << "Filter radius container expression model part and filter model part mismatch."
        << "\n\tFilter = " << *this
        << "\n\tContainerExpression = " << rContainerExpression;

    mpFilterRadiusContainer = rContainerExpression.Clone();
}

template<class TContainerType>
void ExplicitVertexMorphingFilter<TContainerType>::ExplicitVertexMorphingFilter::Update()
{
    KRATOS_TRY

    BuiltinTimer timer;

    KRATOS_INFO("OptimizationApplication") << "Creating search tree to perform mapping..." << std::endl;

    const auto& r_container = ExplicitVertexMorphingFilterHelperUtilities::GetContainer<TContainerType>(mrModelPart);

    if (mEntityPointVector.size() != r_container.size()) {
        mEntityPointVector.resize(r_container.size());
    }

    // now fill the points vector
    IndexPartition<IndexType>(r_container.size()).for_each([&](const IndexType Index) {
        mEntityPointVector[Index] = Kratos::make_shared<EntityPoint<EntityType>>(*(r_container.begin() + Index), Index);
    });

    mpSearchTree =  Kratos::make_shared<ExplicitVertexMorphingFilter::KDTree>(mEntityPointVector.begin(), mEntityPointVector.end(), mBucketSize);

    KRATOS_INFO("OptimizationApplication") << "Search tree created in: " << timer.ElapsedSeconds() << " s" << std::endl;

    KRATOS_CATCH("");
}

template<class TContainerType>
ContainerExpression<TContainerType> ExplicitVertexMorphingFilter<TContainerType>::FilterField(const ContainerExpression<TContainerType>& rContainerExpression) const
{
    KRATOS_TRY

    KRATOS_ERROR_IF(mpFilterRadiusContainer.get() == nullptr)
        << "The filter radius container expression not set. "
        << "Please set it using SetFilterRadius method.\n\t Filter = "
        << *this;

    KRATOS_ERROR_IF_NOT(rContainerExpression.HasExpression())
        << "Uninitialized container expression given. "
        << rContainerExpression;

    KRATOS_ERROR_IF_NOT(&rContainerExpression.GetModelPart() == &mrModelPart)
        << "Filter radius container expression model part and filter model part mismatch."
        << "\n\tFilter = " << *this
        << "\n\tContainerExpression = " << rContainerExpression;

    const IndexType stride = rContainerExpression.GetItemComponentCount();
    const auto& r_origin_expression = rContainerExpression.GetExpression();
    const auto& r_container = rContainerExpression.GetContainer();
    const auto& r_filter_radius_expression = mpFilterRadiusContainer->GetExpression();

    ContainerExpression<TContainerType> result(*rContainerExpression.pGetModelPart());
    auto p_expression = LiteralFlatExpression<double>::Create(result.GetContainer().size(), rContainerExpression.GetItemShape());
    result.SetExpression(p_expression);

    struct TLS
    {
        explicit TLS(const IndexType MaxNumberOfNeighbors)
        {
            mNeighbourEntityPoints.resize(MaxNumberOfNeighbors);
            mResultingSquaredDistances.resize(MaxNumberOfNeighbors);
        }

        EntityPointVector mNeighbourEntityPoints;
        std::vector<double> mResultingSquaredDistances;
    };

    IndexPartition<IndexType>(r_container.size()).for_each(TLS(mMaxNumberOfNeighbors), [&](const IndexType Index, TLS& rTLS){
        const double vertex_morphing_radius = r_filter_radius_expression.Evaluate(Index, Index, 0);

        EntityPoint<EntityType> entity_point(*(r_container.begin() + Index), Index);
        const auto number_of_neighbors = mpSearchTree->SearchInRadius(
                                            entity_point,
                                            vertex_morphing_radius,
                                            rTLS.mNeighbourEntityPoints.begin(),
                                            rTLS.mResultingSquaredDistances.begin(),
                                            mMaxNumberOfNeighbors);

        std::vector<double> list_of_weights(number_of_neighbors, 0.0);
        double sum_of_weights = 0.0;
        ExplicitVertexMorphingFilterHelperUtilities::ComputeWeightForAllNeighbors(
            sum_of_weights, list_of_weights, *mpKernelFunction, vertex_morphing_radius,
            entity_point, rTLS.mNeighbourEntityPoints, number_of_neighbors);

        const IndexType current_data_begin = Index * stride;

        for (IndexType j = 0; j < stride; ++j) {
            double& current_index_value = *(p_expression->begin() + current_data_begin + j);
            current_index_value = 0.0;

            for(IndexType neighbour_index = 0 ; neighbour_index < number_of_neighbors; ++neighbour_index) {
                const IndexType neighbour_id = rTLS.mNeighbourEntityPoints[neighbour_index]->Id();
                const double weight = list_of_weights[neighbour_index] / sum_of_weights;
                const double origin_value = r_origin_expression.Evaluate(neighbour_id, neighbour_id * stride, j);
                current_index_value +=  weight * origin_value;
            }
        }
    });

    return result;

    KRATOS_CATCH("");
}

template<class TContainerType>
ContainerExpression<TContainerType> ExplicitVertexMorphingFilter<TContainerType>::UnFilterField(const ContainerExpression<TContainerType>& rContainerExpression) const
{
    KRATOS_TRY

    KRATOS_ERROR_IF(mpFilterRadiusContainer.get() == nullptr)
        << "The filter radius container expression not set. "
        << "Please set it using SetFilterRadius method.\n\t Filter = "
        << *this;

    KRATOS_ERROR_IF_NOT(rContainerExpression.HasExpression())
        << "Uninitialized container expression given. "
        << rContainerExpression;

    KRATOS_ERROR_IF_NOT(&rContainerExpression.GetModelPart() == &mrModelPart)
        << "Filter radius container expression model part and filter model part mismatch."
        << "\n\tFilter = " << *this
        << "\n\tContainerExpression = " << rContainerExpression;

    const auto& r_origin_expression = rContainerExpression.GetExpression();
    const auto& r_container = rContainerExpression.GetContainer();
    const IndexType stride = rContainerExpression.GetItemComponentCount();
    const auto& r_filter_radius_expression = mpFilterRadiusContainer->GetExpression();

    EntityPointVector neighbour_nodes(mMaxNumberOfNeighbors);
    std::vector<double> resulting_squared_distances(mMaxNumberOfNeighbors);

    ContainerExpression<TContainerType> result(*rContainerExpression.pGetModelPart());
    auto p_expression = LiteralFlatExpression<double>::Create(result.GetContainer().size(), rContainerExpression.GetItemShape());
    result.SetExpression(p_expression);

    IndexPartition<IndexType>(r_container.size() * stride).for_each([&](const IndexType Index) {
        *(p_expression->begin() + Index) = 0.0;
    });

    struct TLS
    {
        explicit TLS(const IndexType MaxNumberOfNeighbors)
        {
            mNeighbourEntityPoints.resize(MaxNumberOfNeighbors);
            mResultingSquaredDistances.resize(MaxNumberOfNeighbors);
        }

        EntityPointVector mNeighbourEntityPoints;
        std::vector<double> mResultingSquaredDistances;
    };

    IndexPartition<IndexType>(r_container.size()).for_each(TLS(mMaxNumberOfNeighbors), [&](const IndexType Index, TLS& rTLS){
        const double vertex_morphing_radius = r_filter_radius_expression.Evaluate(Index, Index, 0);

        EntityPoint<EntityType> entity_point(*(r_container.begin() + Index), Index);
        const auto number_of_neighbors = mpSearchTree->SearchInRadius(
                                            entity_point,
                                            vertex_morphing_radius,
                                            rTLS.mNeighbourEntityPoints.begin(),
                                            rTLS.mResultingSquaredDistances.begin(),
                                            mMaxNumberOfNeighbors);

        std::vector<double> list_of_weights(number_of_neighbors, 0.0);
        double sum_of_weights = 0.0;
        ExplicitVertexMorphingFilterHelperUtilities::ComputeWeightForAllNeighbors(
            sum_of_weights, list_of_weights, *mpKernelFunction, vertex_morphing_radius,
            entity_point, rTLS.mNeighbourEntityPoints, number_of_neighbors);

        for (IndexType j = 0; j < stride; ++j) {
            const double origin_value = r_origin_expression.Evaluate(Index, Index * stride, j);

            for(IndexType neighbour_index = 0 ; neighbour_index < number_of_neighbors; ++neighbour_index) {
                const IndexType neighbour_id = rTLS.mNeighbourEntityPoints[neighbour_index]->Id();
                double weight = list_of_weights[neighbour_index] / sum_of_weights;
                AtomicAdd<double>(*(p_expression->begin() + neighbour_id * stride + j), weight * origin_value);
            }
        }
    });

    return result;

    KRATOS_CATCH("");
}

template<class TContainerType>
std::string ExplicitVertexMorphingFilter<TContainerType>::Info() const
{
    std::stringstream msg;

    msg << "ExplicitVertexMorphingFilter: ";

    if constexpr(std::is_same_v<TContainerType, ModelPart::NodesContainerType>) {
        msg << "Nodal";
    } else if constexpr(std::is_same_v<TContainerType, ModelPart::ConditionsContainerType>) {
        msg << "Condition";
    } else {
        msg << "Element";
    }

    msg << "Container in " <<  mrModelPart.FullName();
    return msg.str();
}

// template instantiations
template class ExplicitVertexMorphingFilter<ModelPart::NodesContainerType>;
template class ExplicitVertexMorphingFilter<ModelPart::ConditionsContainerType>;
template class ExplicitVertexMorphingFilter<ModelPart::ElementsContainerType>;

} // namespace Kratos