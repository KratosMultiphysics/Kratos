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

template<class TContainerType>
Expression::ConstPointer GetNodalDomainSizeExpression(
    const TContainerType& rContainer,
    const ModelPart::NodesContainerType& rNodes)
{
    const IndexType number_of_nodes = rNodes.size();

    auto p_expression = LiteralFlatExpression<double>::Create(number_of_nodes, {});

    // reset the existing values in the expression
    IndexPartition<IndexType>(number_of_nodes).for_each([&p_expression](const IndexType Index) {
        *(p_expression->begin() + Index) = 0.0;
    });

    IndexPartition<IndexType>(rContainer.size()).for_each([&rNodes, &rContainer, &p_expression](const IndexType Index) {
        const auto& r_geometry = (rContainer.begin() + Index)->GetGeometry();
        const double nodal_domain_size = r_geometry.DomainSize() / r_geometry.size();
        for (const auto& r_node : r_geometry) {
            const auto node_ptr = rNodes.find(r_node.Id());
            if (node_ptr != rNodes.end()) {
                const IndexType node_index = std::distance(rNodes.begin(), node_ptr);
                AtomicAdd<double>(*(p_expression->begin() + node_index), nodal_domain_size);
            } else {
                KRATOS_ERROR << "Node with id " << r_node.Id() << " not found.";
            }
        }
    });

    return p_expression;
}

template<class TEntityType>
double GetDomainSize(
    const EntityPoint<TEntityType>& rPoint,
    Expression const * const pExpression)
{
    return rPoint.GetEntity().GetGeometry().DomainSize();
}

template<>
double GetDomainSize(
    const EntityPoint<ModelPart::NodeType>& rPoint,
    Expression const * const pExpression)
{
    return pExpression->Evaluate(rPoint.Id(), rPoint.Id(), 0);
}

struct NonIntegratedWeight
{
    static double Compute(
        const double FilterFunctionWeight,
        const double DomainSize)
    {
        return FilterFunctionWeight;
    }
};

struct IntegratedWeight
{
    static double Compute(
        const double FilterFunctionWeight,
        const double DomainSize)
    {
        return FilterFunctionWeight * DomainSize;
    }
};

template <class TEntityType, class TWeightComputationType>
void ComputeWeightForAllNeighbors(
    double& rSumOfWeights,
    std::vector<double>& rListOfWeights,
    const FilterFunction& rFilterFunction,
    const double VertexMorphingRadius,
    const EntityPoint<TEntityType>& rDesignPoint,
    const std::vector<typename EntityPoint<TEntityType>::Pointer>& rNeighbourNodes,
    const IndexType NumberOfNeighbours,
    Expression const * const pExpression)
{
    for (IndexType i = 0; i < NumberOfNeighbours; ++i) {
        const double domain_size = GetDomainSize(*rNeighbourNodes[i], pExpression);
        const double filter_weight = rFilterFunction.ComputeWeight(rDesignPoint.Coordinates(), rNeighbourNodes[i]->Coordinates(), VertexMorphingRadius);

        rListOfWeights[i] = TWeightComputationType::Compute(filter_weight, domain_size);
        rSumOfWeights += filter_weight * domain_size;
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

    KRATOS_INFO("ExplicitVertexMorphingFilter") << "Creating search tree to perform mapping..." << std::endl;

    const auto& r_container = ExplicitVertexMorphingFilterHelperUtilities::GetContainer<TContainerType>(mrModelPart);

    if (mEntityPointVector.size() != r_container.size()) {
        mEntityPointVector.resize(r_container.size());
    }

    // now fill the points vector
    IndexPartition<IndexType>(r_container.size()).for_each([&](const IndexType Index) {
        mEntityPointVector[Index] = Kratos::make_shared<EntityPoint<EntityType>>(*(r_container.begin() + Index), Index);
    });

    mpSearchTree =  Kratos::make_shared<ExplicitVertexMorphingFilter::KDTree>(mEntityPointVector.begin(), mEntityPointVector.end(), mBucketSize);

    if constexpr(std::is_same_v<TContainerType, ModelPart::NodesContainerType>) {
        const auto& r_nodes = mrModelPart.Nodes();
        const auto& r_conditions = mrModelPart.Conditions();
        const IndexType number_of_conditions = r_conditions.size();

        const auto& r_elements = mrModelPart.Elements();
        const IndexType number_of_elements = r_elements.size();

        KRATOS_ERROR_IF(number_of_conditions > 0 && number_of_elements > 0) << mrModelPart.FullName()
            << " has both elements and conditions. Nodal mapping model parts should only have either conditions or elements.\n";

        if (number_of_conditions > 0) {
            mpNodalDomainSizeExpression = ExplicitVertexMorphingFilterHelperUtilities::GetNodalDomainSizeExpression(r_conditions, r_nodes);
        } else if (number_of_elements > 0) {
            mpNodalDomainSizeExpression = ExplicitVertexMorphingFilterHelperUtilities::GetNodalDomainSizeExpression(r_elements, r_nodes);
        } else {
            KRATOS_ERROR << "Nodal mapping requires atleast either conditions or elements to be present in "
                         << mrModelPart.FullName() << ".\n";
        }
    }

    KRATOS_INFO("ExplicitVertexMorphingFilter") << "Search tree created in: " << timer.ElapsedSeconds() << " s" << std::endl;

    KRATOS_CATCH("");
}

template<class TContainerType>
template<class TWeightIntegrationType>
ContainerExpression<TContainerType> ExplicitVertexMorphingFilter<TContainerType>::GenericFilterField(const ContainerExpression<TContainerType>& rContainerExpression) const
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
        ExplicitVertexMorphingFilterHelperUtilities::ComputeWeightForAllNeighbors<EntityType, TWeightIntegrationType>(
            sum_of_weights, list_of_weights, *mpKernelFunction, vertex_morphing_radius,
            entity_point, rTLS.mNeighbourEntityPoints, number_of_neighbors, this->mpNodalDomainSizeExpression.get());

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
ContainerExpression<TContainerType> ExplicitVertexMorphingFilter<TContainerType>::FilterField(const ContainerExpression<TContainerType>& rContainerExpression) const
{
    return GenericFilterField<ExplicitVertexMorphingFilterHelperUtilities::IntegratedWeight>(rContainerExpression);
}

template<class TContainerType>
ContainerExpression<TContainerType> ExplicitVertexMorphingFilter<TContainerType>::FilterIntegratedField(const ContainerExpression<TContainerType>& rContainerExpression) const
{
    return GenericFilterField<ExplicitVertexMorphingFilterHelperUtilities::NonIntegratedWeight>(rContainerExpression);
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
#define KRATOS_INSTANTIATE_EXPLICIT_VERTEX_MORPHING_FILTER_METHODS(CONTAINER_TYPE)                                                                                                                                                      \
    template class ExplicitVertexMorphingFilter<CONTAINER_TYPE>;                                                                                                                                                                        \
    template ContainerExpression<CONTAINER_TYPE> ExplicitVertexMorphingFilter<CONTAINER_TYPE>::GenericFilterField<ExplicitVertexMorphingFilterHelperUtilities::IntegratedWeight>(const ContainerExpression<CONTAINER_TYPE>&) const;     \
    template ContainerExpression<CONTAINER_TYPE> ExplicitVertexMorphingFilter<CONTAINER_TYPE>::GenericFilterField<ExplicitVertexMorphingFilterHelperUtilities::NonIntegratedWeight>(const ContainerExpression<CONTAINER_TYPE>&) const;

KRATOS_INSTANTIATE_EXPLICIT_VERTEX_MORPHING_FILTER_METHODS(ModelPart::NodesContainerType)
KRATOS_INSTANTIATE_EXPLICIT_VERTEX_MORPHING_FILTER_METHODS(ModelPart::ConditionsContainerType)
KRATOS_INSTANTIATE_EXPLICIT_VERTEX_MORPHING_FILTER_METHODS(ModelPart::ElementsContainerType)

#undef KRATOS_INSTANTIATE_EXPLICIT_VERTEX_MORPHING_FILTER_METHODS

} // namespace Kratos