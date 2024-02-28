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
#include <unordered_map>

// Project includes
#include "expression/literal_expression.h"
#include "expression/literal_flat_expression.h"
#include "utilities/atomic_utilities.h"
#include "utilities/builtin_timer.h"
#include "utilities/parallel_utilities.h"

// Application includes

// Include base h
#include "explicit_filter.h"

namespace Kratos {

namespace ExplicitFilterHelperUtilities
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
    const double Radius,
    const EntityPoint<TEntityType>& rDesignPoint,
    const std::vector<typename EntityPoint<TEntityType>::Pointer>& rNeighbourNodes,
    const IndexType NumberOfNeighbours,
    Expression const * const pExpression)
{
    for (IndexType i = 0; i < NumberOfNeighbours; ++i) {
        const double domain_size = GetDomainSize(*rNeighbourNodes[i], pExpression);
        const double filter_weight = rFilterFunction.ComputeWeight(rDesignPoint.Coordinates(), rNeighbourNodes[i]->Coordinates(), Radius);
        rListOfWeights[i] = TWeightComputationType::Compute(filter_weight, domain_size);
        rSumOfWeights += filter_weight * domain_size;
    }
}

template<class TEntityType>
struct TLS
{
    explicit TLS(const IndexType MaxNumberOfNeighbors)
    {
        mNeighbourEntityPoints.resize(MaxNumberOfNeighbors);
        mResultingSquaredDistances.resize(MaxNumberOfNeighbors);
    }

    std::vector<TEntityType> mNeighbourEntityPoints;
    std::vector<double> mResultingSquaredDistances;
};

}; // namespace ExplicitFilterHelperUtilities

template<class TContainerType>
ExplicitFilter<TContainerType>::ExplicitFilter(
    ModelPart& rModelPart,
    Parameters Settings)
    : mrModelPart(rModelPart)
{
    Parameters default_parameters = Parameters(R"(
    {
        "filter_type"               : "explicit_vertex_morphing",
        "filter_function_type"      : "linear",
        "radius"                    : 0.000000000001,
        "max_nodes_in_filter_radius": 1000,
        "damping_function_type"     : "sigmoidal",
        "damped_sub_model_parts"    : {}
    })" );
    Settings.ValidateAndAssignDefaults(default_parameters);

    mpKernelFunction = Kratos::make_unique<FilterFunction>(Settings["filter_function_type"].GetString());
    mMaxNumberOfNeighbors = Settings["max_nodes_in_filter_radius"].GetInt();

    mpFilterRadiusContainer = Kratos::make_shared<ContainerExpression<TContainerType>>(rModelPart);
    auto p_expression = LiteralExpression<double>::Create(Settings["radius"].GetDouble(), mpFilterRadiusContainer->GetContainer().size());
    mpFilterRadiusContainer->SetExpression(p_expression);

    mpDampingFunction = Kratos::make_unique<FilterFunction>(Settings["damping_function_type"].GetString());
    const auto damped_model_part_settings = Settings["damped_sub_model_parts"];
    for (auto it = damped_model_part_settings.begin(); it != damped_model_part_settings.end(); ++it) {
        std::vector<bool> damped_components;
        for (const auto& r_value : *it) {
            damped_components.push_back(r_value.GetBool());
        }
        mDampingComponentIndices.push_back(damped_components);
        mDampingModelParts.push_back(rModelPart.pGetSubModelPart(it.name()));
    }

    mNumberOfComponents = (mDampingComponentIndices.size() > 0 ? mDampingComponentIndices.front().size() : 0);
    const bool all_damping_valid = std::all_of(mDampingComponentIndices.begin(), mDampingComponentIndices.end(), [&](const auto& rValue) { return rValue.size() == mNumberOfComponents; });
    KRATOS_ERROR_IF_NOT(all_damping_valid) << "Mismatching damping components found.";

    Update();
}

template<class TContainerType>
void ExplicitFilter<TContainerType>::SetFilterRadius(const ContainerExpression<TContainerType>& rContainerExpression)
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
void ExplicitFilter<TContainerType>::ExplicitFilter::Update()
{
    KRATOS_TRY

    BuiltinTimer timer;

    KRATOS_INFO("ExplicitFilter") << "Creating search tree to perform mapping..." << std::endl;

    const auto& r_container = ExplicitFilterHelperUtilities::GetContainer<TContainerType>(mrModelPart);

    // now fill the points vector
    if (mEntityPointVector.size() != r_container.size()) {
        mEntityPointVector.resize(r_container.size());
    }

    std::unordered_map<IndexType, typename EntityPoint<EntityType>::Pointer> id_entity_map;

    for (IndexType i = 0; i < r_container.size(); ++i) {
        const auto& r_entity = *(r_container.begin() + i);
        auto entity_point = Kratos::make_shared<EntityPoint<EntityType>>(r_entity, i);
        mEntityPointVector[i] = entity_point;
        id_entity_map[r_entity.Id()] = entity_point;
        entity_point->SetDampedComponents(std::vector<bool>(mNumberOfComponents, false));
    };

    mpSearchTree =  Kratos::make_shared<ExplicitFilter::KDTree>(mEntityPointVector.begin(), mEntityPointVector.end(), mBucketSize);

    if constexpr(std::is_same_v<TContainerType, ModelPart::NodesContainerType>) {
        const auto& r_nodes = mrModelPart.Nodes();
        const auto& r_conditions = mrModelPart.Conditions();
        const IndexType number_of_conditions = r_conditions.size();

        const auto& r_elements = mrModelPart.Elements();
        const IndexType number_of_elements = r_elements.size();

        if (number_of_elements > 0) {
            mpNodalDomainSizeExpression = ExplicitFilterHelperUtilities::GetNodalDomainSizeExpression(r_elements, r_nodes);
        } else if (number_of_conditions > 0) {
            mpNodalDomainSizeExpression = ExplicitFilterHelperUtilities::GetNodalDomainSizeExpression(r_conditions, r_nodes);
        } else {
            KRATOS_ERROR << "Nodal mapping requires atleast either conditions or elements to be present in "
                         << mrModelPart.FullName() << ".\n";
        }
    }

    for (IndexType i = 0; i < mDampingModelParts.size(); ++i) {
        for (const auto& r_entity : ExplicitFilterHelperUtilities::GetContainer<TContainerType>(*mDampingModelParts[i])) {
            auto p_itr = id_entity_map.find(r_entity.Id());

            KRATOS_ERROR_IF(p_itr == id_entity_map.end())
                << "Entity with id " << r_entity.Id() << " belongs to "
                << mDampingModelParts[i]->FullName()
                << " is not found in main model part " << mrModelPart.FullName() << ".\n";

            p_itr->second->SetDampedComponents(mDampingComponentIndices[i]);
        }
    }

    KRATOS_INFO("ExplicitFilter") << "Search tree created in: " << timer.ElapsedSeconds() << " s" << std::endl;

    KRATOS_CATCH("");
}

template<class TContainerType>
template<class TWeightIntegrationType>
ContainerExpression<TContainerType> ExplicitFilter<TContainerType>::GenericFilterField(const ContainerExpression<TContainerType>& rContainerExpression) const
{
    KRATOS_TRY

    using tls = ExplicitFilterHelperUtilities::TLS<typename EntityPoint<EntityType>::Pointer>;

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

    KRATOS_ERROR_IF(mNumberOfComponents > 0 && mNumberOfComponents != stride)
        << "Damping number of components and expression stride mismatch [ number of damping components = "
        << mNumberOfComponents << ", stride = " << stride << " ].\n";

    const auto& r_origin_expression = rContainerExpression.GetExpression();
    const auto& r_container = rContainerExpression.GetContainer();
    const auto& r_filter_radius_expression = mpFilterRadiusContainer->GetExpression();

    ContainerExpression<TContainerType> result(*rContainerExpression.pGetModelPart());
    auto p_expression = LiteralFlatExpression<double>::Create(result.GetContainer().size(), rContainerExpression.GetItemShape());
    result.SetExpression(p_expression);

    IndexPartition<IndexType>(r_container.size()).for_each(tls(mMaxNumberOfNeighbors), [&](const IndexType Index, tls& rTLS){
        const double radius = r_filter_radius_expression.Evaluate(Index, Index, 0);

        EntityPoint<EntityType> entity_point(*(r_container.begin() + Index), Index);
        const auto number_of_neighbors = mpSearchTree->SearchInRadius(
                                            entity_point,
                                            radius,
                                            rTLS.mNeighbourEntityPoints.begin(),
                                            rTLS.mResultingSquaredDistances.begin(),
                                            mMaxNumberOfNeighbors);

        KRATOS_ERROR_IF(number_of_neighbors >= mMaxNumberOfNeighbors) << "Maximum number of neighbours reached.";

        std::vector<double> list_of_weights(number_of_neighbors, 0.0);
        double sum_of_weights = 0.0;
        ExplicitFilterHelperUtilities::ComputeWeightForAllNeighbors<EntityType, TWeightIntegrationType>(
            sum_of_weights, list_of_weights, *mpKernelFunction, radius,
            entity_point, rTLS.mNeighbourEntityPoints, number_of_neighbors,
            this->mpNodalDomainSizeExpression.get());

        const IndexType current_data_begin = Index * stride;

        for (IndexType j = 0; j < stride; ++j) {
            double& current_index_value = *(p_expression->begin() + current_data_begin + j);
            current_index_value = 0.0;
            double damping_weight = 1.0;
            for(IndexType neighbour_index = 0 ; neighbour_index < number_of_neighbors; ++neighbour_index) {
                const IndexType neighbour_id = rTLS.mNeighbourEntityPoints[neighbour_index]->Id();
                const double weight = list_of_weights[neighbour_index] / sum_of_weights;
                const double origin_value = r_origin_expression.Evaluate(neighbour_id, neighbour_id * stride, j);
                current_index_value +=  weight * origin_value;

                if (mNumberOfComponents > 0 && rTLS.mNeighbourEntityPoints[neighbour_index]->GetDampedComponents()[j]) {
                    damping_weight = std::min(damping_weight, 1 - mpDampingFunction->ComputeWeight(entity_point.Coordinates(), rTLS.mNeighbourEntityPoints[neighbour_index]->Coordinates(), radius));
                }
            }

            // apply damping
            current_index_value *= damping_weight;
        }
    });

    return result;

    KRATOS_CATCH("");
}

template<class TContainerType>
ContainerExpression<TContainerType> ExplicitFilter<TContainerType>::FilterField(const ContainerExpression<TContainerType>& rContainerExpression) const
{
    return GenericFilterField<ExplicitFilterHelperUtilities::IntegratedWeight>(rContainerExpression);
}

template<class TContainerType>
ContainerExpression<TContainerType> ExplicitFilter<TContainerType>::FilterIntegratedField(const ContainerExpression<TContainerType>& rContainerExpression) const
{
    return GenericFilterField<ExplicitFilterHelperUtilities::NonIntegratedWeight>(rContainerExpression);
}

template<class TContainerType>
void ExplicitFilter<TContainerType>::GetIntegrationWeights(ContainerExpression<TContainerType>& rContainerExpression) const
{
    KRATOS_ERROR_IF_NOT(&rContainerExpression.GetModelPart() == &mrModelPart)
        << "The given container expression model part and filter model part mismatch.";

    const IndexType stride = rContainerExpression.GetItemComponentCount();
    const auto& r_container = rContainerExpression.GetContainer();
    auto p_expression = LiteralFlatExpression<double>::Create(r_container.GetContainer().size(), rContainerExpression.GetItemShape());
    rContainerExpression.SetExpression(p_expression);

    IndexPartition<IndexType>(r_container.size()).for_each([&](const IndexType Index){
        const EntityPoint<EntityType> entity(*(r_container.begin() + Index), Index);
        const auto integration_weight = ExplicitFilterHelperUtilities::GetDomainSize(entity, this->mpNodalDomainSizeExpression.get());
        const IndexType current_data_begin = Index * stride;
        for (IndexType j = 0; j < stride; ++j) {
            double& current_index_value = *(p_expression->begin() + current_data_begin + j);
            current_index_value = integration_weight;
        }
    });
}

template<class TContainerType>
std::string ExplicitFilter<TContainerType>::Info() const
{
    std::stringstream msg;

    msg << "ExplicitFilter: ";

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
#define KRATOS_INSTANTIATE_EXPLICIT_FILTER_METHODS(CONTAINER_TYPE)                                                                                                                                                      \
    template class ExplicitFilter<CONTAINER_TYPE>;                                                                                                                                                                        \
    template ContainerExpression<CONTAINER_TYPE> ExplicitFilter<CONTAINER_TYPE>::GenericFilterField<ExplicitFilterHelperUtilities::IntegratedWeight>(const ContainerExpression<CONTAINER_TYPE>&) const;     \
    template ContainerExpression<CONTAINER_TYPE> ExplicitFilter<CONTAINER_TYPE>::GenericFilterField<ExplicitFilterHelperUtilities::NonIntegratedWeight>(const ContainerExpression<CONTAINER_TYPE>&) const;

KRATOS_INSTANTIATE_EXPLICIT_FILTER_METHODS(ModelPart::NodesContainerType)
KRATOS_INSTANTIATE_EXPLICIT_FILTER_METHODS(ModelPart::ConditionsContainerType)
KRATOS_INSTANTIATE_EXPLICIT_FILTER_METHODS(ModelPart::ElementsContainerType)

#undef KRATOS_INSTANTIATE_EXPLICIT_FILTER_METHODS

} // namespace Kratos