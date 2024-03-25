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
#include "utilities/atomic_utilities.h"
#include "utilities/builtin_timer.h"
#include "utilities/model_part_utils.h"
#include "utilities/parallel_utilities.h"

// Application includes
#include "custom_utilities/optimization_utils.h"

// Include base h
#include "explicit_filter_utils.h"

namespace Kratos {

namespace ExplicitFilterUtilsHelperUtilities
{
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

struct MeshIndependentType
{
    static double Compute(
        const double Value,
        const double DomainSize)
    {
        return Value;
    }
};

struct MeshDependentType
{
    static double Compute(
        const double Value,
        const double DomainSize)
    {
        return Value / DomainSize;
    }
};

template <class TEntityType>
void ComputeWeightForAllNeighbors(
    double& rSumOfWeights,
    std::vector<double>& rListOfWeights,
    const FilterFunction& rFilterFunction,
    const double Radius,
    const EntityPoint<TEntityType>& rDesignPoint,
    const std::vector<typename EntityPoint<TEntityType>::Pointer>& rNeighbourNodes,
    const std::vector<double>& rSquaredDistances,
    const IndexType NumberOfNeighbours,
    Expression const * const pExpression)
{
    for (IndexType i = 0; i < NumberOfNeighbours; ++i) {
        const double domain_size = GetDomainSize(*rNeighbourNodes[i], pExpression);
        const double filter_weight = rFilterFunction.ComputeWeight(Radius, std::sqrt(rSquaredDistances[i])) * domain_size;
        rListOfWeights[i] = filter_weight;
        rSumOfWeights += filter_weight;
    }
}

}; // namespace ExplicitFilterUtilsHelperUtilities

template<class TContainerType>
ExplicitFilterUtils<TContainerType>::ExplicitFilterUtils(
    const ModelPart& rModelPart,
    const std::string& rKernelFunctionType,
    const IndexType MaxNumberOfNeighbours,
    const IndexType EchoLevel)
    : mrModelPart(rModelPart),
      mMaxNumberOfNeighbors(MaxNumberOfNeighbours),
      mEchoLevel(EchoLevel)
{
    mpKernelFunction = Kratos::make_unique<FilterFunction>(rKernelFunctionType);
}

template<class TContainerType>
void ExplicitFilterUtils<TContainerType>::SetRadius(const ContainerExpression<TContainerType>& rContainerExpression)
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
ContainerExpression<TContainerType> ExplicitFilterUtils<TContainerType>::GetRadius() const
{
    return *mpFilterRadiusContainer;
}

template<class TContainerType>
void ExplicitFilterUtils<TContainerType>::SetDamping(typename ExplicitDamping<TContainerType>::Pointer pExplicitDamping)
{
    mpDamping = pExplicitDamping;
}

template<class TContainerType>
void ExplicitFilterUtils<TContainerType>::ExplicitFilterUtils::Update()
{
    KRATOS_TRY

    BuiltinTimer timer;

    const auto& r_container = ModelPartUtils::GetContainer<TContainerType>(mrModelPart);

    if (mEntityPointVector.size() != r_container.size()) {
        mEntityPointVector.resize(r_container.size());
    }

    // now fill the points vector
    IndexPartition<IndexType>(r_container.size()).for_each([&](const IndexType Index) {
        mEntityPointVector[Index] = Kratos::make_shared<EntityPoint<EntityType>>(*(r_container.begin() + Index), Index);
    });

    mpSearchTree =  Kratos::make_shared<ExplicitFilterUtils::KDTree>(mEntityPointVector.begin(), mEntityPointVector.end(), mBucketSize);

    if constexpr(std::is_same_v<TContainerType, ModelPart::NodesContainerType>) {
        const auto& r_nodes = mrModelPart.Nodes();
        const auto& r_conditions = mrModelPart.Conditions();
        const IndexType number_of_conditions = r_conditions.size();

        const auto& r_elements = mrModelPart.Elements();
        const IndexType number_of_elements = r_elements.size();

        if (number_of_elements > 0) {
            mpNodalDomainSizeExpression = ExplicitFilterUtilsHelperUtilities::GetNodalDomainSizeExpression(r_elements, r_nodes);
        } else if (number_of_conditions > 0) {
            mpNodalDomainSizeExpression = ExplicitFilterUtilsHelperUtilities::GetNodalDomainSizeExpression(r_conditions, r_nodes);
        } else {
            KRATOS_ERROR << "Nodal mapping requires atleast either conditions or elements to be present in "
                         << mrModelPart.FullName() << ".\n";
        }
    }

    KRATOS_INFO_IF("ExplicitFilterUtils", mEchoLevel > 0) << "Updated filter in: " << timer.ElapsedSeconds() << " s" << std::endl;

    KRATOS_CATCH("");
}

template<class TContainerType>
void ExplicitFilterUtils<TContainerType>::CheckField(const ContainerExpression<TContainerType>& rContainerExpression) const
{
    KRATOS_ERROR_IF(mpFilterRadiusContainer.get() == nullptr)
        << "The filter radius container expression not set. "
        << "Please set it using SetRadius method.\n\t Filter = "
        << *this;

    KRATOS_ERROR_IF_NOT(rContainerExpression.HasExpression())
        << "Uninitialized container expression given. "
        << rContainerExpression;

    KRATOS_ERROR_IF_NOT(&rContainerExpression.GetModelPart() == &mrModelPart)
        << "Filter radius container expression model part and filter model part mismatch."
        << "\n\tFilter = " << *this
        << "\n\tContainerExpression = " << rContainerExpression;

    KRATOS_ERROR_IF(mpDamping.get() == nullptr)
        << "Damping method is not set. "
        << "Please set it using SetDamping method.\n\t Filter = "
        << *this;
}

template<class TContainerType>
ContainerExpression<TContainerType> ExplicitFilterUtils<TContainerType>::ForwardFilterField(const ContainerExpression<TContainerType>& rContainerExpression) const
{
    KRATOS_TRY

    using tls = OptimizationUtils::KDTreeThreadLocalStorage<typename EntityPointType::Pointer>;

    CheckField(rContainerExpression);

    const IndexType stride = rContainerExpression.GetItemComponentCount();
    const auto& r_origin_expression = rContainerExpression.GetExpression();
    const auto& r_container = rContainerExpression.GetContainer();
    const auto& r_filter_radius_expression = mpFilterRadiusContainer->GetExpression();

    KRATOS_ERROR_IF_NOT(stride == mpDamping->GetStride())
        << "Damping stride and expression stride mismatch. [ damping stride = "
        << mpDamping->GetStride() << ", expression stride = " << stride << " ].\n";

    ContainerExpression<TContainerType> result(*rContainerExpression.pGetModelPart());
    auto p_expression = LiteralFlatExpression<double>::Create(result.GetContainer().size(), rContainerExpression.GetItemShape());
    result.SetExpression(p_expression);

    IndexPartition<IndexType>(r_container.size()).for_each(tls(mMaxNumberOfNeighbors, stride), [&](const IndexType Index, tls& rTLS){
        const double radius = r_filter_radius_expression.Evaluate(Index, Index, 0);

        EntityPoint<EntityType> entity_point(*(r_container.begin() + Index), Index);
        const auto number_of_neighbors = mpSearchTree->SearchInRadius(
                                            entity_point,
                                            radius,
                                            rTLS.mNeighbourEntityPoints.begin(),
                                            rTLS.mResultingSquaredDistances.begin(),
                                            mMaxNumberOfNeighbors);

        KRATOS_ERROR_IF(number_of_neighbors >= mMaxNumberOfNeighbors)
            << "Maximum number of allowed neighbours reached when searching for neighbours in "
            << mrModelPart.FullName() << " with radii = " << radius << " [ max number of allowed neighbours = "
            << mMaxNumberOfNeighbors << " ].\n";

        double sum_of_weights = 0.0;
        ExplicitFilterUtilsHelperUtilities::ComputeWeightForAllNeighbors(
            sum_of_weights, rTLS.mListOfWeights, *mpKernelFunction, radius,
            entity_point, rTLS.mNeighbourEntityPoints, rTLS.mResultingSquaredDistances, number_of_neighbors, this->mpNodalDomainSizeExpression.get());

        mpDamping->Apply(rTLS.mListOfDampedWeights, rTLS.mListOfWeights, Index, number_of_neighbors, rTLS.mNeighbourEntityPoints);

        for (IndexType j = 0; j < stride; ++j) {
            const auto& r_damped_weights = rTLS.mListOfDampedWeights[j];
            double& current_index_value = *(p_expression->begin() + Index * stride + j);
            current_index_value = 0.0;
            for(IndexType neighbour_index = 0 ; neighbour_index < number_of_neighbors; ++neighbour_index) {
                const IndexType neighbour_id = rTLS.mNeighbourEntityPoints[neighbour_index]->Id();
                const double weight = r_damped_weights[neighbour_index] / sum_of_weights;
                const double origin_value = r_origin_expression.Evaluate(neighbour_id, neighbour_id * stride, j);
                current_index_value +=  weight * origin_value;
            }
        }
    });

    KRATOS_INFO_IF("ExplicitFilterUtils", mEchoLevel > 1) << "Computed forward filter field." << std::endl;

    return result;

    KRATOS_CATCH("");
}
template<class TContainerType>
template<class TMeshDependencyType>
ContainerExpression<TContainerType> ExplicitFilterUtils<TContainerType>::GenericBackwardFilterField(const ContainerExpression<TContainerType>& rContainerExpression) const
{
    KRATOS_TRY

    using tls = OptimizationUtils::KDTreeThreadLocalStorage<typename EntityPointType::Pointer>;

    CheckField(rContainerExpression);

    const IndexType stride = rContainerExpression.GetItemComponentCount();
    const auto& r_origin_expression = rContainerExpression.GetExpression();
    const auto& r_container = rContainerExpression.GetContainer();
    const auto& r_filter_radius_expression = mpFilterRadiusContainer->GetExpression();

    KRATOS_ERROR_IF_NOT(stride == mpDamping->GetStride())
        << "Damping stride and expression stride mismatch. [ damping stride = "
        << mpDamping->GetStride() << ", expression stride = " << stride << " ].\n";

    ContainerExpression<TContainerType> result(*rContainerExpression.pGetModelPart());
    auto p_expression = LiteralFlatExpression<double>::Create(result.GetContainer().size(), rContainerExpression.GetItemShape());
    result.SetExpression(p_expression);

    IndexPartition<IndexType>(result.GetContainer().size() * stride).for_each([&p_expression](const auto Index) {
        *(p_expression->begin() + Index) = 0.0;
    });

    IndexPartition<IndexType>(r_container.size()).for_each(tls(mMaxNumberOfNeighbors, stride), [&](const IndexType Index, tls& rTLS){
        const double radius = r_filter_radius_expression.Evaluate(Index, Index, 0);

        EntityPoint<EntityType> entity_point(*(r_container.begin() + Index), Index);
        const auto number_of_neighbors = mpSearchTree->SearchInRadius(
                                            entity_point,
                                            radius,
                                            rTLS.mNeighbourEntityPoints.begin(),
                                            rTLS.mResultingSquaredDistances.begin(),
                                            mMaxNumberOfNeighbors);

        KRATOS_ERROR_IF(number_of_neighbors >= mMaxNumberOfNeighbors)
            << "Maximum number of allowed neighbours reached when searching for neighbours in "
            << mrModelPart.FullName() << " with radii = " << radius << " [ max number of allowed neighbours = "
            << mMaxNumberOfNeighbors << " ].\n";

        double sum_of_weights = 0.0;
        ExplicitFilterUtilsHelperUtilities::ComputeWeightForAllNeighbors(
            sum_of_weights, rTLS.mListOfWeights, *mpKernelFunction, radius,
            entity_point, rTLS.mNeighbourEntityPoints, rTLS.mResultingSquaredDistances, number_of_neighbors, this->mpNodalDomainSizeExpression.get());

        mpDamping->Apply(rTLS.mListOfDampedWeights, rTLS.mListOfWeights, Index, number_of_neighbors, rTLS.mNeighbourEntityPoints);

        const IndexType current_data_begin = Index * stride;
        const double domain_size = ExplicitFilterUtilsHelperUtilities::GetDomainSize(entity_point, mpNodalDomainSizeExpression.get());

        for (IndexType j = 0; j < stride; ++j) {
            const auto& r_damped_weights = rTLS.mListOfDampedWeights[j];
            const double origin_value = TMeshDependencyType::Compute(r_origin_expression.Evaluate(Index, current_data_begin, j), domain_size);

            for(IndexType neighbour_index = 0; neighbour_index < number_of_neighbors; ++neighbour_index) {
                const double weight = r_damped_weights[neighbour_index] / sum_of_weights;

                const IndexType neighbour_id = rTLS.mNeighbourEntityPoints[neighbour_index]->Id();
                const IndexType neighbour_data_begin_index = neighbour_id * stride;

                AtomicAdd<double>(*(p_expression->begin() + neighbour_data_begin_index + j), origin_value * weight);
            }
        }
    });

    KRATOS_INFO_IF("ExplicitFilterUtils", mEchoLevel > 1) << "Computed backward filter field." << std::endl;

    return result;

    KRATOS_CATCH("");
}

template<class TContainerType>
ContainerExpression<TContainerType> ExplicitFilterUtils<TContainerType>::BackwardFilterField(const ContainerExpression<TContainerType>& rContainerExpression) const
{
    return GenericBackwardFilterField<ExplicitFilterUtilsHelperUtilities::MeshIndependentType>(rContainerExpression);
}

template<class TContainerType>
ContainerExpression<TContainerType> ExplicitFilterUtils<TContainerType>::BackwardFilterIntegratedField(const ContainerExpression<TContainerType>& rContainerExpression) const
{
    return GenericBackwardFilterField<ExplicitFilterUtilsHelperUtilities::MeshDependentType>(rContainerExpression);
}

template<class TContainerType>
void ExplicitFilterUtils<TContainerType>::GetIntegrationWeights(ContainerExpression<TContainerType>& rContainerExpression) const
{

    KRATOS_ERROR_IF_NOT(&rContainerExpression.GetModelPart() == &mrModelPart)
        << "The given container expression model part and filter model part mismatch.";

    const IndexType stride = rContainerExpression.GetItemComponentCount();
    const auto& r_container = rContainerExpression.GetContainer();
    auto p_expression = LiteralFlatExpression<double>::Create(r_container.GetContainer().size(), rContainerExpression.GetItemShape());
    rContainerExpression.SetExpression(p_expression);

    IndexPartition<IndexType>(r_container.size()).for_each([&](const IndexType Index){
        const EntityPoint<EntityType> entity(*(r_container.begin() + Index), Index);
        const auto integration_weight = ExplicitFilterUtilsHelperUtilities::GetDomainSize(entity, this->mpNodalDomainSizeExpression.get());
        const IndexType current_data_begin = Index * stride;
        for (IndexType j = 0; j < stride; ++j) {
            double& current_index_value = *(p_expression->begin() + current_data_begin + j);
            current_index_value = integration_weight;
        }
    });
}

template<class TContainerType>
std::string ExplicitFilterUtils<TContainerType>::Info() const
{
    std::stringstream msg;

    msg << "ExplicitFilterUtils: ";

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

template<class TContainerType>
void ExplicitFilterUtils<TContainerType>::CalculateMatrix(Matrix& rOutput) const
{
    KRATOS_TRY

    using tls = OptimizationUtils::KDTreeThreadLocalStorage<typename EntityPointType::Pointer>;

    const auto number_of_entities = mEntityPointVector.size();

    const auto& r_filter_radius_expression = mpFilterRadiusContainer->GetExpression();

    if (rOutput.size1 () != number_of_entities || rOutput.size2() != number_of_entities) {
        rOutput.resize(number_of_entities, number_of_entities, false);
    }

    noalias(rOutput) = ZeroMatrix(number_of_entities, number_of_entities);

    IndexPartition<IndexType>(number_of_entities).for_each(tls(mMaxNumberOfNeighbors, 1), [&](const auto Index, auto& rTLS) {
        const double radius = r_filter_radius_expression.Evaluate(Index, Index, 0);

        const auto number_of_neighbors = mpSearchTree->SearchInRadius(
                                            *mEntityPointVector[Index],
                                            radius,
                                            rTLS.mNeighbourEntityPoints.begin(),
                                            rTLS.mResultingSquaredDistances.begin(),
                                            mMaxNumberOfNeighbors);

        KRATOS_ERROR_IF(number_of_neighbors >= mMaxNumberOfNeighbors)
            << "Maximum number of allowed neighbours reached when searching for neighbours in "
            << mrModelPart.FullName() << " with radii = " << radius << " [ max number of allowed neighbours = "
            << mMaxNumberOfNeighbors << " ].\n";

        std::vector<double> list_of_weights(number_of_neighbors, 0.0);
        double sum_of_weights = 0.0;
        ExplicitFilterUtilsHelperUtilities::ComputeWeightForAllNeighbors(
            sum_of_weights, list_of_weights, *mpKernelFunction, radius,
            *mEntityPointVector[Index], rTLS.mNeighbourEntityPoints, rTLS.mResultingSquaredDistances, number_of_neighbors, this->mpNodalDomainSizeExpression.get());

        double* data_begin = (rOutput.data().begin() + Index * number_of_entities);

        for (IndexType neighbour_index = 0; neighbour_index < number_of_neighbors; ++neighbour_index) {
            const IndexType neighbour_id = rTLS.mNeighbourEntityPoints[neighbour_index]->Id();
            *(data_begin + neighbour_id) = list_of_weights[neighbour_index] / sum_of_weights;
        }
    });

    KRATOS_CATCH("");
}

// template instantiations
#define KRATOS_INSTANTIATE_EXPLICIT_FILTER_METHODS(CONTAINER_TYPE)                                                                                                                                                      \
    template class ExplicitFilterUtils<CONTAINER_TYPE>;                                                                                                                                                                        \
    template ContainerExpression<CONTAINER_TYPE> ExplicitFilterUtils<CONTAINER_TYPE>::GenericBackwardFilterField<ExplicitFilterUtilsHelperUtilities::MeshIndependentType>(const ContainerExpression<CONTAINER_TYPE>&) const;     \
    template ContainerExpression<CONTAINER_TYPE> ExplicitFilterUtils<CONTAINER_TYPE>::GenericBackwardFilterField<ExplicitFilterUtilsHelperUtilities::MeshDependentType>(const ContainerExpression<CONTAINER_TYPE>&) const;

KRATOS_INSTANTIATE_EXPLICIT_FILTER_METHODS(ModelPart::NodesContainerType)
KRATOS_INSTANTIATE_EXPLICIT_FILTER_METHODS(ModelPart::ConditionsContainerType)
KRATOS_INSTANTIATE_EXPLICIT_FILTER_METHODS(ModelPart::ElementsContainerType)

#undef KRATOS_INSTANTIATE_EXPLICIT_FILTER_METHODS

} // namespace Kratos