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
void GetNodalDomainSizes(
    std::vector<double>& rNodalDomainSizes,
    const TContainerType& rContainer,
    const ModelPart::NodesContainerType& rNodes)
{
    const IndexType number_of_nodes = rNodes.size();

    rNodalDomainSizes.resize(number_of_nodes);
    std::fill(rNodalDomainSizes.begin(), rNodalDomainSizes.end(), 0.0);

    IndexPartition<IndexType>(rContainer.size()).for_each([&rNodes, &rContainer, &rNodalDomainSizes](const IndexType Index) {
        const auto& r_geometry = (rContainer.begin() + Index)->GetGeometry();
        const double nodal_domain_size = r_geometry.DomainSize() / r_geometry.size();
        for (const auto& r_node : r_geometry) {
            const auto node_ptr = rNodes.find(r_node.Id());
            if (node_ptr != rNodes.end()) {
                const IndexType node_index = std::distance(rNodes.begin(), node_ptr);
                AtomicAdd<double>(rNodalDomainSizes[node_index], nodal_domain_size);
            } else {
                KRATOS_ERROR << "Node with id " << r_node.Id() << " not found.";
            }
        }
    });
}

struct NodeCloudsType
{
    template<class TEntityType>
    static inline double GetDomainSize(
        const EntityPoint<TEntityType>& rPoint,
        const std::vector<double>& rNodalDomainSizes)
    {
        if constexpr(!std::is_same_v<TEntityType, ModelPart::NodeType>) {
            KRATOS_ERROR << "Point / node clouds should only be used with nodes.";
        }
        return 1.0;
    }

    static inline double Compute(
        const double Value,
        const double DomainSize)
    {
        return Value;
    }
};

struct MeshIndependentType
{
    template<class TEntityType>
    static inline double GetDomainSize(
        const EntityPoint<TEntityType>& rPoint,
        const std::vector<double>& rNodalDomainSizes)
    {
        if constexpr(std::is_same_v<TEntityType, ModelPart::NodeType>) {
            return rNodalDomainSizes[rPoint.Id()];
        } else {
            return rPoint.GetEntity().GetGeometry().DomainSize();
        }
    }

    static inline double Compute(
        const double Value,
        const double DomainSize)
    {
        return Value;
    }
};

struct MeshDependentType
{
    template<class TEntityType>
    static inline double GetDomainSize(
        const EntityPoint<TEntityType>& rPoint,
        const std::vector<double>& rNodalDomainSizes)
    {
        if constexpr(std::is_same_v<TEntityType, ModelPart::NodeType>) {
            return rNodalDomainSizes[rPoint.Id()];
        } else {
            return rPoint.GetEntity().GetGeometry().DomainSize();
        }
    }

    static inline double Compute(
        const double Value,
        const double DomainSize)
    {
        return Value / DomainSize;
    }
};

template <class TMeshDependencyType, class TEntityType>
void ComputeWeightForAllNeighbors(
    double& rSumOfWeights,
    std::vector<double>& rListOfWeights,
    const FilterFunction& rFilterFunction,
    const double Radius,
    const EntityPoint<TEntityType>& rDesignPoint,
    const std::vector<typename EntityPoint<TEntityType>::Pointer>& rNeighbourNodes,
    const std::vector<double>& rSquaredDistances,
    const IndexType NumberOfNeighbours,
    const std::vector<double>& rNodalDomainSizes)
{
    for (IndexType i = 0; i < NumberOfNeighbours; ++i) {
        const double domain_size = TMeshDependencyType::GetDomainSize(*rNeighbourNodes[i], rNodalDomainSizes);
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
    const IndexType EchoLevel,
    const bool NodeCloudMesh)
    : mrModelPart(rModelPart),
      mMaxNumberOfNeighbors(MaxNumberOfNeighbours),
      mEchoLevel(EchoLevel),
      mNodeCloudMesh(NodeCloudMesh)
{
    mpKernelFunction = Kratos::make_unique<FilterFunction>(rKernelFunctionType);
}

template<class TContainerType>
void ExplicitFilterUtils<TContainerType>::SetRadius(TensorAdaptor<double>::Pointer pTensorAdaptor)
{
    KRATOS_ERROR_IF_NOT(pTensorAdaptor->Shape().size() == 1)
        << "Only scalar values are allowed for the filter radius tensor adaptor."
        << "Provided tensor adaptor = " << *pTensorAdaptor << ".\n";

    if (std::holds_alternative<ModelPart::NodesContainerType::Pointer>(pTensorAdaptor->GetContainer())) {
        const auto& r_container = *(std::get<ModelPart::NodesContainerType::Pointer>(pTensorAdaptor->GetContainer()));
        KRATOS_ERROR_IF_NOT(&r_container == &mrModelPart.Nodes() || &r_container == &mrModelPart.GetCommunicator().LocalMesh().Nodes())
            << "Filter radius tensor adaptor container and filter model part mismatch."
            << "\n\tFilter = " << *this
            << "\n\tTensorAdaptor = " << *pTensorAdaptor;
    } else if (std::holds_alternative<ModelPart::ConditionsContainerType::Pointer>(pTensorAdaptor->GetContainer())) {
        const auto& r_container = *(std::get<ModelPart::ConditionsContainerType::Pointer>(pTensorAdaptor->GetContainer()));
        KRATOS_ERROR_IF_NOT(&r_container == &mrModelPart.Conditions())
            << "Filter radius tensor adaptor container and filter model part mismatch."
            << "\n\tFilter = " << *this
            << "\n\tTensorAdaptor = " << *pTensorAdaptor;
    } else if (std::holds_alternative<ModelPart::ElementsContainerType::Pointer>(pTensorAdaptor->GetContainer())) {
        const auto& r_container = *(std::get<ModelPart::ElementsContainerType::Pointer>(pTensorAdaptor->GetContainer()));
        KRATOS_ERROR_IF_NOT(&r_container == &mrModelPart.Elements())
            << "Filter radius tensor adaptor container and filter model part mismatch."
            << "\n\tFilter = " << *this
            << "\n\tTensorAdaptor = " << *pTensorAdaptor;
    } else {
        KRATOS_ERROR << "Unsupported container type is found in the tensor adaptor. Only supports nodal, condition or elemental tensor adaptors [ tensor adaptor = "
                     << *pTensorAdaptor << " ].\n";
    }

    mpFilterRadiusTensorAdaptor = pTensorAdaptor;
}

template<class TContainerType>
TensorAdaptor<double>::Pointer ExplicitFilterUtils<TContainerType>::GetRadius() const
{
    return mpFilterRadiusTensorAdaptor;
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
        if (!mNodeCloudMesh) {
            const auto& r_nodes = mrModelPart.Nodes();
            const auto& r_conditions = mrModelPart.Conditions();
            const IndexType number_of_conditions = r_conditions.size();

            const auto& r_elements = mrModelPart.Elements();
            const IndexType number_of_elements = r_elements.size();

            if (number_of_elements > 0) {
                ExplicitFilterUtilsHelperUtilities::GetNodalDomainSizes(mNodalDomainSizes, r_elements, r_nodes);
            } else if (number_of_conditions > 0) {
                ExplicitFilterUtilsHelperUtilities::GetNodalDomainSizes(mNodalDomainSizes, r_conditions, r_nodes);
            } else {
                KRATOS_ERROR << "Nodal mapping requires atleast either conditions or elements to be present in "
                            << mrModelPart.FullName() << ".\n";
            }
        }
    }

    KRATOS_INFO_IF("ExplicitFilterUtils", mEchoLevel > 0) << "Updated filter in: " << timer.ElapsedSeconds() << " s" << std::endl;

    KRATOS_CATCH("");
}

template<class TContainerType>
void ExplicitFilterUtils<TContainerType>::CheckField(const TensorAdaptor<double>& rTensorAdaptor) const
{
    KRATOS_ERROR_IF(mpFilterRadiusTensorAdaptor.get() == nullptr)
        << "The filter radius tensor adaptor not set. "
        << "Please set it using SetRadius method.\n\t Filter = "
        << *this;

    if (std::holds_alternative<ModelPart::NodesContainerType::Pointer>(rTensorAdaptor.GetContainer())) {
        const auto& r_container = *(std::get<ModelPart::NodesContainerType::Pointer>(rTensorAdaptor.GetContainer()));
        KRATOS_ERROR_IF_NOT(&r_container == &mrModelPart.Nodes() || &r_container == &mrModelPart.GetCommunicator().LocalMesh().Nodes())
            << "Filter radius tensor adaptor container and filter model part mismatch."
            << "\n\tFilter = " << *this
            << "\n\tTensorAdaptor = " << rTensorAdaptor;
    } else if (std::holds_alternative<ModelPart::ConditionsContainerType::Pointer>(rTensorAdaptor.GetContainer())) {
        const auto& r_container = *(std::get<ModelPart::ConditionsContainerType::Pointer>(rTensorAdaptor.GetContainer()));
        KRATOS_ERROR_IF_NOT(&r_container == &mrModelPart.Conditions())
            << "Filter radius tensor adaptor container and filter model part mismatch."
            << "\n\tFilter = " << *this
            << "\n\tTensorAdaptor = " << rTensorAdaptor;
    } else if (std::holds_alternative<ModelPart::ElementsContainerType::Pointer>(rTensorAdaptor.GetContainer())) {
        const auto& r_container = *(std::get<ModelPart::ElementsContainerType::Pointer>(rTensorAdaptor.GetContainer()));
        KRATOS_ERROR_IF_NOT(&r_container == &mrModelPart.Elements())
            << "Filter radius tensor adaptor container and filter model part mismatch."
            << "\n\tFilter = " << *this
            << "\n\tTensorAdaptor = " << rTensorAdaptor;
    } else {
        KRATOS_ERROR << "Unsupported container type is found in the tensor adaptor. Only supports nodal, condition or elemental tensor adaptors [ tensor adaptor = "
                     << rTensorAdaptor << " ].\n";
    }

    KRATOS_ERROR_IF(mpDamping.get() == nullptr)
        << "Damping method is not set. "
        << "Please set it using SetDamping method.\n\t Filter = "
        << *this;
}

template<class TContainerType>
template<class TMeshDependencyType>
TensorAdaptor<double>::Pointer ExplicitFilterUtils<TContainerType>::GenericForwardFilterField(const TensorAdaptor<double>& rTensorAdaptor) const
{
    KRATOS_TRY

    using tls = OptimizationUtils::KDTreeThreadLocalStorage<typename EntityPointType::Pointer>;

    CheckField(rTensorAdaptor);

    const auto& r_origin_data_view = rTensorAdaptor.ViewData();
    const auto& r_container = *(std::get<typename TContainerType::Pointer>(rTensorAdaptor.GetContainer()));
    const IndexType stride = rTensorAdaptor.Size() / r_container.size();
    const auto& r_filter_radius_data_view = mpFilterRadiusTensorAdaptor->ViewData();

    KRATOS_ERROR_IF_NOT(stride == mpDamping->GetStride())
        << "Damping stride and tensor adaptor stride mismatch. [ damping stride = "
        << mpDamping->GetStride() << ", tensor adaptor stride = " << stride << " ].\n";

    auto p_result_tensor_adaptor = Kratos::make_shared<TensorAdaptor<double>>(rTensorAdaptor);
    auto result_data_view = p_result_tensor_adaptor->ViewData();

    IndexPartition<IndexType>(r_container.size()).for_each(tls(mMaxNumberOfNeighbors, stride), [&](const IndexType Index, tls& rTLS){
        const double radius = r_filter_radius_data_view[Index];

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
        ExplicitFilterUtilsHelperUtilities::ComputeWeightForAllNeighbors<TMeshDependencyType>(
            sum_of_weights, rTLS.mListOfWeights, *mpKernelFunction, radius,
            entity_point, rTLS.mNeighbourEntityPoints, rTLS.mResultingSquaredDistances, number_of_neighbors, this->mNodalDomainSizes);

        mpDamping->Apply(rTLS.mListOfDampedWeights, rTLS.mListOfWeights, Index, number_of_neighbors, rTLS.mNeighbourEntityPoints);

        for (IndexType j = 0; j < stride; ++j) {
            const auto& r_damped_weights = rTLS.mListOfDampedWeights[j];
            double& current_index_value = result_data_view[Index * stride + j];
            current_index_value = 0.0;
            for(IndexType neighbour_index = 0 ; neighbour_index < number_of_neighbors; ++neighbour_index) {
                const IndexType neighbour_id = rTLS.mNeighbourEntityPoints[neighbour_index]->Id();
                const double weight = r_damped_weights[neighbour_index] / sum_of_weights;
                const double origin_value = r_origin_data_view[neighbour_id * stride + j];
                current_index_value +=  weight * origin_value;
            }
        }
    });

    KRATOS_INFO_IF("ExplicitFilterUtils", mEchoLevel > 1) << "Computed forward filter field." << std::endl;

    return p_result_tensor_adaptor;

    KRATOS_CATCH("");
}

template<class TContainerType>
template<class TMeshDependencyType>
TensorAdaptor<double>::Pointer ExplicitFilterUtils<TContainerType>::GenericBackwardFilterField(const TensorAdaptor<double>& rTensorAdaptor) const
{
    KRATOS_TRY

    using tls = OptimizationUtils::KDTreeThreadLocalStorage<typename EntityPointType::Pointer>;

    CheckField(rTensorAdaptor);

    const auto& r_origin_data_view = rTensorAdaptor.ViewData();
    const auto& r_container = *(std::get<typename TContainerType::Pointer>(rTensorAdaptor.GetContainer()));
    const IndexType stride = r_origin_data_view.size() / r_container.size();
    const auto& r_filter_radius_data_view = mpFilterRadiusTensorAdaptor->ViewData();

    KRATOS_ERROR_IF_NOT(stride == mpDamping->GetStride())
        << "Damping stride and tensor adaptor stride mismatch. [ damping stride = "
        << mpDamping->GetStride() << ", tensor adaptor stride = " << stride << " ].\n";

    auto p_result_tensor_adaptor = Kratos::make_shared<TensorAdaptor<double>>(rTensorAdaptor);
    auto result_data_view = p_result_tensor_adaptor->ViewData();

    IndexPartition<IndexType>(result_data_view.size()).for_each([&result_data_view](const auto Index) {
        result_data_view[Index] = 0.0;
    });

    IndexPartition<IndexType>(r_container.size()).for_each(tls(mMaxNumberOfNeighbors, stride), [&](const IndexType Index, tls& rTLS){
        const double radius = r_filter_radius_data_view[Index];

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
        ExplicitFilterUtilsHelperUtilities::ComputeWeightForAllNeighbors<TMeshDependencyType>(
            sum_of_weights, rTLS.mListOfWeights, *mpKernelFunction, radius,
            entity_point, rTLS.mNeighbourEntityPoints, rTLS.mResultingSquaredDistances, number_of_neighbors, this->mNodalDomainSizes);

        mpDamping->Apply(rTLS.mListOfDampedWeights, rTLS.mListOfWeights, Index, number_of_neighbors, rTLS.mNeighbourEntityPoints);

        const IndexType current_data_begin = Index * stride;
        const double domain_size = TMeshDependencyType::GetDomainSize(entity_point, mNodalDomainSizes);

        for (IndexType j = 0; j < stride; ++j) {
            const auto& r_damped_weights = rTLS.mListOfDampedWeights[j];
            const double origin_value = TMeshDependencyType::Compute(r_origin_data_view[current_data_begin + j], domain_size);
            // KRATOS_WATCH(origin_value)

            for(IndexType neighbour_index = 0; neighbour_index < number_of_neighbors; ++neighbour_index) {
                const double weight = r_damped_weights[neighbour_index] / sum_of_weights;

                const IndexType neighbour_id = rTLS.mNeighbourEntityPoints[neighbour_index]->Id();
                const IndexType neighbour_data_begin_index = neighbour_id * stride;

                AtomicAdd<double>(result_data_view[neighbour_data_begin_index + j], origin_value * weight);
            }
        }
    });

    KRATOS_INFO_IF("ExplicitFilterUtils", mEchoLevel > 1) << "Computed backward filter field." << std::endl;

    return p_result_tensor_adaptor;

    KRATOS_CATCH("");
}

template<class TContainerType>
template<class TMeshDependencyType>
void ExplicitFilterUtils<TContainerType>::GenericGetIntegrationWeights(TensorAdaptor<double>& rTensorAdaptor) const
{
    CheckField(rTensorAdaptor);

    const auto& r_container = *(std::get<typename TContainerType::Pointer>(rTensorAdaptor.GetContainer()));
    const IndexType stride = rTensorAdaptor.Size() / r_container.size();
    auto r_data_view = rTensorAdaptor.ViewData();

    IndexPartition<IndexType>(r_container.size()).for_each([&](const IndexType Index){
        const EntityPoint<EntityType> entity(*(r_container.begin() + Index), Index);
        const auto integration_weight = TMeshDependencyType::GetDomainSize(entity, this->mNodalDomainSizes);
        const IndexType current_data_begin = Index * stride;
        for (IndexType j = 0; j < stride; ++j) {
            r_data_view[current_data_begin + j] = integration_weight;
        }
    });
}

template<class TContainerType>
TensorAdaptor<double>::Pointer ExplicitFilterUtils<TContainerType>::ForwardFilterField(const TensorAdaptor<double>& rTensorAdaptor) const
{
    if (mNodeCloudMesh) {
        return GenericForwardFilterField<ExplicitFilterUtilsHelperUtilities::NodeCloudsType>(rTensorAdaptor);
    } else {
        return GenericForwardFilterField<ExplicitFilterUtilsHelperUtilities::MeshIndependentType>(rTensorAdaptor);
    }
}

template<class TContainerType>
TensorAdaptor<double>::Pointer ExplicitFilterUtils<TContainerType>::BackwardFilterField(const TensorAdaptor<double>& rTensorAdaptor) const
{
    if (mNodeCloudMesh) {
        return GenericBackwardFilterField<ExplicitFilterUtilsHelperUtilities::NodeCloudsType>(rTensorAdaptor);
    } else {
        return GenericBackwardFilterField<ExplicitFilterUtilsHelperUtilities::MeshIndependentType>(rTensorAdaptor);
    }
}

template<class TContainerType>
TensorAdaptor<double>::Pointer ExplicitFilterUtils<TContainerType>::BackwardFilterIntegratedField(const TensorAdaptor<double>& rTensorAdaptor) const
{
    if (mNodeCloudMesh) {
        return GenericBackwardFilterField<ExplicitFilterUtilsHelperUtilities::NodeCloudsType>(rTensorAdaptor);
    } else {
        return GenericBackwardFilterField<ExplicitFilterUtilsHelperUtilities::MeshDependentType>(rTensorAdaptor);
    }
}

template<class TContainerType>
void ExplicitFilterUtils<TContainerType>::GetIntegrationWeights(TensorAdaptor<double>& rTensorAdaptor) const
{
    if (mNodeCloudMesh) {
        GenericGetIntegrationWeights<ExplicitFilterUtilsHelperUtilities::NodeCloudsType>(rTensorAdaptor);
    } else {
        GenericGetIntegrationWeights<ExplicitFilterUtilsHelperUtilities::MeshIndependentType>(rTensorAdaptor);
    }
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

    const auto r_filter_radius_data_view = mpFilterRadiusTensorAdaptor->ViewData();

    if (rOutput.size1 () != number_of_entities || rOutput.size2() != number_of_entities) {
        rOutput.resize(number_of_entities, number_of_entities, false);
    }

    noalias(rOutput) = ZeroMatrix(number_of_entities, number_of_entities);

    IndexPartition<IndexType>(number_of_entities).for_each(tls(mMaxNumberOfNeighbors, 1), [&](const auto Index, auto& rTLS) {
        const double radius = r_filter_radius_data_view[Index];

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
        ExplicitFilterUtilsHelperUtilities::ComputeWeightForAllNeighbors<ExplicitFilterUtilsHelperUtilities::MeshIndependentType>(
            sum_of_weights, list_of_weights, *mpKernelFunction, radius,
            *mEntityPointVector[Index], rTLS.mNeighbourEntityPoints, rTLS.mResultingSquaredDistances, number_of_neighbors, this->mNodalDomainSizes);

        double* data_begin = (rOutput.data().begin() + Index * number_of_entities);

        for (IndexType neighbour_index = 0; neighbour_index < number_of_neighbors; ++neighbour_index) {
            const IndexType neighbour_id = rTLS.mNeighbourEntityPoints[neighbour_index]->Id();
            *(data_begin + neighbour_id) = list_of_weights[neighbour_index] / sum_of_weights;
        }
    });

    KRATOS_CATCH("");
}

// template instantiations
template class ExplicitFilterUtils<ModelPart::NodesContainerType>;
template class ExplicitFilterUtils<ModelPart::ConditionsContainerType>;
template class ExplicitFilterUtils<ModelPart::ElementsContainerType>;

} // namespace Kratos