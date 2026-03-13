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

#pragma once

// System includes
#include <string>
#include <variant>

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "spatial_containers/spatial_containers.h"
#include "tensor_adaptors/tensor_adaptor.h"

// Application includes
#include "entity_point.h"
#include "filter_function.h"
#include "explicit_damping.h"

namespace Kratos {

///@name Kratos Classes
///@{

template<class TContainerType>
class KRATOS_API(OPTIMIZATION_APPLICATION) ExplicitFilterUtils
{
public:
    ///@name Type definitions
    ///@{

    using EntityType = typename TContainerType::value_type;

    using EntityPointType = EntityPoint<EntityType>;

    using EntityPointVector = std::vector<typename EntityPointType::Pointer>;

    // Type definitions for tree-search
    using BucketType = Bucket<3, EntityPoint<EntityType>, EntityPointVector>;

    using KDTree = Tree<KDTreePartition<BucketType>>;

    /// Pointer definition of ContainerMapper
    KRATOS_CLASS_POINTER_DEFINITION(ExplicitFilterUtils);

    ///@}
    ///@name LifeCycle
    ///@{

    ExplicitFilterUtils(
        const ModelPart& rModelPart,
        const std::string& rKernelFunctionType,
        const IndexType MaxNumberOfNeighbours,
        const IndexType EchoLevel,
        const bool NodeCloudMesh);

    ///@}
    ///@name Public operations

    void SetRadius(TensorAdaptor<double>::Pointer pTensorAdaptor);

    TensorAdaptor<double>::Pointer GetRadius() const;

    void SetDamping(typename ExplicitDamping<TContainerType>::Pointer pExplicitDamping);

    /**
     * @brief Updates the internal KD trees or searching neighbours
     *
     */
    void Update();

    /**
     * @brief Filters the given input.
     * @details This method filters the given control space mesh-independent input field using the following formulae where
     *          \f$\tilde{\underline{s}}\f$ represents the control space variables,
     *          \f$\underline{s}\f$ represents the physical space variables, \f$\mathbf{D}\f$ is
     *          the damping matrix and \f$\mathbf{A}\f$ is the smoothing matrix.
     *          \f[
     *              \Delta \underline{s} = \mathbf{D}\mathbf{A}\Delta \tilde{\underline{s}}
     *          \f]
     *
     * @param rTensorAdaptor  mesh-independent update field in control space.
     * @return TensorAdaptor<double>::Pointer Filtered/Smoothened mesh-independent update field in physical space
     */
    TensorAdaptor<double>::Pointer ForwardFilterField(const TensorAdaptor<double>& rTensorAdaptor) const;

    /**
     * @brief Filters the given mesh-independent physical space gradients to mesh-independent control space gradients.
     * @details This method transforms physical space gradients to control space gradients
     *          by using the transpose of the @ref ForwardFilterField method.
     *
     * @param rTensorAdaptor  Mesh-independent physical space gradient.
     * @return TensorAdaptor<double>::Pointer Mesh-independent control space gradient.
     */
    TensorAdaptor<double>::Pointer BackwardFilterField(const TensorAdaptor<double>& rTensorAdaptor) const;

    /**
     * @brief Filters the given mesh-dependent physical space gradients to mesh-independent control space gradients.
     * @details This method transforms physical space gradients to control space gradients
     *          by using the transpose of the @ref ForwardFilterField method.
     *
     * @param rTensorAdaptor  Mesh-dependent physical space gradient.
     * @return TensorAdaptor<double>::Pointer Mesh-independent control space gradient.
     */
    TensorAdaptor<double>::Pointer BackwardFilterIntegratedField(const TensorAdaptor<double>& rTensorAdaptor) const;

    /**
     * @brief Get the Integration Weights object
     */
    void GetIntegrationWeights(TensorAdaptor<double>& rTensorAdaptor) const;

    /**
     * @brief Calculates the filtering matrix used in this filter.
     * @details This method only calculate the filtering matrix and damping is not applied at all.
     *          If one wishes to calculate the filtering matrix with damping, then they can use the
     *          Damping coefficient matrix along with this filtering matrix to do so.
     *
     * @param rOutput           Output filtering matrix
     */
    void CalculateMatrix(Matrix& rOutput) const;

    /**
     * @brief Prints info about the filtering.
     */
    std::string Info() const;

    ///@}
private:
    ///@name Private member variables
    ///@{

    const ModelPart& mrModelPart;

    FilterFunction::UniquePointer mpKernelFunction;

    TensorAdaptor<double>::Pointer mpFilterRadiusTensorAdaptor;

    typename ExplicitDamping<TContainerType>::Pointer mpDamping;

    std::vector<double> mNodalDomainSizes;

    EntityPointVector mEntityPointVector;

    IndexType mBucketSize = 100;

    IndexType mMaxNumberOfNeighbors;

    IndexType mEchoLevel;

    typename KDTree::Pointer mpSearchTree;

    bool mNodeCloudMesh;

    ///@}
    ///@name Private operations
    ///@{

    void CheckField(const TensorAdaptor<double>& rTensorAdaptor) const;

    template<class TMeshDependencyType>
    void GenericGetIntegrationWeights(TensorAdaptor<double>& rTensorAdapto) const;

    template<class TMeshDependencyType>
    TensorAdaptor<double>::Pointer GenericForwardFilterField(const TensorAdaptor<double>& rTensorAdaptor) const;

    template<class TMeshDependencyType>
    TensorAdaptor<double>::Pointer GenericBackwardFilterField(const TensorAdaptor<double>& rTensorAdaptor) const;

    ///@}
};

///@name Input and output
///@{

/// output stream function
template<class TContainerType>
inline std::ostream& operator<<(
    std::ostream& rOStream,
    const ExplicitFilterUtils<TContainerType>& rThis)
{
    return rOStream << rThis.Info();
}

///@}

///@}
} // namespace Kratos