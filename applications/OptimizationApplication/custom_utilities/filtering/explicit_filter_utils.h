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

// External includes
#include "nanoflann/include/nanoflann.hpp"

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/ublas_interface.h"
#include "tensor_adaptors/tensor_adaptor.h"

// Application includes
#include "filter_function.h"
#include "explicit_damping.h"
#include "filter_utils.h"

namespace Kratos {

///@name Kratos Classes
///@{

template<class TContainerType>
class KRATOS_API(OPTIMIZATION_APPLICATION) ExplicitFilterUtils
{
public:
    ///@name Type definitions
    ///@{

    /// Pointer definition of ContainerMapper
    KRATOS_CLASS_POINTER_DEFINITION(ExplicitFilterUtils);

    ///@}
    ///@name Nanoflann KD Tree type definitions
    ///@{

    using PositionAdapter = NanoFlannSingleContainerPositionAdapter<TContainerType>;

    using ResultVectorType = typename PositionAdapter::ResultVectorType;

    using PointerVectorType = typename PositionAdapter::PointerVectorType;

    using DistanceMetricType = typename nanoflann::metric_L2_Simple::traits<double, PositionAdapter>::distance_t;

    using KDTreeIndexType = nanoflann::KDTreeSingleIndexAdaptor<DistanceMetricType, PositionAdapter, 3>;

    using KDTreeThreadLocalStorage = NanoFlannKDTreeThreadLocalStorage<PointerVectorType>;

    ///@}
    ///@name LifeCycle
    ///@{

    ExplicitFilterUtils(
        const ModelPart& rModelPart,
        const std::string& rKernelFunctionType,
        const IndexType MaxLeafSize = 10,
        const IndexType EchoLevel = 0,
        const bool NodeCloudMesh = false,
        const bool StoreFilterMatrix = false);

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

    IndexType mLeafMaxSize;

    IndexType mEchoLevel;

    std::unique_ptr<PositionAdapter> mpAdapter;

    std::unique_ptr<KDTreeIndexType> mpKDTreeIndex;

    bool mNodeCloudMesh;

    bool mStoreFilteringMatrix;

    // The storage of the filtering matrix mFilteringMatrix.
    //      number of rows of first dense vector -> number of items in the container of interest in mrModelPart
    //      number of rows of second dense vector -> number of neighbours for each item of each row of the first dense vector
    //      std::get<0>(unsigned int) -> neighbour index
    //      std::get<1>(DenseVector)  -> filter matrix coefficient corresponding to the component in stride.
    DenseVector<DenseVector<std::pair<unsigned int, DenseVector<double>>>> mFilteringMatrix;

    ///@}
    ///@name Private operations
    ///@{

    void CheckField(const TensorAdaptor<double>& rTensorAdaptor) const;

    template<class TMeshDependencyType>
    void ComputeForwardFilteringMatrix();

    template<class TMeshDependencyType>
    void GenericGetIntegrationWeights(TensorAdaptor<double>& rTensorAdapto) const;

    template<class TMeshDependencyType, bool TUseFilterMatrix = false>
    TensorAdaptor<double>::Pointer GenericForwardFilterField(const TensorAdaptor<double>& rTensorAdaptor) const;

    template<class TMeshDependencyType, bool TUseFilterMatrix = false>
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