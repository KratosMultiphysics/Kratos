//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: SystemIdentificationApplication/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//

#pragma once

// System includes
#include <vector>

// External includes

// Project includes
#include "includes/define.h"
#include "tensor_adaptors/tensor_adaptor.h"

// Application includes

namespace Kratos {
///@name Kratos Classes
///@{

class KRATOS_API(SYSTEM_IDENTIFICATION_APPLICATION) MaskUtils
{
public:
    ///@name Type definitions
    ///@{

    using IndexType = std::size_t;

    ///@}
    ///@name Public static operations
    ///@{

    /**
     * @brief Get the Mask Size satisfying the minimum redundancy specified.
     *
     * This method returns number of entities of which the value is greater
     * than or equal to RequiredMinimumRedundancy.
     *
     * @param rMask                                     Input mask.
     * @param RequiredMinimumRedundancy                 Required minimum redundancy.
     * @return IndexType                                Number of entities above or equal to minimum redundancy.
     */
    static IndexType GetMaskSize(
        const TensorAdaptor<double>& rMask,
        const IndexType RequiredMinimumRedundancy = 1);

    /**
     * @brief Get the Mask for the given scalar tensor adaptor
     *
     * The mask is computed such that the inner product between
     * the mask and the scalar tensor adaptor is maximized.
     *
     * @param rScalarTensorAdaptor                         Input scalar tensor adaptor
     * @return TensorAdaptor<double>::Pointer              TensorAdaptor with values either 0 or 1.
     */
    static TensorAdaptor<double>::Pointer GetMask(
        const TensorAdaptor<double>& rScalarTensorAdaptor);

    /**
     * @brief Get the Mask Threshold
     * @details The mask threshold is computed such that the given scalar tensor adaptor
     *          is aligned to the mask with minimum cosine distance
     *
     * @param rScalarTensorAdaptor                         Input scalar tensor adaptor
     * @return double                                      Mask threshold
     */
    static double GetMaskThreshold(const TensorAdaptor<double>& rScalarTensorAdaptor);

    /**
     * @brief Get the Mask for the given scalar tensor adaptor using the threshold.
     *
     * The mask is computed such that all the entity values in the scalar tensor adaptor
     * which is higher than the threshold will have 1.0, others will have 0.0.
     *
     * @param rScalarTensorAdaptor                         Input scalar tensor adaptor
     * @param Threshold                                    Threshold for mask 0 or 1 detection.
     * @return TensorAdaptor<double>::Pointer              TensorAdaptor with values either 0 or 1.
     */
    static TensorAdaptor<double>::Pointer GetMask(
        const TensorAdaptor<double>& rScalarTensorAdaptor,
        const double Threshold);

    /**
     * @brief Get the union of two masks.
     *
     * This method returns a mask with the value of RequiredMinimumRedundancy for entities
     * which has either in rMask1 or rMask2 a value equal or grater than RequiredMinimumRedundancy.
     *
     * @param rMask1                                    Mask 1 tensor adaptor
     * @param rMask2                                    Mask 2 tensor adaptor
     * @param RequiredMinimumRedundancy                 Required minimum redundancy.
     * @return TensorAdaptor<double>::Pointer              Union mask
     */
    static TensorAdaptor<double>::Pointer Union(
        const TensorAdaptor<double>& rMask1,
        const TensorAdaptor<double>& rMask2,
        const IndexType RequiredMinimumRedundancy = 1);

    /**
     * @brief Get the intersection of two masks.
     *
     * This method returns a mask with the value of RequiredMinimumRedundancy for entities
     * which has equal or greater than RequiredMinimumRedundancy in both rMask1 and rMask2.
     *
     * @param rMask1                                    Mask 1 tensor adaptor
     * @param rMask2                                    Mask 2 tensor adaptor
     * @param RequiredMinimumRedundancy                 Required minimum redundancy.
     * @return TensorAdaptor<double>::Pointer              Intersection mask
     */
    static TensorAdaptor<double>::Pointer Intersect(
        const TensorAdaptor<double>& rMask1,
        const TensorAdaptor<double>& rMask2,
        const IndexType RequiredMinimumRedundancy = 1);

    /**
     * @brief Get the subtraction of two masks.
     *
     * @param rMask1                                    Mask 1 tensor adaptor
     * @param rMask2                                    Mask 2 tensor adaptor
     * @param RequiredMinimumRedundancy                 Required minimum redundancy.
     * @return TensorAdaptor<double>::Pointer           Subtraction mask which subtracts rMask2 from rMask1.
     */
    static TensorAdaptor<double>::Pointer Subtract(
        const TensorAdaptor<double>& rMask1,
        const TensorAdaptor<double>& rMask2,
        const IndexType RequiredMinimumRedundancy = 1);


    /**
     * @brief Scale the scalar tensor adaptor with a mask.
     *
     * This method scales the rScalarTensorAdaptor with a mask of [0, 1] where mask will
     * have 1 if that corresponding entity in rMask has a value equal or greater than the
     * RequiredMinimumRedundancy otherwise 0.
     *
     * @param rScalarTensorAdaptor                      Input scalar tensor adaptor
     * @param rMask                                     Mask to be used.
     * @param RequiredMinimumRedundancy                 Required minimum redundancy.
     * @return TensorAdaptor<double>::Pointer           Scaled scalar tensor adaptor with the mask.
     */
    static TensorAdaptor<double>::Pointer Scale(
        const TensorAdaptor<double>& rScalarTensorAdaptor,
        const TensorAdaptor<double>& rMask,
        const IndexType RequiredMinimumRedundancy = 1);

    /**
     * @brief Returns list of indices list each list for each cluster
     *
     * This method returns a list of tuples having indices list and a mask. Each indices list is a cluster
     * The mask for that specific cluster is stored in the 2nd element of the tuple.
     *
     * @param rMasksList                                List of masks.
     * @param RequiredMinimumRedundancy                 Required minimum redundancy.
     * @return std::vector<std::tuple<std::vector<IndexType>, TensorAdaptor<double>::Pointer>> Tuple having list of mask indices and masks for each cluster.
     */
    static std::vector<std::tuple<std::vector<IndexType>, TensorAdaptor<double>::Pointer>> ClusterMasks(
        const std::vector<TensorAdaptor<double>::Pointer>& rMasksList,
        const IndexType RequiredMinimumRedundancy = 1);

    /**
     * @brief Get the Masks Dividing Reference Mask.
     *
     * This method returns indices of masks which can divide the given rReferenceMask.
     *
     * @param rReferenceMask                            Reference mask.
     * @param rMasksList                                List of masks.
     * @param RequiredMinimumRedundancy                 Required minimum redundancy.
     * @return std::vector<IndexType>                   List of indices.
     */
    static std::vector<IndexType> GetMasksDividingReferenceMask(
        const TensorAdaptor<double>& rReferenceMask,
        const std::vector<TensorAdaptor<double>::Pointer>& rMasksList,
        const IndexType RequiredMinimumRedundancy = 1);

    ///@}

private:
    ///@name Static private operations
    ///@{

    static void CheckCompatibility(
        const TensorAdaptor<double>& rMask1,
        const TensorAdaptor<double>& rMask2);

    ///@}
};
} // namespace Kratos