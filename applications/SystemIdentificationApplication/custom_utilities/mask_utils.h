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
#include "expression/container_expression.h"

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
     * @tparam TContainerType
     * @param rMask                                     Input mask.
     * @param RequiredMinimumRedundancy                 Required minimum redundancy.
     * @return IndexType                                Number of entities above or equal to minimum redundancy.
     */
    template<class TContainerType>
    static IndexType GetMaskSize(
        const ContainerExpression<TContainerType>& rMask,
        const IndexType RequiredMinimumRedundancy = 1);

    /**
     * @brief Get the Mask for the given scalar expression
     *
     * The mask is computed such that the inner product between
     * the mask and the scalar expression is maximized.
     *
     * @tparam TContainerType                           Container type.
     * @param rScalarExpression                         Input scalar expression.
     * @return ContainerExpression<TContainerType>      Expression with values either 0 or 1.
     */
    template<class TContainerType>
    static ContainerExpression<TContainerType> GetMask(
        const ContainerExpression<TContainerType>& rScalarExpression);

    /**
     * @brief Get the Mask Threshold
     * @details The mask threshold is computed such that the given scalar expression
     *          is aligned to the mask with minimum cosine distance
     *
     * @tparam TContainerType                           Container type.
     * @param rScalarExpression                         Input scalar expression.
     * @return double                                   Mask threshold
     */
    template<class TContainerType>
    static double GetMaskThreshold(const ContainerExpression<TContainerType>& rScalarExpression);

    /**
     * @brief Get the Mask for the given scalar expression using the threshold.
     *
     * The mask is computed such that all the entity values in the scalar expression
     * which is higher than the threshold will have 1.0, others will have 0.0.
     *
     * @tparam TContainerType                           Container type.
     * @param rScalarExpression                         Input scalar expression.
     * @param Threshold                                 Threshold for mask 0 or 1 detection.
     * @return ContainerExpression<TContainerType>      Expression with values either 0 or 1.
     */
    template<class TContainerType>
    static ContainerExpression<TContainerType> GetMask(
        const ContainerExpression<TContainerType>& rScalarExpression,
        const double Threshold);

    /**
     * @brief Get the union of two masks.
     *
     * This method returns a mask with the value of RequiredMinimumRedundancy for entities
     * which has either in rMask1 or rMask2 a value equal or grater than RequiredMinimumRedundancy.
     *
     * @tparam TContainerType                           Container type.
     * @param rMask1                                    Mask 1 expression
     * @param rMask2                                    Mask 2 expression
     * @param RequiredMinimumRedundancy                 Required minimum redundancy.
     * @return ContainerExpression<TContainerType>      Union mask
     */
    template<class TContainerType>
    static ContainerExpression<TContainerType> Union(
        const ContainerExpression<TContainerType>& rMask1,
        const ContainerExpression<TContainerType>& rMask2,
        const IndexType RequiredMinimumRedundancy = 1);

    /**
     * @brief Get the intersection of two masks.
     *
     * This method returns a mask with the value of RequiredMinimumRedundancy for entities
     * which has equal or greater than RequiredMinimumRedundancy in both rMask1 and rMask2.
     *
     * @tparam TContainerType                           Container type.
     * @param rMask1                                    Mask 1 expression
     * @param rMask2                                    Mask 2 expression
     * @param RequiredMinimumRedundancy                 Required minimum redundancy.
     * @return ContainerExpression<TContainerType>      Intersection mask
     */
    template<class TContainerType>
    static ContainerExpression<TContainerType> Intersect(
        const ContainerExpression<TContainerType>& rMask1,
        const ContainerExpression<TContainerType>& rMask2,
        const IndexType RequiredMinimumRedundancy = 1);

    /**
     * @brief Get the subtraction of two masks.
     *
     * @tparam TContainerType                           Container type.
     * @param rMask1                                    Mask 1 expression
     * @param rMask2                                    Mask 2 expression
     * @param RequiredMinimumRedundancy                 Required minimum redundancy.
     * @return ContainerExpression<TContainerType>      Subtraction mask which subtracts rMask2 from rMask1.
     */
    template<class TContainerType>
    static ContainerExpression<TContainerType> Subtract(
        const ContainerExpression<TContainerType>& rMask1,
        const ContainerExpression<TContainerType>& rMask2,
        const IndexType RequiredMinimumRedundancy = 1);


    /**
     * @brief Scale the scalar expression with a mask.
     *
     * This method scales the rScalarExpression with a mask of [0, 1] where mask will
     * have 1 if that corresponding entity in rMask has a value equal or greater than the
     * RequiredMinimumRedundancy otherwise 0.
     *
     * @tparam TContainerType                           Container type.
     * @param rScalarExpression                         Input scalar expression.
     * @param rMask                                     Mask to be used.
     * @param RequiredMinimumRedundancy                 Required minimum redundancy.
     * @return ContainerExpression<TContainerType>      Scaled scalar expression with the mask.
     */
    template<class TContainerType>
    static ContainerExpression<TContainerType> Scale(
        const ContainerExpression<TContainerType>& rScalarExpression,
        const ContainerExpression<TContainerType>& rMask,
        const IndexType RequiredMinimumRedundancy = 1);

    /**
     * @brief Returns list of indices list each list for each cluster
     *
     * This method returns a list of tuples having indices list and a mask. Each indices list is a cluster
     * The mask for that specific cluster is stored in the 2nd element of the tuple.
     *
     * @tparam TContainerType                           Container type.
     * @param rMasksList                                List of masks.
     * @param RequiredMinimumRedundancy                 Required minimum redundancy.
     * @return std::vector<std::tuple<std::vector<IndexType>, typename ContainerExpression<TContainerType>::Pointer>> Tuple having list of mask indices and masks for each cluster.
     */
    template<class TContainerType>
    static std::vector<std::tuple<std::vector<IndexType>, typename ContainerExpression<TContainerType>::Pointer>> ClusterMasks(
        const std::vector<ContainerExpression<TContainerType>>& rMasksList,
        const IndexType RequiredMinimumRedundancy = 1);

    /**
     * @brief Fills the given model part with the entities where the pClusterMask has minimum value of RequiredMinimumRedundancy.
     * @details This method will fill the given model part with the elements where the pClusterMask has the value for
     *          corresponding entity greater than or equal to the RequiredMinimumRedundancy.
     *
     *          In the case of TContainerType = ModelPart::ElementsContainerType or ModelPart::ConditionsContainerType,
     *          then corresponding nodes of those valid elements or conditions will also be added to the model part.
     *
     * @throws if the model part is not empty.
     * @tparam TContainerType       Container type of the cluster expression.
     * @param rModelPart            Model part to be filled with the cluster entities.
     * @param pClusterMask          Cluster mask.
     */
    template<class TContainerType>
    static void FillModelPartUsingClusterMask(
        ModelPart& rModelPart,
        const ContainerExpression<TContainerType>& rClusterMask,
        const IndexType RequiredMinimumRedundancy = 1);

    /**
     * @brief Get the Masks Dividing Reference Mask.
     *
     * This method returns indices of masks which can divide the given rReferenceMask.
     *
     * @tparam TContainerType                           Container type.
     * @param rReferenceMask                            Reference mask.
     * @param rMasksList                                List of masks.
     * @param RequiredMinimumRedundancy                 Required minimum redundancy.
     * @return std::vector<IndexType>                   List of indices.
     */
    template<class TContainerType>
    static std::vector<IndexType> GetMasksDividingReferenceMask(
        const ContainerExpression<TContainerType>& rReferenceMask,
        const std::vector<typename ContainerExpression<TContainerType>::Pointer>& rMasksList,
        const IndexType RequiredMinimumRedundancy = 1);

    ///@}

private:
    ///@name Static private operations
    ///@{

    template<class TContainerType>
    static void CheckCompatibility(
        const ContainerExpression<TContainerType>& rMask1,
        const ContainerExpression<TContainerType>& rMask2);

    ///@}
};
} // namespace Kratos