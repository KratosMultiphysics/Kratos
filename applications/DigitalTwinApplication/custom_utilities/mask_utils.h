//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: DigitalTwinApplication/license.txt
//
//  Main authors:    Suneth Wranakulasuriya
//

#pragma once

// System includes

// External includes

// Project includes
#include "expression/container_expression.h"

// Application includes

namespace Kratos {
///@name Kratos Classes
///@{

class MaskUtils
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
     * @brief Get the substraction of two masks.
     *
     * @tparam TContainerType                           Container type.
     * @param rMask1                                    Mask 1 expression
     * @param rMask2                                    Mask 2 expression
     * @param RequiredMinimumRedundancy                 Required minimum redundancy.
     * @return ContainerExpression<TContainerType>      Substraction mask which substracts rMask2 from rMask1.
     */
    template<class TContainerType>
    static ContainerExpression<TContainerType> Substract(
        const ContainerExpression<TContainerType>& rMask1,
        const ContainerExpression<TContainerType>& rMask2,
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