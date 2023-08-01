// KRATOS    ______            __             __  _____ __                  __                   __
//          / ____/___  ____  / /_____ ______/ /_/ ___// /________  _______/ /___  ___________ _/ /
//         / /   / __ \/ __ \/ __/ __ `/ ___/ __/\__ \/ __/ ___/ / / / ___/ __/ / / / ___/ __ `/ / 
//        / /___/ /_/ / / / / /_/ /_/ / /__/ /_ ___/ / /_/ /  / /_/ / /__/ /_/ /_/ / /  / /_/ / /  
//        \____/\____/_/ /_/\__/\__,_/\___/\__//____/\__/_/   \__,_/\___/\__/\__,_/_/   \__,_/_/  MECHANICS
//
//  License:         BSD License
//                   license: ContactStructuralMechanicsApplication/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

#pragma once

// System includes

// External includes

// Project includes
#include "includes/model_part.h"

namespace Kratos
{
///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

///@}
///@name  Enum's
///@{

///@}
///@name  Functions
///@{

/**
 * @namespace ActiveSetUtilities
 * @ingroup ContactStructuralMechanicsApplication
 * @brief This namespace includes some utilities used for contact active set computations
 * @author Vicente Mataix Ferrandiz
 */
namespace ActiveSetUtilities
{
    ///@name Type Definitions
    ///@{

    // Some geometrical definitions
    typedef Node                                              NodeType;
    typedef Point::CoordinatesArrayType              CoordinatesArrayType;

    /// Definition of geometries
    typedef Geometry<NodeType>                               GeometryType;

    /// Index type definition
    typedef std::size_t                                         IndexType;

    /// Size type definition
    typedef std::size_t                                          SizeType;

    ///@}
    ///@name  Functions
    ///@{

    /**
     * @brief This function computes the active set for penalty frictionless cases
     * @param rThisModelPart The modelpart to compute
     */
    std::size_t KRATOS_API(CONTACT_STRUCTURAL_MECHANICS_APPLICATION) ComputePenaltyFrictionlessActiveSet(ModelPart& rModelPart);

    /**
     * @brief This function computes the active set for penalty frictional cases
     * @param rThisModelPart The modelpart to compute
     * @param PureSlip If we are considering pure slip case
     * @param EchoLevel The echo level considered
     */
    array_1d<std::size_t, 2> KRATOS_API(CONTACT_STRUCTURAL_MECHANICS_APPLICATION) ComputePenaltyFrictionalActiveSet(
        ModelPart& rModelPart,
        const bool PureSlip = false,
        const SizeType EchoLevel = 0
        );

    /**
     * @brief This function computes the active set for penalty frictionless cases
     * @param rThisModelPart The modelpart to compute
     */
    std::size_t KRATOS_API(CONTACT_STRUCTURAL_MECHANICS_APPLICATION) ComputeALMFrictionlessActiveSet(ModelPart& rModelPart);

    /**
     * @brief This function computes the active set for penalty frictionless cases
     * @param rThisModelPart The modelpart to compute
     */
    std::size_t KRATOS_API(CONTACT_STRUCTURAL_MECHANICS_APPLICATION) ComputeALMFrictionlessComponentsActiveSet(ModelPart& rModelPart);

    /**
     * @brief This function computes the active set for penalty frictional cases
     * @param rThisModelPart The modelpart to compute
     * @param PureSlip If we are considering pure slip case
     * @param EchoLevel The echo level considered
     */
    array_1d<std::size_t, 2> KRATOS_API(CONTACT_STRUCTURAL_MECHANICS_APPLICATION) ComputeALMFrictionalActiveSet(
        ModelPart& rModelPart,
        const bool PureSlip = false,
        const SizeType EchoLevel = 0
        );

};// namespace ActiveSetUtilities

} // namespace Kratos