// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:             BSD License
//                                       license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

#if !defined(KRATOS_ACTIVE_SET_UTILITIES)
#define KRATOS_ACTIVE_SET_UTILITIES

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
    typedef Node<3>                                              NodeType;
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

#endif /* KRATOS_ACTIVE_SET_UTILITIES defined */
