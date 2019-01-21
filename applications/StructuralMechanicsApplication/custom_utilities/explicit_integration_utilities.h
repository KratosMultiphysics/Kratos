// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Klaus B Sautter
//                   Vicente Mataix Ferrandiz
//

#if !defined(KRATOS_EXPLICIT_LAW_UTILITIES)
#define KRATOS_EXPLICIT_LAW_UTILITIES

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

///@}
///@name Kratos Classes
///@{
/**
 * @namespace ExplicitIntegrationUtilities
 * @ingroup StructuralMechanicsApplication
 * @brief This namespace includes several utilities necessaries for the computation of the explicit integration
 * @author Klaus B Sautter
 * @author Vicente Mataix Ferrandiz
 */
namespace ExplicitIntegrationUtilities
{
    /// The size type definition
    typedef std::size_t SizeType;

    /// The index type definition
    typedef std::size_t IndexType;

    /// The arrays of elements and nodes
    typedef ModelPart::ElementsContainerType ElementsArrayType;
    typedef ModelPart::NodesContainerType NodesArrayType;

    /**
     * @brief This method computes the necessry delta time to avoid numerical instabilities
     * @param rModelPart The model of the problem to solve
     * @param PredictionLevel The prediction level
     * @param MaximumDeltaTime The maximum delta time to be considered
     * @param SafetyFactor The factor to not consider exactly the theoretical value
     * @return The critical delta time
     */
    double CalculateDeltaTime(
        ModelPart& rModelPart,
        const double PredictionLevel = 2.0,
        const double MaximumDeltaTime = 1.0e-3,
        const double SafetyFactor = 0.5
        );

}; // namespace ExplicitIntegrationUtilities
}  // namespace Kratos
#endif /* KRATOS_EXPLICIT_LAW_UTILITIES defined */
