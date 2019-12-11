// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//


#if !defined(KRATOS_SELF_CONTACT_UTILITIES)
#define KRATOS_SELF_CONTACT_UTILITIES

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
 * @namespace SelfContactUtilities
 * @ingroup ContactStructuralMechanicsApplication
 * @brief This namespace includes several utilities necessaries for the computation of self-contact
 * @author Vicente Mataix Ferrandiz
 */
namespace SelfContactUtilities
{
    /**
     * @brief This method computes the pairing for self-contact
     * @param rModelPart Reference to the ModelPart containing the problem
     * @param EchoLevel The verbosity level
     */
    void KRATOS_API(CONTACT_STRUCTURAL_MECHANICS_APPLICATION) ComputeSelfContactPairing(
        ModelPart& rModelPart,
        const std::size_t EchoLevel = 0
        );

    /**
     * @brief This method assign all the potential pairs (brute-force)
     * @param rModelPart Reference to the ModelPart containing the problem
     */
    void KRATOS_API(CONTACT_STRUCTURAL_MECHANICS_APPLICATION) FullAssignmentOfPairs(ModelPart& rModelPart);

    /**
     * @brief This method is used in case of not predefined master/slave we assign the master/slave nodes and conditions
     * @param rModelPart The model part to assign the flags
     */
    void KRATOS_API(CONTACT_STRUCTURAL_MECHANICS_APPLICATION) NotPredefinedMasterSlave(ModelPart& rModelPart);

}; // namespace SelfContactUtilities
}  // namespace Kratos
#endif /* KRATOS_SELF_CONTACT_UTILITIES defined */
