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

#if !defined(KRATOS_RAYLEIGH_DAMPING_COEFFICIENTS_UTILITIES)
#define KRATOS_RAYLEIGH_DAMPING_COEFFICIENTS_UTILITIES

// System includes

// External includes

// Project includes
#include "includes/kratos_parameters.h"

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
 * @namespace RayleighDampingCoefficientsUtilities
 * @ingroup StructuralMechanicsApplication
 * @brief This utility computes the two first eigen values of the system and estimates the alpha and beta Rayleigh damping coefficients
 * @details It uses the well stablished formulation: Wilson, E. L. (2004). Static and Dynamic Analysis of Structures (4th ed.). Berkeley, CA: Computers and Structures, Inc.
 * @author Vicente Mataix Ferrandiz
*/
namespace RayleighDampingCoefficientsUtilities
{
    /**
    * @brief This utility computes the two first eigen values of the system and estimates the alpha and beta Rayleigh damping coefficients
    * @details It uses the well stablished formulation: Wilson, E. L. (2004). Static and Dynamic Analysis of Structures (4th ed.). Berkeley, CA: Computers and Structures, Inc.
    * @param ThisParameters The configuration parameters
    */
    Vector KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) ComputeDampingCoefficients(Parameters ThisParameters);

} /// namespace RayleighDampingCoefficientsUtilities
} /// namespace Kratos
#endif /* KRATOS_RAYLEIGH_DAMPING_COEFFICIENTS_UTILITIES defined */
