//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher
//

#if !defined( KRATOS_STRUCTURAL_MECHANICS_ELEMENT_UTILITIES_H_INCLUDED )
#define  KRATOS_STRUCTURAL_MECHANICS_ELEMENT_UTILITIES_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/element.h"

namespace Kratos {
namespace StructuralMechanicsElementUtilities {

/**
 * @brief Method to specify if the lumped or the consistent mass-matrix should be computed
 * @param rProperites The Properties where it is specified
 * @param rCurrentProcessInfo The ProcessInfo where it is specified
 * @return whether to compute the lumped mass-matrix
 */
bool ComputeLumpedMassMatrix(
    const Properties& rProperites,
    const ProcessInfo& rCurrentProcessInfo);

/**
 * @brief Method to specify if rayligh-damping is specified
 * @param rProperites The Properties where it is specified
 * @param rCurrentProcessInfo The ProcessInfo where it is specified
 * @return whether rayleigh-damping was specified
 */
bool HasRayleighDamping(
    const Properties& rProperites,
    const ProcessInfo& rCurrentProcessInfo);

/**
 * @brief Method to get the rayleigh-alpha parameter
 * @param rProperites The Properties where it is specified
 * @param rCurrentProcessInfo The ProcessInfo where it is specified
 * @return rayleigh-alpha
 */
double GetRayleighAlpha(
    const Properties& rProperites,
    const ProcessInfo& rCurrentProcessInfo);

/**
 * @brief Method to get the rayleigh-beta parameter
 * @param rProperites The Properties where it is specified
 * @param rCurrentProcessInfo The ProcessInfo where it is specified
 * @return rayleigh-beta
 */
double GetRayleighBeta(
    const Properties& rProperites,
    const ProcessInfo& rCurrentProcessInfo);

/**
 * @brief Method to claculate the rayleigh damping-matrix
 * @param rElement The Element for which the damping-matrix should be computed
 * @param rDampingMatrix The current ProcessInfo
 * @param rCurrentProcessInfo The ProcessInfo where it is specified
 * @param MatrixSize The size of the damping-matrix
 */
void CalculateRayleighDampingMatrix(
    Element& rElement,
    Element::MatrixType& rDampingMatrix,
    /*const*/ ProcessInfo& rCurrentProcessInfo,
    const std::size_t MatrixSize);

} // namespace StructuralMechanicsElementUtilities.
}  // namespace Kratos.

#endif // KRATOS_STRUCTURAL_MECHANICS_ELEMENT_UTILITIES_H_INCLUDED  defined
