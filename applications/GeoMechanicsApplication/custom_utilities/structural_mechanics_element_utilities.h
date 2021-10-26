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
//                   Vahid Galavi

#if !defined( KRATOS_GEO_STRUCTURAL_MECHANICS_ELEMENT_UTILITIES_H_INCLUDED )
#define  KRATOS_GEO_STRUCTURAL_MECHANICS_ELEMENT_UTILITIES_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/element.h"

namespace Kratos {
class GeoStructuralMechanicsElementUtilities {

public:
/**
 * @brief Method to specify if the lumped or the consistent mass-matrix should be computed
 * @param rProperites The Properties where it is specified
 * @param rCurrentProcessInfo The ProcessInfo where it is specified
 * @return whether to compute the lumped mass-matrix
 */
static bool ComputeLumpedMassMatrix(
    const Properties& rProperites,
    const ProcessInfo& rCurrentProcessInfo);

/**
 * @brief Method to specify if rayligh-damping is specified
 * @param rProperites The Properties where it is specified
 * @param rCurrentProcessInfo The ProcessInfo where it is specified
 * @return whether rayleigh-damping was specified
 */
static bool HasRayleighDamping(
    const Properties& rProperites,
    const ProcessInfo& rCurrentProcessInfo);

/**
 * @brief Method to get the rayleigh-alpha parameter
 * @param rProperites The Properties where it is specified
 * @param rCurrentProcessInfo The ProcessInfo where it is specified
 * @return rayleigh-alpha
 */
static double GetRayleighAlpha(
    const Properties& rProperites,
    const ProcessInfo& rCurrentProcessInfo);

/**
 * @brief Method to get the rayleigh-beta parameter
 * @param rProperites The Properties where it is specified
 * @param rCurrentProcessInfo The ProcessInfo where it is specified
 * @return rayleigh-beta
 */
static double GetRayleighBeta(
    const Properties& rProperites,
    const ProcessInfo& rCurrentProcessInfo);

/**
 * @brief Method to returns the density to be consider for the mass-matrix computation
 * @param rElement The Element for which the mass-matrix should be computed
 * @return The density after apply the mass factor to the element
 */
static double GetDensityForMassMatrixComputation(const Element& rElement);

/**
 * @brief Method to calculate the rayleigh damping-matrix
 * @param rElement The Element for which the damping-matrix should be computed
 * @param rDampingMatrix The damping-matrix of the element
 * @param rCurrentProcessInfo The ProcessInfo where it is specified
 * @param MatrixSize The size of the damping-matrix
 */
static void CalculateRayleighDampingMatrix(
    Element& rElement,
    Element::MatrixType& rDampingMatrix,
    const ProcessInfo& rCurrentProcessInfo,
    const std::size_t MatrixSize);

/**
 * @brief This function calculates the reference length for 2D2N elements
 * @param rElement The Element for which the reference length should be computed
 * @return reference length
 */
static double CalculateReferenceLength2D2N(const Element& rElement);

/**
 * @brief This function calculates the current length for 2D2N elements
 * @param rElement The Element for which the current length should be computed
 * @return current length
 */
static double CalculateCurrentLength2D2N(const Element& rElement);

/**
 * @brief This function calculates the reference length for 3D2N elements
 * @param rElement The Element for which the reference length should be computed
 * @return reference length
 */
static double CalculateReferenceLength3D2N(const Element& rElement);

/**
 * @brief This function calculates the current length for 3D2N elements
 * @param rElement The Element for which the current length should be computed
 * @return current length
 */
static double CalculateCurrentLength3D2N(const Element& rElement);

}; // class StructuralMechanicsElementUtilities.
}  // namespace Kratos.

#endif // KRATOS_GEO_STRUCTURAL_MECHANICS_ELEMENT_UTILITIES_H_INCLUDED  defined
