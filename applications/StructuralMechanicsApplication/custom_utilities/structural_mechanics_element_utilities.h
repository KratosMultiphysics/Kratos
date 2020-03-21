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
//                   Vicente Mataix Ferrandiz
//                   Riccardo Rossi
//                   Ruben Zorrilla
//

#if !defined( KRATOS_STRUCTURAL_MECHANICS_ELEMENT_UTILITIES_H_INCLUDED )
#define  KRATOS_STRUCTURAL_MECHANICS_ELEMENT_UTILITIES_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/element.h"

namespace Kratos {
/**
 * @namespace ElementUtilities
 * @ingroup StructuralMechanicsApplication
 * @brief This class includes several utilities necessaries for the computation of the different elements
 * @author Philipp Bucher
 * @author Vicente Mataix Ferrandiz
 * @author Riccardo Rossi
 * @author Ruben Zorrilla
 */
namespace StructuralMechanicsElementUtilities {

/// The size type definition
typedef std::size_t SizeType;

/// The index type definition
typedef std::size_t IndexType;

/// Node type definition
typedef Node<3> NodeType;

/// Geometry definitions
typedef Geometry<NodeType> GeometryType;

/**
 * @brief This method performs commons checks on the solid elements
 * @param rElement Reference to the element
 * @param rCurrentProcessInfo The current process info instance
 * @param rConstitutiveLaws The vector containing CL
 */
int SolidElementCheck(
    const Element& rElement,
    const ProcessInfo& rCurrentProcessInfo,
    std::vector<ConstitutiveLaw::Pointer>& rConstitutiveLaws
    );

/**
 * @brief This method computes the deformation gradient F (for small deformation solid elements)
 * @param rElement Reference to the element
 * @param rStrainTensor The strain tensor
 * @param rF The deformation gradient F
 * @tparam TVectorType The vector type
 * @tparam TMatrixType The matrix type
 */
template<class TVectorType, class TMatrixType>
void ComputeEquivalentF(
    const Element& rElement,
    const TVectorType& rStrainTensor,
    TMatrixType& rF
    )
{
    const auto& r_geometry = rElement.GetGeometry();
    const SizeType dimension = r_geometry.WorkingSpaceDimension();

    if(dimension == 2) {
        rF(0,0) = 1.0+rStrainTensor(0);
        rF(0,1) = 0.5*rStrainTensor(2);
        rF(1,0) = 0.5*rStrainTensor(2);
        rF(1,1) = 1.0+rStrainTensor(1);
    } else {
        rF(0,0) = 1.0+rStrainTensor(0);
        rF(0,1) = 0.5*rStrainTensor(3);
        rF(0,2) = 0.5*rStrainTensor(5);
        rF(1,0) = 0.5*rStrainTensor(3);
        rF(1,1) = 1.0+rStrainTensor(1);
        rF(1,2) = 0.5*rStrainTensor(4);
        rF(2,0) = 0.5*rStrainTensor(5);
        rF(2,1) = 0.5*rStrainTensor(4);
        rF(2,2) = 1.0+rStrainTensor(2);
    }
}

/**
 * @brief This method computes the deformation tensor B (for small deformation solid elements)
 * @param rElement Reference to the element
 * @param rDN_DX The shape function derivatives
 * @param rB The deformation tensor B
 * @tparam TMatrixType1 The first matrix type
 * @tparam TMatrixType2 The second matrix type
 */
template<class TMatrixType1, class TMatrixType2>
void CalculateB(
    const Element& rElement,
    const TMatrixType1& rDN_DX,
    TMatrixType2& rB
    )
{
    const auto& r_geometry = rElement.GetGeometry();
    const SizeType number_of_nodes = r_geometry.PointsNumber();
    const SizeType dimension = r_geometry.WorkingSpaceDimension();

    rB.clear();

    if(dimension == 2) {
        for ( IndexType i = 0; i < number_of_nodes; ++i ) {
            const IndexType initial_index = i*2;
            rB(0, initial_index    ) = rDN_DX(i, 0);
            rB(1, initial_index + 1) = rDN_DX(i, 1);
            rB(2, initial_index    ) = rDN_DX(i, 1);
            rB(2, initial_index + 1) = rDN_DX(i, 0);
        }
    } else if(dimension == 3) {
        for ( IndexType i = 0; i < number_of_nodes; ++i ) {
            const IndexType initial_index = i*3;
            rB(0, initial_index    ) = rDN_DX(i, 0);
            rB(1, initial_index + 1) = rDN_DX(i, 1);
            rB(2, initial_index + 2) = rDN_DX(i, 2);
            rB(3, initial_index    ) = rDN_DX(i, 1);
            rB(3, initial_index + 1) = rDN_DX(i, 0);
            rB(4, initial_index + 1) = rDN_DX(i, 2);
            rB(4, initial_index + 2) = rDN_DX(i, 1);
            rB(5, initial_index    ) = rDN_DX(i, 2);
            rB(5, initial_index + 2) = rDN_DX(i, 0);
        }
    }
}

/**
 * @brief This method returns the computed the computed body force
 * @param rElement Reference to the element
 * @param rIntegrationPoints The integrations points
 * @param PointNumber The integration point number
 */
array_1d<double, 3> GetBodyForce(
    const Element& rElement,
    const GeometryType::IntegrationPointsArrayType& rIntegrationPoints,
    const IndexType PointNumber
    );

/**
 * @brief Method to specify if the lumped or the consistent mass-matrix should be computed
 * @param rProperties The Properties where it is specified
 * @param rCurrentProcessInfo The ProcessInfo where it is specified
 * @return whether to compute the lumped mass-matrix
 */
bool ComputeLumpedMassMatrix(
    const Properties& rProperties,
    const ProcessInfo& rCurrentProcessInfo);

/**
 * @brief Method to specify if rayligh-damping is specified
 * @param rProperties The Properties where it is specified
 * @param rCurrentProcessInfo The ProcessInfo where it is specified
 * @return whether rayleigh-damping was specified
 */
bool HasRayleighDamping(
    const Properties& rProperties,
    const ProcessInfo& rCurrentProcessInfo);

/**
 * @brief Method to get the rayleigh-alpha parameter
 * @param rProperties The Properties where it is specified
 * @param rCurrentProcessInfo The ProcessInfo where it is specified
 * @return rayleigh-alpha
 */
double GetRayleighAlpha(
    const Properties& rProperties,
    const ProcessInfo& rCurrentProcessInfo);

/**
 * @brief Method to get the rayleigh-beta parameter
 * @param rProperties The Properties where it is specified
 * @param rCurrentProcessInfo The ProcessInfo where it is specified
 * @return rayleigh-beta
 */
double GetRayleighBeta(
    const Properties& rProperties,
    const ProcessInfo& rCurrentProcessInfo);

/**
 * @brief Method to returns the density to be consider for the mass-matrix computation
 * @param rElement The Element for which the mass-matrix should be computed
 * @return The density after apply the mass factor to the element
 */
double GetDensityForMassMatrixComputation(const Element& rElement);

/**
 * @brief Method to calculate the rayleigh damping-matrix
 * @param rElement The Element for which the damping-matrix should be computed
 * @param rDampingMatrix The damping-matrix of the element
 * @param rCurrentProcessInfo The ProcessInfo where it is specified
 * @param MatrixSize The size of the damping-matrix
 */
void CalculateRayleighDampingMatrix(
    Element& rElement,
    Element::MatrixType& rDampingMatrix,
    /*const*/ ProcessInfo& rCurrentProcessInfo,
    const std::size_t MatrixSize);

/**
 * @brief This function calculates the reference length for 2D2N elements
 * @param rElement The Element for which the reference length should be computed
 * @return reference length
 */
double CalculateReferenceLength2D2N(const Element& rElement);

/**
 * @brief This function calculates the current length for 2D2N elements
 * @param rElement The Element for which the current length should be computed
 * @return current length
 */
double CalculateCurrentLength2D2N(const Element& rElement);

/**
 * @brief This function calculates the reference length for 3D2N elements
 * @param rElement The Element for which the reference length should be computed
 * @return reference length
 */
double CalculateReferenceLength3D2N(const Element& rElement);

/**
 * @brief This function calculates the current length for 3D2N elements
 * @param rElement The Element for which the current length should be computed
 * @return current length
 */
double CalculateCurrentLength3D2N(const Element& rElement);

} // namespace StructuralMechanicsElementUtilities.
}  // namespace Kratos.

#endif // KRATOS_STRUCTURAL_MECHANICS_ELEMENT_UTILITIES_H_INCLUDED  defined
