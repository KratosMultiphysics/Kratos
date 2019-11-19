// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//                   Riccardo Rossi
//                   Ruben Zorrilla
//

#if !defined(KRATOS_ELEMENT_UTILITIES_INCLUDE_H)
#define KRATOS_ELEMENT_UTILITIES_INCLUDE_H

// System includes

// External includes

// Project includes
#include "includes/element.h"

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
 * @class ElementUtilities
 * @ingroup StructuralMechanicsApplication
 * @brief This class includes several utilities necessaries for the computation of the different elements
 * @details The methods are static, so it can be called without constructing the class
 * @author Vicente Mataix Ferrandiz
 * @author Riccardo Rossi
 * @author Ruben Zorrilla
 */
class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) ElementUtilities
{
  public:
    ///@name Type definitions
    ///@{

    /// The size type definition
    typedef std::size_t SizeType;

    /// The index type definition
    typedef std::size_t IndexType;

    /// Node type definition
    typedef Node<3> NodeType;

    /// Geometry definitions
    typedef Geometry<NodeType> GeometryType;

    /// The zero tolerance
    static constexpr double tolerance = std::numeric_limits<double>::epsilon();

    ///@}
    ///@name  Enum's
    ///@{

    ///@}
    ///@name Life Cycle
    ///@{

    ///@}
    ///@name Operations
    ///@{

    /**
     * @brief This method performs commons checks on the elements
     * @param pElement Pointer to the element
     * @param rCurrentProcessInfo The current process info instance
     */
    static int BaseElementCheck(
        const Element* pElement,
        const ProcessInfo& rCurrentProcessInfo
        );

    /**
     * @brief This method computes the deformation gradient F (for small deformation solid elements)
     * @param pElement Pointer to the element
     * @param rF The deformation gradient F
     * @param rStrainTensor The strain tensor
     */
    static void ComputeEquivalentF(
        const Element* pElement,
        Matrix& rF,
        const Vector& rStrainTensor
        );

    /**
     * @brief This method computes the deformation tensor B (for small deformation solid elements)
     * @param pElement Pointer to the element
     * @param rB The deformation tensor B
     * @param rDN_DX The shape function derivatives
     */
    static void CalculateB(
        const Element* pElement,
        Matrix& rB,
        const Matrix& rDN_DX
        );

    /**
     * @brief This method returns the computed the computed body force
     * @param pElement Pointer to the element
     * @param rIntegrationPoints The integrations points
     * @param PointNumber The integration point number
     */
    static array_1d<double, 3> GetBodyForce(
        const Element* pElement,
        const GeometryType::IntegrationPointsArrayType& rIntegrationPoints,
        const IndexType PointNumber
        );

  private:

}; // class ElementUtilities
} // namespace Kratos
#endif /* KRATOS_ELEMENT_UTILITIES_INCLUDE_H defined */
