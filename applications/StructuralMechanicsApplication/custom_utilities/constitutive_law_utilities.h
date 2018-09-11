// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:		 BSD License
//					 license: structural_mechanics_application/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//                   Alejandro Cornejo
//

#if !defined(KRATOS_CONSTITUTIVE_LAW_UTILITIES)
#define KRATOS_CONSTITUTIVE_LAW_UTILITIES

// System includes

// External includes

// Project includes

namespace Kratos
{
///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

    /// The size type definition
    typedef std::size_t SizeType;

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
 * @class ConstitutiveLawUtilities
 * @ingroup StructuralMechanicsApplication
 * @brief This class includes several utilities necessaries for the computation of the constitutive law
 * @details The methods are static, so it can be called without constructing the class
 * @tparam TVoigtSize The number of components on the Voigt notation
 * @author Vicente Mataix Ferrandiz
 * @author Alejandro Cornejo
 * @todo Adapt for 2D dimension
 */
template <SizeType TVoigtSize = 6>
class ConstitutiveLawUtilities
{
  public:
    ///@name Type definitions
    ///@{

    /// The index type definition
    typedef std::size_t IndexType;

    /// We define the dimension
    static constexpr SizeType Dimension = TVoigtSize == 6 ? 3 : 2;

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
     * @brief This method computes the first invariant from a given stress vector
     * @param StressVector The stress vector on Voigt notation
     * @param rI1 The first invariant
     * @todo Adapt for 2D dimension
     */
    static void CalculateI1Invariant(
        const Vector& StressVector,
        double& rI1
        );

    /**
     * @brief This method computes the second invariant from a given stress vector
     * @param StressVector The stress vector on Voigt notation
     * @param rI2 The second invariant
     * @todo Adapt for 2D dimension
     */
    static void CalculateI2Invariant(
        const Vector& StressVector,
        double& rI2
        );

    /**
     * @brief This method computes the third invariant from a given stress vector
     * @param StressVector The stress vector on Voigt notation
     * @param rI3 The third invariant
     * @todo Adapt for 2D dimension
     */
    static void CalculateI3Invariant(
        const Vector& StressVector,
        double& rI3
        );

    /**
     * @brief This method computes the second invariant of J
     * @param StressVector The stress vector on Voigt notation
     * @param I1 The first invariant
     * @param Deviator The deviator of the stress
     * @param rJ2 The second invariant of J
     * @todo Adapt for 2D dimension
     */
    static void CalculateJ2Invariant(
        const Vector& StressVector,
        const double I1,
        Vector& rDeviator,
        double& rJ2
        );

    /**
     * @brief This method computes the third invariant of J
     * @param Deviator The deviator of the stress
     * @param rJ3 The third invariant of J
     */
    static void CalculateJ3Invariant(
        const Vector& Deviator,
        double& rJ3
        );

    /**
     * @brief This method computes the first vector
     * @param FirstVector The first vector
     */
    static void CalculateFirstVector(Vector& FirstVector);

    /**
     * @brief This method computes the second vector
     * @param Deviator The deviator of the stress
     * @param J2 The resultant J2 stress
     * @param SecondVector The second vector
     */
    static void CalculateSecondVector(
        const Vector& Deviator,
        const double J2,
        Vector& SecondVector
        );

    /**
     * @brief This method computes the third vector
     * @param Deviator The deviator of the stress
     * @param J2 The resultant J2 stress
     * @param ThirdVector The third vector
     * @todo Adapt for 2D dimension
     */
    static void CalculateThirdVector(
        const Vector& Deviator,
        const double J2,
        Vector& ThirdVector
        );

    /**
     * @brief This method computes the lode angle
     * @param J2 The resultant J2 stress
     * @param J3 The resultant J3 stress
     * @param LodeAngle The lode angle
     * @todo Adapt for 2D dimension
     */
    static void CalculateLodeAngle(
        const double J2,
        const double J3,
        double& LodeAngle
        );

    /**
     * @brief This method computes the principal stresses vector
     * @param rPrincipalStressVector The vector of principal stresses
     * @param rStressVector The vector of stresses
     * @todo Adapt for 2D dimension
     */
    static void CalculatePrincipalStresses(
        Vector& rPrincipalStressVector,
        const Vector& rStressVector
        );

  private:
}; // class ConstitutiveLawUtilities
} // namespace Kratos
#endif /* KRATOS_CONSTITUTIVE_LAW_UTILITIES defined */
