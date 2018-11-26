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
#include "includes/ublas_interface.h"

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
 * @todo Finish adapt for 2D dimension
 */
template <SizeType TVoigtSize = 6>
class KRATOS_API(STRUCTURAL_MECHANICS_APPLICATION) ConstitutiveLawUtilities
{
  public:
    ///@name Type definitions
    ///@{

    /// The index type definition
    typedef std::size_t IndexType;

    /// We define the dimension
    static constexpr SizeType Dimension = TVoigtSize == 6 ? 3 : 2;

    /// We define the Voigt size
    static constexpr SizeType VoigtSize = TVoigtSize;

    /// The matrix type definition
    typedef Matrix MatrixType;

    /// the vector type definition
    typedef Vector VectorType;

    /// The definition of the bounded matrix type
    typedef array_1d<double, VoigtSize> BoundedVectorType;

    /// The definition of the bounded matrix type
    typedef BoundedMatrix<double, Dimension, Dimension> BoundedMatrixType;

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
     * @param rStressVector The stress vector on Voigt notation
     * @param rI1 The first invariant
     */
    static void CalculateI1Invariant(
        const BoundedVectorType& rStressVector,
        double& rI1
        );

    /**
     * @brief This method computes the second invariant from a given stress vector
     * @param rStressVector The stress vector on Voigt notation
     * @param rI2 The second invariant
     * @todo Adapt for 2D dimension
     */
    static void CalculateI2Invariant(
        const BoundedVectorType& rStressVector,
        double& rI2
        );

    /**
     * @brief This method computes the third invariant from a given stress vector
     * @param rStressVector The stress vector on Voigt notation
     * @param rI3 The third invariant
     * @todo Adapt for 2D dimension
     */
    static void CalculateI3Invariant(
        const BoundedVectorType& rStressVector,
        double& rI3
        );

    /**
     * @brief This method computes the second invariant of J
     * @param rStressVector The stress vector on Voigt notation
     * @param I1 The first invariant
     * @param rDeviator The deviator of the stress
     * @param rJ2 The second invariant of J
     */
    static void CalculateJ2Invariant(
        const BoundedVectorType& rStressVector,
        const double I1,
        BoundedVectorType& rDeviator,
        double& rJ2
        );

    /**
     * @brief This method computes the third invariant of J
     * @param rDeviator The deviator of the stress
     * @param rJ3 The third invariant of J
     */
    static void CalculateJ3Invariant(
        const BoundedVectorType& rDeviator,
        double& rJ3
        );

    /**
     * @brief This method computes the first vector
     * @param rFirstVector The first vector
     */
    static void CalculateFirstVector(BoundedVectorType& rFirstVector);

    /**
     * @brief This method computes the second vector
     * @param rDeviator The deviator of the stress
     * @param J2 The resultant J2 stress
     * @param rSecondVector The second vector
     */
    static void CalculateSecondVector(
        const BoundedVectorType& rDeviator,
        const double J2,
        BoundedVectorType& rSecondVector
        );

    /**
     * @brief This method computes the third vector
     * @param rDeviator The deviator of the stress
     * @param J2 The resultant J2 stress
     * @param rThirdVector The third vector
     * @todo Adapt for 2D dimension
     */
    static void CalculateThirdVector(
        const BoundedVectorType& rDeviator,
        const double J2,
        BoundedVectorType& rThirdVector
        );

    /**
     * @brief This method computes the lode angle
     * @param J2 The resultant J2 stress
     * @param J3 The resultant J3 stress
     * @param rLodeAngle The lode angle
     */
    static void CalculateLodeAngle(
        const double J2,
        const double J3,
        double& rLodeAngle
        );

    /**
     * @brief Calculation of the Green-Lagrange strain vector
     * @details See https://en.wikipedia.org/wiki/Finite_strain_theory#Seth%E2%80%93Hill_family_of_generalized_strain_tensors
     * @param rCauchyTensor The right Cauchy tensor
     * @param rStrainVector The Green-Lagrange strain vector
     */
    static void CalculateGreenLagrangianStrain(
        const MatrixType& rCauchyTensor,
        VectorType& rStrainVector
        );

    /**
     * @brief Calculation of the Almansi strain vector
     * @details See https://en.wikipedia.org/wiki/Finite_strain_theory#Seth%E2%80%93Hill_family_of_generalized_strain_tensors
     * @param rLeftCauchyTensor The left Cauchy tensor
     * @param rStrainVector The Almansi strain vector
     */
    static void CalculateAlmansiStrain(
        const MatrixType& rLeftCauchyTensor,
        VectorType& rStrainVector
        );

    /**
     * @brief Calculation of the Hencky strain vector (true strain, natural strain, logarithmic strain)
     * @details See https://en.wikipedia.org/wiki/Finite_strain_theory#Seth%E2%80%93Hill_family_of_generalized_strain_tensors
     * @param rCauchyTensor The right Cauchy tensor
     * @param rStrainVector The Hencky strain vector
     */
    static void CalculateHenckyStrain(
        const MatrixType& rCauchyTensor,
        VectorType& rStrainVector
        );

    /**
     * @brief Calculation of the Biot strain vector
     * @details See https://en.wikipedia.org/wiki/Finite_strain_theory#Seth%E2%80%93Hill_family_of_generalized_strain_tensors
     * @param rCauchyTensor The right Cauchy tensor
     * @param rStrainVector The Biot strain vector
     */
    static void CalculateBiotStrain(
        const MatrixType& rCauchyTensor,
        VectorType& rStrainVector
        );
    
    /**
     * @brief The deformation gradient F, like any invertible second-order tensor, can be decomposed, using the polar decomposition theorem, into a product of two second-order tensors (Truesdell and Noll, 1965): an orthogonal tensor and a positive definite symmetric tensor, i.e F = R U
     * @details See https://en.wikipedia.org/wiki/Finite_strain_theory#Polar_decomposition_of_the_deformation_gradient_tensor
     * @param rFDeformationGradient The deformation gradient
     * @param rRMatrix The rotation component
     * @param rUMatrix The pure displacement component
     */
    static void PolarDecomposition(
        const MatrixType& rFDeformationGradient,
        MatrixType& rRMatrix,
        MatrixType& rUMatrix
        );

    /**
     * @brief This method computes the principal stresses vector
     * @details http://www.continuummechanics.org/principalstress.html
     * @param rPrincipalStressVector The vector of principal stresses
     * @param rStressVector The vector of stresses
     * @todo Adapt for 2D dimension
     */
    static void CalculatePrincipalStresses(
        array_1d<double, Dimension>& rPrincipalStressVector,
        const BoundedVectorType& rStressVector
        );

    /**
     * @brief This method computes the principal stresses vector
     * @details Using Cardano formula and renormalizing (TODO)
     * @param rPrincipalStressVector The vector of principal stresses
     * @param rStressVector The vector of stresses
     * @todo Adapt for 2D dimension
     */
    static void CalculatePrincipalStressesWithCardano(
        array_1d<double, Dimension>& rPrincipalStressVector,
        const BoundedVectorType& rStressVector
        );

    /**
     * @brief This method calculates the projection operator
     * and calculates the Projection Operator
     * @details see "An energy-Equivalent" d+/d- Damage model with Enhanced
     * Microcrack Closure/Reopening Capabilities for Cohesive-Frictional
     * Materials" - M. Cervera and C. Tesei.
     * @param rStrainVector The Strain Vector
     * @param rProjectionOperator The projection operator
     */
    static void CalculateProjectionOperator(
        const Vector& rStrainVector,
        MatrixType& rProjectionOperator
        );

    /**
     * @brief This method performs Spectral Decomposition of the Stress Vector/Tensor
     * @details see "An energy-Equivalent" d+/d- Damage model with Enhanced
     * Microcrack Closure/Reopening Capabilities for Cohesive-Frictional
     * Materials" - M. Cervera and C. Tesei.
     * @param rStressVector The Stress Vector
     * @param rStressVectorTension The Stress Vector
     * @param rStressVectorCompression The Stress Vector
     * @param rMatrixTension The Stress Vector
     * @param rMatrixCompression The Stress Vector
     */
    static void SpectralDecomposition(
        const BoundedVectorType& rStressVector,
        BoundedVectorType& rStressVectorTension,
        BoundedVectorType& rStressVectorCompression
        );

    /**
     * @brief This computes the linear plastic deformation gradient increment
     * @param rPlasticPotentialDerivative The derivative of the plastic potential
     * @param PlasticConsistencyFactorIncrement The incremenetal of plastic flow
     */
    static MatrixType CalculateLinearPlasticDeformationGradientIncrement(
        const BoundedVectorType& rPlasticPotentialDerivative,
        const double PlasticConsistencyFactorIncrement
        );

    /**
     * @brief This computes the exponential plastic deformation gradient increment
     * @param rPlasticPotentialDerivative The derivative of the plastic potential
     * @param PlasticConsistencyFactorIncrement The incremenetal of plastic flow
     * @param rRe The rotation decomposition of the elastic eformation
     */
    static MatrixType CalculateExponentialPlasticDeformationGradientIncrement(
        const BoundedVectorType& rPlasticPotentialDerivative,
        const double PlasticConsistencyFactorIncrement,
        const MatrixType& rRe
        );

  private:

}; // class ConstitutiveLawUtilities
} // namespace Kratos
#endif /* KRATOS_CONSTITUTIVE_LAW_UTILITIES defined */
