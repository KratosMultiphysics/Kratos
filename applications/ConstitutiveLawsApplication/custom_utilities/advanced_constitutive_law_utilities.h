// KRATOS ___                _   _ _         _   _             __                       _
//       / __\___  _ __  ___| |_(_) |_ _   _| |_(_)_   _____  / /  __ ___      _____   /_\  _ __  _ __
//      / /  / _ \| '_ \/ __| __| | __| | | | __| \ \ / / _ \/ /  / _` \ \ /\ / / __| //_\\| '_ \| '_  |
//     / /__| (_) | | | \__ \ |_| | |_| |_| | |_| |\ V /  __/ /__| (_| |\ V  V /\__ \/  _  \ |_) | |_) |
//     \____/\___/|_| |_|___/\__|_|\__|\__,_|\__|_| \_/ \___\____/\__,_| \_/\_/ |___/\_/ \_/ .__/| .__/
//                                                                                         |_|   |_|
//
//  License:         BSD License
//                   license: structural_mechanics_application/license.txt
//
//  Main authors:    Alejandro Cornejo
//

#pragma once

// System includes

// External includes

// Project includes
#include "includes/ublas_interface.h"
#include "includes/node.h"
#include "includes/constitutive_law.h"
#include "geometries/geometry.h"

namespace Kratos
{
///@name Kratos Globals
///@{

///@}
///@name Type Definitions
///@{

    /// The size type definition
    using SizeType = std::size_t;

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
 * @class AdvancedConstitutiveLawUtilities
 * @ingroup StructuralMechanicsApplication
 * @brief This class includes several utilities necessaries for the computation of the constitutive law
 * @details The methods are static, so it can be called without constructing the class
 * @tparam TVoigtSize The number of components on the Voigt notation
 * @author Alejandro Cornejo
 */
template <SizeType TVoigtSize = 6>
class KRATOS_API(CONSTITUTIVE_LAWS_APPLICATION) AdvancedConstitutiveLawUtilities
{
  public:
    ///@name Type definitions
    ///@{

    /// The index type definition
    using IndexType = std::size_t;

    /// We define the dimension
    static constexpr SizeType Dimension = TVoigtSize == 6 ? 3 : 2;

    /// We define the Voigt size
    static constexpr SizeType VoigtSize = TVoigtSize;

    /// The matrix type definition
    using MatrixType = Matrix;

    /// the vector type definition
    using VectorType = Vector;

    /// The definition of the bounded vector type
    using BoundedVectorType = array_1d<double, VoigtSize>;

    /// The definition of the bounded matrix type
    using BoundedMatrixType = BoundedMatrix<double, Dimension, Dimension>;

    /// The definition of the bounded matrix type
    using BoundedMatrixVoigtType = BoundedMatrix<double, VoigtSize, VoigtSize>;

    /// Node type definition
    using NodeType = Node;

    /// Geometry definitions
    using GeometryType = Geometry<NodeType>;

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
     * @tparam TVector The themplate for the vector class
     */
    template<class TVector>
    static void CalculateI1Invariant(
        const TVector& rStressVector,
        double& rI1
        )
    {
        rI1 = rStressVector[0];
        for (IndexType i = 1; i < Dimension; ++i)
            rI1 += rStressVector[i];
    }

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
     * @tparam TVector The themplate for the vector class
     */
    template<class TVector>
    static void CalculateJ2Invariant(
        const TVector& rStressVector,
        const double I1,
        BoundedVectorType& rDeviator,
        double& rJ2
        )
    {
        noalias(rDeviator) = rStressVector;
        const double p_mean = I1 / 3.0;
        if constexpr (Dimension == 3) {
            for (IndexType i = 0; i < Dimension; ++i)
                rDeviator[i] -= p_mean;
            rJ2 = 0.0;
            for (IndexType i = 0; i < Dimension; ++i)
                rJ2 += 0.5 * std::pow(rDeviator[i], 2);
            for (IndexType i = Dimension; i < 6; ++i)
                rJ2 += std::pow(rDeviator[i], 2);
        } else {
            for (IndexType i = 0; i < Dimension; ++i)
                rDeviator[i] -= p_mean;
            rJ2 = 0.5 * (std::pow(rDeviator[0], 2.0) + std::pow(rDeviator[1], 2.0) + std::pow(p_mean, 2.0)) + std::pow(rDeviator[2], 2.0);
        }
    }

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
     * to be used in the derivative of the yield surface
     * @param rFirstVector The first vector is the
     * d(I1)/d(Stress)
     */
    static void CalculateFirstVector(BoundedVectorType& rFirstVector);

    /**
     * @brief This method computes the first vector
     * to be used in the derivative of the yield surface
     * @param rDeviator The deviator of the stress
     * @param J2 The resultant J2 stress
     * @param rSecondVector The second vector is
     * 1/(2*sqrt(J2))*d(J2)/d(Stress)
     */
    static void CalculateSecondVector(
        const BoundedVectorType& rDeviator,
        const double J2,
        BoundedVectorType& rSecondVector
        );

    /**
     * @brief This method computes the third vector
     * to be used in the derivative of the yield surface
     * @param rDeviator The deviator of the stress
     * @param J2 The resultant J2 stress
     * @param rThirdVector The third vector is
     * d(J3)/d(Stress)
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
     * @brief Calculates the maximal distance between corner node of a geometry and its center
     * @param rGeometry The geometry to compute
     * @return The characteristic length
     */
    static double CalculateCharacteristicLength(const GeometryType& rGeometry);

    /**
     * @brief Calculates the maximal distance between corner node of a geometry and its center (on reference configuration)
     * @param rGeometry The geometry to compute
     * @return The characteristic length
     */
    static double CalculateCharacteristicLengthOnReferenceConfiguration(const GeometryType& rGeometry);

    /**
     * @brief This method computes the equivalent deformation gradient for the elements which provide the deformation gradient as input
     * @param rStrainVector The strain vector
     */
    static Matrix ComputeEquivalentSmallDeformationDeformationGradient(const Vector& rStrainVector);

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
     * @brief This computes the elastic def gradient tensor
     * Fe = F*(Fp)^-1
     * @param rF The total def gradient tensor
     * @param rFp The plastic def gradient tensor
     */
    static Matrix CalculateElasticDeformationGradient(
        const MatrixType& rF,
        const MatrixType& rFp
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
     * @brief This computes the Fp from F and Fe
     * @param rF The total def gradient tensor
     * @param rFp The plastic def gradient tensor
     */
    static Matrix CalculatePlasticDeformationGradientFromElastic(
        const MatrixType& rF,
        const MatrixType& rFp
        );

    /**
     * @brief This updates the exponential elastic deformation gradient
     * @param rElasticTrial The elastic trial deformation gradient
     * @param rPlasticPotentialDerivative The derivative of the plastic potential
     * @param PlasticConsistencyFactorIncrement The incremenetal of plastic flow
     * @param rRe The rotation decomposition of the elastic eformation
     * @note Formula 14.75 in Souza book
     */
    static MatrixType CalculateExponentialElasticDeformationGradient(
        const MatrixType& rTrialFe,
        const BoundedVectorType& rPlasticPotentialDerivative,
        const double PlasticConsistencyFactorIncrement,
        const MatrixType& rRe
        );

    /**
     * @brief This computes the plastic strain from Fp
     * @param rFp The plastic deformation gradient
     * @param rPlasticStrainVector The plastic strain vector
     */
    static void CalculatePlasticStrainFromFp(
        const MatrixType& rFp,
        Vector& rPlasticStrainVector
        );

    /**
     * @brief This computes the elastic deformation gradient
     * @param rElasticTrial The elastic trial deformation gradient
     * @param rPlasticPotentialDerivative The derivative of the plastic potential
     * @param PlasticConsistencyFactorIncrement The incremenetal of plastic flow
     * @param rRe The rotation decomposition of the elastic eformation
     */
    static MatrixType CalculateDirectElasticDeformationGradient(
        const MatrixType& rElasticTrial,
        const BoundedVectorType& rPlasticPotentialDerivative,
        const double PlasticConsistencyFactorIncrement,
        const MatrixType& rRe
        );

    /**
     * @brief This computes the exponential plastic deformation gradient increment
     * @param rPlasticPotentialDerivative The derivative of the plastic potential
     * @param PlasticConsistencyFactorIncrement The incremenetal of plastic flow
     * @param rRe The rotation decomposition of the elastic eformation
     * @note Formula 14.73 in Souza book
     */
    static MatrixType CalculateExponentialPlasticDeformationGradientIncrement(
        const BoundedVectorType& rPlasticPotentialDerivative,
        const double PlasticConsistencyFactorIncrement,
        const MatrixType& rRe
        );

    /**
     * @brief This computes the exponential plastic deformation gradient increment
     * @param rPlasticPotentialDerivative The derivative of the plastic potential
     * @param PlasticConsistencyFactorIncrement The incremenetal of plastic flow
     * @param rRe The rotation decomposition of the elastic eformation
     * @note Formula 14.74 in Souza book
     */
    static MatrixType CalculateDirectPlasticDeformationGradientIncrement(
        const BoundedVectorType& rPlasticPotentialDerivative,
        const double PlasticConsistencyFactorIncrement,
        const MatrixType& rRe
        );

    /**
     * @brief This computes the rotation matrix for the 1st Euler angle
     * http://mathworld.wolfram.com/EulerAngles.html
     */
    static void CalculateRotationOperatorEuler1(
        const double EulerAngle1,
        BoundedMatrix<double, 3, 3>& rRotationOperator
        );

    /**
     * @brief This computes the rotation matrix for the 2nd Euler angle
     * http://mathworld.wolfram.com/EulerAngles.html
     */
    static void CalculateRotationOperatorEuler2(
        const double EulerAngle2,
        BoundedMatrix<double, 3, 3>& rRotationOperator
        );

    /**
     * @brief This computes the rotation matrix for the 3rd Euler angle
     * http://mathworld.wolfram.com/EulerAngles.html
     */
    static void CalculateRotationOperatorEuler3(
        const double EulerAngle3,
        BoundedMatrix<double, 3, 3>& rRotationOperator
        );

    /**
     * @brief This computes the total rotation matrix
     * rotates from global to local coordinates.
     * The so-called "x convention" is used.
     * Order of the rotations:
     *    1. The first rotation PHI around the Z-axis
     *    2. The second rotation THETA around the X'-axis
     *    3. The third rotation HI around the former Z'-axis
     * more info: http://mathworld.wolfram.com/EulerAngles.html
     */
    static void CalculateRotationOperator(
        const double EulerAngle1, // phi
        const double EulerAngle2, // theta
        const double EulerAngle3, // hi
        BoundedMatrix<double, 3, 3>& rRotationOperator
        );

     /**
     * @brief This computes the MacaullyBrackets of a double
     */
    static double MacaullyBrackets(const double Number);

    /**
     * @brief This substracts the thermal strain contribution to a vector
     */
    static void SubstractThermalStrain(
        ConstitutiveLaw::StrainVectorType &rStrainVector,
        const double ReferenceTemperature,
        ConstitutiveLaw::Parameters &rParameters,
        const bool IsPlaneStrain = false
        );

    /**
     * @brief This retrieves an interpolated nodal variable to a GP
     */
    static double CalculateInGaussPoint(
        const Variable<double> &rVariableInput,
        ConstitutiveLaw::Parameters &rParameters,
        unsigned int step = 0
        );

    /**
     * @brief This retrieves a double type variable checking the accessor
     */
    static double GetMaterialPropertyThroughAccessor(
        const Variable<double>& rVariable,
        ConstitutiveLaw::Parameters &rValues
        );

    /**
     * @brief This retrieves a double type variable from a table if exists, assumes TEMPERATURE to be the independent variable
     */
    static double GetPropertyFromTemperatureTable(
        const Variable<double>& rVariable,
        ConstitutiveLaw::Parameters &rValues,
        const double Temperature
        );

}; // class AdvancedConstitutiveLawUtilities
} // namespace Kratos