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
#include "includes/node.h"
#include "geometries/geometry.h"

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
 * @author Alejandro Cornejo
 * @author Vicente Mataix Ferrandiz
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

    /// The definition of the bounded vector type
    typedef array_1d<double, VoigtSize> BoundedVectorType;

    /// The definition of the bounded matrix type
    typedef BoundedMatrix<double, Dimension, Dimension> BoundedMatrixType;

    /// The definition of the bounded matrix type
    typedef BoundedMatrix<double, VoigtSize, VoigtSize> BoundedMatrixVoigtType;

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
        if (Dimension == 3) {
            rDeviator = rStressVector;
            const double p_mean = I1 / 3.0;
            for (IndexType i = 0; i < Dimension; ++i)
                rDeviator[i] -= p_mean;
            rJ2 = 0.0;
            for (IndexType i = 0; i < Dimension; ++i)
                rJ2 += 0.5 * std::pow(rDeviator[i], 2);
            for (IndexType i = Dimension; i < 6; ++i)
                rJ2 += std::pow(rDeviator[i], 2);
        } else {
            rDeviator = rStressVector;
            const double p_mean = I1 / 3.0;
            for (IndexType i = 0; i < Dimension; ++i)
                rDeviator[i] -= p_mean;
            rJ2 = 0.5 * (std::pow(rDeviator[0], 2.0) + std::pow(rDeviator[1], 2.0) + std::pow(p_mean, 2.0)) + std::pow(rDeviator[2], 2.0);
        }
    }

    /**
     * @brief This method the uniaxial equivalent stress for Von Mises
     * @param rStressVector The stress vector S = C:E
     * @return The VM equivalent stress
     * @tparam TVector The themplate for the vector class
     */
    template<class TVector>
    static double CalculateVonMisesEquivalentStress(const TVector& rStressVector)
    {
        double I1, J2;
        array_1d<double, VoigtSize> deviator = ZeroVector(VoigtSize);

        ConstitutiveLawUtilities<VoigtSize>::CalculateI1Invariant(rStressVector, I1);
        ConstitutiveLawUtilities<VoigtSize>::CalculateJ2Invariant(rStressVector, I1, deviator, J2);

        return std::sqrt(3.0 * J2);
    }

private:

}; // class ConstitutiveLawUtilities
} // namespace Kratos
#endif /* KRATOS_CONSTITUTIVE_LAW_UTILITIES defined */
