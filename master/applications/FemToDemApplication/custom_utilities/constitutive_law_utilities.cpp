//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \.
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics FemDem Application
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Alejandro Cornejo Velazquez
//                   Vicente Mataix Ferrandiz
//

// System includes

// External includes

// Project includes
#include "includes/global_variables.h"
#include "utilities/math_utils.h"
#include "custom_utilities/constitutive_law_utilities.h"

namespace Kratos
{

template<>
void ConstitutiveLawUtilities<6>::CalculateElasticMatrix(
    Matrix& rConstitutiveMatrix,
    const double E,
    const double nu
    )
{
    if (rConstitutiveMatrix.size1() != VoigtSize || rConstitutiveMatrix.size2() != VoigtSize)
        rConstitutiveMatrix.resize(VoigtSize, VoigtSize, false);

    const double c1 = E / ((1.00 + nu) * (1 - 2 * nu));
    const double c2 = c1 * (1 - nu);
    const double c3 = c1 * nu;
    const double c4 = c1 * 0.5 * (1 - 2 * nu);

    rConstitutiveMatrix(0, 0) = c2;
    rConstitutiveMatrix(0, 1) = c3;
    rConstitutiveMatrix(0, 2) = c3;
    rConstitutiveMatrix(1, 0) = c3;
    rConstitutiveMatrix(1, 1) = c2;
    rConstitutiveMatrix(1, 2) = c3;
    rConstitutiveMatrix(2, 0) = c3;
    rConstitutiveMatrix(2, 1) = c3;
    rConstitutiveMatrix(2, 2) = c2;
    rConstitutiveMatrix(3, 3) = c4;
    rConstitutiveMatrix(4, 4) = c4;
    rConstitutiveMatrix(5, 5) = c4;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void ConstitutiveLawUtilities<3>::CalculateElasticMatrix(
    Matrix& rConstitutiveMatrix,
    const double E,
    const double nu
    )
{
    // Assuming plane strain TODO
    if (rConstitutiveMatrix.size1() != VoigtSize || rConstitutiveMatrix.size2() != VoigtSize)
        rConstitutiveMatrix.resize(VoigtSize, VoigtSize, false);

    const double c0 = E / ((1.00 + nu) * (1 - 2 * nu));
    const double c1 = (1.00 - nu) * c0;
    const double c2 = c0 * nu;
    const double c3 = (0.5 - nu) * c0;

    rConstitutiveMatrix(0, 0) = c1;
    rConstitutiveMatrix(0, 1) = c2;
    rConstitutiveMatrix(1, 0) = c2;
    rConstitutiveMatrix(1, 1) = c1;
    rConstitutiveMatrix(2, 2) = c3;
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TVoigtSize>
void ConstitutiveLawUtilities<TVoigtSize>::CalculateDeviatoricStrainVector(
    const Vector& rStrainVector,
    const Vector& rVolumetricStrainVector,
    Vector &rDeviatoricStrainVector
    )
{
    if (rDeviatoricStrainVector.size() != VoigtSize)
        rDeviatoricStrainVector.resize(VoigtSize);

    noalias(rDeviatoricStrainVector) = rStrainVector - rVolumetricStrainVector;
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TVoigtSize>
void ConstitutiveLawUtilities<TVoigtSize>::CalculateVolumetricStrainVector(
    const Vector& rStrainVector,
    Vector& rVolumetricStrainVector
    )
{
    if (rVolumetricStrainVector.size() != VoigtSize)
        rVolumetricStrainVector.resize(VoigtSize);

    BoundedVectorType r_identity_vector = ZeroVector(VoigtSize);
    CalculateIdentityVector(r_identity_vector);
    double strain_trace = 0.0;

    // Compute the trace of strain vector
    for (IndexType i = 0; i < Dimension; ++i)
        strain_trace += rStrainVector[i];

    noalias(rVolumetricStrainVector) = 1.0 / 3.0 * (strain_trace)*r_identity_vector;
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TVoigtSize>
double ConstitutiveLawUtilities<TVoigtSize>::CalculateShearModulus(
    const double YoungModulus,
    const double PoissonRatio
    )
{
    return YoungModulus / (3.0 * (1.0 - 2.0 * PoissonRatio));
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TVoigtSize>
double ConstitutiveLawUtilities<TVoigtSize>::CalculateBulkModulus(
    const double YoungModulus,
    const double PoissonRatio
    )
{
    return YoungModulus / (2.0 * (1.0 + PoissonRatio));
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TVoigtSize>
void ConstitutiveLawUtilities<TVoigtSize>::CalculateI1Invariant(
    const BoundedVectorType& rStressVector,
    double& rI1
    )
{
    rI1 = rStressVector[0];
    for (IndexType i = 1; i < Dimension; ++i)
        rI1 += rStressVector[i];
}


/***********************************************************************************/
/***********************************************************************************/

template<>
void ConstitutiveLawUtilities<6>::CalculateI2Invariant(
    const BoundedVectorType& rStressVector,
    double& rI2
    )
{
    rI2 = (rStressVector[0] + rStressVector[2]) * rStressVector[1] + rStressVector[0] * rStressVector[2] +
            -rStressVector[3] * rStressVector[3] - rStressVector[4] * rStressVector[4] - rStressVector[5] * rStressVector[5];
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void ConstitutiveLawUtilities<3>::CalculateI2Invariant(
    const BoundedVectorType& rStressVector,
    double& rI2
    )
{
    rI2 = rStressVector[0] * rStressVector[1] - std::pow(rStressVector[2], 2);
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void ConstitutiveLawUtilities<6>::CalculateI3Invariant(
    const BoundedVectorType& rStressVector,
    double& rI3
    )
{
    rI3 = (rStressVector[1] * rStressVector[2] - rStressVector[4] * rStressVector[4]) * rStressVector[0] -
            rStressVector[1] * rStressVector[5] * rStressVector[5] - rStressVector[2] * rStressVector[3] * rStressVector[3] +
            2.0 * rStressVector[3] * rStressVector[4] * rStressVector[5];
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void ConstitutiveLawUtilities<3>::CalculateI3Invariant(
    const BoundedVectorType& rStressVector,
    double& rI3
    )
{
    rI3 = rStressVector[0] * rStressVector[1] - std::pow(rStressVector[2], 2);
}
/***********************************************************************************/
/***********************************************************************************/

template<>
void ConstitutiveLawUtilities<6>::CalculateJ2Invariant(
    const BoundedVectorType& rStressVector,
    const double I1,
    BoundedVectorType& rDeviator,
    double& rJ2
    )
{
    rDeviator = rStressVector;
    const double p_mean = I1 / static_cast<double>(Dimension);

    for (IndexType i = 0; i < Dimension; ++i)
        rDeviator[i] -= p_mean;

    rJ2 = 0.0;
    for (IndexType i = 0; i < Dimension; ++i)
        rJ2 += 0.5 * std::pow(rDeviator[i], 2);
    for (IndexType i = Dimension; i < 6; ++i)
        rJ2 += std::pow(rDeviator[i], 2);
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void ConstitutiveLawUtilities<3>::CalculateJ2Invariant(
    const BoundedVectorType& rStressVector,
    const double I1,
    BoundedVectorType& rDeviator,
    double& rJ2
    )
{
    rDeviator = rStressVector;
    const double p_mean = I1 / 3.0;

    for (IndexType i = 0; i < Dimension; ++i)
        rDeviator[i] -= p_mean;

    rJ2 = 0.5 * (std::pow(rDeviator[0], 2.0) + std::pow(rDeviator[1], 2.0) + std::pow(p_mean, 2.0)) +
          std::pow(rDeviator[2], 2.0);
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void ConstitutiveLawUtilities<6>::CalculateJ3Invariant(
    const BoundedVectorType& rDeviator,
    double& rJ3
    )
{
    rJ3 = rDeviator[0] * (rDeviator[1] * rDeviator[2] - rDeviator[4] * rDeviator[4]) +
            rDeviator[3] * (-rDeviator[3] * rDeviator[2] + rDeviator[5] * rDeviator[4]) +
            rDeviator[5] * (rDeviator[3] * rDeviator[4] - rDeviator[5] * rDeviator[1]);
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void ConstitutiveLawUtilities<3>::CalculateJ3Invariant(
    const BoundedVectorType& rDeviator,
    double& rJ3
    )
{
    rJ3 = rDeviator[0] * rDeviator[1] - std::pow(rDeviator[2], 2);
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void ConstitutiveLawUtilities<6>::CalculateFirstVector(BoundedVectorType& rFirstVector)
{
    rFirstVector[0] = 1.0;
    rFirstVector[1] = 1.0;
    rFirstVector[2] = 1.0;
    rFirstVector[3] = 0.0;
    rFirstVector[4] = 0.0;
    rFirstVector[5] = 0.0;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void ConstitutiveLawUtilities<3>::CalculateFirstVector(BoundedVectorType& rFirstVector)
{
    rFirstVector[0] = 1.0;
    rFirstVector[1] = 1.0;
    rFirstVector[2] = 0.0;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void ConstitutiveLawUtilities<6>::CalculateSecondVector(
    const BoundedVectorType& rDeviator,
    const double J2,
    BoundedVectorType& rSecondVector
    )
{
    const double twosqrtJ2 = 2.0 * std::sqrt(J2);

    if (twosqrtJ2 > tolerance) {
        for (IndexType i = 0; i < 6; ++i) {
            rSecondVector[i] = rDeviator[i] / (twosqrtJ2);
        }

        for (IndexType i = Dimension; i < 6; ++i)
            rSecondVector[i] *= 2.0;
    } else {
        noalias(rSecondVector) = ZeroVector(VoigtSize);
    }

}

/***********************************************************************************/
/***********************************************************************************/

template<>
void ConstitutiveLawUtilities<3>::CalculateSecondVector(
    const BoundedVectorType& rDeviator,
    const double J2,
    BoundedVectorType& rSecondVector
    )
{
    const double twosqrtJ2 = 2.0 * std::sqrt(J2);
    for (IndexType i = 0; i < 3; ++i) {
        rSecondVector[i] = rDeviator[i] / (twosqrtJ2);
    }
    rSecondVector[2] *= 2.0;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void ConstitutiveLawUtilities<6>::CalculateThirdVector(
    const BoundedVectorType& rDeviator,
    const double J2,
    BoundedVectorType& rThirdVector
    )
{
    const double J2thirds = J2 / 3.0;

    rThirdVector[0] = rDeviator[1] * rDeviator[2] - rDeviator[4] * rDeviator[4] + J2thirds;
    rThirdVector[1] = rDeviator[0] * rDeviator[2] - rDeviator[5] * rDeviator[5] + J2thirds;
    rThirdVector[2] = rDeviator[0] * rDeviator[1] - rDeviator[3] * rDeviator[3] + J2thirds;
    rThirdVector[3] = 2.0 * (rDeviator[4] * rDeviator[5] - rDeviator[3] * rDeviator[2]);
    rThirdVector[4] = 2.0 * (rDeviator[3] * rDeviator[4] - rDeviator[1] * rDeviator[5]);
    rThirdVector[5] = 2.0 * (rDeviator[5] * rDeviator[3] - rDeviator[0] * rDeviator[4]);
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void ConstitutiveLawUtilities<3>::CalculateThirdVector(
    const BoundedVectorType& rDeviator,
    const double J2,
    BoundedVectorType& rThirdVector
    )
{
    const double J2thirds = J2 / 3.0;

    rThirdVector[0] = rDeviator[1] * rDeviator[2] + J2thirds;
    rThirdVector[1] = rDeviator[0] * rDeviator[2] + J2thirds;
    rThirdVector[2] = rDeviator[0] * rDeviator[1] - std::pow(rDeviator[3], 2) + J2thirds;
    rThirdVector[3] = -2.0 * rDeviator[3] * rDeviator[2];
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TVoigtSize>
void ConstitutiveLawUtilities<TVoigtSize>::CalculateLodeAngle(
    const double J2,
    const double J3,
    double& rLodeAngle
    )
{
    if (J2 > tolerance) {
        double sint3 = (-3.0 * std::sqrt(3.0) * J3) / (2.0 * J2 * std::sqrt(J2));
        if (sint3 < -0.95)
            sint3 = -1.0;
        else if (sint3 > 0.95)
            sint3 = 1.0;
        rLodeAngle = std::asin(sint3) / 3.0;
    } else {
        rLodeAngle = 0.0;
    }
}


/***********************************************************************************/
/***********************************************************************************/

template<>
void ConstitutiveLawUtilities<6>::CalculatePrincipalStresses(
    array_1d<double, Dimension>& rPrincipalStressVector,
    const BoundedVectorType& rStressVector
    )
{
    double I1, I2, I3;
    BoundedMatrix<double, Dimension, Dimension> tensor = MathUtils<double>::VectorToSymmetricTensor<BoundedVectorType, BoundedMatrix<double, Dimension, Dimension>>(rStressVector);
    double norm = norm_frobenius(tensor);
    if (norm < tolerance) norm = 1.0;
    const BoundedVectorType norm_stress_vector = rStressVector/norm;
    CalculateI1Invariant(norm_stress_vector, I1);
    CalculateI2Invariant(norm_stress_vector, I2);
    CalculateI3Invariant(norm_stress_vector, I3);

    const double II1 = std::pow(I1, 2);

    const double R = (2.0 * II1 * I1 - 9.0 * I2 * I1 + 27.0 * I3) / 54.0;
    const double Q = (3.0 * I2 - II1) / 9.0;

    if (std::abs(Q) > tolerance) {
        double cos_phi = R / (std::sqrt(-std::pow(Q, 3)));
        if (cos_phi >= 1.0)
            cos_phi = 1.0;
        else if (cos_phi <= -1.0)
            cos_phi = -1.0;
        const double phi = std::acos(cos_phi);
        const double phi_3 = phi / 3.0;

        const double aux1 = 2.0 * std::sqrt(-Q);
        const double aux2 = I1 / 3.0;
        const double deg_120 = 2.0 / 3.0 * Globals::Pi;

        for (IndexType i = 0; i < 3; ++i) {
            rPrincipalStressVector[i] = norm * (aux2 + aux1 * std::cos(phi_3 + deg_120 * i));
        }
    } else {
        for (IndexType i = 0; i < Dimension; ++i) {
            rPrincipalStressVector[i] = rStressVector[i];
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void ConstitutiveLawUtilities<3>::CalculatePrincipalStresses(
    array_1d<double, Dimension>& rPrincipalStressVector,
    const BoundedVectorType& rStressVector
    )
{
    rPrincipalStressVector[0] = 0.5 * (rStressVector[0] + rStressVector[1]) +
        std::sqrt(std::pow(0.5 * (rStressVector[0] - rStressVector[1]), 2)  +
        std::pow(rStressVector[2], 2));

    rPrincipalStressVector[1] = 0.5 * (rStressVector[0] + rStressVector[1]) -
        std::sqrt(std::pow(0.5 * (rStressVector[0] - rStressVector[1]), 2)  +
        std::pow(rStressVector[2], 2));
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TVoigtSize>
void ConstitutiveLawUtilities<TVoigtSize>::CalculatePrincipalStressesWithCardano(
    array_1d<double, Dimension>& rPrincipalStressVector,
    const BoundedVectorType& rStressVector
    )
{
    double a, b, c;
    BoundedMatrix<double, Dimension, Dimension> tensor = MathUtils<double>::VectorToSymmetricTensor<BoundedVectorType, BoundedMatrix<double, Dimension, Dimension>>(rStressVector);
    double norm = norm_frobenius(tensor);
    if (norm < tolerance) norm = 1.0;
    const BoundedVectorType norm_stress_vector = rStressVector/norm;
    CalculateI1Invariant(norm_stress_vector, a);
    CalculateI2Invariant(norm_stress_vector, b);
    CalculateI3Invariant(norm_stress_vector, c);

    const double p = b - std::pow(a, 2)/3.0;
    const double q = 2.0 * std::pow(a, 3)/27.0 - (a * b)/3.0 + c;
    const double discriminant = std::pow(q, 2) + 4.0/27.0 * std::pow(p, 3);

    if (std::abs(p) > tolerance) {
        if (discriminant > tolerance) { // This is bad news (complex numbers)
            KRATOS_ERROR << "Complex conjugated solutions" << std::endl;
        } else if (discriminant < - tolerance) {
            const double aux = 2.0 * std::sqrt(-p/3.0);
            const double base_sol = a / 3.0;
            const double phi_3 = 1.0/3.0 * std::acos(-3.0*q/(p * 2.0) * std::sqrt(-3.0/p));
            for (IndexType i = 0; i < 3; ++i) {
                rPrincipalStressVector[i] = (base_sol + aux * std::cos(phi_3 - 2.0/3.0 * Globals::Pi * i)) * norm;
            }
        } else { // Equal to zero
            rPrincipalStressVector[0] = 3.0 * q/p;
            rPrincipalStressVector[1] = -3.0/2.0 * q/p;
            rPrincipalStressVector[2] = rPrincipalStressVector[1];
        }
    } else {
        for (IndexType i = 0; i < Dimension; ++i) {
            rPrincipalStressVector[i] = rStressVector[i];
        }
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TVoigtSize>
void ConstitutiveLawUtilities<TVoigtSize>::CalculateHenckyStrain(
    const MatrixType& rCauchyTensor,
    Vector& rStrainVector
    )
{
    // Doing resize in case is needed
    if (rStrainVector.size() != VoigtSize)
        rStrainVector.resize(VoigtSize, false);

    // Declare the different matrix
    BoundedMatrixType eigen_values_matrix, eigen_vectors_matrix;

    // Decompose matrix
    MathUtils<double>::GaussSeidelEigenSystem(rCauchyTensor, eigen_vectors_matrix, eigen_values_matrix, 1.0e-16, 20);

    // Calculate the eigenvalues of the E matrix
    for (IndexType i = 0; i < Dimension; ++i) {
        eigen_values_matrix(i, i) = 0.5 * std::log(eigen_values_matrix(i, i));
    }

    // Calculate E matrix
    BoundedMatrixType E_matrix;
    MathUtils<double>::BDBtProductOperation(E_matrix, eigen_values_matrix, eigen_vectors_matrix);

    // Hencky Strain Calculation
    rStrainVector = MathUtils<double>::StrainTensorToVector(E_matrix, TVoigtSize);
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TVoigtSize>
void ConstitutiveLawUtilities<TVoigtSize>::CalculateBiotStrain(
    const MatrixType& rCauchyTensor,
    Vector& rStrainVector
    )
{
    // Doing resize in case is needed
    if (rStrainVector.size() != VoigtSize)
        rStrainVector.resize(VoigtSize, false);

    // Declare the different matrix
    BoundedMatrixType eigen_values_matrix, eigen_vectors_matrix;

    // Decompose matrix
    MathUtils<double>::GaussSeidelEigenSystem(rCauchyTensor, eigen_vectors_matrix, eigen_values_matrix, 1.0e-16, 20);

    // Calculate the eigenvalues of the E matrix
    for (IndexType i = 0; i < Dimension; ++i) {
        eigen_values_matrix(i, i) = std::sqrt(eigen_values_matrix(i, i));
    }

    // Calculate E matrix
    BoundedMatrixType E_matrix;
    MathUtils<double>::BDBtProductOperation(E_matrix, eigen_values_matrix, eigen_vectors_matrix);

    // Biot Strain Calculation
    rStrainVector = MathUtils<double>::StrainTensorToVector(E_matrix, TVoigtSize);
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TVoigtSize>
void ConstitutiveLawUtilities<TVoigtSize>::CalculateAlmansiStrain(
    const MatrixType& rLeftCauchyTensor,
    Vector& rStrainVector
    )
{
    // Doing resize in case is needed
    if (rStrainVector.size() != VoigtSize)
        rStrainVector.resize(VoigtSize, false);

    // Identity matrix
    MatrixType identity_matrix(Dimension, Dimension);
    for (IndexType i = 0; i < Dimension; ++i) {
        for (IndexType j = 0; j < Dimension; ++j) {
            if (i == j) identity_matrix(i, j) = 1.0;
            else identity_matrix(i, j) = 0.0;
        }
    }

    // Calculating the inverse of the left Cauchy tensor
    MatrixType inverse_B_tensor ( Dimension, Dimension );
    double aux_det_b = 0;
    MathUtils<double>::InvertMatrix( rLeftCauchyTensor, inverse_B_tensor, aux_det_b);

    // Calculate E matrix
    const BoundedMatrixType E_matrix = 0.5 * (identity_matrix - inverse_B_tensor);

    // Almansi Strain Calculation
    rStrainVector = MathUtils<double>::StrainTensorToVector(E_matrix, TVoigtSize);
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TVoigtSize>
void ConstitutiveLawUtilities<TVoigtSize>::CalculateGreenLagrangianStrain(
    const MatrixType& rCauchyTensor,
    Vector& rStrainVector
    )
{
    // Doing resize in case is needed
    if (rStrainVector.size() != VoigtSize)
        rStrainVector.resize(VoigtSize, false);

    // Identity matrix
    MatrixType identity_matrix(Dimension, Dimension);
    for (IndexType i = 0; i < Dimension; ++i) {
        for (IndexType j = 0; j < Dimension; ++j) {
            if (i == j) identity_matrix(i, j) = 1.0;
            else identity_matrix(i, j) = 0.0;
        }
    }

    // Calculate E matrix
    const BoundedMatrixType E_matrix = 0.5 * (rCauchyTensor - identity_matrix);

    // Green-Lagrangian Strain Calculation
    rStrainVector = MathUtils<double>::StrainTensorToVector(E_matrix, TVoigtSize);
}

/***********************************************************************************/
/***********************************************************************************/

template class ConstitutiveLawUtilities<3>;
template class ConstitutiveLawUtilities<6>;

} // namespace Kratos
