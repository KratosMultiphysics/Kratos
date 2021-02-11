// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                     license: structural_mechanics_application/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//                   Alejandro Cornejo
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

template<SizeType TVoigtSize>
double ConstitutiveLawUtilities<TVoigtSize>::CalculateCharacteristicLength(const GeometryType& rGeometry)
{
    double radius = 0.0;
    const Point& r_center = rGeometry.Center();

    for(IndexType i_node = 0; i_node < rGeometry.PointsNumber(); ++i_node)  {
        const array_1d<double, 3>& aux_vector = r_center.Coordinates() - rGeometry[i_node].Coordinates();
        double aux_value = inner_prod(aux_vector, aux_vector);
        if(aux_value > radius) radius = aux_value;
    }

    return std::sqrt(radius);
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TVoigtSize>
double ConstitutiveLawUtilities<TVoigtSize>::CalculateCharacteristicLengthOnReferenceConfiguration(const GeometryType& rGeometry)
{
    double radius = 0.0;

    const SizeType points_number = rGeometry.size();
	array_1d<double, 3> center = rGeometry[0].GetInitialPosition().Coordinates();
    for ( IndexType i_node = 1 ; i_node < points_number ; ++i_node ) {
        center += rGeometry[i_node].GetInitialPosition().Coordinates();
    }
    center /= static_cast<double>( points_number );

    for(IndexType i_node = 0; i_node < points_number; ++i_node)  {
        const array_1d<double, 3>& aux_vector = center - rGeometry[i_node].GetInitialPosition().Coordinates();
        double aux_value = inner_prod(aux_vector, aux_vector);
        if(aux_value > radius) radius = aux_value;
    }

    return std::sqrt(radius);
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TVoigtSize>
Matrix ConstitutiveLawUtilities<TVoigtSize>::ComputeEquivalentSmallDeformationDeformationGradient(const Vector& rStrainVector)
{
    // We update the deformation gradient
    Matrix equivalent_F(Dimension, Dimension);

    if(Dimension == 2) {
        equivalent_F(0,0) = 1.0 + rStrainVector[0];
        equivalent_F(0,1) = 0.5 * rStrainVector[2];
        equivalent_F(1,0) = 0.5 * rStrainVector[2];
        equivalent_F(1,1) = 1.0 + rStrainVector[1];
    } else {
        equivalent_F(0,0) = 1.0 + rStrainVector[0];
        equivalent_F(0,1) = 0.5 * rStrainVector[3];
        equivalent_F(0,2) = 0.5 * rStrainVector[5];
        equivalent_F(1,0) = 0.5 * rStrainVector[3];
        equivalent_F(1,1) = 1.0 + rStrainVector[1];
        equivalent_F(1,2) = 0.5 * rStrainVector[4];
        equivalent_F(2,0) = 0.5 * rStrainVector[5];
        equivalent_F(2,1) = 0.5 * rStrainVector[4];
        equivalent_F(2,2) = 1.0 + rStrainVector[2];
    }

    return equivalent_F;
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

    // Compute square root matrix
    BoundedMatrixType E_matrix;
    MathUtils<double>::MatrixSquareRoot(rCauchyTensor, E_matrix, 1.0e-16, 20);

    // Biot Strain Calculation
    rStrainVector = MathUtils<double>::StrainTensorToVector(E_matrix, TVoigtSize);
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TVoigtSize>
void ConstitutiveLawUtilities<TVoigtSize>::PolarDecomposition(
    const MatrixType& rFDeformationGradient,
    MatrixType& rRMatrix,
    MatrixType& rUMatrix
    )
{
    // Doing resize in case is needed
    if (rRMatrix.size1() != Dimension || rRMatrix.size2() != Dimension)
        rRMatrix.resize(Dimension, Dimension, false);
    if (rUMatrix.size1() != Dimension || rUMatrix.size2() != Dimension)
        rUMatrix.resize(Dimension, Dimension, false);

    // We compute Right Cauchy tensor
    const MatrixType C = prod(trans(rFDeformationGradient), rFDeformationGradient);

    // Decompose matrix C
    BoundedMatrix<double, Dimension, Dimension> eigen_vector_matrix, eigen_values_matrix;
    MathUtils<double>::GaussSeidelEigenSystem(C, eigen_vector_matrix, eigen_values_matrix, 1.0e-16, 200);

    for (IndexType i = 0; i < Dimension; ++i)
        eigen_values_matrix(i, i) = std::sqrt(eigen_values_matrix(i, i));

    noalias(rUMatrix) = prod(eigen_values_matrix, trans(eigen_vector_matrix));

    double aux_det;
    MatrixType invU(Dimension, Dimension);
    MathUtils<double>::InvertMatrix(rUMatrix, invU, aux_det);
    noalias(rRMatrix) = prod(rFDeformationGradient, invU);
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

template<>
void ConstitutiveLawUtilities<6>::CalculateProjectionOperator(
    const Vector& rStrainVector,
    MatrixType& rProjectionOperatorTensor
    )
{
    BoundedMatrix<double, Dimension, Dimension> strain_tensor;
    strain_tensor = MathUtils<double>::StrainVectorToTensor(rStrainVector);
    BoundedMatrix<double, Dimension, Dimension> eigen_vectors_matrix;
    BoundedMatrix<double, Dimension, Dimension> eigen_values_matrix;

    MathUtils<double>::GaussSeidelEigenSystem(strain_tensor, eigen_vectors_matrix, eigen_values_matrix, 1.0e-16, 20);

    std::vector<Vector> eigen_vectors_container;

    Vector auxiliar_vector = ZeroVector(Dimension);
    for (IndexType i = 0; i < Dimension; ++i) {
        auxiliar_vector[0] = eigen_vectors_matrix(i, 0);
        auxiliar_vector[1] = eigen_vectors_matrix(i, 1);
        auxiliar_vector[2] = eigen_vectors_matrix(i, 2);
        eigen_vectors_container.push_back(auxiliar_vector);
    }

    Vector sigma_tension_vector;
    Matrix sigma_tension_tensor;
    for (IndexType i = 0; i < Dimension; ++i) {
        if (eigen_values_matrix(i, i) > 0.0) {
            sigma_tension_tensor = outer_prod(eigen_vectors_container[i], eigen_vectors_container[i]); // p_i x p_i
            sigma_tension_vector = MathUtils<double>::StressTensorToVector(sigma_tension_tensor);
            rProjectionOperatorTensor += outer_prod(sigma_tension_vector, sigma_tension_vector);
        }
    }

    Matrix indexes_ij;
    indexes_ij.resize(3, 2, false);
    indexes_ij(0, 0) = 0;
    indexes_ij(0, 1) = 1;
    indexes_ij(1, 0) = 1;
    indexes_ij(1, 1) = 2;
    indexes_ij(2, 0) = 0;
    indexes_ij(2, 1) = 2;

    IndexType i, j;
    double h_i = 0.0, h_j = 0.0;
    Matrix cross_p_ij_tensor;
    Vector cross_p_ij_vector;

    for (IndexType index = 0; index < Dimension; ++index) {
        i = indexes_ij(index, 0);
        j = indexes_ij(index, 1);

        if (eigen_values_matrix(i, i) > 0)
            h_i = 1.0;
        if (eigen_values_matrix(j, j) > 0)
            h_j = 1.0;

        cross_p_ij_tensor = 0.5 * (outer_prod(eigen_vectors_container[i], eigen_vectors_container[j]) +
                                   outer_prod(eigen_vectors_container[j], eigen_vectors_container[i]));
        cross_p_ij_vector = MathUtils<double>::StressTensorToVector(cross_p_ij_tensor);
        rProjectionOperatorTensor += (h_i + h_j) * (outer_prod(cross_p_ij_vector, cross_p_ij_vector));

        h_i = 0.0;
        h_j = 0.0;
    }
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void ConstitutiveLawUtilities<3>::CalculateProjectionOperator(
    const Vector& rStrainVector,
    MatrixType& rProjectionOperatorTensor
    )
{
    BoundedMatrix<double, Dimension, Dimension> strain_tensor;
    strain_tensor = MathUtils<double>::StrainVectorToTensor(rStrainVector);
    BoundedMatrix<double, Dimension, Dimension> eigen_vectors_matrix;
    BoundedMatrix<double, Dimension, Dimension> eigen_values_matrix;

    MathUtils<double>::GaussSeidelEigenSystem(strain_tensor, eigen_vectors_matrix, eigen_values_matrix, 1.0e-16, 20);

    std::vector<Vector> eigen_vectors_container;

    Vector auxiliar_vector = ZeroVector(Dimension);
    for (IndexType i = 0; i < Dimension; ++i) {
        auxiliar_vector[0] = eigen_vectors_matrix(i, 0);
        auxiliar_vector[1] = eigen_vectors_matrix(i, 1);
        eigen_vectors_container.push_back(auxiliar_vector);
    }

    Vector sigma_tension_vector;
    Matrix sigma_tension_tensor;
    for (IndexType i = 0; i < Dimension; ++i) {
        if (eigen_values_matrix(i, i) > 0.0) {
            sigma_tension_tensor = outer_prod(eigen_vectors_container[i], eigen_vectors_container[i]); // p_i x p_i
            sigma_tension_vector = MathUtils<double>::StressTensorToVector(sigma_tension_tensor);
            rProjectionOperatorTensor += outer_prod(sigma_tension_vector, sigma_tension_vector);
        }
    }

    double h_i = 0.0, h_j = 0.0;
    Matrix cross_p_ij_tensor;
    Vector cross_p_ij_vector;

    if (eigen_values_matrix(0, 0) > 0)
        h_i = 1.0;
    if (eigen_values_matrix(1, 1) > 0)
        h_j = 1.0;

    cross_p_ij_tensor = 0.5 * (outer_prod(eigen_vectors_container[0], eigen_vectors_container[1]) +
                               outer_prod(eigen_vectors_container[1], eigen_vectors_container[0]));
    cross_p_ij_vector = MathUtils<double>::StressTensorToVector(cross_p_ij_tensor);

    rProjectionOperatorTensor += (h_i + h_j) * (outer_prod(cross_p_ij_vector, cross_p_ij_vector));
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TVoigtSize>
void ConstitutiveLawUtilities<TVoigtSize>::SpectralDecomposition(
    const BoundedVectorType& rStressVector,
    BoundedVectorType& rStressVectorTension,
    BoundedVectorType& rStressVectorCompression
    )
{
    rStressVectorTension     = ZeroVector(TVoigtSize);
    rStressVectorCompression = ZeroVector(TVoigtSize);

    BoundedMatrix<double, Dimension, Dimension> stress_tensor;
    stress_tensor = MathUtils<double>::StressVectorToTensor(rStressVector);
    BoundedMatrix<double, Dimension, Dimension> eigen_vectors_matrix;
    BoundedMatrix<double, Dimension, Dimension> eigen_values_matrix;

    MathUtils<double>::GaussSeidelEigenSystem(stress_tensor, eigen_vectors_matrix, eigen_values_matrix, 1.0e-16, 20);

    std::vector<Vector> eigen_vectors_container;
    Vector auxiliar_vector = ZeroVector(Dimension);
    for (IndexType i = 0; i < Dimension; ++i) {
        for (IndexType j = 0; j < Dimension; ++j) {
            auxiliar_vector[j] = eigen_vectors_matrix(j, i);
        }
        eigen_vectors_container.push_back(auxiliar_vector);
    }

    Vector sigma_tension_vector;
    Matrix sigma_tension_tensor;
    for (IndexType i = 0; i < Dimension; ++i) {
        if (eigen_values_matrix(i, i) > 0.0) {
            sigma_tension_tensor = eigen_values_matrix(i, i) * outer_prod(eigen_vectors_container[i], eigen_vectors_container[i]); // p_i x p_i
            rStressVectorTension += MathUtils<double>::StressTensorToVector(sigma_tension_tensor);
        }
    }
    rStressVectorCompression = rStressVector - rStressVectorTension;
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TVoigtSize>
Matrix ConstitutiveLawUtilities<TVoigtSize>::CalculateLinearPlasticDeformationGradientIncrement(
    const BoundedVectorType& rPlasticPotentialDerivative,
    const double PlasticConsistencyFactorIncrement
    )
{
    const MatrixType plastic_deformation_gradient_increment  = IdentityMatrix(Dimension) + PlasticConsistencyFactorIncrement * MathUtils<double>::StrainVectorToTensor<BoundedVectorType, MatrixType>(rPlasticPotentialDerivative);

    return plastic_deformation_gradient_increment;
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TVoigtSize>
Matrix ConstitutiveLawUtilities<TVoigtSize>::CalculateExponentialElasticDeformationGradient(
    const MatrixType& rElasticTrial,
    const BoundedVectorType& rPlasticPotentialDerivative,
    const double PlasticConsistencyFactorIncrement,
    const MatrixType& rRe
    )
{
    // Define elastic deformation gradient
    MatrixType elastic_deformation_gradient(Dimension, Dimension);

    // Define plastic flow
    const BoundedMatrixType plastic_flow = PlasticConsistencyFactorIncrement * MathUtils<double>::StrainVectorToTensor<BoundedVectorType, MatrixType>(rPlasticPotentialDerivative);

    // Declare the different eigen decomposition matrices
    BoundedMatrixType eigen_values_matrix, eigen_vectors_matrix;

    // We compute the exponential matrix
    // Decompose matrix
    MathUtils<double>::GaussSeidelEigenSystem(plastic_flow, eigen_vectors_matrix, eigen_values_matrix, 1.0e-16, 20);

    // Calculate the eigenvalues of the E matrix
    for (std::size_t i = 0; i < Dimension; ++i) {
        eigen_values_matrix(i, i) = std::exp(-eigen_values_matrix(i, i));
    }

    // Calculate exponential matrix
    MathUtils<double>::BDBtProductOperation(elastic_deformation_gradient, eigen_values_matrix, eigen_vectors_matrix);

    // Pre and post multiply by Re
    elastic_deformation_gradient = prod(elastic_deformation_gradient, rRe);
    elastic_deformation_gradient = prod(trans(rRe), elastic_deformation_gradient);
    elastic_deformation_gradient = prod(rElasticTrial, elastic_deformation_gradient);

    return elastic_deformation_gradient;
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TVoigtSize>
Matrix ConstitutiveLawUtilities<TVoigtSize>::CalculateDirectElasticDeformationGradient(
    const MatrixType& rElasticTrial,
    const BoundedVectorType& rPlasticPotentialDerivative,
    const double PlasticConsistencyFactorIncrement,
    const MatrixType& rRe
    )
{
    // Define elastic deformation gradient
    MatrixType elastic_deformation_gradient(Dimension, Dimension);
    MatrixType auxiliar_deformation_gradient_increment(Dimension, Dimension);
    MatrixType inverse_plastic_deformation_gradient_increment(Dimension, Dimension);

    // Define plastic flow
    const BoundedMatrixType inverse_plastic_flow = - PlasticConsistencyFactorIncrement * MathUtils<double>::StrainVectorToTensor<BoundedVectorType, MatrixType>(rPlasticPotentialDerivative);

    // Pre and post multiply by Re
    noalias(auxiliar_deformation_gradient_increment) = prod(inverse_plastic_flow, rRe);
    noalias(auxiliar_deformation_gradient_increment) = prod(trans(rRe), inverse_plastic_flow);

    auxiliar_deformation_gradient_increment = IdentityMatrix(Dimension) - auxiliar_deformation_gradient_increment;

    double aux_det;
    MathUtils<double>::InvertMatrix(auxiliar_deformation_gradient_increment, inverse_plastic_deformation_gradient_increment, aux_det);

    // Pre and post multiply by Re
    noalias(elastic_deformation_gradient) = prod(rElasticTrial, inverse_plastic_deformation_gradient_increment);

    return elastic_deformation_gradient;
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TVoigtSize>
Matrix ConstitutiveLawUtilities<TVoigtSize>::CalculateExponentialPlasticDeformationGradientIncrement(
    const BoundedVectorType& rPlasticPotentialDerivative,
    const double PlasticConsistencyFactorIncrement,
    const MatrixType& rRe
    )
{
    // Define DeltaFp
    MatrixType plastic_deformation_gradient_increment(Dimension, Dimension);

    // Define plastic flow
    const BoundedMatrixType plastic_flow = PlasticConsistencyFactorIncrement * MathUtils<double>::StrainVectorToTensor<BoundedVectorType, MatrixType>(rPlasticPotentialDerivative);

    // Declare the different eigen decomposition matrices
    BoundedMatrixType eigen_values_matrix, eigen_vectors_matrix;

    // We compute the exponential matrix
    // Decompose matrix
    MathUtils<double>::GaussSeidelEigenSystem(plastic_flow, eigen_vectors_matrix, eigen_values_matrix, 1.0e-16, 20);

    // Calculate the eigenvalues of the E matrix
    for (std::size_t i = 0; i < Dimension; ++i) {
        eigen_values_matrix(i, i) = std::exp(eigen_values_matrix(i, i));
    }

    // Calculate exponential matrix
    MathUtils<double>::BDBtProductOperation(plastic_deformation_gradient_increment, eigen_values_matrix, eigen_vectors_matrix);

    // Pre and post multiply by Re
    plastic_deformation_gradient_increment = prod(plastic_deformation_gradient_increment, rRe);
    plastic_deformation_gradient_increment = prod(trans(rRe), plastic_deformation_gradient_increment);

    return plastic_deformation_gradient_increment;
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TVoigtSize>
Matrix ConstitutiveLawUtilities<TVoigtSize>::CalculateDirectPlasticDeformationGradientIncrement(
    const BoundedVectorType& rPlasticPotentialDerivative,
    const double PlasticConsistencyFactorIncrement,
    const MatrixType& rRe
    )
{
    // Define DeltaFp
    MatrixType auxiliar_deformation_gradient_increment(Dimension, Dimension);
    MatrixType plastic_deformation_gradient_increment(Dimension, Dimension);

    // Define plastic flow
    const BoundedMatrixType plastic_flow = PlasticConsistencyFactorIncrement * MathUtils<double>::StrainVectorToTensor<BoundedVectorType, MatrixType>(rPlasticPotentialDerivative);

    // Pre and post multiply by Re
    auxiliar_deformation_gradient_increment = prod(plastic_flow, rRe);
    auxiliar_deformation_gradient_increment = prod(trans(rRe), plastic_flow);

    auxiliar_deformation_gradient_increment = IdentityMatrix(Dimension) - auxiliar_deformation_gradient_increment;

    double aux_det;
    MathUtils<double>::InvertMatrix(auxiliar_deformation_gradient_increment, plastic_deformation_gradient_increment, aux_det);

    return plastic_deformation_gradient_increment;
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TVoigtSize>
void ConstitutiveLawUtilities<TVoigtSize>::CalculateRotationOperatorEuler1(
    const double EulerAngle1,
    BoundedMatrix<double, 3, 3>& rRotationOperator
)
{
    noalias(rRotationOperator) = ZeroMatrix(Dimension, Dimension);

    const double cos_angle = std::cos(EulerAngle1 * Globals::Pi / 180.0);
    const double sin_angle = std::sin(EulerAngle1 * Globals::Pi / 180.0);

    rRotationOperator(0, 0) = cos_angle;
    rRotationOperator(0, 1) = sin_angle;
    rRotationOperator(1, 0) = -sin_angle;
    rRotationOperator(1, 1) = cos_angle;
    rRotationOperator(2, 2) = 1.0;
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TVoigtSize>
void ConstitutiveLawUtilities<TVoigtSize>::CalculateRotationOperatorEuler2(
    const double EulerAngle2,
    BoundedMatrix<double, 3, 3>& rRotationOperator
)
{
    noalias(rRotationOperator) = ZeroMatrix(Dimension, Dimension);

    const double cos_angle = std::cos(EulerAngle2 * Globals::Pi / 180.0);
    const double sin_angle = std::sin(EulerAngle2 * Globals::Pi / 180.0);

    rRotationOperator(0, 0) = 1.0;
    rRotationOperator(1, 1) = cos_angle;
    rRotationOperator(1, 2) = sin_angle;
    rRotationOperator(2, 1) = -sin_angle;
    rRotationOperator(2, 2) = cos_angle;
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TVoigtSize>
void ConstitutiveLawUtilities<TVoigtSize>::CalculateRotationOperatorEuler3(
    const double EulerAngle3,
    BoundedMatrix<double, 3, 3>& rRotationOperator
)
{
    ConstitutiveLawUtilities<TVoigtSize>::CalculateRotationOperatorEuler1(EulerAngle3, rRotationOperator);
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TVoigtSize>
void ConstitutiveLawUtilities<TVoigtSize>::CalculateRotationOperator(
    const double EulerAngle1, // phi
    const double EulerAngle2, // theta
    const double EulerAngle3, // hi
    BoundedMatrix<double, 3, 3>& rRotationOperator // global to local coordinates
)
{
    noalias(rRotationOperator) = ZeroMatrix(Dimension, Dimension);

    const double pi_over_180 = Globals::Pi / 180.0;
    const double cos1 = std::cos(EulerAngle1 * pi_over_180);
    const double sin1 = std::sin(EulerAngle1 * pi_over_180);
    const double cos2 = std::cos(EulerAngle2 * pi_over_180);
    const double sin2 = std::sin(EulerAngle2 * pi_over_180);
    const double cos3 = std::cos(EulerAngle3 * pi_over_180);
    const double sin3 = std::sin(EulerAngle3 * pi_over_180);

    rRotationOperator(0, 0) = cos1 * cos3 - sin1 * cos2 * sin3;
    rRotationOperator(0, 1) = sin1 * cos3 + cos1 * cos2 * sin3;
    rRotationOperator(0, 2) = sin2 * sin3;
    rRotationOperator(1, 0) = -cos1 * sin3 - sin1 * cos2 * cos3;
    rRotationOperator(1, 1) = -sin1 * sin3 + cos1 * cos2 * cos3;
    rRotationOperator(1, 2) = sin2 * cos3;
    rRotationOperator(2, 0) = sin1 * sin2;
    rRotationOperator(2, 1) = -cos1 * sin2;
    rRotationOperator(2, 2) = cos2;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void ConstitutiveLawUtilities<6>::CalculateRotationOperatorVoigt(
    const BoundedMatrixType& rEulerOperator,
    BoundedMatrixVoigtType& rVoigtOperator
    )
{
    rVoigtOperator.clear();

    const double l1 = rEulerOperator(0, 0);
    const double l2 = rEulerOperator(1, 0);
    const double l3 = rEulerOperator(2, 0);
    const double m1 = rEulerOperator(0, 1);
    const double m2 = rEulerOperator(1, 1);
    const double m3 = rEulerOperator(2, 1);
    const double n1 = rEulerOperator(0, 2);
    const double n2 = rEulerOperator(1, 2);
    const double n3 = rEulerOperator(2, 2);

    rVoigtOperator(0, 0) = std::pow(l1, 2);
    rVoigtOperator(0, 1) = std::pow(m1, 2);
    rVoigtOperator(0, 2) = std::pow(n1, 2);
    rVoigtOperator(0, 3) = l1 * m1;
    rVoigtOperator(0, 4) = m1 * n1;
    rVoigtOperator(0, 5) = n1 * l1;

    rVoigtOperator(1, 0) = std::pow(l2, 2);
    rVoigtOperator(1, 1) = std::pow(m2, 2);
    rVoigtOperator(1, 2) = std::pow(n2, 2);
    rVoigtOperator(1, 3) = l2 * m2;
    rVoigtOperator(1, 4) = m2 * n2;
    rVoigtOperator(1, 5) = n2 * l2;

    rVoigtOperator(2, 0) = std::pow(l3, 2);
    rVoigtOperator(2, 1) = std::pow(m3, 2);
    rVoigtOperator(2, 2) = std::pow(n3, 2);
    rVoigtOperator(2, 3) = l3 * m3;
    rVoigtOperator(2, 4) = m3 * n3;
    rVoigtOperator(2, 5) = n3 * l3;

    rVoigtOperator(3, 0) = 2.0 * l1 * l2;
    rVoigtOperator(3, 1) = 2.0 * m1 * m2;
    rVoigtOperator(3, 2) = 2.0 * n1 * n2;
    rVoigtOperator(3, 3) = l1 * m2 + l2 * m1;
    rVoigtOperator(3, 4) = m1 * n2 + m2 * n1;
    rVoigtOperator(3, 5) = n1 * l2 + n2 * l1;

    rVoigtOperator(4, 0) = 2.0 * l2 * l3;
    rVoigtOperator(4, 1) = 2.0 * m2 * m3;
    rVoigtOperator(4, 2) = 2.0 * n2 * n3;
    rVoigtOperator(4, 3) = l2 * m3 + l3 * m2;
    rVoigtOperator(4, 4) = m2 * n3 + m3 * n2;
    rVoigtOperator(4, 5) = n2 * l3 + n3 * l2;

    rVoigtOperator(5, 0) = 2.0 * l3 * l1;
    rVoigtOperator(5, 1) = 2.0 * m3 * m1;
    rVoigtOperator(5, 2) = 2.0 * n3 * n1;
    rVoigtOperator(5, 3) = l3 * m1 + l1 * m3;
    rVoigtOperator(5, 4) = m3 * n1 + m1 * n3;
    rVoigtOperator(5, 5) = n3 * l1 + n1 * l3;
}

/***********************************************************************************/
/***********************************************************************************/

template<>
void ConstitutiveLawUtilities<3>::CalculateRotationOperatorVoigt(
    const BoundedMatrixType& rEulerOperator,
    BoundedMatrixVoigtType& rVoigtOperator
    )
{
    const double c = rEulerOperator(0, 0);
    const double s = rEulerOperator(0, 1);

    rVoigtOperator(0, 0) = std::pow(c, 2);
    rVoigtOperator(0, 1) = std::pow(s, 2);
    rVoigtOperator(0, 2) = c * s;

    rVoigtOperator(1, 0) = std::pow(s, 2);
    rVoigtOperator(1, 1) = std::pow(c, 2);
    rVoigtOperator(1, 2) = -c * s;

    rVoigtOperator(2, 0) = -2.0 * c * s;
    rVoigtOperator(2, 1) = 2.0 * c * s;
    rVoigtOperator(2, 2) = std::pow(c, 2) - std::pow(s, 2);
}

/***********************************************************************************/
/***********************************************************************************/

template class ConstitutiveLawUtilities<3>;
template class ConstitutiveLawUtilities<6>;

} // namespace Kratos
