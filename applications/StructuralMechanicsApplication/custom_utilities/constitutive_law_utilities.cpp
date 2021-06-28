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

template class ConstitutiveLawUtilities<3>;
template class ConstitutiveLawUtilities<6>;

} // namespace Kratos
