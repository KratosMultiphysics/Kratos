// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: StructuralMechanicsApplication/license.txt
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

    Vector auxiliary_vector = ZeroVector(Dimension);
    for (IndexType i = 0; i < Dimension; ++i) {
        auxiliary_vector[0] = eigen_vectors_matrix(i, 0);
        auxiliary_vector[1] = eigen_vectors_matrix(i, 1);
        auxiliary_vector[2] = eigen_vectors_matrix(i, 2);
        eigen_vectors_container.push_back(auxiliary_vector);
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

    Vector auxiliary_vector = ZeroVector(Dimension);
    for (IndexType i = 0; i < Dimension; ++i) {
        auxiliary_vector[0] = eigen_vectors_matrix(i, 0);
        auxiliary_vector[1] = eigen_vectors_matrix(i, 1);
        eigen_vectors_container.push_back(auxiliary_vector);
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
void ConstitutiveLawUtilities<TVoigtSize>::CalculateRotationOperatorVoigt(
    const BoundedMatrixType& rEulerOperator,
    BoundedMatrixVoigtType& rVoigtOperator
    )
{
    if constexpr (Dimension == 2) {
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
    } else {
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
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TVoigtSize>
void ConstitutiveLawUtilities<TVoigtSize>::CalculateElasticMatrixPlaneStress(
    MatrixType& rC,
    const double YoungModulus,
    const double PoissonRatio)
{
    if (rC.size1() != VoigtSize || rC.size2() != VoigtSize)
        rC.resize(VoigtSize, VoigtSize, false);
    rC.clear();

    const double c1 = YoungModulus / (1.0 - PoissonRatio * PoissonRatio);
    const double c2 = c1 * PoissonRatio;
    const double c3 = 0.5 * YoungModulus / (1.0 + PoissonRatio);

    rC(0, 0) = c1;
    rC(0, 1) = c2;
    rC(1, 0) = c2;
    rC(1, 1) = c1;
    rC(2, 2) = c3;
}

/***********************************************************************************/
/***********************************************************************************/
template<SizeType TVoigtSize>
void ConstitutiveLawUtilities<TVoigtSize>::CalculateElasticMatrixPlaneStrain(
    MatrixType& rC,
    const double YoungModulus,
    const double PoissonRatio)
{
    if (rC.size1() != VoigtSize || rC.size2() != VoigtSize)
        rC.resize(VoigtSize, VoigtSize, false);
    rC.clear();

    const double c0 = YoungModulus / ((1.0 + PoissonRatio) * (1.0 - 2.0 * PoissonRatio));
    const double c1 = (1.0 - PoissonRatio) * c0;
    const double c2 = c0 * PoissonRatio;
    const double c3 = (0.5 - PoissonRatio) * c0;

    rC(0, 0) = c1;
    rC(0, 1) = c2;
    rC(1, 0) = c2;
    rC(1, 1) = c1;
    rC(2, 2) = c3;
}

/***********************************************************************************/
/***********************************************************************************/
template<SizeType TVoigtSize>
void ConstitutiveLawUtilities<TVoigtSize>::CalculateElasticMatrix(
    MatrixType& rC,
    const double YoungModulus,
    const double PoissonRatio)
{
    if (rC.size1() != VoigtSize || rC.size2() != VoigtSize)
        rC.resize(VoigtSize, VoigtSize, false);
    rC.clear();

    const double c1 = YoungModulus / ((1.0 + PoissonRatio) * (1.0 - 2.0 * PoissonRatio));
    const double c2 = c1 * (1.0 - PoissonRatio);
    const double c3 = c1 * PoissonRatio;
    const double c4 = c1 * 0.5 * (1.0 - 2.0 * PoissonRatio);

    rC(0, 0) = c2;
    rC(0, 1) = c3;
    rC(0, 2) = c3;
    rC(1, 0) = c3;
    rC(1, 1) = c2;
    rC(1, 2) = c3;
    rC(2, 0) = c3;
    rC(2, 1) = c3;
    rC(2, 2) = c2;
    rC(3, 3) = c4;
    rC(4, 4) = c4;
    rC(5, 5) = c4;
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TVoigtSize>
void ConstitutiveLawUtilities<TVoigtSize>::CalculateCauchyGreenStrain(
    ConstitutiveLaw::Parameters& rValues,
    Vector& rStrainVector
    )
{
    if constexpr (Dimension == 2) {
        //1.-Compute total deformation gradient
        const ConstitutiveLaw::DeformationGradientMatrixType& F = rValues.GetDeformationGradientF();

        // for shells/membranes in case the DeformationGradient is of size 3x3
        BoundedMatrix<double,2,2> F2x2;
        for (unsigned int i = 0; i<2; ++i)
            for (unsigned int j = 0; j<2; ++j)
                F2x2(i, j) = F(i, j);

        BoundedMatrix<double,2,2> E_tensor = prod(trans(F2x2), F2x2);

        for (unsigned int i = 0; i<2; ++i)
            E_tensor(i, i) -= 1.0;

        E_tensor *= 0.5;

        noalias(rStrainVector) = MathUtils<double>::StrainTensorToVector(E_tensor);
    } else {
        const SizeType space_dimension = rValues.GetElementGeometry().WorkingSpaceDimension();
        //1.-Compute total deformation gradient
        const ConstitutiveLaw::DeformationGradientMatrixType& F = rValues.GetDeformationGradientF();
        KRATOS_DEBUG_ERROR_IF(F.size1()!= space_dimension || F.size2() != space_dimension)
            << "expected size of F " << space_dimension << "x" << space_dimension << ", got " << F.size1() << "x" << F.size2() << std::endl;

        ConstitutiveLaw::DeformationGradientMatrixType E_tensor = prod(trans(F),F);
        for(unsigned int i = 0; i < space_dimension; ++i)
            E_tensor(i,i) -= 1.0;
        E_tensor *= 0.5;

        noalias(rStrainVector) = MathUtils<double>::StrainTensorToVector(E_tensor);
    }

}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TVoigtSize>
void ConstitutiveLawUtilities<TVoigtSize>::CalculatePK2StressFromStrain(
    ConstitutiveLaw::StrainVectorType& rStressVector,
    const ConstitutiveLaw::StressVectorType &rStrainVector,
    const double E,
    const double NU
    )
{
    const double c1 = E / ((1.0 + NU) * (1.0 - 2.0 * NU));
    const double c2 = c1 * (1 - NU);
    const double c3 = c1 * NU;
    const double c4 = c1 * 0.5 * (1 - 2 * NU);

    rStressVector[0] = c2 * rStrainVector[0] + c3 * rStrainVector[1] + c3 * rStrainVector[2];
    rStressVector[1] = c3 * rStrainVector[0] + c2 * rStrainVector[1] + c3 * rStrainVector[2];
    rStressVector[2] = c3 * rStrainVector[0] + c3 * rStrainVector[1] + c2 * rStrainVector[2];
    rStressVector[3] = c4 * rStrainVector[3];
    rStressVector[4] = c4 * rStrainVector[4];
    rStressVector[5] = c4 * rStrainVector[5];
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TVoigtSize>
void ConstitutiveLawUtilities<TVoigtSize>::CalculatePK2StressFromStrainPlaneStress(
    ConstitutiveLaw::StrainVectorType& rStressVector,
    const ConstitutiveLaw::StressVectorType &rStrainVector,
    const double E,
    const double NU
    )
{
    const double c1 = E / (1.0 - NU * NU);
    const double c2 = c1 * NU;
    const double c3 = 0.5 * E / (1.0 + NU);

    rStressVector[0] = c1 * rStrainVector[0] + c2 * rStrainVector[1];
    rStressVector[1] = c2 * rStrainVector[0] + c1 * rStrainVector[1];
    rStressVector[2] = c3 * rStrainVector[2];
}

/***********************************************************************************/
/***********************************************************************************/

template<SizeType TVoigtSize>
void ConstitutiveLawUtilities<TVoigtSize>::CalculatePK2StressFromStrainPlaneStrain(
    ConstitutiveLaw::StrainVectorType& rStressVector,
    const ConstitutiveLaw::StressVectorType &rStrainVector,
    const double E,
    const double NU
    )
{
    const double c0 = E / ((1.0 + NU)*(1.0 - 2.0 * NU));
    const double c1 = (1.0 - NU)*c0;
    const double c2 = c0 * NU;
    const double c3 = (0.5 - NU)*c0;

    rStressVector[0] = c1 * rStrainVector[0] + c2 * rStrainVector[1];
    rStressVector[1] = c2 * rStrainVector[0] + c1 * rStrainVector[1];
    rStressVector[2] = c3 * rStrainVector[2];
}

/***********************************************************************************/
/***********************************************************************************/

template class ConstitutiveLawUtilities<3>;
template class ConstitutiveLawUtilities<6>;

} // namespace Kratos
