// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Richard Faasse,
//                   Gennady Markelov
//

#include "stress_strain_utilities.h"
#include "geo_mechanics_application_constants.h"
#include "utilities/math_utils.h"
#include <cmath>

namespace Kratos
{

double StressStrainUtilities::CalculateVonMisesStress(const Vector& rStressVector)
{
    const Matrix LocalStressTensor =
        MathUtils<double>::StressVectorToTensor(rStressVector); // reduced dimension stress tensor

    Matrix StressTensor(3, 3); // 3D stress tensor
    noalias(StressTensor) = ZeroMatrix(3, 3);
    for (std::size_t i = 0; i < LocalStressTensor.size1(); ++i) {
        for (std::size_t j = 0; j < LocalStressTensor.size2(); ++j) {
            StressTensor(i, j) = LocalStressTensor(i, j);
        }
    }

    const double SigmaEquivalent =
        0.5 * ((StressTensor(0, 0) - StressTensor(1, 1)) * (StressTensor(0, 0) - StressTensor(1, 1)) +
               (StressTensor(1, 1) - StressTensor(2, 2)) * (StressTensor(1, 1) - StressTensor(2, 2)) +
               (StressTensor(2, 2) - StressTensor(0, 0)) * (StressTensor(2, 2) - StressTensor(0, 0)) +
               6.0 * (StressTensor(0, 1) * StressTensor(1, 0) + StressTensor(1, 2) * StressTensor(2, 1) +
                      StressTensor(2, 0) * StressTensor(0, 2)));

    return std::sqrt(std::max(SigmaEquivalent, 0.));
}

double StressStrainUtilities::CalculateTrace(const Vector& rStressVector)
{
    const Matrix StressTensor =
        MathUtils<double>::StressVectorToTensor(rStressVector); // reduced dimension stress tensor

    double trace = 0.0;
    for (std::size_t i = 0; i < StressTensor.size1(); ++i) {
        trace += StressTensor(i, i);
    }

    return trace;
}

double StressStrainUtilities::CalculateMeanStress(const Vector& rStressVector)
{
    return CalculateTrace(rStressVector) / (rStressVector.size() == 3 ? 2.0 : 3.0);
}

double StressStrainUtilities::CalculateLodeAngle(const Vector& rStressVector)
{
    KRATOS_ERROR_IF(rStressVector.size() < 4);

    const double p                   = CalculateMeanStress(rStressVector);
    const double q                   = CalculateVonMisesStress(rStressVector);
    const Matrix local_stress_tensor = MathUtils<double>::StressVectorToTensor(rStressVector);
    Matrix       sigma_princi;
    Matrix       eigen_vectors;
    MathUtils<double>::GaussSeidelEigenSystem(local_stress_tensor, eigen_vectors, sigma_princi, 1.0e-16, 20);
    const double numerator = (sigma_princi(0, 0) - p) * (sigma_princi(1, 1) - p) * (sigma_princi(2, 2) - p);
    if (std::abs(numerator) < 1.0E-12) return 0.;
    return std::asin((-27. / 2.) * numerator / (q * q * q)) / 3.0;
}

double StressStrainUtilities::CalculateMohrCoulombShearCapacity(const Vector& rStressVector, double C, double Phi)
{
    KRATOS_ERROR_IF(rStressVector.size() < 4);

    const double q_mc = CalculateQMohrCoulomb(rStressVector, C, Phi);
    const double q    = CalculateVonMisesStress(rStressVector);

    return q / q_mc;
}

double StressStrainUtilities::CalculateQMohrCoulomb(const Vector& rStressVector, double C, double Phi)
{
    const double denominator = CalculateDenominator(rStressVector, Phi);
    const double p           = -CalculateMeanStress(rStressVector);
    return 3. * (p * std::sin(Phi) + C * std::cos(Phi)) / denominator;
}

double StressStrainUtilities::CalculateDenominator(const Vector& rStressVector, double Phi)
{
    const double lode_angle = CalculateLodeAngle(rStressVector);
    return std::sqrt(3.) * std::cos(lode_angle) - std::sin(lode_angle) * std::sin(Phi);
}

double StressStrainUtilities::CalculateMohrCoulombPressureCapacity(const Vector& rStressVector, double C, double Phi)
{
    KRATOS_ERROR_IF(rStressVector.size() < 4);

    const double denominator = CalculateDenominator(rStressVector, Phi);
    const double q_mc        = CalculateQMohrCoulomb(rStressVector, C, Phi);
    const double q           = CalculateVonMisesStress(rStressVector);

    return 3. * std::sin(Phi) * (q_mc - q) / denominator;
}

double StressStrainUtilities::CalculateVonMisesStrain(const Vector& rStrainVector)
{
    return (2.0 / 3.0) * CalculateVonMisesStress(rStrainVector);
}

Vector StressStrainUtilities::CalculateHenckyStrain(const Matrix& rDeformationGradient, size_t VoigtSize)
{
    // right Cauchy Green deformation tensor C
    Matrix C = prod(trans(rDeformationGradient), rDeformationGradient);
    // Eigenvalues of C matrix, so principal right Cauchy Green deformation tensor C
    Matrix EigenValuesMatrix;
    Matrix EigenVectorsMatrix;
    MathUtils<double>::GaussSeidelEigenSystem(C, EigenVectorsMatrix, EigenValuesMatrix, 1.0e-16, 20);
    // Compute natural strain == Logarithmic strain == Hencky strain from principal strains
    for (std::size_t i = 0; i < rDeformationGradient.size1(); ++i) {
        EigenValuesMatrix(i, i) = 0.5 * std::log(EigenValuesMatrix(i, i));
    }

    // Rotate from principal strains back to the used coordinate system
    Matrix ETensor;
    MathUtils<double>::BDBtProductOperation(ETensor, EigenValuesMatrix, EigenVectorsMatrix);

    // From tensor to vector
    if (rDeformationGradient.size1() == 2 && VoigtSize == 4) {
        // Plane strain
        Vector StrainVector2D;
        StrainVector2D = MathUtils<double>::StrainTensorToVector(ETensor, 3);
        Vector StrainVector(4);
        StrainVector[INDEX_2D_PLANE_STRAIN_XX] = StrainVector2D[0];
        StrainVector[INDEX_2D_PLANE_STRAIN_YY] = StrainVector2D[1];
        StrainVector[INDEX_2D_PLANE_STRAIN_ZZ] = 0.0;
        StrainVector[INDEX_2D_PLANE_STRAIN_XY] = StrainVector2D[2];
        return StrainVector;
    } else {
        return MathUtils<double>::StrainTensorToVector(ETensor, VoigtSize);
    }
}

Matrix StressStrainUtilities::CalculateGreenLagrangeStrainTensor(const Matrix& rDeformationGradient)
{
    return 0.5 * (prod(trans(rDeformationGradient), rDeformationGradient) -
                  IdentityMatrix(rDeformationGradient.size1()));
}

Vector StressStrainUtilities::CalculateCauchyStrain(const Matrix& rB, const Vector& rDisplacements)
{
    return prod(rB, rDisplacements);
}

std::vector<Vector> StressStrainUtilities::CalculateStrains(const std::vector<Matrix>& rDeformationGradients,
                                                            const std::vector<Matrix>& rBs,
                                                            const Vector& rDisplacements,
                                                            bool          UseHenckyStrain,
                                                            std::size_t   VoigtSize)
{
    std::vector<Vector> result;
    result.reserve(rDeformationGradients.size());
    std::transform(
        rDeformationGradients.begin(), rDeformationGradients.end(), rBs.begin(), std::back_inserter(result),
        [&rDisplacements, UseHenckyStrain, VoigtSize](const auto& rDeformationGradient, const auto& rB) {
        return UseHenckyStrain ? CalculateHenckyStrain(rDeformationGradient, VoigtSize)
                               : CalculateCauchyStrain(rB, rDisplacements);
    });

    return result;
}

void StressStrainUtilities::CalculatePrincipalStresses(const Vector& rCauchyStressVector,
                                                       Vector&       rPrincipalStressVector,
                                                       Matrix&       rEigenVectorsMatrix)
{
    auto   stress_tensor = MathUtils<double>::StressVectorToTensor(rCauchyStressVector);
    Matrix PrincipalStressMatrix;
    Matrix EigenVectorsMatrix;
    MathUtils<double>::GaussSeidelEigenSystem(stress_tensor, rEigenVectorsMatrix,
                                              PrincipalStressMatrix, 1.0e-16, 20);
    rPrincipalStressVector = ZeroVector(3);
    for (int i = 0; i < 3; ++i) {
        rPrincipalStressVector(i) = PrincipalStressMatrix(i, i);
    }
    KRATOS_INFO("CalculatePrincipalStresses")
        << "rPrincipalStressVector (1): " << rPrincipalStressVector << std::endl;
    KRATOS_INFO("CalculatePrincipalStresses")
        << "rEigenVectorsMatrix (1): " << rEigenVectorsMatrix << std::endl;
    std::vector<std::size_t> indices(rPrincipalStressVector.size());
    std::iota(indices.begin(), indices.end(), 0);
    KRATOS_INFO("CalculatePrincipalStresses") << "indices (1): " << indices << std::endl;
    std::sort(indices.begin(), indices.end(), [&rPrincipalStressVector](const std::size_t i, const std::size_t j) {
        return rPrincipalStressVector[j] < rPrincipalStressVector[i];
    });
    KRATOS_INFO("CalculatePrincipalStresses")
        << "rPrincipalStressVector (2): " << rPrincipalStressVector << std::endl;
    KRATOS_INFO("CalculatePrincipalStresses") << "indices (2): " << indices << std::endl;

    std::vector<double> tmp_vector;
    tmp_vector.reserve(rPrincipalStressVector.size());
    for (const auto index : indices) {
        tmp_vector.push_back(rPrincipalStressVector[index]);
    }

    std::copy(tmp_vector.begin(), tmp_vector.end(), rPrincipalStressVector.begin());
    KRATOS_INFO("CalculatePrincipalStresses")
        << "rPrincipalStressVector (3): " << rPrincipalStressVector << std::endl;

    Matrix tmp_matrix(rEigenVectorsMatrix.size1(), rEigenVectorsMatrix.size2());
    for (auto i = std::size_t{0}; i < rEigenVectorsMatrix.size1(); ++i) {
        for (auto j = std::size_t{0}; j < rEigenVectorsMatrix.size2(); ++j) {
            tmp_matrix(i, j) = rEigenVectorsMatrix(i, indices[j]);
        }
    }

    noalias(rEigenVectorsMatrix) = tmp_matrix;
    KRATOS_INFO("CalculatePrincipalStresses")
        << "rEigenVectorsMatrix (3): " << rEigenVectorsMatrix << std::endl;
}

} // namespace Kratos
