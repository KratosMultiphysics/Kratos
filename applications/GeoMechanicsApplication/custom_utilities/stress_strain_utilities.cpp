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
#include "custom_constitutive/geo_principal_stresses.hpp"
#include "custom_constitutive/geo_sigma_tau.h"
#include "custom_utilities/generic_utilities.hpp"
#include "custom_utilities/math_utilities.hpp"
#include "custom_utilities/ublas_utilities.h"
#include "geo_mechanics_application_constants.h"
#include "includes/process_info.h"
#include "utilities/math_utils.h"

namespace Kratos
{

double StressStrainUtilities::CalculateVonMisesStress(const Vector& rStressVector)
{
    const auto LocalStressTensor =
        MathUtils<double>::StressVectorToTensor(rStressVector); // reduced dimension stress tensor

    Matrix StressTensor(3, 3); // 3D stress tensor
    noalias(StressTensor) = ZeroMatrix(3, 3);
    for (std::size_t i = 0; i < LocalStressTensor.size1(); ++i) {
        for (std::size_t j = 0; j < LocalStressTensor.size2(); ++j) {
            StressTensor(i, j) = LocalStressTensor(i, j);
        }
    }

    const auto SigmaEquivalent =
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

    const auto q = CalculateVonMisesStress(rStressVector);
    if (constexpr auto tolerance = 1.0E-12; q < tolerance) return 0.0;

    Matrix sigma_princi;
    Matrix eigen_vectors;
    MathUtils<>::GaussSeidelEigenSystem(MathUtils<>::StressVectorToTensor(rStressVector),
                                        eigen_vectors, sigma_princi, 1.0e-16, 20);
    const auto p = CalculateMeanStress(rStressVector);
    const auto numerator = (sigma_princi(0, 0) - p) * (sigma_princi(1, 1) - p) * (sigma_princi(2, 2) - p);
    // Avoid a domain error when computing the arc sine (which results in a NaN value)
    const auto arg_to_asin = std::clamp((-27. / 2.) * numerator / (q * q * q), -1.0, 1.0);

    return std::asin(arg_to_asin) / 3.0;
}

double StressStrainUtilities::CalculateMohrCoulombShearCapacity(const Vector& rStressVector, double C, double PhiInRadians)
{
    KRATOS_ERROR_IF(rStressVector.size() < 4);
    KRATOS_ERROR_IF(PhiInRadians < 0.0 || PhiInRadians > Globals::Pi / 2.0)
        << "Friction angle must be in the range [0, 90] (degrees) : " << PhiInRadians * 180.0 / Globals::Pi
        << std::endl;
    const auto q_mc = CalculateQMohrCoulomb(rStressVector, C, PhiInRadians);
    if (q_mc < 1e-10) return 1.0;
    const auto q = CalculateVonMisesStress(rStressVector);

    return q / q_mc;
}

double StressStrainUtilities::CalculateQMohrCoulomb(const Vector& rStressVector, double C, double PhiInRadians)
{
    const auto denominator = CalculateDenominator(rStressVector, PhiInRadians);
    const auto p           = -CalculateMeanStress(rStressVector);
    return 3. * (p * std::sin(PhiInRadians) + C * std::cos(PhiInRadians)) / denominator;
}

double StressStrainUtilities::CalculateDenominator(const Vector& rStressVector, double PhiInRadians)
{
    const auto lode_angle = CalculateLodeAngle(rStressVector);
    return std::sqrt(3.) * std::cos(lode_angle) - std::sin(lode_angle) * std::sin(PhiInRadians);
}

double StressStrainUtilities::CalculateMohrCoulombPressureCapacity(const Vector& rStressVector, double C, double PhiInRadians)
{
    KRATOS_ERROR_IF(rStressVector.size() < 4);

    const auto denominator = CalculateDenominator(rStressVector, PhiInRadians);
    const auto q_mc        = CalculateQMohrCoulomb(rStressVector, C, PhiInRadians);
    const auto q           = CalculateVonMisesStress(rStressVector);

    return 3. * std::sin(PhiInRadians) * (q_mc - q) / denominator;
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

Geo::PrincipalStresses StressStrainUtilities::CalculatePrincipalStresses(const Vector& rStressVector)
{
    auto principal_stresses = Vector{};
    auto eigen_vectors      = Matrix{};
    CalculatePrincipalStresses(rStressVector, principal_stresses, eigen_vectors);

    return Geo::PrincipalStresses{principal_stresses};
}

void StressStrainUtilities::CalculatePrincipalStresses(const Vector& rCauchyStressVector,
                                                       Vector&       rPrincipalStressVector,
                                                       Matrix&       rEigenVectorsMatrix)
{
    Matrix principal_stress_matrix;
    MathUtils<>::GaussSeidelEigenSystem(MathUtils<>::StressVectorToTensor(rCauchyStressVector),
                                        rEigenVectorsMatrix, principal_stress_matrix, 1.0e-16, 20);
    rPrincipalStressVector = GeoMechanicsMathUtilities::DiagonalOfMatrixToVector(principal_stress_matrix);
    ReorderEigenValuesAndVectors(rPrincipalStressVector, rEigenVectorsMatrix);
}

void StressStrainUtilities::ReorderEigenValuesAndVectors(Vector& rPrincipalStressVector, Matrix& rEigenVectorsMatrix)
{
    std::vector<std::size_t> indices(rPrincipalStressVector.size());
    std::iota(indices.begin(), indices.end(), std::size_t{0});
    std::sort(indices.begin(), indices.end(), [&rPrincipalStressVector](auto i, auto j) {
        return rPrincipalStressVector[j] < rPrincipalStressVector[i];
    });

    rPrincipalStressVector = GenericUtilities::PermutedVector(rPrincipalStressVector, indices);
    rEigenVectorsMatrix = GenericUtilities::MatrixWithPermutedColumns(rEigenVectorsMatrix, indices);
}

Vector StressStrainUtilities::RotatePrincipalStresses(const Vector& rPrincipalStressVector,
                                                      const Matrix& rRotationMatrix,
                                                      std::size_t   StressVectorSize)
{
    const auto principal_stress_matrix =
        GeoMechanicsMathUtilities::VectorToDiagonalMatrix(rPrincipalStressVector);
    const auto rotated_stress_matrix =
        GeoMechanicsMathUtilities::RotateSecondOrderTensor(principal_stress_matrix, rRotationMatrix);
    return MathUtils<>::StressTensorToVector(rotated_stress_matrix, StressVectorSize);
}

Vector StressStrainUtilities::TransformPrincipalStressesToSigmaTau(const Vector& rPrincipalStresses)
{
    auto result = Vector{2};
    result[0]   = 0.5 * (rPrincipalStresses[0] + rPrincipalStresses[2]);
    result[1]   = 0.5 * (rPrincipalStresses[0] - rPrincipalStresses[2]);
    return result;
}

Geo::SigmaTau StressStrainUtilities::TransformPrincipalStressesToSigmaTau(const Geo::PrincipalStresses& rPrincipalStresses)
{
    return Geo::SigmaTau{TransformPrincipalStressesToSigmaTau(rPrincipalStresses.CopyTo<Vector>())};
}

Vector StressStrainUtilities::TransformSigmaTauToPrincipalStresses(const Vector& rSigmaTau,
                                                                   const Vector& rPrincipalStresses)
{
    auto result = Vector{3};
    result[0]   = rSigmaTau[0] + rSigmaTau[1];
    result[1]   = rPrincipalStresses[1];
    result[2]   = rSigmaTau[0] - rSigmaTau[1];
    return result;
}

Vector StressStrainUtilities::TransformPrincipalStressesToPandQ(const Vector& rPrincipalStresses)
{
    auto stress_vector = Vector(6, 0.0);
    std::ranges::copy(rPrincipalStresses, stress_vector.begin());
    return UblasUtilities::CreateVector(
        {CalculateMeanStress(stress_vector), CalculateVonMisesStress(stress_vector)});
}

std::vector<Vector> StressStrainUtilities::CalculateStressVectorsFromStrainVectors(
    const std::vector<Vector>&                   rStrains,
    const ProcessInfo&                           rProcessInfo,
    const Properties&                            rProperties,
    const std::vector<ConstitutiveLaw::Pointer>& rConstitutiveLaws)
{
    // We have to make a copy of each strain vector, since setting it at the
    // constitutive law parameters requires a reference to a _mutable_ object!
    auto calculate_stress = [&rProperties, &rProcessInfo](auto Strain, auto& rpLaw) {
        auto law_parameters = ConstitutiveLaw::Parameters{};
        law_parameters.SetStrainVector(Strain);
        auto result = Vector{};
        result.resize(rpLaw->GetStrainSize());
        law_parameters.SetStressVector(result);
        law_parameters.SetMaterialProperties(rProperties);
        law_parameters.Set(ConstitutiveLaw::COMPUTE_STRESS);
        law_parameters.SetProcessInfo(rProcessInfo);
        rpLaw->CalculateMaterialResponseCauchy(law_parameters);
        return result;
    };
    auto result = std::vector<Vector>{};
    result.reserve(rStrains.size());
    std::ranges::transform(rStrains, rConstitutiveLaws, std::back_inserter(result), calculate_stress);

    return result;
}

} // namespace Kratos
