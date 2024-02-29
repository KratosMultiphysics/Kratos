// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Richard Faasse
//

#include "stress_strain_utilities.h"
#include "custom_utilities/math_utilities.hpp"
#include "geo_mechanics_application_constants.h"

namespace Kratos
{

double StressStrainUtilities::CalculateVonMisesStress(const Vector& StressVector)
{
    Matrix LocalStressTensor =
        MathUtils<double>::StressVectorToTensor(StressVector); // reduced dimension stress tensor

    Matrix StressTensor(3, 3); // 3D stress tensor
    noalias(StressTensor) = ZeroMatrix(3, 3);
    for (std::size_t i = 0; i < LocalStressTensor.size1(); ++i) {
        for (std::size_t j = 0; j < LocalStressTensor.size2(); ++j) {
            StressTensor(i, j) = LocalStressTensor(i, j);
        }
    }

    double SigmaEquivalent =
        0.5 * ((StressTensor(0, 0) - StressTensor(1, 1)) * (StressTensor(0, 0) - StressTensor(1, 1)) +
               (StressTensor(1, 1) - StressTensor(2, 2)) * (StressTensor(1, 1) - StressTensor(2, 2)) +
               (StressTensor(2, 2) - StressTensor(0, 0)) * (StressTensor(2, 2) - StressTensor(0, 0)) +
               6.0 * (StressTensor(0, 1) * StressTensor(1, 0) + StressTensor(1, 2) * StressTensor(2, 1) +
                      StressTensor(2, 0) * StressTensor(0, 2)));

    return std::sqrt(std::max(SigmaEquivalent, 0.));
}

double StressStrainUtilities::CalculateTrace(const Vector& StressVector)
{
    Matrix StressTensor = MathUtils<double>::StressVectorToTensor(StressVector); // reduced dimension stress tensor

    double trace = 0.0;
    for (std::size_t i = 0; i < StressTensor.size1(); ++i) {
        trace += StressTensor(i, i);
    }

    return trace;
}

double StressStrainUtilities::CalculateMeanStress(const Vector& StressVector)
{
    return CalculateTrace(StressVector) / (StressVector.size() == 3 ? 2.0 : 3.0);
}

double StressStrainUtilities::CalculateVonMisesStrain(const Vector& StrainVector)
{
    return (2.0 / 3.0) * CalculateVonMisesStress(StrainVector);
}

Vector StressStrainUtilities::CalculateHenckyStrain(const Matrix& DeformationGradient, size_t VoigtSize)
{
    // right Cauchy Green deformation tensor C
    Matrix C = prod(trans(DeformationGradient), DeformationGradient);
    // Eigenvalues of C matrix, so principal right Cauchy Green deformation tensor C
    Matrix EigenValuesMatrix;
    Matrix EigenVectorsMatrix;
    MathUtils<double>::GaussSeidelEigenSystem(C, EigenVectorsMatrix, EigenValuesMatrix, 1.0e-16, 20);
    // Compute natural strain == Logarithmic strain == Hencky strain from principal strains
    for (std::size_t i = 0; i < DeformationGradient.size1(); ++i) {
        EigenValuesMatrix(i, i) = 0.5 * std::log(EigenValuesMatrix(i, i));
    }

    // Rotate from principal strains back to the used coordinate system
    Matrix ETensor;
    MathUtils<double>::BDBtProductOperation(ETensor, EigenValuesMatrix, EigenVectorsMatrix);

    // From tensor to vector
    if (DeformationGradient.size1() == 2 && VoigtSize == 4) {
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

} // namespace Kratos
