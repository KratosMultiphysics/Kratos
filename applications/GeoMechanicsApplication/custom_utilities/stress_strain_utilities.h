// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    JMCarbonell,
//                   Vahid Galavi,
//                   Richard Faasse,
//                   Gennady Markelov
//

#pragma once

/* Project includes */
#include "geo_mechanics_application_constants.h"
#include "includes/define.h"
#include "includes/ublas_interface.h"

namespace Kratos
{

class KRATOS_API(GEO_MECHANICS_APPLICATION) StressStrainUtilities
{
public:
    static double CalculateVonMisesStress(const Vector& rStressVector);
    static double CalculateTrace(const Vector& rStressVector);
    static double CalculateMeanStress(const Vector& rStressVector);
    static double CalculateLodeAngle(const Vector& rStressVector);
    static double CalculateMohrCoulombShearCapacity(const Vector& rStressVector, double C, double Phi);
    static double CalculateMohrCoulombPressureCapacity(const Vector& rStressVector, double C, double Phi);
    static double CalculateVonMisesStrain(const Vector& rStrainVector);
    static Vector CalculateHenckyStrain(const Matrix& rDeformationGradient, size_t VoigtSize);
    static Matrix CalculateGreenLagrangeStrainTensor(const Matrix& rDeformationGradient);
    static Vector CalculateCauchyStrain(const Matrix& rB, const Vector& rDisplacements);
    static std::vector<Vector> CalculateStrains(const std::vector<Matrix>& rDeformationGradients,
                                                const std::vector<Matrix>& rBs,
                                                const Vector&              rDisplacements,
                                                bool                       UseHenckyStrain,
                                                std::size_t                VoigtSize);

    static constexpr std::size_t GetVoigtSize(std::size_t Dimension)
    {
        return Dimension == N_DIM_3D ? VOIGT_SIZE_3D : VOIGT_SIZE_2D_PLANE_STRAIN;
    }

    static constexpr std::size_t GetStressTensorSize(std::size_t Dimension)
    {
        return Dimension == N_DIM_3D ? STRESS_TENSOR_SIZE_3D : STRESS_TENSOR_SIZE_2D;
    }

    static const Vector& GetVoigtVector(std::size_t Dimension);

private:
    static double       CalculateQMohrCoulomb(const Vector& rStressVector, double C, double Phi);
    static double       CalculateDenominator(const Vector& rStressVector, double Phi);
    static Vector       DefineVoigtVector(std::size_t Dimension);
    static const Vector VoigtVector2D;
    static const Vector VoigtVector3D;
};

} // namespace Kratos