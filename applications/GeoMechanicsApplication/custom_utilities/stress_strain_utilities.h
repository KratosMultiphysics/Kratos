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
#include "custom_constitutive/geo_principal_stresses.hpp"
#include "includes/constitutive_law.h"
#include "includes/define.h"
#include "includes/ublas_interface.h"

namespace Kratos
{

namespace Geo
{
struct PrincipalStresses;
struct SigmaTau;
} // namespace Geo

class KRATOS_API(GEO_MECHANICS_APPLICATION) StressStrainUtilities
{
public:
    static double CalculateVonMisesStress(const Vector& rStressVector);
    static double CalculateTrace(const Vector& rStressVector);
    static double CalculateMeanStress(const Vector& rStressVector);
    static double CalculateLodeAngle(const Vector& rStressVector);
    static double CalculateMohrCoulombShearCapacity(const Vector& rStressVector, double C, double PhiInRadians);
    static double CalculateMohrCoulombPressureCapacity(const Vector& rStressVector, double C, double PhiInRadians);
    static double CalculateVonMisesStrain(const Vector& rStrainVector);
    static Vector CalculateHenckyStrain(const Matrix& rDeformationGradient, size_t VoigtSize);
    static Matrix CalculateGreenLagrangeStrainTensor(const Matrix& rDeformationGradient);
    static Vector CalculateCauchyStrain(const Matrix& rB, const Vector& rDisplacements);
    static std::vector<Vector> CalculateStrains(const std::vector<Matrix>& rDeformationGradients,
                                                const std::vector<Matrix>& rBs,
                                                const Vector&              rDisplacements,
                                                bool                       UseHenckyStrain,
                                                std::size_t                VoigtSize);
    static std::pair<Geo::PrincipalStresses, Matrix> CalculatePrincipalStresses(const Vector& rStressVector);
    static void CalculatePrincipalStresses(const Vector& rCauchyStressVector,
                                           Vector&       rPrincipalStressVector,
                                           Matrix&       rEigenVectorsMatrix);
    static void ReorderEigenValuesAndVectors(Vector& rPrincipalStressVector, Matrix& rEigenVectorsMatrix);
    static Vector RotatePrincipalStresses(const Vector& rPrincipalStressVector,
                                          const Matrix& rRotationMatrix,
                                          std::size_t   StressVectorSize);

    static Vector TransformPrincipalStressesToSigmaTau(const Vector& rPrincipalStresses);
    static Geo::SigmaTau TransformPrincipalStressesToSigmaTau(const Geo::PrincipalStresses& rPrincipalStresses);
    static Vector TransformSigmaTauToPrincipalStresses(const Vector& rSigmaTau, const Vector& rPrincipalStresses);
    static Geo::PrincipalStresses TransformSigmaTauToPrincipalStresses(const Geo::SigmaTau& rSigmaTau,
                                                                       const Geo::PrincipalStresses& rPrincipalStresses);
    static Vector TransformPrincipalStressesToPandQ(const Vector& rPrincipalStresses);

    /// @brief This function calculates stresses from strains using the constitutive laws.
    /// However, it can also be used to calculate tractions from relative displacements.
    /// In that case, the relative displacements should be input as strains, and the return
    /// is a vector with tractions
    static std::vector<Vector> CalculateStressVectorsFromStrainVectors(
        const std::vector<Vector>&                   rStrains,
        const ProcessInfo&                           rProcessInfo,
        const Properties&                            rProperties,
        const std::vector<ConstitutiveLaw::Pointer>& rConstitutiveLaws);

private:
    static double CalculateQMohrCoulomb(const Vector& rStressVector, double C, double PhiInRadians);
    static double CalculateDenominator(const Vector& rStressVector, double PhiInRadians);
};

} // namespace Kratos