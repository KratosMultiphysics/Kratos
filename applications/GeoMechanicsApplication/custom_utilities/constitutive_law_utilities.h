// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Anne van de Graaf,
//                   Wijtze Pieter Kikstra
//

#pragma once

#include "containers/variable.h"
#include "geo_mechanics_application_constants.h"
#include "includes/constitutive_law.h"

namespace Kratos
{

class KRATOS_API(GEO_MECHANICS_APPLICATION) ConstitutiveLawUtilities
{
public:
    static int GetStateVariableIndex(const Variable<double>& rThisVariable);

    static void SetConstitutiveParameters(ConstitutiveLaw::Parameters& rConstitutiveParameters,
                                          Vector&                      rStrainVector,
                                          Matrix&                      rConstitutiveMatrix,
                                          const Vector&                rN,
                                          const Matrix&                rGradNpT,
                                          const Matrix&                rF,
                                          double                       detF);

    static double GetCohesion(const Properties& rProperties);
    static bool   HasFrictionAngle(const Properties& rProperties);
    static void   ValidateFrictionAngle(const Properties& rProperties, IndexType ElementId);
    static double GetFrictionAngleInDegrees(const Properties& rProperties);
    static double GetFrictionAngleInRadians(const Properties& rProperties);

    static std::pair<double, double> GetOrCalculateElasticProperties(const Properties& rProperties,
                                                                     bool Undrained = false);

    static Matrix MakeInterfaceElasticConstitutiveTensor(double      NormalStiffness,
                                                         double      ShearStiffness,
                                                         std::size_t TractionSize,
                                                         std::size_t NumberOfNormalComponents);

    static void CheckStrainSize(const Properties& rProperties, std::size_t ExpectedSize, std::size_t ElementId);

    static void CheckHasStrainMeasure_Infinitesimal(const Properties& rProperties, std::size_t ElementId);

    [[nodiscard]] static double CalculateK0NCFromFrictionAngleInRadians(double FrictionAngleInRadians);

    [[nodiscard]] static double CalculateUndrainedYoungsModulus(const Properties& rProperties,
                                                                double UndrainedPoissonsRatio);

    [[nodiscard]] static double CalculateUndrainedPoissonsRatio(const Properties& rProperties);

    [[nodiscard]] static double GetOrCalculateUndrainedPoissonsRatio(const Properties& rProperties);

    [[nodiscard]] static double GetOrCalculateSkemptonB(const Properties& rProperties);

    [[nodiscard]] static Matrix MakeContinuumElasticConstitutiveTensor(double      YoungsModulus,
                                                                       double      PoissonsRatio,
                                                                       std::size_t StrainSize,
                                                                       std::size_t NumberOfNormalComponents);

    [[nodiscard]] static DrainageType StringToDrainageType(const std::string& rDrainageTypeName);
    [[nodiscard]] static bool         IsUndrained(const Properties& rProperties);
    [[nodiscard]] static bool         IsConstantWaterPressure(const Properties& rProperties);
    [[nodiscard]] static std::size_t  GetNumberOfNormalStrainComponents(const Properties& rProperties);
    [[nodiscard]] static double       CalculateVolumetricStrain(const Vector&         rStrainVector,
                                                                         const Properties&     rProperties);
    static void                       ReplaceIgnoreUndrainedByDrainageType(Properties& rProperties);

    [[nodiscard]] static double CalculateExcessPorePressureIncrement(const Properties& rProperties,
                                                                     double VolumetricStrainIncrement);

    [[nodiscard]] static Vector CalculateExcessPorePressureForce(const Properties&     rProperties,
                                                                 const Vector&         rStrainVector,
                                                                 const Matrix&         rB,
                                                                 const Vector&         rVoigtVector,
                                                                 double                IntegrationCoefficient,
                                                                 std::size_t           IntegrationPoint,
                                                                 const Vector&         rExcessPorePressurePrevious);

    static void AssembleExcessPorePressureForces(Vector&                   rResultVector,
                                                 const Properties&         rProperties,
                                                 const std::vector<Vector>& rStrainVectors,
                                                 const std::vector<Matrix>& rBMatrices,
                                                 const Vector&             rVoigtVector,
                                                 const std::vector<double>& rIntegrationCoefficients,
                                                 const Vector&             rExcessPorePressurePrevious);
}; /* Class ConstitutiveLawUtilities*/

} // namespace Kratos
