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

#include <optional>

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

    static Matrix MakeInterfaceConstitutiveMatrix(double      NormalStiffness,
                                                  double      ShearStiffness,
                                                  std::size_t TractionSize,
                                                  std::size_t NumberOfNormalComponents);

    static void CheckStrainSize(const Properties& rProperties, std::size_t ExpectedSize, std::size_t ElementId);

    static void CheckHasStrainMeasure_Infinitesimal(const Properties& rProperties, std::size_t ElementId);

    [[nodiscard]] static double CalculateK0NCFromFrictionAngleInRadians(double FrictionAngleInRadians);
    [[nodiscard]] static DrainageType StringToDrainageType(const std::string& rDrainageTypeName);
    [[nodiscard]] static bool         IsConstantWaterPressure(const Properties& rProperties);
    static void                       ReplaceIgnoreUndrainedByDrainageType(Properties* pProperties);
}; /* Class ConstitutiveLawUtilities*/

} // namespace Kratos
