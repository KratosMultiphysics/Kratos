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
    static double GetFrictionAngleInDegrees(const Properties& rProperties);
    static double GetFrictionAngleInRadians(const Properties& rProperties);

    static Matrix MakeInterfaceConstitutiveMatrix(double      NormalStiffness,
                                                  double      ShearStiffness,
                                                  std::size_t TractionSize,
                                                  std::size_t NumberOfNormalComponents);

    static void CheckStrainSize(const Properties& rProperties, std::size_t ExpectedSize, std::size_t ElementId);

    static void CheckHasStrainMeasure_Infinitesimal(const Properties& rProperties, std::size_t ElementId);

    [[nodiscard]] static double CalculateK0NCFromFrictionAngleInDegrees(double FrictionAngleInDegrees);
}; /* Class ConstitutiveLawUtilities*/

} // namespace Kratos
