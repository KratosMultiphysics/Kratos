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

    static void CheckProperty(const Properties&       rMaterialProperties,
                              const Variable<double>& rVariable,
                              std::optional<double>   MaxValue = std::nullopt);

    static Matrix ConstitutiveLawUtilities::MakeInterfaceConstitutiveMatrix(double NormalStiffness,
                                                                            double ShearStiffness,
                                                                            std::size_t TractionSize);
}; /* Class ConstitutiveLawUtilities*/

} // namespace Kratos
