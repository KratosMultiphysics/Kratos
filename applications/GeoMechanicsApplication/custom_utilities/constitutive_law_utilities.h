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

// Project includes
#include "containers/variable.h"
#include "custom_constitutive/coulomb_yield_surface.h"
#include "custom_constitutive/tension_cutoff.h"
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
    static double GetFrictionAngleInDegrees(const Properties& rProperties);

    static Vector MapStressesInMorhCoulomb(const Properties&          rProperties,
                                           Vector&                    rSigmaTau,
                                           const CoulombYieldSurface& rCoulombYieldSurface,
                                           const TensionCutoff&       rTensionCutOff);
    static bool   IsAdmissiblePrincipalStressState(const Vector&              rSigmaTau,
                                                   const CoulombYieldSurface& rCoulombYieldSurface,
                                                   const TensionCutoff&       rTensionCutOff);

    static double GetFrictionAngleInRadians(const Properties& rProperties);
}; /* Class ConstitutiveLawUtilities*/

} // namespace Kratos
