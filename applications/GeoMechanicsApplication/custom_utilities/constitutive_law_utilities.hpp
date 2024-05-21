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
#include "includes/constitutive_law.h"

namespace Kratos
{

class ConstitutiveLawUtilities
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

}; /* Class ConstitutiveLawUtilities*/

} // namespace Kratos