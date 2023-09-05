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

namespace Kratos
{

class ConstitutiveLawUtilities
{

public:

    static int GetStateVariableIndex(const Variable<double>& rThisVariable);

}; /* Class ConstitutiveLawUtilities*/

}