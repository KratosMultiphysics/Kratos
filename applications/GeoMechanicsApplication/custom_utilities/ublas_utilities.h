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

#pragma once
#include "includes/ublas_interface.h"

namespace Kratos
{
class KRATOS_API(GEO_MECHANICS_APPLICATION) UblasUtilities
{
public:
    static Vector CreateVector(const std::initializer_list<double>& rInitializerList);
    static Matrix CreateMatrix(const std::initializer_list<std::initializer_list<double>>& rRows);
};

} // namespace Kratos