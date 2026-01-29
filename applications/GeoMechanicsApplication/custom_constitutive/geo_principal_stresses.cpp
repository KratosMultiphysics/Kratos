// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//
//  Main authors:    Anne van de Graaf
//

#include "custom_constitutive/geo_principal_stresses.hpp"

namespace Kratos::Geo
{

PrincipalStresses::PrincipalStresses(const std::initializer_list<double>& rValues)
    : PrincipalStresses{rValues.begin(), rValues.end()}
{
}

} // namespace Kratos::Geo
