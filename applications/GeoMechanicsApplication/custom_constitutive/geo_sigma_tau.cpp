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

#include "custom_constitutive/geo_sigma_tau.hpp"

namespace Kratos::Geo
{

SigmaTau::SigmaTau() : SigmaTau{{}} {}

SigmaTau::SigmaTau(const std::initializer_list<double>& rValues)
    : SigmaTau{rValues.begin(), rValues.end()}
{
}

SigmaTau::SigmaTau(const SigmaTau& rOther) : sigma{values[0]}, tau{values[1]} { *this = rOther; }

SigmaTau& SigmaTau::operator=(const SigmaTau& rOther)
{
    if (&rOther != this) std::ranges::copy(rOther.values, values.begin());

    return *this;
}

} // namespace Kratos::Geo
