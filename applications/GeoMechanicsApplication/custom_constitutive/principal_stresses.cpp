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

#include "principal_stresses.hpp"

namespace Kratos::Geo
{
PrincipalStresses::PrincipalStresses(const std::initializer_list<double>& rValues)
    : PrincipalStresses{std::begin(rValues), std::end(rValues)}
{
}

const PrincipalStresses::InternalVectorType& PrincipalStresses::Values() const { return mValues; }

PrincipalStresses::InternalVectorType& PrincipalStresses::Values() { return mValues; }

} // namespace Kratos::Geo