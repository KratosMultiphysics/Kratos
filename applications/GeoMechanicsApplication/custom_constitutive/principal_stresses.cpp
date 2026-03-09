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
PrincipalStresses::PrincipalStresses(double Sigma1, double Sigma2, double Sigma3)
{
    mValues[0] = Sigma1;
    mValues[1] = Sigma2;
    mValues[2] = Sigma3;
}

const PrincipalStresses::InternalVectorType& PrincipalStresses::Values() const { return mValues; }

PrincipalStresses::InternalVectorType& PrincipalStresses::Values() { return mValues; }

PrincipalStresses& PrincipalStresses::operator+=(const PrincipalStresses& rRhs)
{
    std::ranges::transform(mValues, rRhs.mValues, mValues.begin(), std::plus{});
    return *this;
}

PrincipalStresses operator+(PrincipalStresses Lhs, const PrincipalStresses& rRhs)
{
    Lhs += rRhs;
    return Lhs;
}

} // namespace Kratos::Geo