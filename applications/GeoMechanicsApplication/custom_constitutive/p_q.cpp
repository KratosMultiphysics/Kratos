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

#include "custom_constitutive/p_q.hpp"

namespace Kratos::Geo
{

PQ::PQ(double P, double Q)
{
    mValues[0] = P;
    mValues[1] = Q;
}

const PQ::InternalVectorType& PQ::Values() const noexcept { return mValues; }

double PQ::P() const noexcept { return mValues[0]; }

double& PQ::P() noexcept { return mValues[0]; }

double PQ::Q() const noexcept { return mValues[1]; }

double& PQ::Q() noexcept { return mValues[1]; }

PQ& PQ::operator+=(const PQ& rRhs)
{
    std::ranges::transform(mValues, rRhs.mValues, mValues.begin(), std::plus{});
    return *this;
}

PQ operator+(PQ Lhs, const PQ& rRhs)
{
    Lhs += rRhs;
    return Lhs;
}

} // namespace Kratos::Geo