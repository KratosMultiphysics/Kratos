// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Anne van de Graaf
//

#include "custom_constitutive/p_q.hpp"

namespace Kratos::Geo
{

PQTheta::PQTheta(double P, double Q, double Theta)
{
    mValues[0] = P;
    mValues[1] = Q;
    mValues[2] = Theta;
}

const PQTheta::InternalVectorType& PQTheta::Values() const noexcept { return mValues; }

double PQTheta::P() const noexcept { return mValues[0]; }

double& PQTheta::P() noexcept { return mValues[0]; }

double PQTheta::Q() const noexcept { return mValues[1]; }

double& PQTheta::Q() noexcept { return mValues[1]; }

double PQTheta::Theta() const noexcept { return mValues[2]; }

double& PQTheta::Theta() noexcept { return mValues[2]; }

PQTheta& PQTheta::operator+=(const PQTheta& rRhs)
{
    std::ranges::transform(mValues, rRhs.mValues, mValues.begin(), std::plus{});
    return *this;
}

PQTheta operator+(PQTheta Lhs, const PQTheta& rRhs)
{
    Lhs += rRhs;
    return Lhs;
}

} // namespace Kratos::Geo