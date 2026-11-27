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
//                   Wijtze Pieter Kikstra
//

#include <cmath> // for std::abs
#include <functional>
#include <numbers>

#include "custom_constitutive/p_q_theta.hpp"

namespace Kratos::Geo
{

PQTheta::PQTheta(double P, double Q, double Theta) : PQTheta{MakeVector(P, Q, Theta)} {}

PQTheta::PQTheta(const InternalVectorType& rValues) : mValues{rValues} { CheckValues(); }

const PQTheta::InternalVectorType& PQTheta::Values() const noexcept { return mValues; }

double PQTheta::P() const noexcept { return mValues[0]; }

double PQTheta::Q() const noexcept { return mValues[1]; }

double PQTheta::Theta() const noexcept { return mValues[2]; }

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

PQTheta::InternalVectorType PQTheta::MakeVector(double P, double Q, double Theta)
{
    InternalVectorType result;
    result[0] = P;
    result[1] = Q;
    result[2] = Theta;
    return result;
}

void PQTheta::CheckValues() const
{
    KRATOS_ERROR_IF(mValues[1] < 0.0) << "Invalid PQTheta: Q = " << mValues[1] << " is negative.\n";
    KRATOS_ERROR_IF(std::abs(mValues[2]) > std::numbers::pi / 6.0)
        << "Invalid PQTheta: Theta = " << mValues[2] << " is out of the range [-pi/6, pi/6].\n";
}

} // namespace Kratos::Geo