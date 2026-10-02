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

#include "custom_constitutive/sigma_tau.hpp"

namespace Kratos::Geo
{

SigmaTau::SigmaTau(double Sigma, double Tau)
{
    mValues[0] = Sigma;
    mValues[1] = Tau;
}

const SigmaTau::InternalVectorType& SigmaTau::Values() const noexcept { return mValues; }

double SigmaTau::Sigma() const noexcept { return mValues[0]; }

double& SigmaTau::Sigma() noexcept { return mValues[0]; }

double SigmaTau::Tau() const noexcept { return mValues[1]; }

double& SigmaTau::Tau() noexcept { return mValues[1]; }

SigmaTau& SigmaTau::operator+=(const SigmaTau& rRhs)
{
    std::ranges::transform(mValues, rRhs.mValues, mValues.begin(), std::plus{});
    return *this;
}

SigmaTau operator+(SigmaTau Lhs, const SigmaTau& rRhs)
{
    Lhs += rRhs;
    return Lhs;
}

} // namespace Kratos::Geo