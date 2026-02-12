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

#include "custom_constitutive/sigma_tau.hpp"

namespace Kratos::Geo
{

SigmaTau::SigmaTau(const std::initializer_list<double>& rValues)
    : SigmaTau{rValues.begin(), rValues.end()}
{
}

const SigmaTau::InternalVectorType& SigmaTau::Values() const { return mValues; }

double SigmaTau::Sigma() const { return mValues[0]; }

double& SigmaTau::Sigma() { return mValues[0]; }

double SigmaTau::Tau() const { return mValues[1]; }

double& SigmaTau::Tau() { return mValues[1]; }

} // namespace Kratos::Geo