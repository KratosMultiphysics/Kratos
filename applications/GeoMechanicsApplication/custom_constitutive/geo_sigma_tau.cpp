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

const SigmaTau::InternalVectorType& SigmaTau::Values() const { return mValues; }

double SigmaTau::Sigma() const { return mValues[0]; }

double SigmaTau::Tau() const { return mValues[1]; }

} // namespace Kratos::Geo