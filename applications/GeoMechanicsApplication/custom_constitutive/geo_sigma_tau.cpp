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

SigmaTau::SigmaTau(const std::initializer_list<double>& rValues)
{
    KRATOS_DEBUG_ERROR_IF(rValues.size() != msVectorSize && rValues.size() != 0)
        << "Cannot construct a SigmaTau instance: the given initializer list is nonempty and has "
        << rValues.size() << " entry/ies, but expected " << msVectorSize << "\n";
    std::ranges::copy(rValues, mValues.begin());
}

const SigmaTau::InternalVectorType& SigmaTau::Values() const { return mValues; }

double SigmaTau::Sigma() const { return mValues[0]; }

double SigmaTau::Tau() const { return mValues[1]; }

} // namespace Kratos::Geo