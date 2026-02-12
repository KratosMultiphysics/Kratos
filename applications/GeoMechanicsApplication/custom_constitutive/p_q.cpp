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

PQ::PQ(const std::initializer_list<double>& rValues) : PQ{rValues.begin(), rValues.end()} {}

const PQ::InternalArrayType& PQ::Values() const { return mValues; }

double PQ::P() const { return mValues[0]; }

double& PQ::P() { return mValues[0]; }

double PQ::Q() const { return mValues[1]; }

double& PQ::Q() { return mValues[1]; }

} // namespace Kratos::Geo