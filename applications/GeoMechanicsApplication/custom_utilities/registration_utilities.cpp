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
#include "registration_utilities.hpp"

namespace Kratos
{

ScopedSerializerRegistration::~ScopedSerializerRegistration()
{
    for (const auto& r_name : mNames) {
        Serializer::Deregister(r_name);
    }
}
} // namespace Kratos