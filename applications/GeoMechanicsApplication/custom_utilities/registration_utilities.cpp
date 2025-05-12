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
#include "registration_utilities.h"

namespace Kratos
{

ScopedSerializerRegistration::~ScopedSerializerRegistration() { Serializer::Deregister(mName); }

} // namespace Kratos