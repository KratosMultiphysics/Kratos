// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Richard Faasse
//

#include "principal_stresses.h"

namespace Kratos::Geo
{

PrincipalStresses::InternalVectorType PrincipalStresses::Values() const { return mValues; }
} // namespace Kratos::Geo