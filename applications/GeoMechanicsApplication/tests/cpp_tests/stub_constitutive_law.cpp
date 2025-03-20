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
#include "stub_constitutive_law.h"

namespace Kratos::Testing
{

bool StubConstitutiveLaw::RequiresInitializeMaterialResponse() { return false; }

ConstitutiveLaw::Pointer StubConstitutiveLaw::Clone() const
{
    return std::make_shared<StubConstitutiveLaw>();
}

} // namespace Kratos::Testing
