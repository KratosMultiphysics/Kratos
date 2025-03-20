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
#pragma once

#include "includes/constitutive_law.h"

namespace Kratos::Testing
{

class StubConstitutiveLaw : public ConstitutiveLaw
{
public:
    bool                                   RequiresInitializeMaterialResponse() override;
    [[nodiscard]] ConstitutiveLaw::Pointer Clone() const override;
};

} // namespace Kratos::Testing
