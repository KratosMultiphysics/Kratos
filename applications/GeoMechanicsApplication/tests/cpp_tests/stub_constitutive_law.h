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

#include <gmock/gmock.h>

namespace Kratos::Testing
{

ConstitutiveLaw::Pointer MakeStubConstitutiveLaw();

class StubConstitutiveLaw : public ConstitutiveLaw
{
public:
    StubConstitutiveLaw()
    {
        ON_CALL(*this, RequiresInitializeMaterialResponse).WillByDefault(testing::Return(false));
        ON_CALL(*this, Clone).WillByDefault([]() { return MakeStubConstitutiveLaw(); });
    }

    MOCK_METHOD(bool, RequiresInitializeMaterialResponse, (), (override));
    MOCK_METHOD(ConstitutiveLaw::Pointer, Clone, (), (const, override));
};

} // namespace Kratos::Testing
