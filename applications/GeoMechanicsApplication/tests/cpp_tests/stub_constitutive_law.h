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

class StubConstitutiveLaw : public ConstitutiveLaw
{
public:
    StubConstitutiveLaw()
    {
        ON_CALL(*this, RequiresInitializeMaterialResponse).WillByDefault(testing::Return(false));
        ON_CALL(*this, Clone).WillByDefault([]() { return std::make_shared<StubConstitutiveLaw>(); });
        ON_CALL(*this, GetStrainSize).WillByDefault(testing::Return(4));
    }

    MOCK_METHOD(bool, RequiresInitializeMaterialResponse, (), (override));
    MOCK_METHOD(ConstitutiveLaw::Pointer, Clone, (), (const, override));
    MOCK_METHOD(SizeType, GetStrainSize, (), (const, override));

    void GetLawFeatures(ConstitutiveLaw::Features& rFeatures) override
    {
        rFeatures.mStrainMeasures.clear();
        rFeatures.mStrainMeasures.push_back(StrainMeasure_Infinitesimal);
    }
};

} // namespace Kratos::Testing
