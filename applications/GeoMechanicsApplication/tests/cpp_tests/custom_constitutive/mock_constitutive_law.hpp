// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Wijtze Pieter Kikstra
//                   Gennady Markelov
//

#pragma once

#include "geo_mechanics_application_variables.h"
#include "includes/constitutive_law.h"
#include <gmock/gmock.h>

namespace Kratos::Testing
{

class MockConstitutiveLaw : public ConstitutiveLaw
{
public:
    [[nodiscard]] ConstitutiveLaw::Pointer Clone() const override
    {
        return std::make_shared<MockConstitutiveLaw>();
    }

    [[nodiscard]] SizeType GetStrainSize() const override { return 4; }

    void AddStrainMeasure_Infinitesimal(bool Add_StrainMeasure_Infinitesimal)
    {
        mAdd_StrainMeasure_Infinitesimal = Add_StrainMeasure_Infinitesimal;
    }

    void GetLawFeatures(Features& rFeatures) override
    {
        if (mAdd_StrainMeasure_Infinitesimal)
            rFeatures.mStrainMeasures.push_back(StrainMeasure_Infinitesimal);
        rFeatures.mStrainMeasures.push_back(StrainMeasure_Deformation_Gradient);
    }

    void SetValue(const Variable<Vector>& rVariable, const Vector& rValue, const ProcessInfo& rCurrentProcessInfo) override
    {
        if (rVariable == STATE_VARIABLES) {
            mStateVariables = rValue;
        }
    }

    using ConstitutiveLaw::SetValue;

    MOCK_METHOD(bool, RequiresInitializeMaterialResponse, (), (override));
    MOCK_METHOD(void, CalculateMaterialResponseCauchy, (Parameters&), (override));
    MOCK_METHOD(void, FinalizeMaterialResponseCauchy, (Parameters&), (override));

    Vector& GetValue(const Variable<Vector>& rVariable, Vector& rValue) override
    {
        if (rVariable == STATE_VARIABLES) {
            rValue = mStateVariables;
        }
        return rValue;
    }

    using ConstitutiveLaw::GetValue;

    bool Has(const Variable<Vector>& rThisVariable) override
    {
        return rThisVariable == STATE_VARIABLES;
    }

    using ConstitutiveLaw::Has;

    Vector mStateVariables;
    bool   mAdd_StrainMeasure_Infinitesimal = false;
};
} // namespace Kratos::Testing
