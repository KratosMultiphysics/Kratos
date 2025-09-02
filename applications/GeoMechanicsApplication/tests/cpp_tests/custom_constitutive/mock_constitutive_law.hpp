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

#include "geo_mechanics_application_variables.h"
#include "includes/constitutive_law.h"

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

    void GetLawFeatures(Features& rFeatures) override
    {
        mNumberOfCalls++;
        if (mNumberOfCalls == 2) rFeatures.mStrainMeasures.push_back(StrainMeasure_Infinitesimal);
        if (mNumberOfCalls == 1)
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

private:
    Vector mStateVariables;
    int    mNumberOfCalls = 0;
};
} // namespace Kratos::Testing
