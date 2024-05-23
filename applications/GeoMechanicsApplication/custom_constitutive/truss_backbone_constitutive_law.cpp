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
//  Main authors:    Wijtze Pieter Kikstra
//

// Project includes
#include "truss_backbone_constitutive_law.h"
#include "geo_mechanics_application_variables.h"
#include "includes/properties.h"

namespace Kratos
{

TrussBackboneConstitutiveLaw::TrussBackboneConstitutiveLaw() : ConstitutiveLaw() {}

TrussBackboneConstitutiveLaw::TrussBackboneConstitutiveLaw(const TrussBackboneConstitutiveLaw& rOther)
    : ConstitutiveLaw(rOther)
{
}

ConstitutiveLaw::Pointer TrussBackboneConstitutiveLaw::Clone() const
{
    return Kratos::make_shared<TrussBackboneConstitutiveLaw>(*this);
}

TrussBackboneConstitutiveLaw::~TrussBackboneConstitutiveLaw() = default;

void TrussBackboneConstitutiveLaw::GetLawFeatures(Features& rFeatures)
{
    // Set the strain size
    rFeatures.mStrainSize = 1;

    // Set the space dimension
    rFeatures.mSpaceDimension = 3;
}

double& TrussBackboneConstitutiveLaw::GetValue(const Variable<double>& rThisVariable, double& rValue)
{
    if (rThisVariable == ACCUMULATED_STRAIN) {
        rValue = mAccumulatedStrain;
    } else {
        KRATOS_ERROR << "Can't get the specified value for " << rThisVariable << std::endl;
    }
    return rValue;
}

array_1d<double, 3>& TrussBackboneConstitutiveLaw::GetValue(const Variable<array_1d<double, 3>>& rThisVariable,
                                                            array_1d<double, 3>& rValue)
{
    KRATOS_ERROR << "Can't get the specified value for " << rThisVariable << std::endl;
    return rValue;
}

void TrussBackboneConstitutiveLaw::SetValue(const Variable<double>& rThisVariable,
                                            const double&           rValue,
                                            const ProcessInfo&      rCurrentProcessInfo)
{
    if (rThisVariable == ACCUMULATED_STRAIN) {
        mAccumulatedStrain = rValue;
    } else {
        KRATOS_ERROR << "Can't set the specified value for " << rThisVariable << std::endl;
    }
}

double& TrussBackboneConstitutiveLaw::CalculateValue(ConstitutiveLaw::Parameters& rParameterValues,
                                                     const Variable<double>&      rThisVariable,
                                                     double&                      rValue)
{
    if (rThisVariable == TANGENT_MODULUS) {
        const auto axial_strain(rParameterValues.GetStrainVector()[0]);
        if (TrussBackboneConstitutiveLaw::IsWithinUnReLoading(axial_strain, rParameterValues)) {
            rValue = rParameterValues.GetMaterialProperties()[YOUNG_MODULUS];
        } else {
            rValue = TrussBackboneConstitutiveLaw::BackboneStiffness(axial_strain);
        }
    } else {
        KRATOS_ERROR << "Can't calculate the specified value" << std::endl;
    }
    return rValue;
}

void TrussBackboneConstitutiveLaw::CalculateMaterialResponsePK2(Parameters& rValues)
{
    // backbone loading, un- and reloading laws
    const auto axial_strain(rValues.GetStrainVector()[0]);
    const auto axial_strain_increment = axial_strain - mPreviousAxialStrain;
    auto&      axial_stress_vector    = rValues.GetStressVector();
    if (TrussBackboneConstitutiveLaw::IsWithinUnReLoading(axial_strain, rValues)) {
        // Un- or reloading
        axial_stress_vector[0] =
            rValues.GetMaterialProperties()[YOUNG_MODULUS] * (axial_strain - mUnReLoadCenter);
    } else {
        // backbone
        axial_stress_vector[0] = TrussBackboneConstitutiveLaw::BackboneStress(
            mAccumulatedStrain + (std::abs(axial_strain - mUnReLoadCenter) -
                                  (TrussBackboneConstitutiveLaw::CalculateUnReLoadAmplitude(rValues) / 2.)));
        axial_stress_vector[0] = (axial_strain - mUnReLoadCenter) >= 0. ? axial_stress_vector[0]
                                                                        : -axial_stress_vector[0];
    }
    rValues.SetStressVector(axial_stress_vector);
}

void TrussBackboneConstitutiveLaw::FinalizeMaterialResponsePK2(Parameters& rValues)
{
    CalculateMaterialResponsePK2(rValues);
    const auto axial_strain(rValues.GetStrainVector()[0]);
    if (!TrussBackboneConstitutiveLaw::IsWithinUnReLoading(axial_strain, rValues)) {
        // update accumulated backbone strain and un- reload center strain
        auto pos_side = (axial_strain - mUnReLoadCenter) > 0.;
        mAccumulatedStrain += std::abs(axial_strain - mUnReLoadCenter) -
                              (TrussBackboneConstitutiveLaw::CalculateUnReLoadAmplitude(rValues) / 2.);
        mUnReLoadCenter =
            pos_side
                ? axial_strain - TrussBackboneConstitutiveLaw::CalculateUnReLoadAmplitude(rValues) / 2.
                : axial_strain + TrussBackboneConstitutiveLaw::CalculateUnReLoadAmplitude(rValues) / 2.;
    }
    mPreviousAxialStrain = axial_strain;
}

int TrussBackboneConstitutiveLaw::Check(const Properties&   rMaterialProperties,
                                        const GeometryType& rElementGeometry,
                                        const ProcessInfo&  rCurrentProcessInfo) const
{
    KRATOS_CHECK(rMaterialProperties.Has(YOUNG_MODULUS))
    // andere input voor backbone nalopen
    return 0;
}

double TrussBackboneConstitutiveLaw::BackboneStress(const double Strain) const
{
    double young_modulus = 200.;
    return 0.5 * young_modulus * Strain;
}

double TrussBackboneConstitutiveLaw::BackboneStiffness([[maybe_unused]] double Strain) const
{
    double young_modulus = 200.;
    return 0.5 * young_modulus;
}

double TrussBackboneConstitutiveLaw::CalculateUnReLoadCenter(const double Strain,
                                                             const double Stress,
                                                             ConstitutiveLaw::Parameters& rParameterValues) const
{
    return Stress >= 0
               ? Strain - TrussBackboneConstitutiveLaw::CalculateUnReLoadAmplitude(rParameterValues)
               : Strain - TrussBackboneConstitutiveLaw::CalculateUnReLoadAmplitude(rParameterValues);
}

double TrussBackboneConstitutiveLaw::CalculateUnReLoadAmplitude(ConstitutiveLaw::Parameters& rParameterValues) const
{
    return 2. * TrussBackboneConstitutiveLaw::BackboneStress(mAccumulatedStrain) /
           rParameterValues.GetMaterialProperties()[YOUNG_MODULUS];
}

bool TrussBackboneConstitutiveLaw::IsWithinUnReLoading(const double Strain,
                                                       ConstitutiveLaw::Parameters& rParameterValues) const
{
    return std::abs(Strain - mUnReLoadCenter) <
           (TrussBackboneConstitutiveLaw::CalculateUnReLoadAmplitude(rParameterValues) / 2.);
}

} // Namespace Kratos
