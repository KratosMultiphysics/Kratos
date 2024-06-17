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
ConstitutiveLaw::Pointer TrussBackboneConstitutiveLaw::Clone() const
{
    return Kratos::make_shared<TrussBackboneConstitutiveLaw>(*this);
}

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
        const auto youngs_modulus = rParameterValues.GetMaterialProperties()[YOUNG_MODULUS];
        const auto axial_strain(rParameterValues.GetStrainVector()[0]);
        if (IsWithinUnReLoading(axial_strain, youngs_modulus)) {
            rValue = youngs_modulus;
        } else {
            rValue = BackboneStiffness(axial_strain);
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
    auto&      axial_stress_vector = rValues.GetStressVector();
    const auto youngs_modulus      = rValues.GetMaterialProperties()[YOUNG_MODULUS];
    if (IsWithinUnReLoading(axial_strain, youngs_modulus)) {
        // Un- or reloading
        axial_stress_vector[0] = youngs_modulus * (axial_strain - mUnReLoadCenter);
    } else {
        // backbone
        axial_stress_vector[0] =
            BackboneStress(mAccumulatedStrain + (std::abs(axial_strain - mUnReLoadCenter) -
                                                 (CalculateUnReLoadAmplitude(youngs_modulus) / 2.)));
        axial_stress_vector[0] = (axial_strain - mUnReLoadCenter) >= 0. ? axial_stress_vector[0]
                                                                        : -axial_stress_vector[0];
    }
    rValues.SetStressVector(axial_stress_vector);
}

void TrussBackboneConstitutiveLaw::FinalizeMaterialResponsePK2(Parameters& rValues)
{
    CalculateMaterialResponsePK2(rValues);
    const auto axial_strain(rValues.GetStrainVector()[0]);
    if (const auto youngs_modulus = rValues.GetMaterialProperties()[YOUNG_MODULUS];
        !IsWithinUnReLoading(axial_strain, youngs_modulus)) {
        // update accumulated backbone strain and un- reload center strain
        auto pos_side = (axial_strain - mUnReLoadCenter) > 0.;
        mAccumulatedStrain += std::abs(axial_strain - mUnReLoadCenter) -
                              (CalculateUnReLoadAmplitude(youngs_modulus) / 2.);
        mUnReLoadCenter = pos_side ? axial_strain - CalculateUnReLoadAmplitude(youngs_modulus) / 2.
                                   : axial_strain + CalculateUnReLoadAmplitude(youngs_modulus) / 2.;
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

double TrussBackboneConstitutiveLaw::BackboneStress(double Strain) const
{
    double young_modulus = 200.;
    return 0.5 * young_modulus * Strain;
}

double TrussBackboneConstitutiveLaw::BackboneStiffness([[maybe_unused]] double Strain) const
{
    double young_modulus = 200.;
    return 0.5 * young_modulus;
}

double TrussBackboneConstitutiveLaw::CalculateUnReLoadAmplitude(double YoungsModulus) const
{
    return 2. * BackboneStress(mAccumulatedStrain) / YoungsModulus;
}

bool TrussBackboneConstitutiveLaw::IsWithinUnReLoading(double Strain, double YoungsModulus) const
{
    return std::abs(Strain - mUnReLoadCenter) < (CalculateUnReLoadAmplitude(YoungsModulus) / 2.);
}

bool TrussBackboneConstitutiveLaw::RequiresInitializeMaterialResponse() { return false; }

SizeType TrussBackboneConstitutiveLaw::GetStrainSize() const { return 1; }

void TrussBackboneConstitutiveLaw::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, BaseType)
    rSerializer.save("AccumulatedStrain", mAccumulatedStrain);
    rSerializer.save("PreviousAxialStrain", mPreviousAxialStrain);
    rSerializer.save("UnReload", mUnReLoadCenter);
}

void TrussBackboneConstitutiveLaw::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, BaseType)
    rSerializer.load("AccumulatedStrain", mAccumulatedStrain);
    rSerializer.load("PreviousAxialStrain", mPreviousAxialStrain);
    rSerializer.load("UnReload", mUnReLoadCenter);
}
} // Namespace Kratos
