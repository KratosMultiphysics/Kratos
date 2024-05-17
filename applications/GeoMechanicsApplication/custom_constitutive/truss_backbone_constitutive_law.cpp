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
        KRATOS_ERROR << "Can't get the specified value" << std::endl;
    }
    return rValue;
}

array_1d<double, 3>& TrussBackboneConstitutiveLaw::GetValue(const Variable<array_1d<double, 3>>& rThisVariable,
                                                            array_1d<double, 3>& rValue)
{
    KRATOS_ERROR << "Can't get the specified value" << std::endl;
    return rValue;
}

void TrussBackboneConstitutiveLaw::SetValue(const Variable<double>& rThisVariable,
                                            const double&           rValue,
                                            const ProcessInfo&      rCurrentProcessInfo)
{
    if (rThisVariable == ACCUMULATED_STRAIN) {
        mAccumulatedStrain = rValue;
    } else {
        KRATOS_ERROR << "Can't set the specified value" << std::endl;
    }
}

double& TrussBackboneConstitutiveLaw::CalculateValue(ConstitutiveLaw::Parameters& rParameterValues,
                                                     const Variable<double>&      rThisVariable,
                                                     double&                      rValue)
{
    if (rThisVariable == TANGENT_MODULUS) {
        const auto axial_strain(rParameterValues.GetStrainVector()[0]);
        const auto un_re_loading_strain_amplitude =
            2. * TrussBackboneConstitutiveLaw::BackboneStress(mAccumulatedStrain);
        if (axial_strain >= mAccumulatedStrain ||
            axial_strain <= mAccumulatedStrain - un_re_loading_strain_amplitude) {
            // backbone
            if (axial_strain > mAccumulatedStrain){
                rValue = TrussBackboneConstitutiveLaw::BackboneStiffness(axial_strain);
            } else {
                // on backbone opposite side
                auto accumulated_strain_increment = mAccumulatedStrain - un_re_loading_strain_amplitude - axial_strain;
                rValue = TrussBackboneConstitutiveLaw::BackboneStiffness(mAccumulatedStrain+accumulated_strain_increment);
            }
        } else {
            // unloading reloading
            rValue = rParameterValues.GetMaterialProperties()[YOUNG_MODULUS];
        }
    } else {
        KRATOS_ERROR << "Can't calculate the specified value" << std::endl;
    }
    KRATOS_INFO("TrussBackboneConstitutiveLaw::CalculateValue") << "stijfheid = " << rValue << std::endl;
    return rValue;
}

void TrussBackboneConstitutiveLaw::CalculateMaterialResponsePK2(Parameters& rValues)
{
    // backbone loading, un- and reloading laws
    const auto   axial_strain(rValues.GetStrainVector()[0]);
    auto&        axial_stress_vector = rValues.GetStressVector();
    const auto   un_re_loading_strain_amplitude =
        2. * TrussBackboneConstitutiveLaw::BackboneStress(mAccumulatedStrain);

    bool on_back_bone = axial_strain >= mAccumulatedStrain ||
                        axial_strain <= mAccumulatedStrain - un_re_loading_strain_amplitude;
    if (on_back_bone) {
        if (axial_strain >= mAccumulatedStrain){
            axial_stress_vector[0] = TrussBackboneConstitutiveLaw::BackboneStress(axial_strain);
        } else {
            // on backbone opposite side
            auto accumulated_strain_increment = mAccumulatedStrain - un_re_loading_strain_amplitude - axial_strain;
            axial_stress_vector[0] = -TrussBackboneConstitutiveLaw::BackboneStress(mAccumulatedStrain+accumulated_strain_increment);
        }
    } else {
        // unloading reloading
        axial_stress_vector[0] =
            TrussBackboneConstitutiveLaw::BackboneStress(mAccumulatedStrain) -
            rValues.GetMaterialProperties()[YOUNG_MODULUS] * (mAccumulatedStrain - axial_strain);
    }
    rValues.SetStressVector(axial_stress_vector);

/*
    if (axial_strain > mAccumulatedStrain) {
        // loading
        axial_stress_vector[0] = TrussBackboneConstitutiveLaw::BackboneStress(axial_strain);
        rValues.SetStressVector(axial_stress_vector);
    } else if (axial_strain >= mAccumulatedStrain - un_re_loading_strain_amplitude) {
        // unloading reloading
        axial_stress_vector[0] =
            TrussBackboneConstitutiveLaw::BackboneStress(mAccumulatedStrain) -
            rValues.GetMaterialProperties()[YOUNG_MODULUS] * (mAccumulatedStrain - axial_strain);
        rValues.SetStressVector(axial_stress_vector);
    } else {
        // loading opposite direction
        /// nog niet goed!!
        axial_stress_vector[0] = TrussBackboneConstitutiveLaw::BackboneStress(mAccumulatedStrain) -
                                 rValues.GetMaterialProperties()[YOUNG_MODULUS] *
                                     (mAccumulatedStrain - un_re_loading_strain_amplitude);
        rValues.SetStressVector(axial_stress_vector);
    }
*/
}

void TrussBackboneConstitutiveLaw::FinalizeMaterialResponsePK2(Parameters& rValues)
{
    CalculateMaterialResponsePK2(rValues);
    const auto axial_strain(rValues.GetStrainVector()[0]);
    const auto un_re_loading_strain_amplitude =
        2. * TrussBackboneConstitutiveLaw::BackboneStress(mAccumulatedStrain);

    if (axial_strain > mAccumulatedStrain) {
        mAccumulatedStrain = std::max(mAccumulatedStrain, axial_strain);
    } else if (axial_strain < mAccumulatedStrain - un_re_loading_strain_amplitude) {
        auto accumulated_strain_increment = std::abs(mAccumulatedStrain - un_re_loading_strain_amplitude - axial_strain);
        mAccumulatedStrain += accumulated_strain_increment;
    }
    /* mAccumulatedStrain = std::max(mAccumulatedStrain, axial_strain); */
    // mPreviousStrain = rValues.GetStrainVector()[0];
}

int TrussBackboneConstitutiveLaw::Check(const Properties&   rMaterialProperties,
                                        const GeometryType& rElementGeometry,
                                        const ProcessInfo&  rCurrentProcessInfo) const
{
    KRATOS_CHECK(rMaterialProperties.Has(YOUNG_MODULUS));
    // andere input voor backbone nalopen
    return 0;
}

double TrussBackboneConstitutiveLaw::BackboneStress(const double Strain) const
{
    double young_modulus = 200.;
    return 0.5 * young_modulus * Strain;
}

double TrussBackboneConstitutiveLaw::BackboneStiffness(const double Strain) const
{
    double young_modulus = 200.;
    return 0.5 * young_modulus;
}

} // Namespace Kratos
