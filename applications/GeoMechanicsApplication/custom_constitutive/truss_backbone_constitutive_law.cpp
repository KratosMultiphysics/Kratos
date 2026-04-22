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
using namespace std::string_literals;

ConstitutiveLaw::Pointer TrussBackboneConstitutiveLaw::Clone() const
{
    return make_shared<TrussBackboneConstitutiveLaw>(*this);
}

void TrussBackboneConstitutiveLaw::GetLawFeatures(Features& rFeatures)
{
    rFeatures.mStrainSize     = 1;
    rFeatures.mSpaceDimension = 3;
}

double& TrussBackboneConstitutiveLaw::GetValue(const Variable<double>& rThisVariable, double& rValue)
{
    if (rThisVariable == ACCUMULATED_STRAIN) {
        rValue = mAccumulatedStrain;
    } else {
        rValue = BaseType::GetValue(rThisVariable, rValue);
    }

    return rValue;
}

void TrussBackboneConstitutiveLaw::SetValue(const Variable<double>& rThisVariable,
                                            const double&           rValue,
                                            const ProcessInfo&      rCurrentProcessInfo)
{
    if (rThisVariable == ACCUMULATED_STRAIN) {
        mAccumulatedStrain = rValue;
    } else {
        BaseType::SetValue(rThisVariable, rValue, rCurrentProcessInfo);
    }
}

double& TrussBackboneConstitutiveLaw::CalculateValue(ConstitutiveLaw::Parameters& rParameterValues,
                                                     const Variable<double>&      rThisVariable,
                                                     double&                      rValue)
{
    if (rThisVariable == TANGENT_MODULUS) {
        const auto youngs_modulus = rParameterValues.GetMaterialProperties()[YOUNG_MODULUS];
        const auto axial_strain   = rParameterValues.GetStrainVector()[0];
        if (IsWithinUnReLoading(axial_strain, youngs_modulus)) {
            rValue = youngs_modulus;
        } else {
            rValue = BackboneStiffness(axial_strain);
        }
    } else {
        rValue = BaseType::CalculateValue(rParameterValues, rThisVariable, rValue);
    }

    return rValue;
}

void TrussBackboneConstitutiveLaw::CalculateMaterialResponsePK2(Parameters& rValues)
{
    // backbone loading, un- and reloading laws
    const auto axial_strain        = rValues.GetStrainVector()[0];
    auto&      axial_stress_vector = rValues.GetStressVector();
    const auto youngs_modulus      = rValues.GetMaterialProperties()[YOUNG_MODULUS];
    if (IsWithinUnReLoading(axial_strain, youngs_modulus)) {
        axial_stress_vector[0] = youngs_modulus * (axial_strain - mUnReLoadCenter);
        return;
    }

    // backbone
    axial_stress_vector[0] =
        BackboneStress(mAccumulatedStrain + (std::abs(axial_strain - mUnReLoadCenter) -
                                             CalculateUnReLoadAmplitude(youngs_modulus)));
    if (axial_strain - mUnReLoadCenter < 0.0) axial_stress_vector[0] *= -1.0;
}

void TrussBackboneConstitutiveLaw::FinalizeMaterialResponsePK2(Parameters& rValues)
{
    CalculateMaterialResponsePK2(rValues);
    const auto axial_strain = rValues.GetStrainVector()[0];
    if (const auto youngs_modulus = rValues.GetMaterialProperties()[YOUNG_MODULUS];
        !IsWithinUnReLoading(axial_strain, youngs_modulus)) {
        // update accumulated backbone strain and un- reload center strain
        const auto difference = axial_strain - mUnReLoadCenter;
        mAccumulatedStrain += std::abs(difference) - CalculateUnReLoadAmplitude(youngs_modulus);
        mUnReLoadCenter = difference > 0.0 ? axial_strain - CalculateUnReLoadAmplitude(youngs_modulus)
                                           : axial_strain + CalculateUnReLoadAmplitude(youngs_modulus);
    }
    mPreviousAxialStrain = axial_strain;
}

int TrussBackboneConstitutiveLaw::Check(const Properties&   rMaterialProperties,
                                        const GeometryType& rElementGeometry,
                                        const ProcessInfo&  rCurrentProcessInfo) const
{
    if (const auto exit_code = BaseType::Check(rMaterialProperties, rElementGeometry, rCurrentProcessInfo);
        exit_code != 0)
        return exit_code;

    KRATOS_ERROR_IF_NOT(rMaterialProperties.Has(YOUNG_MODULUS)) << "No YOUNGS_MODULUS found" << std::endl;

    CheckStressStrainDiagram(rMaterialProperties);

    return 0;
}

double TrussBackboneConstitutiveLaw::BackboneStress(double Strain) const
{
    return mStressStrainTable.GetValue(Strain);
}

double TrussBackboneConstitutiveLaw::BackboneStiffness(double Strain) const
{
    return mStressStrainTable.GetDerivative(Strain);
}

double TrussBackboneConstitutiveLaw::CalculateUnReLoadAmplitude(double YoungsModulus) const
{
    return BackboneStress(mAccumulatedStrain) / YoungsModulus;
}

bool TrussBackboneConstitutiveLaw::IsWithinUnReLoading(double Strain, double YoungsModulus) const
{
    return std::abs(Strain - mUnReLoadCenter) < CalculateUnReLoadAmplitude(YoungsModulus);
}

bool TrussBackboneConstitutiveLaw::RequiresInitializeMaterialResponse() { return true; }

SizeType TrussBackboneConstitutiveLaw::GetStrainSize() const { return 1; }

void TrussBackboneConstitutiveLaw::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, BaseType)
    rSerializer.save("AccumulatedStrain", mAccumulatedStrain);
    rSerializer.save("PreviousAxialStrain", mPreviousAxialStrain);
    rSerializer.save("UnReload", mUnReLoadCenter);
    rSerializer.save("StressStrainTable", mStressStrainTable);
}

void TrussBackboneConstitutiveLaw::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, BaseType)
    rSerializer.load("AccumulatedStrain", mAccumulatedStrain);
    rSerializer.load("PreviousAxialStrain", mPreviousAxialStrain);
    rSerializer.load("UnReload", mUnReLoadCenter);
    rSerializer.load("StressStrainTable", mStressStrainTable);
}

void TrussBackboneConstitutiveLaw::InitializeMaterial(const Properties& rMaterialProperties,
                                                      const ConstitutiveLaw::GeometryType& rElementGeometry,
                                                      const Vector& rShapeFunctionsValues)
{
    BaseType::InitializeMaterial(rMaterialProperties, rElementGeometry, rShapeFunctionsValues);

    mStressStrainTable.Clear();

    const auto& r_strains_of_piecewise_linear_law = rMaterialProperties[STRAINS_OF_PIECEWISE_LINEAR_LAW];
    const auto& r_stresses_of_piecewise_linear_law = rMaterialProperties[STRESSES_OF_PIECEWISE_LINEAR_LAW];
    for (auto i = std::size_t{0}; i < r_strains_of_piecewise_linear_law.size(); ++i) {
        mStressStrainTable.PushBack(r_strains_of_piecewise_linear_law[i],
                                    r_stresses_of_piecewise_linear_law[i]);
    }
}

std::string TrussBackboneConstitutiveLaw::Info() const { return "TrussBackboneConstitutiveLaw"s; }

void TrussBackboneConstitutiveLaw::CheckStressStrainDiagram(const Properties& rMaterialProperties)
{
    KRATOS_ERROR_IF_NOT(rMaterialProperties.Has(STRAINS_OF_PIECEWISE_LINEAR_LAW))
        << "No STRAINS_OF_PIECEWISE_LINEAR_LAW found" << std::endl;
    KRATOS_ERROR_IF_NOT(rMaterialProperties.Has(STRESSES_OF_PIECEWISE_LINEAR_LAW))
        << "No STRESSES_OF_PIECEWISE_LINEAR_LAW found" << std::endl;

    KRATOS_ERROR_IF_NOT(rMaterialProperties[STRAINS_OF_PIECEWISE_LINEAR_LAW].size() ==
                        rMaterialProperties[STRESSES_OF_PIECEWISE_LINEAR_LAW].size())
        << "The number of strain components does not match the number of stress components" << std::endl;
    KRATOS_ERROR_IF(rMaterialProperties[STRAINS_OF_PIECEWISE_LINEAR_LAW].empty())
        << "STRAINS_OF_PIECEWISE_LINEAR_LAW is empty";

    CheckStrainValuesAreAscending(rMaterialProperties[STRAINS_OF_PIECEWISE_LINEAR_LAW]);

    CheckBackboneStiffnessesDontExceedYoungsModulus(
        rMaterialProperties[STRAINS_OF_PIECEWISE_LINEAR_LAW],
        rMaterialProperties[STRESSES_OF_PIECEWISE_LINEAR_LAW], rMaterialProperties[YOUNG_MODULUS]);
}

void TrussBackboneConstitutiveLaw::CheckStrainValuesAreAscending(const Vector& rStrains)
{
    auto first_ge_second = [](const auto& First, const auto& Second) { return First >= Second; };
    auto pos             = std::adjacent_find(rStrains.cbegin(), rStrains.cend(), first_ge_second);
    KRATOS_ERROR_IF(pos != rStrains.cend())
        << "Values in STRAINS_OF_PIECEWISE_LINEAR_LAW are not ascending: " << *(pos + 1)
        << " (at index " << std::distance(rStrains.begin(), pos) + 2 << ") does not exceed " << *pos
        << " (at index " << std::distance(rStrains.begin(), pos) + 1 << ")" << std::endl;
}

void TrussBackboneConstitutiveLaw::CheckBackboneStiffnessesDontExceedYoungsModulus(const Vector& rStrains,
                                                                                   const Vector& rStresses,
                                                                                   double YoungsModulus)
{
    if (rStrains.size() > 1) { // only check for potentially non-constant laws
        for (auto i = std::size_t{0}; i < rStrains.size() - 1; ++i) {
            const auto backbone_stiffness =
                (rStresses[i + 1] - rStresses[i]) / (rStrains[i + 1] - rStrains[i]);
            KRATOS_ERROR_IF(YoungsModulus < backbone_stiffness)
                << "YOUNGS_MODULUS (" << YoungsModulus << ") is smaller than the backbone stiffness of line segment "
                << i + 1 << " (" << backbone_stiffness << ")" << std::endl;
        }
    }
}

} // Namespace Kratos
