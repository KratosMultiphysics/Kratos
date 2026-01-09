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
//                   Anne van de Graaf
//

#include "incremental_linear_elastic_interface_law.h"
#include "constitutive_law_dimension.h"
#include "custom_utilities/check_utilities.h"
#include "geo_mechanics_application_variables.h"

namespace Kratos
{

GeoIncrementalLinearElasticInterfaceLaw::~GeoIncrementalLinearElasticInterfaceLaw() = default;

GeoIncrementalLinearElasticInterfaceLaw::GeoIncrementalLinearElasticInterfaceLaw(std::unique_ptr<ConstitutiveLawDimension> pConstitutiveLawDimension)
    : mpConstitutiveLawDimension(std::move(pConstitutiveLawDimension))
{
}

ConstitutiveLaw::Pointer GeoIncrementalLinearElasticInterfaceLaw::Clone() const
{
    return std::make_shared<GeoIncrementalLinearElasticInterfaceLaw>(mpConstitutiveLawDimension->Clone());
}

ConstitutiveLaw::SizeType GeoIncrementalLinearElasticInterfaceLaw::WorkingSpaceDimension()
{
    return mpConstitutiveLawDimension->GetDimension();
}

ConstitutiveLaw::SizeType GeoIncrementalLinearElasticInterfaceLaw::GetStrainSize() const
{
    return mpConstitutiveLawDimension->GetStrainSize();
}

Vector& GeoIncrementalLinearElasticInterfaceLaw::GetValue(const Variable<Vector>& rThisVariable, Vector& rValue)
{
    if (rThisVariable == GEO_RELATIVE_DISPLACEMENT_VECTOR) {
        rValue = mPreviousRelativeDisplacement;
    } else if (rThisVariable == GEO_EFFECTIVE_TRACTION_VECTOR) {
        rValue = mPreviousTraction;
    } else {
        KRATOS_ERROR << "Can't get value of " << rThisVariable.Name() << ": unsupported variable\n";
    }

    return rValue;
}

Matrix& GeoIncrementalLinearElasticInterfaceLaw::CalculateValue(ConstitutiveLaw::Parameters& rParameterValues,
                                                                const Variable<Matrix>& rThisVariable,
                                                                Matrix& rValue)
{
    if (rThisVariable == CONSTITUTIVE_MATRIX) {
        rValue = mpConstitutiveLawDimension->CalculateElasticMatrix(rParameterValues.GetMaterialProperties());
    } else {
        KRATOS_ERROR << "Can't calculate value of " << rThisVariable.Name() << ": unsupported variable\n";
    }

    return rValue;
}

ConstitutiveLaw::StressMeasure GeoIncrementalLinearElasticInterfaceLaw::GetStressMeasure()
{
    return ConstitutiveLaw::StressMeasure_Cauchy;
}

bool GeoIncrementalLinearElasticInterfaceLaw::IsIncremental() { return true; }

void GeoIncrementalLinearElasticInterfaceLaw::InitializeMaterial(const Properties&,
                                                                 const ConstitutiveLaw::GeometryType&,
                                                                 const Vector&)
{
    mPreviousRelativeDisplacement =
        HasInitialState() ? GetInitialState().GetInitialStrainVector() : ZeroVector{GetStrainSize()};
    mPreviousTraction =
        HasInitialState() ? GetInitialState().GetInitialStressVector() : ZeroVector{GetStrainSize()};
}

void GeoIncrementalLinearElasticInterfaceLaw::CalculateMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues)
{
    rValues.GetStressVector() =
        mPreviousTraction +
        prod(mpConstitutiveLawDimension->CalculateElasticMatrix(rValues.GetMaterialProperties()),
             rValues.GetStrainVector() - mPreviousRelativeDisplacement);
}

bool GeoIncrementalLinearElasticInterfaceLaw::RequiresInitializeMaterialResponse() { return false; }

void GeoIncrementalLinearElasticInterfaceLaw::FinalizeMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues)
{
    mPreviousRelativeDisplacement = rValues.GetStrainVector();
    mPreviousTraction             = rValues.GetStressVector();
}

int GeoIncrementalLinearElasticInterfaceLaw::Check(const Properties& rMaterialProperties,
                                                   const ConstitutiveLaw::GeometryType& rElementGeometry,
                                                   const ProcessInfo& rCurrentProcessInfo) const
{
    const auto result = BaseType::Check(rMaterialProperties, rElementGeometry, rCurrentProcessInfo);

    const CheckProperties check_properties(rMaterialProperties, "material properties",
                                           CheckProperties::Bounds::AllExclusive);
    check_properties.Check(INTERFACE_NORMAL_STIFFNESS);
    check_properties.Check(INTERFACE_SHEAR_STIFFNESS);

    return result;
}

void GeoIncrementalLinearElasticInterfaceLaw::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, BaseType)
    rSerializer.save("PreviousRelativeDisplacement", mPreviousRelativeDisplacement);
    rSerializer.save("PreviousTraction", mPreviousTraction);
    rSerializer.save("ConstitutiveLawDimension", mpConstitutiveLawDimension);
}

void GeoIncrementalLinearElasticInterfaceLaw::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, BaseType)
    rSerializer.load("PreviousRelativeDisplacement", mPreviousRelativeDisplacement);
    rSerializer.load("PreviousTraction", mPreviousTraction);
    rSerializer.load("ConstitutiveLawDimension", mpConstitutiveLawDimension);
}

// Instances of this class can not be copied but can be moved. Check that at compile time.
static_assert(!std::is_copy_constructible_v<GeoIncrementalLinearElasticInterfaceLaw>);
static_assert(!std::is_copy_assignable_v<GeoIncrementalLinearElasticInterfaceLaw>);
static_assert(std::is_move_constructible_v<GeoIncrementalLinearElasticInterfaceLaw>);
static_assert(std::is_move_assignable_v<GeoIncrementalLinearElasticInterfaceLaw>);
} // namespace Kratos
