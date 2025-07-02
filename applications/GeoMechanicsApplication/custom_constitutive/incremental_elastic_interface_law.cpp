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

#include "incremental_elastic_interface_law.h"
#include "custom_geometries/line_interface_geometry.h"
#include "custom_utilities/constitutive_law_utilities.h"
#include "geo_mechanics_application_constants.h"
#include "geo_mechanics_application_variables.h"

namespace Kratos
{

ConstitutiveLaw::Pointer GeoIncrementalElasticInterfaceLaw::Clone() const
{
    return std::make_shared<GeoIncrementalElasticInterfaceLaw>(mConstitutiveLawDimension->Clone());
}

ConstitutiveLaw::SizeType GeoIncrementalElasticInterfaceLaw::WorkingSpaceDimension()
{
    return mConstitutiveLawDimension->GetDimension();
}

ConstitutiveLaw::SizeType GeoIncrementalElasticInterfaceLaw::GetStrainSize() const
{
    return mConstitutiveLawDimension->GetStrainSize();
}

Vector& GeoIncrementalElasticInterfaceLaw::GetValue(const Variable<Vector>& rThisVariable, Vector& rValue)
{
    if (rThisVariable == STRAIN) {
        rValue = mPreviousRelativeDisplacement;
    } else if (rThisVariable == CAUCHY_STRESS_VECTOR) {
        rValue = mPreviousTraction;
    } else {
        KRATOS_ERROR << "Can't get value of " << rThisVariable.Name() << ": unsupported variable\n";
    }

    return rValue;
}

Matrix& GeoIncrementalElasticInterfaceLaw::CalculateValue(ConstitutiveLaw::Parameters& rParameterValues,
                                                                const Variable<Matrix>& rThisVariable,
                                                                Matrix& rValue)
{
    if (rThisVariable == CONSTITUTIVE_MATRIX) {
        const auto& r_properties = rParameterValues.GetMaterialProperties();
        rValue                   = ConstitutiveLawUtilities::MakeInterfaceConstitutiveMatrix(
            r_properties[INTERFACE_NORMAL_STIFFNESS], r_properties[INTERFACE_SHEAR_STIFFNESS], GetStrainSize());
    } else {
        KRATOS_ERROR << "Can't calculate value of " << rThisVariable.Name() << ": unsupported variable\n";
    }

    return rValue;
}

ConstitutiveLaw::StressMeasure GeoIncrementalElasticInterfaceLaw::GetStressMeasure()
{
    return ConstitutiveLaw::StressMeasure_Cauchy;
}

bool GeoIncrementalElasticInterfaceLaw::IsIncremental() { return true; }

void GeoIncrementalElasticInterfaceLaw::InitializeMaterial(const Properties&,
                                                                 const ConstitutiveLaw::GeometryType&,
                                                                 const Vector&)
{
    mPreviousRelativeDisplacement =
        HasInitialState() ? GetInitialState().GetInitialStrainVector() : ZeroVector{GetStrainSize()};
    mPreviousTraction =
        HasInitialState() ? GetInitialState().GetInitialStressVector() : ZeroVector{GetStrainSize()};
}

void GeoIncrementalElasticInterfaceLaw::CalculateMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues)
{
    rValues.GetStressVector() =
        mPreviousTraction +
        prod(ConstitutiveLawUtilities::MakeInterfaceConstitutiveMatrix(
                 rValues.GetMaterialProperties()[INTERFACE_NORMAL_STIFFNESS],
                 rValues.GetMaterialProperties()[INTERFACE_SHEAR_STIFFNESS], GetStrainSize()),
             rValues.GetStrainVector() - mPreviousRelativeDisplacement);
}

bool GeoIncrementalElasticInterfaceLaw::RequiresInitializeMaterialResponse() { return false; }

void GeoIncrementalElasticInterfaceLaw::FinalizeMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues)
{
    mPreviousRelativeDisplacement = rValues.GetStrainVector();
    mPreviousTraction             = rValues.GetStressVector();
}

int GeoIncrementalElasticInterfaceLaw::Check(const Properties& rMaterialProperties,
                                                   const ConstitutiveLaw::GeometryType& rElementGeometry,
                                                   const ProcessInfo& rCurrentProcessInfo) const
{
    const auto result = BaseType::Check(rMaterialProperties, rElementGeometry, rCurrentProcessInfo);

    KRATOS_ERROR_IF_NOT(rMaterialProperties.Has(INTERFACE_NORMAL_STIFFNESS))
        << "No interface normal stiffness is defined" << std::endl;

    KRATOS_ERROR_IF_NOT(rMaterialProperties[INTERFACE_NORMAL_STIFFNESS] > 0.0)
        << "Interface normal stiffness must be positive, but got "
        << rMaterialProperties[INTERFACE_NORMAL_STIFFNESS] << std::endl;

    KRATOS_ERROR_IF_NOT(rMaterialProperties.Has(INTERFACE_SHEAR_STIFFNESS))
        << "No interface shear stiffness is defined" << std::endl;

    KRATOS_ERROR_IF_NOT(rMaterialProperties[INTERFACE_SHEAR_STIFFNESS] > 0.0)
        << "Interface shear stiffness must be positive, but got "
        << rMaterialProperties[INTERFACE_SHEAR_STIFFNESS] << std::endl;

    return result;
}

void GeoIncrementalElasticInterfaceLaw::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, BaseType)
    rSerializer.save("PreviousRelativeDisplacement", mPreviousRelativeDisplacement);
    rSerializer.save("PreviousTraction", mPreviousTraction);
}

void GeoIncrementalElasticInterfaceLaw::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, BaseType)
    rSerializer.load("PreviousRelativeDisplacement", mPreviousRelativeDisplacement);
    rSerializer.load("PreviousTraction", mPreviousTraction);
}

} // namespace Kratos
