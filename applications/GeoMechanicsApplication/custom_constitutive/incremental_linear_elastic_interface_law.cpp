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
#include "custom_geometries/line_interface_geometry.h"
#include "geo_mechanics_application_constants.h"
#include "geo_mechanics_application_variables.h"

namespace Kratos
{

ConstitutiveLaw::Pointer GeoIncrementalLinearElasticInterfaceLaw::Clone() const
{
    return std::make_shared<GeoIncrementalLinearElasticInterfaceLaw>(*this);
}

ConstitutiveLaw::SizeType GeoIncrementalLinearElasticInterfaceLaw::WorkingSpaceDimension()
{
    return 2;
}

ConstitutiveLaw::SizeType GeoIncrementalLinearElasticInterfaceLaw::GetStrainSize() const
{
    return VOIGT_SIZE_2D_INTERFACE;
}

Vector& GeoIncrementalLinearElasticInterfaceLaw::GetValue(const Variable<Vector>& rThisVariable, Vector& rValue)
{
    if (rThisVariable == STRAIN) {
        rValue = mPreviousRelativeDisplacement;
    } else if (rThisVariable == CAUCHY_STRESS_VECTOR) {
        rValue = mPreviousTraction;
    }

    return rValue;
}

Matrix& GeoIncrementalLinearElasticInterfaceLaw::CalculateValue(ConstitutiveLaw::Parameters& rParameterValues,
                                                                const Variable<Matrix>& rThisVariable,
                                                                Matrix& rValue)
{
    if (rThisVariable == CONSTITUTIVE_MATRIX) {
        const auto& r_properties = rParameterValues.GetMaterialProperties();
        rValue                   = MakeConstitutiveMatrix(r_properties[INTERFACE_NORMAL_STIFFNESS],
                                                          r_properties[INTERFACE_SHEAR_STIFFNESS]);
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
        prod(MakeConstitutiveMatrix(rValues.GetMaterialProperties()[INTERFACE_NORMAL_STIFFNESS],
                                    rValues.GetMaterialProperties()[INTERFACE_SHEAR_STIFFNESS]),
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

void GeoIncrementalLinearElasticInterfaceLaw::save(Serializer& rSerializer) const
{
    KRATOS_SERIALIZE_SAVE_BASE_CLASS(rSerializer, BaseType)
    rSerializer.save("PreviousRelativeDisplacement", mPreviousRelativeDisplacement);
    rSerializer.save("PreviousTraction", mPreviousTraction);
}

void GeoIncrementalLinearElasticInterfaceLaw::load(Serializer& rSerializer)
{
    KRATOS_SERIALIZE_LOAD_BASE_CLASS(rSerializer, BaseType)
    rSerializer.load("PreviousRelativeDisplacement", mPreviousRelativeDisplacement);
    rSerializer.load("PreviousTraction", mPreviousTraction);
}

Matrix GeoIncrementalLinearElasticInterfaceLaw::MakeConstitutiveMatrix(double NormalStiffness,
                                                                       double ShearStiffness) const
{
    auto result  = Matrix{ZeroMatrix{GetStrainSize(), GetStrainSize()}};
    result(0, 0) = NormalStiffness;
    result(1, 1) = ShearStiffness;
    return result;
}

} // namespace Kratos
