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

ConstitutiveLaw::SizeType GeoIncrementalLinearElasticInterfaceLaw::WorkingSpaceDimension()
{
    return 2;
}

ConstitutiveLaw::SizeType GeoIncrementalLinearElasticInterfaceLaw::GetStrainSize() const
{
    return VOIGT_SIZE_2D_INTERFACE;
}

ConstitutiveLaw::StressMeasure GeoIncrementalLinearElasticInterfaceLaw::GetStressMeasure()
{
    return ConstitutiveLaw::StressMeasure_Cauchy;
}

bool GeoIncrementalLinearElasticInterfaceLaw::IsIncremental() { return true; }

int GeoIncrementalLinearElasticInterfaceLaw::Check(const Properties& rMaterialProperties,
                                                   const ConstitutiveLaw::GeometryType& rElementGeometry,
                                                   const ProcessInfo& rCurrentProcessInfo) const
{
    KRATOS_ERROR_IF_NOT(rMaterialProperties.Has(INTERFACE_NORMAL_STIFFNESS))
        << "No interface normal stiffness defined" << std::endl;

    KRATOS_ERROR_IF_NOT(rMaterialProperties[INTERFACE_NORMAL_STIFFNESS] > 0.0)
        << "Interface normal stiffness must be positive, but got "
        << rMaterialProperties[INTERFACE_NORMAL_STIFFNESS] << std::endl;

    KRATOS_ERROR_IF_NOT(rMaterialProperties.Has(INTERFACE_SHEAR_STIFFNESS))
        << "No interface shear stiffness defined" << std::endl;

    KRATOS_ERROR_IF_NOT(rMaterialProperties[INTERFACE_SHEAR_STIFFNESS] > 0.0)
        << "Interface shear stiffness must be positive, but got "
        << rMaterialProperties[INTERFACE_SHEAR_STIFFNESS] << std::endl;

    KRATOS_ERROR_IF_NOT(dynamic_cast<const LineInterfaceGeometry*>(&rElementGeometry))
        << "Expected a line interface geometry, but got " << rElementGeometry.Info() << std::endl;

    return 0;
}

void GeoIncrementalLinearElasticInterfaceLaw::CalculateMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues)
{
    const auto& relative_displacement = rValues.GetStrainVector();
    auto&       traction              = rValues.GetStressVector();
    auto        constitutive_matrix   = Matrix{ZeroMatrix{GetStrainSize(), GetStrainSize()}};
    constitutive_matrix(0, 0)         = rValues.GetMaterialProperties()[INTERFACE_NORMAL_STIFFNESS];
    constitutive_matrix(1, 1)         = rValues.GetMaterialProperties()[INTERFACE_SHEAR_STIFFNESS];
    traction = mPreviousTraction + prod(constitutive_matrix, relative_displacement - mPreviousRelativeDisplacement);
}

void GeoIncrementalLinearElasticInterfaceLaw::InitializeMaterial(const Properties&,
                                                                 const ConstitutiveLaw::GeometryType&,
                                                                 const Vector&)
{
    mPreviousRelativeDisplacement = ZeroVector{GetStrainSize()};
    mPreviousTraction             = ZeroVector{GetStrainSize()};
}

bool GeoIncrementalLinearElasticInterfaceLaw::RequiresInitializeMaterialResponse()
{
    return false;
}

void GeoIncrementalLinearElasticInterfaceLaw::FinalizeMaterialResponseCauchy(ConstitutiveLaw::Parameters& rValues)
{
    mPreviousRelativeDisplacement = rValues.GetStrainVector();
    mPreviousTraction             = rValues.GetStressVector();
}

} // namespace Kratos
