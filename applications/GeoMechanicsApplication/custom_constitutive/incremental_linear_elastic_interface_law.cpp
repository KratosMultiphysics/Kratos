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

    return 0;
}

} // namespace Kratos
