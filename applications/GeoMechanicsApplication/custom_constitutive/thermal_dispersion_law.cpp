// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Mohamed Nabi
//                   John van Esch
//                   Gennady Markelov
//

#include "custom_constitutive/thermal_dispersion_law.hpp"
#include "custom_retention/retention_law_factory.h"
#include "geo_mechanics_application_variables.h"
#include <iostream>

namespace Kratos {

GeoThermalDispersionLaw::GeoThermalDispersionLaw(std::size_t NumberOfDimensions)
    : mNumberOfDimensions{NumberOfDimensions}
{
}

ConstitutiveLaw::Pointer GeoThermalDispersionLaw::Clone() const
{
    return Kratos::make_shared<GeoThermalDispersionLaw>(*this);
}

SizeType GeoThermalDispersionLaw::WorkingSpaceDimension()
{
    return mNumberOfDimensions;
}

Matrix GeoThermalDispersionLaw::CalculateThermalDispersionMatrix(
    const Properties& rProp, const ProcessInfo& rProcessInfo, const Geometry<Node>& rElementGeometry) const
{
    KRATOS_TRY

    Matrix result = ZeroMatrix(mNumberOfDimensions, mNumberOfDimensions);

    RetentionLaw::Parameters parameters(rElementGeometry, rProp, rProcessInfo);
    auto retention_law = RetentionLawFactory::Clone(rProp);
    const double saturation = retention_law->CalculateSaturation(parameters);
    const double coefficient_water = rProp[POROSITY] * saturation;
    const double coefficient_solid = 1.0 - rProp[POROSITY];

    result(INDEX_THERMAL_FLUX_X, INDEX_THERMAL_FLUX_X) =
        coefficient_solid * rProp[THERMAL_CONDUCTIVITY_SOLID_XX] +
        coefficient_water * rProp[THERMAL_CONDUCTIVITY_WATER];
    result(INDEX_THERMAL_FLUX_X, INDEX_THERMAL_FLUX_Y) =
        coefficient_solid * rProp[THERMAL_CONDUCTIVITY_SOLID_XY];
    result(INDEX_THERMAL_FLUX_Y, INDEX_THERMAL_FLUX_Y) =
        coefficient_solid * rProp[THERMAL_CONDUCTIVITY_SOLID_YY] +
        coefficient_water * rProp[THERMAL_CONDUCTIVITY_WATER];
    result(INDEX_THERMAL_FLUX_Y, INDEX_THERMAL_FLUX_X) =
        result(INDEX_THERMAL_FLUX_X, INDEX_THERMAL_FLUX_Y);

    if (mNumberOfDimensions == 3) {
        result(INDEX_THERMAL_FLUX_Y, INDEX_THERMAL_FLUX_Z) =
            coefficient_solid * rProp[THERMAL_CONDUCTIVITY_SOLID_YZ];
        result(INDEX_THERMAL_FLUX_Z, INDEX_THERMAL_FLUX_X) =
            coefficient_solid * rProp[THERMAL_CONDUCTIVITY_SOLID_XZ];
        result(INDEX_THERMAL_FLUX_Z, INDEX_THERMAL_FLUX_Z) =
            coefficient_solid * rProp[THERMAL_CONDUCTIVITY_SOLID_ZZ] +
            coefficient_water * rProp[THERMAL_CONDUCTIVITY_WATER];
        result(INDEX_THERMAL_FLUX_Z, INDEX_THERMAL_FLUX_Y) =
            result(INDEX_THERMAL_FLUX_Y, INDEX_THERMAL_FLUX_Z);
        result(INDEX_THERMAL_FLUX_X, INDEX_THERMAL_FLUX_Z) =
            result(INDEX_THERMAL_FLUX_Z, INDEX_THERMAL_FLUX_X);
    }

    return result;

    KRATOS_CATCH("")
}

} // Namespace Kratos
