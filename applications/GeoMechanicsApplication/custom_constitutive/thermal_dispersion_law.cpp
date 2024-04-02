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

#include "custom_constitutive/thermal_dispersion_law.h"
#include "custom_retention/retention_law_factory.h"
#include "geo_mechanics_application_variables.h"

namespace Kratos
{

GeoThermalDispersionLaw::GeoThermalDispersionLaw() : mNumberOfDimensions{2} {}

GeoThermalDispersionLaw::GeoThermalDispersionLaw(std::size_t NumberOfDimensions)
    : mNumberOfDimensions{NumberOfDimensions}
{
    KRATOS_ERROR_IF(mNumberOfDimensions < 2 || mNumberOfDimensions > 3)
        << "Got invalid number of dimensions: " << mNumberOfDimensions << std::endl;
}

ConstitutiveLaw::Pointer GeoThermalDispersionLaw::Clone() const
{
    return Kratos::make_shared<GeoThermalDispersionLaw>(*this);
}

SizeType GeoThermalDispersionLaw::WorkingSpaceDimension() { return mNumberOfDimensions; }

Matrix GeoThermalDispersionLaw::CalculateThermalDispersionMatrix(const Properties& rProp,
                                                                 const ProcessInfo& rProcessInfo) const
{
    KRATOS_TRY

    Matrix result = ZeroMatrix(mNumberOfDimensions, mNumberOfDimensions);

    RetentionLaw::Parameters parameters(rProp, rProcessInfo);
    auto                     retention_law  = RetentionLawFactory::Clone(rProp);
    const double             saturation     = retention_law->CalculateSaturation(parameters);
    const double             water_fraction = rProp[POROSITY] * saturation;
    const double             solid_fraction = 1.0 - rProp[POROSITY];

    const auto x = static_cast<int>(indexThermalFlux::X);
    const auto y = static_cast<int>(indexThermalFlux::Y);

    result(x, x) = solid_fraction * rProp[THERMAL_CONDUCTIVITY_SOLID_XX] +
                   water_fraction * rProp[THERMAL_CONDUCTIVITY_WATER];
    result(x, y) = solid_fraction * rProp[THERMAL_CONDUCTIVITY_SOLID_XY];
    result(y, y) = solid_fraction * rProp[THERMAL_CONDUCTIVITY_SOLID_YY] +
                   water_fraction * rProp[THERMAL_CONDUCTIVITY_WATER];
    result(y, x) = result(x, y);

    if (mNumberOfDimensions == 3) {
        const auto z = static_cast<int>(indexThermalFlux::Z);
        result(y, z) = solid_fraction * rProp[THERMAL_CONDUCTIVITY_SOLID_YZ];
        result(z, x) = solid_fraction * rProp[THERMAL_CONDUCTIVITY_SOLID_XZ];
        result(z, z) = solid_fraction * rProp[THERMAL_CONDUCTIVITY_SOLID_ZZ] +
                       water_fraction * rProp[THERMAL_CONDUCTIVITY_WATER];
        result(z, y) = result(y, z);
        result(x, z) = result(z, x);
    }

    return result;

    KRATOS_CATCH("")
}

} // Namespace Kratos
