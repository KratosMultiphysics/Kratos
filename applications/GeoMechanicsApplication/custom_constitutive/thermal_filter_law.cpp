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
//

#include "custom_constitutive/thermal_filter_law.h"
#include "geo_mechanics_application_variables.h"

namespace Kratos
{

GeoThermalFilterLaw::GeoThermalFilterLaw() : mNumberOfDimensions{1} {}

GeoThermalFilterLaw::GeoThermalFilterLaw(std::size_t NumberOfDimensions)
    : mNumberOfDimensions{NumberOfDimensions}
{
    KRATOS_ERROR_IF(mNumberOfDimensions != 1)
        << "Got invalid number of dimensions: " << mNumberOfDimensions << std::endl;
}

ConstitutiveLaw::Pointer GeoThermalFilterLaw::Clone() const
{
    return Kratos::make_shared<GeoThermalFilterLaw>(*this);
}

SizeType GeoThermalFilterLaw::WorkingSpaceDimension() { return 1; }

Matrix GeoThermalFilterLaw::CalculateThermalFilterMatrix(const Properties& rProp, const ProcessInfo& rProcessInfo) const
{
    KRATOS_TRY

    Matrix result = ZeroMatrix(1, 1);
    result(0, 0) = rProp[THERMAL_CONDUCTIVITY_WATER];
    return result;

    KRATOS_CATCH("")
}

} // Namespace Kratos
