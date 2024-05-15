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

GeoThermalFilterLaw::GeoThermalFilterLaw() : GeoThermalLaw() {}

GeoThermalFilterLaw::GeoThermalFilterLaw(SizeType NumberOfDimensions)
    : GeoThermalLaw(NumberOfDimensions)
{
    KRATOS_ERROR_IF(mNumberOfDimensions != 1)
        << "Got invalid number of dimensions. The dimension has to be 1, but got: " << mNumberOfDimensions
        << std::endl;
}

GeoThermalLaw::Pointer GeoThermalFilterLaw::Clone() const
{
    return Kratos::make_shared<GeoThermalFilterLaw>(*this);
}

Matrix GeoThermalFilterLaw::CalculateThermalDispersionMatrix(const Properties& rProp, const ProcessInfo&) const
{
    return ScalarMatrix(1, 1, rProp[THERMAL_CONDUCTIVITY_WATER]);
}

} // Namespace Kratos
