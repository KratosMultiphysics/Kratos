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
#include "includes/global_variables.h"


namespace Kratos
{

Matrix GeoThermalFilterLaw::CalculateThermalDispersionMatrix(const Properties&  rProp,
                                                             const ProcessInfo& rProcessInfo) const
{
    return ScalarMatrix(1, 1, rProp[THERMAL_CONDUCTIVITY_WATER]);
}

Matrix GeoThermalFilterLaw::CalculateThermalDispersionMatrix(const Properties&  rProp,
                                                             const ProcessInfo& rProcessInfo,
                                                             const Vector&      rDischargeVector,
                                                             const double       waterDensity) const
{
    Matrix       result = this->CalculateThermalDispersionMatrix(rProp, rProcessInfo);
    const double equivalent_radius_square = rProp[CROSS_AREA] / Globals::Pi;
    result(0, 0) += waterDensity * rProp[SPECIFIC_HEAT_CAPACITY_WATER] * equivalent_radius_square *
                    rDischargeVector(0) * rDischargeVector(0) / (48.0 * rProp[THERMAL_CONDUCTIVITY_WATER]);
    return result;
}

} // Namespace Kratos
