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

Matrix GeoThermalFilterLaw::CalculateThermalDispersionMatrix(const Properties&  rProp,
                                                             const ProcessInfo& rProcessInfo) const
{
    return ScalarMatrix(1, 1, rProp[THERMAL_CONDUCTIVITY_WATER]);
}

Matrix GeoThermalFilterLaw::CalculateThermalDispersionMatrix(const Properties& rProp,
                                                             const ProcessInfo& rProcessInfo,
                                                             const Vector& rDischargeVector) const
{
    Matrix result = this->CalculateThermalDispersionMatrix(rProp, rProcessInfo);
    // Add discharge here
    return result;
}

} // Namespace Kratos
