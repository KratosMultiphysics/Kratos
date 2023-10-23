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

// System includes
#include <iostream>

// External includes

// Project includes
#include "custom_constitutive/thermal_dispersion_2D_law.hpp"

#include "geo_mechanics_application_variables.h"

namespace Kratos
{
    // ============================================================================================
    // ============================================================================================
    ConstitutiveLaw::Pointer GeoThermalDispersion2DLaw::Clone() const
    {
        return Kratos::make_shared<GeoThermalDispersion2DLaw>(*this);
    }

    // ============================================================================================
    // ============================================================================================
    void GeoThermalDispersion2DLaw::CalculateThermalDispersionMatrix(
        Matrix& C,
        const Properties& rProp)
    {
        KRATOS_TRY

        const double porosity = rProp[POROSITY];
        const double saturation = rProp[SATURATION];
        const double thermalCondWater = rProp[THERMAL_CONDUCTIVITY_WATER];
        const double thermalCondSolidXX = rProp[THERMAL_CONDUCTIVITY_SOLID_XX];
        const double thermalCondSolidXY = rProp[THERMAL_CONDUCTIVITY_SOLID_XY];
        const double thermalCondSolidYX = rProp[THERMAL_CONDUCTIVITY_SOLID_YX];
        const double thermalCondSolidYY = rProp[THERMAL_CONDUCTIVITY_SOLID_YY];

        const double cWater = porosity * saturation ;
        const double cSolid = 1.0 - porosity;
        
        C(INDEX_2D_HEAT_X, INDEX_2D_HEAT_X) = cSolid * thermalCondSolidXX + cWater * thermalCondWater;
        C(INDEX_2D_HEAT_X, INDEX_2D_HEAT_Y) = cSolid * thermalCondSolidXY;
        C(INDEX_2D_HEAT_Y, INDEX_2D_HEAT_X) = cSolid * thermalCondSolidYX;
        C(INDEX_2D_HEAT_Y, INDEX_2D_HEAT_Y) = cSolid * thermalCondSolidYY + cWater * thermalCondWater;
        
        KRATOS_CATCH("")
    }

} // Namespace Kratos
