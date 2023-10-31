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

/*        const double cWater = rProp[POROSITY] * rProp[SATURATION];
        const double cSolid = 1.0 - rProp[POROSITY];

        C(INDEX_2D_THERMAL_FLUX_X, INDEX_2D_THERMAL_FLUX_X) = cSolid * rProp[THERMAL_CONDUCTIVITY_SOLID_XX]
                                                            + cWater * rProp[THERMAL_CONDUCTIVITY_WATER];
        C(INDEX_2D_THERMAL_FLUX_X, INDEX_2D_THERMAL_FLUX_Y) = cSolid * rProp[THERMAL_CONDUCTIVITY_SOLID_XY];
        C(INDEX_2D_THERMAL_FLUX_Y, INDEX_2D_THERMAL_FLUX_Y) = cSolid * rProp[THERMAL_CONDUCTIVITY_SOLID_YY]
                                                            + cWater * rProp[THERMAL_CONDUCTIVITY_WATER];
        C(INDEX_2D_THERMAL_FLUX_Y, INDEX_2D_THERMAL_FLUX_X) = C(INDEX_2D_THERMAL_FLUX_X, INDEX_2D_THERMAL_FLUX_Y);
        
        if (C.size1() == 3 && C.size2() == 3) {
            C(INDEX_2D_THERMAL_FLUX_Y, INDEX_2D_THERMAL_FLUX_Z) = cSolid * rProp[THERMAL_CONDUCTIVITY_SOLID_YZ];
            C(INDEX_2D_THERMAL_FLUX_Z, INDEX_2D_THERMAL_FLUX_X) = cSolid * rProp[THERMAL_CONDUCTIVITY_SOLID_ZX];
            C(INDEX_2D_THERMAL_FLUX_Z, INDEX_2D_THERMAL_FLUX_Z) = cSolid * rProp[THERMAL_CONDUCTIVITY_SOLID_ZZ]
                                                                + cWater * rProp[THERMAL_CONDUCTIVITY_WATER];
            C(INDEX_2D_THERMAL_FLUX_Z, INDEX_2D_THERMAL_FLUX_Y) = C(INDEX_2D_THERMAL_FLUX_Y, INDEX_2D_THERMAL_FLUX_Z);
            C(INDEX_2D_THERMAL_FLUX_X, INDEX_2D_THERMAL_FLUX_Z) = C(INDEX_2D_THERMAL_FLUX_Z, INDEX_2D_THERMAL_FLUX_X);
        }*/

        KRATOS_CATCH("")
    }

} // Namespace Kratos
