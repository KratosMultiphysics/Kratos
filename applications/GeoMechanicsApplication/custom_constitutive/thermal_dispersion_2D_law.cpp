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
        const Properties& rProp,
        const bool isPressureCoupled,
        const Vector& dischargeVector)
    {
        KRATOS_TRY

        const double cWater = rProp[POROSITY] * rProp[SATURATION];
        const double cSolid = 1.0 - rProp[POROSITY];
        
        C(INDEX_2D_THERMAL_FLUX_X, INDEX_2D_THERMAL_FLUX_X) = cSolid * rProp[THERMAL_CONDUCTIVITY_SOLID_XX]
                                                            + cWater * rProp[THERMAL_CONDUCTIVITY_WATER];
        C(INDEX_2D_THERMAL_FLUX_X, INDEX_2D_THERMAL_FLUX_Y) = cSolid * rProp[THERMAL_CONDUCTIVITY_SOLID_XY];
        C(INDEX_2D_THERMAL_FLUX_Y, INDEX_2D_THERMAL_FLUX_X) = cSolid * rProp[THERMAL_CONDUCTIVITY_SOLID_YX];
        C(INDEX_2D_THERMAL_FLUX_Y, INDEX_2D_THERMAL_FLUX_Y) = cSolid * rProp[THERMAL_CONDUCTIVITY_SOLID_YY]
                                                            + cWater * rProp[THERMAL_CONDUCTIVITY_WATER];

        if (isPressureCoupled) {
            double totalDischarge = MathUtils<>::Norm(dischargeVector);
            if (totalDischarge > 0.0) {
                const double c2 = rProp[DENSITY_WATER] * rProp[SPECIFIC_HEAT_CAPACITY_WATER];
                const double c3 = (rProp[LONGITUDINAL_DISPERSIVITY] - rProp[TRANSVERSE_DISPERSIVITY]) / totalDischarge;
                const double c4 = rProp[SOLID_COMPRESSIBILITY] * totalDischarge;
                C(INDEX_2D_THERMAL_FLUX_X, INDEX_2D_THERMAL_FLUX_X) += c2 * (c3 * dischargeVector[0] * dischargeVector[0] + c4);
                C(INDEX_2D_THERMAL_FLUX_X, INDEX_2D_THERMAL_FLUX_Y) += c2 * c3 * dischargeVector[0] * dischargeVector[1];
                C(INDEX_2D_THERMAL_FLUX_Y, INDEX_2D_THERMAL_FLUX_X) += c2 * c3 * dischargeVector[0] * dischargeVector[1];
                C(INDEX_2D_THERMAL_FLUX_Y, INDEX_2D_THERMAL_FLUX_Y) += c2 * (c3 * dischargeVector[1] * dischargeVector[1] + c4);
            }

        }

        KRATOS_CATCH("")
    }

} // Namespace Kratos
