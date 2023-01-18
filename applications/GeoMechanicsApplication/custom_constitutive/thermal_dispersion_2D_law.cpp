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
    GeoThermalDispersion2DLaw::GeoThermalDispersion2DLaw() : LinearPlaneStrainK0Law()
    {
    }

    // ============================================================================================
    // ============================================================================================
    GeoThermalDispersion2DLaw::GeoThermalDispersion2DLaw(const GeoThermalDispersion2DLaw& rOther)
        : LinearPlaneStrainK0Law(rOther) {}

    // ============================================================================================
    // ============================================================================================
    ConstitutiveLaw::Pointer GeoThermalDispersion2DLaw::Clone() const
    {
        return Kratos::make_shared<GeoThermalDispersion2DLaw>(*this);
    }

    // ============================================================================================
    // ============================================================================================
    GeoThermalDispersion2DLaw::~GeoThermalDispersion2DLaw() {}

    // ============================================================================================
    // ============================================================================================
    void GeoThermalDispersion2DLaw::CalculateThermalDispersionMatrix(
        Matrix& C,
        ConstitutiveLaw::Parameters& rValues)
    {
        KRATOS_TRY

        const Properties& r_material_properties = rValues.GetMaterialProperties();
        const double porosity = r_material_properties[POROSITY];
        const double saturation = r_material_properties[SATURATION];
        const double thermalCondWater = r_material_properties[THERMAL_CONDUCTIVITY_WATER];
        const double thermalCondSolidXX = r_material_properties[THERMAL_CONDUCTIVITY_SOLID_XX];
        const double thermalCondSolidXY = r_material_properties[THERMAL_CONDUCTIVITY_SOLID_XY];
        const double thermalCondSolidYX = r_material_properties[THERMAL_CONDUCTIVITY_SOLID_YX];
        const double thermalCondSolidYY = r_material_properties[THERMAL_CONDUCTIVITY_SOLID_YY];

        const double c0 = porosity * saturation * thermalCondWater;
        const double c1 = 1.0 - porosity;
        
        C(INDEX_2D_HEAT_X, INDEX_2D_HEAT_X) = c1 * thermalCondSolidXX + c0;
        C(INDEX_2D_HEAT_X, INDEX_2D_HEAT_Y) = c1 * thermalCondSolidXY;
        C(INDEX_2D_HEAT_Y, INDEX_2D_HEAT_X) = c1 * thermalCondSolidYX;
        C(INDEX_2D_HEAT_Y, INDEX_2D_HEAT_Y) = c1 * thermalCondSolidYY + c0;
        
        KRATOS_CATCH("")
    }

} // Namespace Kratos
