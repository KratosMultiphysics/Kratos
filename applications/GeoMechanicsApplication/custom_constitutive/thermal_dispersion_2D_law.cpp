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
#include "custom_constitutive/thermal_dispersion_2D_law.h"

#include "geo_mechanics_application_variables.h"

namespace Kratos
{

//******************************CONSTRUCTOR*******************************************
/***********************************************************************************/

GeoThermalDispersion2DLaw::GeoThermalDispersion2DLaw() : ConstitutiveLaw()
{
}

//******************************COPY CONSTRUCTOR**************************************
/***********************************************************************************/

GeoThermalDispersion2DLaw::
GeoThermalDispersion2DLaw(const GeoThermalDispersion2DLaw& rOther) : ConstitutiveLaw(rOther)
{
	
}

//********************************CLONE***********************************************
/***********************************************************************************/

ConstitutiveLaw::Pointer GeoThermalDispersion2DLaw::Clone() const
{
    return Kratos::make_shared <GeoThermalDispersion2DLaw> (*this);
}

//*******************************DESTRUCTOR*******************************************
/***********************************************************************************/

GeoThermalDispersion2DLaw::~GeoThermalDispersion2DLaw() {}

/***********************************************************************************/
/***********************************************************************************/

bool& GeoThermalDispersion2DLaw::GetValue(const Variable<bool>& rThisVariable, bool& rValue)
{
    // This Constitutive Law has been checked with Stenberg Stabilization
    if (rThisVariable == STENBERG_SHEAR_STABILIZATION_SUITABLE) {
        rValue = true;
    }

    return rValue;
}

/***********************************************************************************/
/***********************************************************************************/
void GeoThermalDispersion2DLaw::CalculateThermalDispersionMatrix(Matrix& C)
{
    KRATOS_TRY

    const double c0 = POROSITY * SATURATION * THERMAL_CONDUCTIVITY_WATER;
    const double c1 = 1.0 - POROSITY;

    C(0, 0) = c1 * THERMAL_CONDUCTIVITY_SOLID_XX + c0;
    C(0, 1) = c1 * THERMAL_CONDUCTIVITY_SOLID_XY;
    C(1, 0) = c1 * THERMAL_CONDUCTIVITY_SOLID_YX;
    C(1, 1) = c1 * THERMAL_CONDUCTIVITY_SOLID_YY + c0;

    KRATOS_CATCH("")
}


} // Namespace Kratos
