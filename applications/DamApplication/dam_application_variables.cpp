//
//   Project Name:        KratosDamApplication   $   
//   Last Modified by:    $Author:     LGracia   $
//   Date:                $Date:      March 2016 $
//   Revision:            $Revision:         1.0 $
//

#include "dam_application_variables.h"

namespace Kratos
{
//Create Variables //Note that the application variables must not be defined if they already exist in KRATOS

//Bofang, Hidrostatic and uplift variables for evolution changes
KRATOS_CREATE_VARIABLE( std::string, GRAVITY_DIRECTION )
KRATOS_CREATE_VARIABLE( double, COORDINATE_BASE_DAM )
KRATOS_CREATE_VARIABLE( double, SURFACE_TEMP )
KRATOS_CREATE_VARIABLE( double, BOTTOM_TEMP )
KRATOS_CREATE_VARIABLE( double, HEIGHT_DAM )
KRATOS_CREATE_VARIABLE( double, AMPLITUDE )
KRATOS_CREATE_VARIABLE( double, DAY_MAXIMUM )
KRATOS_CREATE_VARIABLE( double, SPECIFIC_WEIGHT )
KRATOS_CREATE_VARIABLE( std::string, UPLIFT_DIRECTION )
KRATOS_CREATE_VARIABLE( double, COORDINATE_BASE_DAM_UPLIFT )
KRATOS_CREATE_VARIABLE( double, BASE_OF_DAM )

// Thermal Variables
KRATOS_CREATE_VARIABLE( Matrix, THERMAL_STRESS_TENSOR )
KRATOS_CREATE_VARIABLE( Matrix, MECHANICAL_STRESS_TENSOR )
KRATOS_CREATE_VARIABLE( Matrix, THERMAL_STRAIN_TENSOR )

KRATOS_CREATE_VARIABLE( Vector, THERMAL_STRESS_VECTOR )
KRATOS_CREATE_VARIABLE( Vector, MECHANICAL_STRESS_VECTOR )
KRATOS_CREATE_VARIABLE( Vector, THERMAL_STRAIN_VECTOR )

}// namespace Kratos.
