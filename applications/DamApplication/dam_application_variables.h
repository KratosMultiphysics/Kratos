//
//   Project Name:        KratosDamApplication $
//   Last Modified by:    $Author:     LGracia $
//   Date:                $Date:    March 2016 $
//   Revision:            $Revision:       1.0 $
//

#if !defined(KRATOS_DAM_APPLICATION_VARIABLES_H_INCLUDED )
#define  KRATOS_DAM_APPLICATION_VARIABLES_H_INCLUDED

// External includes
#include "includes/define.h"
#include "includes/kratos_application.h"
#include "includes/variables.h"
#include "solid_mechanics_application_variables.h"
#include "poromechanics_application_variables.h"   


namespace Kratos
{
//Define Variables

//Bofang, Hidrostatic and uplift variables for evolution changes
KRATOS_DEFINE_VARIABLE( std::string, GRAVITY_DIRECTION )
KRATOS_DEFINE_VARIABLE( double, COORDINATE_BASE_DAM )
KRATOS_DEFINE_VARIABLE( double, SURFACE_TEMP )
KRATOS_DEFINE_VARIABLE( double, BOTTOM_TEMP )
KRATOS_DEFINE_VARIABLE( double, HEIGHT_DAM )
KRATOS_DEFINE_VARIABLE( double, AMPLITUDE )
KRATOS_DEFINE_VARIABLE( double, DAY_MAXIMUM )
KRATOS_DEFINE_VARIABLE( double, SPECIFIC_WEIGHT )
KRATOS_DEFINE_VARIABLE( std::string, UPLIFT_DIRECTION )
KRATOS_DEFINE_VARIABLE( double, COORDINATE_BASE_DAM_UPLIFT )
KRATOS_DEFINE_VARIABLE( double, BASE_OF_DAM )

// Thermal Variables
KRATOS_DEFINE_VARIABLE( Matrix, THERMAL_STRESS_TENSOR )
KRATOS_DEFINE_VARIABLE( Matrix, MECHANICAL_STRESS_TENSOR )
KRATOS_DEFINE_VARIABLE( Matrix, THERMAL_STRAIN_TENSOR )

KRATOS_DEFINE_VARIABLE( Vector, THERMAL_STRESS_VECTOR )
KRATOS_DEFINE_VARIABLE( Vector, MECHANICAL_STRESS_VECTOR )
KRATOS_DEFINE_VARIABLE( Vector, THERMAL_STRAIN_VECTOR )

}  // namespace Kratos.

#endif /* KRATOS_DAM_APPLICATION_VARIABLES_H_INCLUDED  */
