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
    typedef array_1d<double,3> Vector3;
    
    //Define Variables
    KRATOS_DEFINE_VARIABLE( double, THERMAL_EXPANSION )

    // Thermal Variables
    KRATOS_DEFINE_VARIABLE( Matrix, THERMAL_STRESS_TENSOR )
    KRATOS_DEFINE_VARIABLE( Matrix, MECHANICAL_STRESS_TENSOR )
    KRATOS_DEFINE_VARIABLE( Matrix, THERMAL_STRAIN_TENSOR )

    KRATOS_DEFINE_VARIABLE( Vector, THERMAL_STRESS_VECTOR )
    KRATOS_DEFINE_VARIABLE( Vector, MECHANICAL_STRESS_VECTOR )
    KRATOS_DEFINE_VARIABLE( Vector, THERMAL_STRAIN_VECTOR )

    KRATOS_DEFINE_VARIABLE( double, ALPHA_HEAT_SOURCE )   
    KRATOS_DEFINE_VARIABLE( double, TIME_ACTIVATION )    
        
    // Output Variables
    KRATOS_DEFINE_VARIABLE( Vector3, Vi_POSITIVE )
    KRATOS_DEFINE_VARIABLE( Vector3, Viii_POSITIVE )
    KRATOS_DEFINE_VARIABLE( double, NODAL_JOINT_WIDTH )
    KRATOS_DEFINE_VARIABLE( double, NODAL_JOINT_AREA )
    
    // Wave Equation
    KRATOS_DEFINE_VARIABLE( double, Dt_PRESSURE )
    KRATOS_DEFINE_VARIABLE( double, Dt2_PRESSURE )
    KRATOS_DEFINE_VARIABLE( double, VELOCITY_PRESSURE_COEFFICIENT )
    KRATOS_DEFINE_VARIABLE( double, ACCELERATION_PRESSURE_COEFFICIENT )

    // Others
    KRATOS_DEFINE_VARIABLE( double, NODAL_YOUNG_MODULUS )
    KRATOS_DEFINE_VARIABLE( double, ADDED_MASS )    
    

}  // namespace Kratos.

#endif /* KRATOS_DAM_APPLICATION_VARIABLES_H_INCLUDED  */
