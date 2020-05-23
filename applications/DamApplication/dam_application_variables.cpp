//
//   Project Name:        KratosDamApplication   $
//   Last Modified by:    $Author:     LGracia   $
//   Date:                $Date:      March 2016 $
//   Revision:            $Revision:         1.0 $
//

#include "dam_application_variables.h"

namespace Kratos
{
    typedef array_1d<double,3> Vector3;

    //Create Variables //Note that the application variables must not be defined if they already exist in KRATOS
    KRATOS_CREATE_VARIABLE( double, TIME_UNIT_CONVERTER )


    KRATOS_CREATE_VARIABLE( double, THERMAL_EXPANSION )

    // Thermal Variables
    KRATOS_CREATE_VARIABLE( Matrix, THERMAL_STRESS_TENSOR )
    KRATOS_CREATE_VARIABLE( Matrix, MECHANICAL_STRESS_TENSOR )
    KRATOS_CREATE_VARIABLE( Matrix, THERMAL_STRAIN_TENSOR )

    KRATOS_CREATE_VARIABLE( Vector, THERMAL_STRESS_VECTOR )
    KRATOS_CREATE_VARIABLE( Vector, MECHANICAL_STRESS_VECTOR )
    KRATOS_CREATE_VARIABLE( Vector, THERMAL_STRAIN_VECTOR )

    KRATOS_CREATE_VARIABLE( double, ALPHA_HEAT_SOURCE )
    KRATOS_CREATE_VARIABLE( double, TIME_ACTIVATION )

    // Output Variables
    KRATOS_CREATE_VARIABLE( Vector3, Vi_POSITIVE )
    KRATOS_CREATE_VARIABLE( Vector3, Viii_POSITIVE )

    // Wave Equation
    KRATOS_CREATE_VARIABLE( double, Dt_PRESSURE )
    KRATOS_CREATE_VARIABLE( double, Dt2_PRESSURE )
    KRATOS_CREATE_VARIABLE( double, VELOCITY_PRESSURE_COEFFICIENT )
    KRATOS_CREATE_VARIABLE( double, ACCELERATION_PRESSURE_COEFFICIENT )

    // Others
    KRATOS_CREATE_VARIABLE( double, NODAL_YOUNG_MODULUS )
    KRATOS_CREATE_VARIABLE( double, ADDED_MASS )
    KRATOS_CREATE_VARIABLE( double, NODAL_REFERENCE_TEMPERATURE )
    KRATOS_CREATE_VARIABLE( Matrix, NODAL_CAUCHY_STRESS_TENSOR )
    KRATOS_CREATE_VARIABLE( Matrix, INITIAL_NODAL_CAUCHY_STRESS_TENSOR )
    KRATOS_CREATE_VARIABLE( double, PLACEMENT_TEMPERATURE )

    // From Solid
    KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( FORCE_LOAD )
    KRATOS_CREATE_VARIABLE(bool, COMPUTE_CONSISTENT_MASS_MATRIX)

}// namespace Kratos.
