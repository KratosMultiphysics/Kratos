//
//   Project Name:        KratosPfemFluidDynamicsApplication $
//   Created by:          $Author:               JMCarbonell $
//   Last modified by:    $Co-Author:                        $
//   Date:                $Date:               February 2016 $
//   Revision:            $Revision:                     0.0 $
//
//

#include "pfem_fluid_dynamics_application_variables.h"

namespace Kratos
{

  ///@name Type Definitions
  ///@{

  ///@}

  ///@name Kratos Globals
  ///@{

  //Create Variables 

  // KRATOS_CREATE_VARIABLE(double, M_MODULUS )
  // KRATOS_CREATE_VARIABLE(int, PATCH_INDEX )
  // KRATOS_CREATE_VARIABLE(double, NORMVELOCITY )
  KRATOS_CREATE_VARIABLE(bool, FREESURFACE )
  KRATOS_CREATE_VARIABLE(double, INITIAL_DELTA_TIME )
  KRATOS_CREATE_VARIABLE(double, CURRENT_DELTA_TIME )
  KRATOS_CREATE_VARIABLE(bool, TIME_INTERVAL_CHANGED )
  KRATOS_CREATE_VARIABLE(bool, BAD_VELOCITY_CONVERGENCE )
  KRATOS_CREATE_VARIABLE(bool, BAD_PRESSURE_CONVERGENCE )

  //Papanastasiou variables
  KRATOS_CREATE_VARIABLE(double, FLOW_INDEX )
  KRATOS_CREATE_VARIABLE(double, YIELDED )
  KRATOS_CREATE_VARIABLE(double, YIELD_SHEAR )
  KRATOS_CREATE_VARIABLE(double, ADAPTIVE_EXPONENT )

  //mu(I)-rheology variables
  KRATOS_CREATE_VARIABLE( double, STATIC_FRICTION )
  KRATOS_CREATE_VARIABLE( double, DYNAMIC_FRICTION )
  KRATOS_CREATE_VARIABLE( double, INERTIAL_NUMBER_ZERO )
  KRATOS_CREATE_VARIABLE( double, GRAIN_DIAMETER )
  KRATOS_CREATE_VARIABLE( double, GRAIN_DENSITY )
  KRATOS_CREATE_VARIABLE( double, REGULARIZATION_COEFFICIENT )
  KRATOS_CREATE_VARIABLE( double, INFINITE_FRICTION )
  KRATOS_CREATE_VARIABLE( double, INERTIAL_NUMBER_ONE )
  KRATOS_CREATE_VARIABLE( double, ALPHA_PARAMETER )

  KRATOS_CREATE_VARIABLE(double, PRESSURE_VELOCITY )
  KRATOS_CREATE_VARIABLE(double, PRESSURE_ACCELERATION )


  ///@}

}
