//
//   Project Name:        KratosPfemSolidMechanicsApplication $
//   Created by:          $Author:                JMCarbonell $
//   Last modified by:    $Co-Author:                         $
//   Date:                $Date:                February 2016 $
//   Revision:            $Revision:                      0.0 $
//
//

#include "pfem_solid_mechanics_application_variables.h"

namespace Kratos
{
  ///@name Type Definitions
  ///@{
  typedef array_1d<double,3> Vector3;
  typedef array_1d<double,6> Vector6;
  ///@}

  ///@name Kratos Globals
  ///@{

  //Create Variables

  //solution
  KRATOS_CREATE_VARIABLE(double, IMPOSED_WATER_PRESSURE )

  //material
  KRATOS_CREATE_VARIABLE(double, PRE_CONSOLIDATION_STRESS )
  KRATOS_CREATE_VARIABLE(double, OVER_CONSOLIDATION_RATIO )
  KRATOS_CREATE_VARIABLE(double, INITIAL_SHEAR_MODULUS )
  KRATOS_CREATE_VARIABLE(double, WATER_BULK_MODULUS )
  KRATOS_CREATE_VARIABLE(double, PERMEABILITY )
  KRATOS_CREATE_VARIABLE(double, NORMAL_COMPRESSION_SLOPE )
  KRATOS_CREATE_VARIABLE(double, SWELLING_SLOPE )
  KRATOS_CREATE_VARIABLE(double, CRITICAL_STATE_LINE )
  KRATOS_CREATE_VARIABLE(double, ALPHA_SHEAR )
  KRATOS_CREATE_VARIABLE(double, INITIAL_POROSITY )
  KRATOS_CREATE_VARIABLE(double, COHESION )
  KRATOS_CREATE_VARIABLE(double, INTERNAL_DILATANCY_ANGLE )

  //element
  KRATOS_CREATE_VARIABLE(Vector, DARCY_FLOW )
  KRATOS_CREATE_VARIABLE(Matrix, TOTAL_CAUCHY_STRESS )

  // transfer variables and initial
  KRATOS_CREATE_VARIABLE(Matrix, ELASTIC_LEFT_CAUCHY_GREEN_TENSOR )
  KRATOS_CREATE_VARIABLE(Vector, ELASTIC_LEFT_CAUCHY_GREEN_VECTOR )

  KRATOS_CREATE_VARIABLE(Matrix, KIRCHHOFF_STRESS_TENSOR )
  KRATOS_CREATE_VARIABLE(Vector, ELASTIC_LEFT_CAUCHY_FROM_KIRCHHOFF_STRESS )

  //thermal

  //mechanical

  //geometrical
  KRATOS_CREATE_VARIABLE( double, MEAN_RADIUS )

  //domain definition
  KRATOS_CREATE_VARIABLE(double, WALL_TIP_RADIUS )
  KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( WALL_REFERENCE_POINT )
  KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( WALL_VELOCITY )

  // some post process variables + stress invariants
  KRATOS_CREATE_VARIABLE(double, PRECONSOLIDATION )
  KRATOS_CREATE_VARIABLE(double, STRESS_INV_P )
  KRATOS_CREATE_VARIABLE(double, STRESS_INV_J2 )
  KRATOS_CREATE_VARIABLE(double, STRESS_INV_THETA )
  KRATOS_CREATE_VARIABLE(double, VOLUMETRIC_PLASTIC )
  KRATOS_CREATE_VARIABLE(double, INCR_SHEAR_PLASTIC )

  KRATOS_CREATE_VARIABLE(double, M_MODULUS )

  //deprecated
  KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( IMPOSED_DISPLACEMENT )
  KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( IMPOSED_ROTATION )

  ///@}

}
