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

  //scheme

  //solution
  KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( WATER_DISPLACEMENT_REACTION )
  KRATOS_CREATE_VARIABLE(double, WATER_PRESSURE_VELOCITY )


  KRATOS_CREATE_VARIABLE(double, JACOBIAN )
  KRATOS_CREATE_VARIABLE(double, REACTION_JACOBIAN )

  //material
  KRATOS_CREATE_VARIABLE(double, WATER_BULK_MODULUS )
  KRATOS_CREATE_VARIABLE(double, PERMEABILITY )
  KRATOS_CREATE_VARIABLE(bool, KOZENY_CARMAN)
  KRATOS_CREATE_VARIABLE(double, INITIAL_POROSITY )
  KRATOS_CREATE_VARIABLE(double, VOID_RATIO)


  //element
  KRATOS_CREATE_VARIABLE(Matrix, TOTAL_CAUCHY_STRESS )
  KRATOS_CREATE_VARIABLE(Vector, DARCY_FLOW )

  KRATOS_CREATE_VARIABLE( double, STABILIZATION_FACTOR_J )
  KRATOS_CREATE_VARIABLE( double, STABILIZATION_FACTOR_P )
  KRATOS_CREATE_VARIABLE( double, STABILIZATION_FACTOR_WP )

  // transfer variables and initial
  KRATOS_CREATE_VARIABLE(Matrix, ELASTIC_LEFT_CAUCHY_GREEN_TENSOR )
  KRATOS_CREATE_VARIABLE(Vector, ELASTIC_LEFT_CAUCHY_GREEN_VECTOR )


  KRATOS_CREATE_VARIABLE(Matrix, INVERSE_DEFORMATION_GRADIENT )

  //thermal

  //mechanical

  //geometrical

  //domain definition
  KRATOS_CREATE_VARIABLE(double, WALL_TIP_RADIUS )
  KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( WALL_REFERENCE_POINT )
  KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( WALL_VELOCITY )

  // some post process variables + stress invariants
  KRATOS_CREATE_VARIABLE(double, PRECONSOLIDATION )
  KRATOS_CREATE_VARIABLE(double, VOLUMETRIC_PLASTIC )
  KRATOS_CREATE_VARIABLE(double, INCR_SHEAR_PLASTIC )

  KRATOS_CREATE_VARIABLE(double, M_MODULUS )

  //deprecated
  KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( IMPOSED_DISPLACEMENT )
  KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( IMPOSED_ROTATION )

  KRATOS_CREATE_VARIABLE( double, PENALTY_PARAMETER )
  ///@}

}
