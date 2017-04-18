//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:            JMCarbonell $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:            February 2016 $
//   Revision:            $Revision:                  0.0 $
//
//

#include "solid_mechanics_application_variables.h"

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

  //explicit schemes
  KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( MIDDLE_VELOCITY )

  //solution
  KRATOS_CREATE_VARIABLE( int, WRITE_ID )
  KRATOS_CREATE_VARIABLE( double, PREVIOUS_DELTA_TIME )
  KRATOS_CREATE_VARIABLE( double, RAYLEIGH_ALPHA )
  KRATOS_CREATE_VARIABLE( double, RAYLEIGH_BETA )

  //nodal load variables
  KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( POINT_LOAD )
  KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( LINE_LOAD )
  KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( SURFACE_LOAD )
  
  KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( FOLLOWER_POINT_LOAD )
  KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( FOLLOWER_LINE_LOAD )
  KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( FOLLOWER_SURFACE_LOAD )
  
  KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( POINT_MOMENT )
  KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( FOLLOWER_POINT_MOMENT )
  
  //condition load variables
  KRATOS_CREATE_VARIABLE( Vector, POINT_LOADS_VECTOR )
  KRATOS_CREATE_VARIABLE( Vector, LINE_LOADS_VECTOR )
  KRATOS_CREATE_VARIABLE( Vector, SURFACE_LOADS_VECTOR )
  KRATOS_CREATE_VARIABLE( Vector, POSITIVE_FACE_PRESSURES_VECTOR )
  KRATOS_CREATE_VARIABLE( Vector, NEGATIVE_FACE_PRESSURES_VECTOR )
  

  //element
  KRATOS_CREATE_VARIABLE( Matrix, ALMANSI_STRAIN_TENSOR )
  KRATOS_CREATE_VARIABLE( Vector, GREEN_LAGRANGE_STRAIN_VECTOR )
  KRATOS_CREATE_VARIABLE( Vector, ALMANSI_STRAIN_VECTOR )


  KRATOS_CREATE_VARIABLE( double, VON_MISES_STRESS )

  //nodal dofs
  KRATOS_CREATE_VARIABLE( double, PRESSURE_REACTION )


  ///@}

}
