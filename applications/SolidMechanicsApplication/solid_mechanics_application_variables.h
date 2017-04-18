//--------------------------------------------------------------------
//    |  /           |                                               .
//    ' /   __| _` | __|  _ \   __|                                  .
//    . \  |   (   | |   (   |\__ \                                  .
//   _|\_\_|  \__,_|\__|\___/ ____/                                  .
//                 KRATOS  __|   _ \  |   |  _ \                     .
//                       \__ \  (   | |   | | , )                    .      
//                       |___/ \___/ ___|_| ___/ MECHANICS           .            
//			                                             .
//   License:(BSD)	  SolidMechanicsApplication/license.txt      .
//   Main authors:        Josep Maria Carbonell                      .
//                        ..                                         .
//--------------------------------------------------------------------
//
//   Project Name:        KratosSolidMechanicsApplication $
//   Created by:          $Author:            JMCarbonell $
//   Last modified by:    $Co-Author:                     $
//   Date:                $Date:            February 2016 $
//   Revision:            $Revision:                  0.0 $
//
//

#if !defined(KRATOS_SOLID_MECHANICS_APPLICATION_VARIABLES_H_INCLUDED )
#define  KRATOS_SOLID_MECHANICS_APPLICATION_VARIABLES_H_INCLUDED

// System includes

// External includes


// Project includes
#include "includes/define.h"
#include "includes/kratos_application.h"
#include "includes/variables.h"

#include "constitutive_models_application_variables.h"

namespace Kratos
{
  ///@name Type Definitions
  ///@{
  typedef array_1d<double,3> Vector3;
  typedef array_1d<double,6> Vector6;
  ///@}

  ///@name Kratos Globals
  ///@{

  //Define Variables

  //for explicit schemes
  KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS( SOLID_MECHANICS_APPLICATION, MIDDLE_VELOCITY )
  //KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS( SOLID_MECHANICS_APPLICATION, DELTA_DISPLACEMENT )

  //solution
  KRATOS_DEFINE_APPLICATION_VARIABLE( SOLID_MECHANICS_APPLICATION, int, WRITE_ID )
  KRATOS_DEFINE_APPLICATION_VARIABLE( SOLID_MECHANICS_APPLICATION, double, PREVIOUS_DELTA_TIME )
  KRATOS_DEFINE_APPLICATION_VARIABLE( SOLID_MECHANICS_APPLICATION, double, RAYLEIGH_ALPHA )
  KRATOS_DEFINE_APPLICATION_VARIABLE( SOLID_MECHANICS_APPLICATION, double, RAYLEIGH_BETA )

  //nodal load variables
  KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS( SOLID_MECHANICS_APPLICATION, POINT_LOAD )
  KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS( SOLID_MECHANICS_APPLICATION, LINE_LOAD )
  KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS( SOLID_MECHANICS_APPLICATION, SURFACE_LOAD )

  KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS( SOLID_MECHANICS_APPLICATION, FOLLOWER_POINT_LOAD )
  KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS( SOLID_MECHANICS_APPLICATION, FOLLOWER_LINE_LOAD )
  KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS( SOLID_MECHANICS_APPLICATION, FOLLOWER_SURFACE_LOAD )

  KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS( SOLID_MECHANICS_APPLICATION, POINT_MOMENT )
  KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS( SOLID_MECHANICS_APPLICATION, FOLLOWER_POINT_MOMENT )
    
  //condition load variables
  KRATOS_DEFINE_APPLICATION_VARIABLE( SOLID_MECHANICS_APPLICATION, Vector, POINT_LOADS_VECTOR )
  KRATOS_DEFINE_APPLICATION_VARIABLE( SOLID_MECHANICS_APPLICATION, Vector, LINE_LOADS_VECTOR )
  KRATOS_DEFINE_APPLICATION_VARIABLE( SOLID_MECHANICS_APPLICATION, Vector, SURFACE_LOADS_VECTOR )
  KRATOS_DEFINE_APPLICATION_VARIABLE( SOLID_MECHANICS_APPLICATION, Vector, POSITIVE_FACE_PRESSURES_VECTOR )
  KRATOS_DEFINE_APPLICATION_VARIABLE( SOLID_MECHANICS_APPLICATION, Vector, NEGATIVE_FACE_PRESSURES_VECTOR )
    
  //element
  KRATOS_DEFINE_APPLICATION_VARIABLE( SOLID_MECHANICS_APPLICATION, Matrix, ALMANSI_STRAIN_TENSOR )
  KRATOS_DEFINE_APPLICATION_VARIABLE( SOLID_MECHANICS_APPLICATION, Vector, GREEN_LAGRANGE_STRAIN_VECTOR )
  KRATOS_DEFINE_APPLICATION_VARIABLE( SOLID_MECHANICS_APPLICATION, Vector, ALMANSI_STRAIN_VECTOR )

  //KRATOS_DEFINE_APPLICATION_VARIABLE( SOLID_MECHANICS_APPLICATION, Matrix, MATERIAL_STIFFNESS_MATRIX )
  //KRATOS_DEFINE_APPLICATION_VARIABLE( SOLID_MECHANICS_APPLICATION, Matrix, GEOMETRIC_STIFFNESS_MATRIX )

  KRATOS_DEFINE_APPLICATION_VARIABLE( SOLID_MECHANICS_APPLICATION, double, VON_MISES_STRESS )

  //nodal dofs
  KRATOS_DEFINE_APPLICATION_VARIABLE( SOLID_MECHANICS_APPLICATION, double, PRESSURE_REACTION )


  ///@}

}

#endif	/* KRATOS_SOLID_MECHANICS_APPLICATION_VARIABLES_H_INCLUDED */
