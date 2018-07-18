//
//   Project Name:        KratosPfemApplication     $
//   Created by:          $Author:      JMCarbonell $
//   Last modified by:    $Co-Author:               $
//   Date:                $Date:           May 2018 $
//   Revision:            $Revision:            0.0 $
//
//

#include "pfem_application_variables.h"

namespace Kratos
{
  ///@name Type Definitions
  ///@{
  ///@}

  typedef PointerVectorSet<Properties, IndexedObject> PropertiesContainerType;
  typedef typename PropertiesContainerType::Pointer   PropertiesContainerPointerType;


  ///@name Kratos Globals
  ///@{


  //Create Variables
  KRATOS_CREATE_VARIABLE( double, FLUID_PRESSURE )
  KRATOS_CREATE_VARIABLE( double, FLUID_PRESSURE_VELOCITY )
  KRATOS_CREATE_VARIABLE( double, FLUID_PRESSURE_ACCELERATION )
  KRATOS_CREATE_VARIABLE( double, FLUID_PRESSURE_REACTION )

  KRATOS_CREATE_VARIABLE( PropertiesContainerPointerType, PROPERTIES_VECTOR )
  KRATOS_CREATE_VARIABLE( Vector, MATERIAL_PERCENTAGE )

  //Adaptive time step (review needed)
  KRATOS_CREATE_VARIABLE(   bool, TIME_INTERVAL_CHANGED )
  KRATOS_CREATE_VARIABLE(   bool, BAD_VELOCITY_CONVERGENCE )
  KRATOS_CREATE_VARIABLE(   bool, BAD_PRESSURE_CONVERGENCE )
  KRATOS_CREATE_VARIABLE( double, INITIAL_DELTA_TIME )
  KRATOS_CREATE_VARIABLE( double, CURRENT_DELTA_TIME )


  ///@}

}
