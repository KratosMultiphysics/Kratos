//
//   Project Name:        KratosPfemBaseApplication $
//   Created by:          $Author:      JMCarbonell $
//   Last modified by:    $Co-Author:               $
//   Date:                $Date:      February 2016 $
//   Revision:            $Revision:            0.0 $
//
//

#include "pfem_base_application_variables.h"

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

  //geometrical definition
  KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( OFFSET )
  KRATOS_CREATE_VARIABLE(double, SHRINK_FACTOR )

  //domain definition
  KRATOS_CREATE_VARIABLE(bool, INITIALIZED_DOMAINS )
  KRATOS_CREATE_VARIABLE(std::string, MODEL_PART_NAME )
  KRATOS_CREATE_VARIABLE(unsigned int, DOMAIN_LABEL )

  //boundary definition
  KRATOS_CREATE_VARIABLE(int,                               RIGID_WALL )
  KRATOS_CREATE_VARIABLE(Condition::Pointer,          MASTER_CONDITION )
  KRATOS_CREATE_VARIABLE(WeakPointerVector< Element >, MASTER_ELEMENTS )
  KRATOS_CREATE_VARIABLE(WeakPointerVector< Node<3> >,    MASTER_NODES )

  //condition variables
  KRATOS_CREATE_VARIABLE(ConditionContainerType,   CHILDREN_CONDITIONS )
  
  //modeler criteria
  KRATOS_CREATE_VARIABLE(double, MEAN_ERROR )

  ///@}

}
