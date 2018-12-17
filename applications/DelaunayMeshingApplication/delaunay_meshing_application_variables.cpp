//
//   Project Name:        KratosDelaunayMeshingApplication $
//   Created by:          $Author:             JMCarbonell $
//   Last modified by:    $Co-Author:                      $
//   Date:                $Date:                April 2018 $
//   Revision:            $Revision:                   0.0 $
//
//

#include "delaunay_meshing_application_variables.h"

namespace Kratos
{
  ///@name Type Definitions
  ///@{
  typedef array_1d<double,3> Vector3;
  typedef array_1d<double,6> Vector6;

  typedef Node<3>* NodePointerType;
  typedef Element* ElementPointerType;
  typedef Condition* ConditionPointerType;

  typedef std::vector<NodePointerType > NodePointerVectorType;
  typedef std::vector<ElementPointerType > ElementPointerVectorType;
  typedef std::vector<ConditionPointerType > ConditionPointerVectorType;

  ///@}

  ///@name Kratos Globals
  ///@{

  //Create Variables

  //geometrical definition
  KRATOS_CREATE_3D_VARIABLE_WITH_COMPONENTS( OFFSET )
  KRATOS_CREATE_VARIABLE(double, SHRINK_FACTOR )

  //domain definition
  KRATOS_CREATE_VARIABLE(bool, INITIALIZED_DOMAINS )
  KRATOS_CREATE_VARIABLE(double, MESHING_STEP_TIME )
  KRATOS_CREATE_VARIABLE(std::string, MODEL_PART_NAME )

  //boundary definition
  KRATOS_CREATE_VARIABLE(int,                                 RIGID_WALL )

  //custom neighbor and masters
  KRATOS_CREATE_VARIABLE(NodePointerType,                    MASTER_NODE )
  KRATOS_CREATE_VARIABLE(ElementPointerType,              MASTER_ELEMENT )
  KRATOS_CREATE_VARIABLE(ConditionPointerType,          MASTER_CONDITION )

  KRATOS_CREATE_VARIABLE(NodePointerVectorType,             MASTER_NODES )
  KRATOS_CREATE_VARIABLE(ElementPointerVectorType,       MASTER_ELEMENTS )
  KRATOS_CREATE_VARIABLE(ConditionPointerVectorType,   MASTER_CONDITIONS )

  KRATOS_CREATE_VARIABLE(NodePointerVectorType,           NEIGHBOR_NODES )
  KRATOS_CREATE_VARIABLE(ElementPointerVectorType,     NEIGHBOR_ELEMENTS )
  KRATOS_CREATE_VARIABLE(ConditionPointerVectorType, NEIGHBOR_CONDITIONS )

  KRATOS_CREATE_VARIABLE(std::vector<std::string>,      MODEL_PART_NAMES )

  //condition variables
  KRATOS_CREATE_VARIABLE(ConditionContainerType,     CHILDREN_CONDITIONS )

  //mesher criteria
  KRATOS_CREATE_VARIABLE(double, MEAN_ERROR )

  ///@}

}
