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

  typedef Kratos::weak_ptr<Node<3> > NodeWeakPtrType;
  typedef Kratos::weak_ptr<Element> ElementWeakPtrType;
  typedef Kratos::weak_ptr<Condition> ConditionWeakPtrType;

  typedef WeakPointerVector<Node<3> > NodeWeakPtrVectorType;
  typedef WeakPointerVector<Element> ElementWeakPtrVectorType;
  typedef WeakPointerVector<Condition> ConditionWeakPtrVectorType;
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
  KRATOS_CREATE_VARIABLE(std::vector<std::string>, MODEL_PART_NAMES )

  //boundary definition
  KRATOS_CREATE_VARIABLE(int,                                 RIGID_WALL )

  //custom neighbor and masters
  KRATOS_CREATE_VARIABLE(NodeWeakPtrType,                    MASTER_NODE )
  KRATOS_CREATE_VARIABLE(ElementWeakPtrType,              MASTER_ELEMENT )
  KRATOS_CREATE_VARIABLE(ConditionWeakPtrType,          MASTER_CONDITION )

  KRATOS_CREATE_VARIABLE(NodeWeakPtrVectorType,             MASTER_NODES )
  KRATOS_CREATE_VARIABLE(ElementWeakPtrVectorType,       MASTER_ELEMENTS )
  KRATOS_CREATE_VARIABLE(ConditionWeakPtrVectorType,   MASTER_CONDITIONS )

  //condition variables
  KRATOS_CREATE_VARIABLE(ConditionContainerType,     CHILDREN_CONDITIONS )

  //mesher criteria
  KRATOS_CREATE_VARIABLE(double, MEAN_ERROR )

  ///@}

}
