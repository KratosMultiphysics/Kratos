//----------------------------------------------------------
//         ___      _                                      .
//  KRATOS|   \ ___| |__ _ _  _ _ _  __ _ _  _             .
//        | |) / -_| / _` | || | ' \/ _` | || |            .
//        |___/\___|_\__,_|\_,_|_||_\__,_|\_, |MESHING     .
//                                        |__/             .
//                                                         .
//  License:(BSD)   DelaunayMeshingApplication/license.txt .
//  Main authors:   Josep Maria Carbonell                  .
//                        ..                               .
//----------------------------------------------------------
//
//   Project Name:        KratosDelaunayMeshingApplication $
//   Created by:          $Author:             JMCarbonell $
//   Last modified by:    $Co-Author:                      $
//   Date:                $Date:                April 2018 $
//   Revision:            $Revision:                   0.0 $
//
//

#if !defined(KRATOS_DELAUNAY_MESHING_APPLICATION_VARIABLES_H_INCLUDED )
#define  KRATOS_DELAUNAY_MESHING_APPLICATION_VARIABLES_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/variables.h"
#include "includes/kratos_flags.h"
#include "containers/pointer_vector_set.h"
#include "includes/element.h"
#include "includes/condition.h"
#include "utilities/indexed_object.h"
#include "includes/global_pointer_variables.h"

namespace Kratos
{
  ///@name Type Definitions
  ///@{
  typedef PointerVectorSet<Condition, IndexedObject> ConditionContainerType;

  typedef Node<3>::WeakPointer NodeWeakPtrType;
  typedef Element::WeakPointer ElementWeakPtrType;
  typedef Condition::WeakPointer ConditionWeakPtrType;

  typedef GlobalPointersVector<Node<3> > NodeWeakPtrVectorType;
  typedef GlobalPointersVector<Element> ElementWeakPtrVectorType;
  typedef GlobalPointersVector<Condition> ConditionWeakPtrVectorType;
  ///@}

  ///@name Kratos Globals
  ///@{

  //Define Variables


  //nodal dofs
  KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS( DELAUNAY_MESHING_APPLICATION, OFFSET )
  KRATOS_DEFINE_APPLICATION_VARIABLE( DELAUNAY_MESHING_APPLICATION, double, SHRINK_FACTOR )

  //domain definition
  KRATOS_DEFINE_APPLICATION_VARIABLE( DELAUNAY_MESHING_APPLICATION, bool, INITIALIZED_DOMAINS )
  KRATOS_DEFINE_APPLICATION_VARIABLE( DELAUNAY_MESHING_APPLICATION, double, MESHING_STEP_TIME )
  KRATOS_DEFINE_APPLICATION_VARIABLE( DELAUNAY_MESHING_APPLICATION, std::string, MODEL_PART_NAME )
  KRATOS_DEFINE_APPLICATION_VARIABLE( DELAUNAY_MESHING_APPLICATION, std::vector<std::string>, MODEL_PART_NAMES )

  //boundary definition
  KRATOS_DEFINE_APPLICATION_VARIABLE( DELAUNAY_MESHING_APPLICATION, int,                                 RIGID_WALL )
  KRATOS_DEFINE_APPLICATION_VARIABLE( DELAUNAY_MESHING_APPLICATION, int, PROPERTY_ID )


  //custom neighbor and masters
  KRATOS_DEFINE_APPLICATION_VARIABLE( DELAUNAY_MESHING_APPLICATION, NodeWeakPtrType,                    MASTER_NODE )
  KRATOS_DEFINE_APPLICATION_VARIABLE( DELAUNAY_MESHING_APPLICATION, ElementWeakPtrType,              MASTER_ELEMENT )
  KRATOS_DEFINE_APPLICATION_VARIABLE( DELAUNAY_MESHING_APPLICATION, ConditionWeakPtrType,          MASTER_CONDITION )

  KRATOS_DEFINE_APPLICATION_VARIABLE( DELAUNAY_MESHING_APPLICATION, NodeWeakPtrVectorType,             MASTER_NODES )
  KRATOS_DEFINE_APPLICATION_VARIABLE( DELAUNAY_MESHING_APPLICATION, ElementWeakPtrVectorType,       MASTER_ELEMENTS )
  KRATOS_DEFINE_APPLICATION_VARIABLE( DELAUNAY_MESHING_APPLICATION, ConditionWeakPtrVectorType,   MASTER_CONDITIONS )

  //condition variables
  KRATOS_DEFINE_APPLICATION_VARIABLE( DELAUNAY_MESHING_APPLICATION, ConditionContainerType,      CHILDREN_CONDITIONS)

  //mesher criteria
  KRATOS_DEFINE_APPLICATION_VARIABLE( DELAUNAY_MESHING_APPLICATION, double, MEAN_ERROR )


  ///@}

}

#endif	/* KRATOS_DELAUNAY_MESHING_APPLICATION_VARIABLES_H_INCLUDED */
