//-------------------------------------------------------------
//         ___  __           ___          
//  KRATOS| _ \/ _|___ _ __ | _ ) __ _ ___ ___
//        |  _/  _/ -_) '  \| _ \/ _` (_-</ -_)
//        |_| |_| \___|_|_|_|___/\__,_/__/\___|APPLICATION
//                                                                
//  License:(BSD)    PfemFluidMechanicsApplication/license.txt
//
//  Main authors:    Josep Maria Carbonell
//                   ..                  
//-------------------------------------------------------------
//
//   Project Name:        KratosPfemBaseApplication $
//   Created by:          $Author:      JMCarbonell $
//   Last modified by:    $Co-Author:               $
//   Date:                $Date:      February 2016 $
//   Revision:            $Revision:            0.0 $
//
//

#if !defined(KRATOS_PFEM_BASE_APPLICATION_VARIABLES_H_INCLUDED )
#define  KRATOS_PFEM_BASE_APPLICATION_VARIABLES_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/variables.h"
#include "includes/kratos_flags.h"
#include "includes/kratos_application.h"


namespace Kratos
{
  ///@name Type Definitions
  ///@{
  typedef array_1d<double,3> Vector3;
  typedef array_1d<double,6> Vector6;
  typedef PointerVectorSet<Condition, IndexedObject> ConditionContainerType;
  ///@}

  ///@name Kratos Globals
  ///@{

  //Define Variables

  //nodal dofs
  KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS( PFEM_BASE_APPLICATION, OFFSET )
  KRATOS_DEFINE_APPLICATION_VARIABLE( PFEM_BASE_APPLICATION, double, SHRINK_FACTOR )



  //domain definition
  KRATOS_DEFINE_APPLICATION_VARIABLE( PFEM_BASE_APPLICATION, bool, INITIALIZED_DOMAINS )
  KRATOS_DEFINE_APPLICATION_VARIABLE( PFEM_BASE_APPLICATION, std::string, MODEL_PART_NAME )
  KRATOS_DEFINE_APPLICATION_VARIABLE( PFEM_BASE_APPLICATION, unsigned int, DOMAIN_LABEL )


  //boundary definition
  KRATOS_DEFINE_APPLICATION_VARIABLE( PFEM_BASE_APPLICATION, int,                               RIGID_WALL )
  KRATOS_DEFINE_APPLICATION_VARIABLE( PFEM_BASE_APPLICATION, Condition::Pointer,          MASTER_CONDITION )
  KRATOS_DEFINE_APPLICATION_VARIABLE( PFEM_BASE_APPLICATION, WeakPointerVector< Element >, MASTER_ELEMENTS )
  KRATOS_DEFINE_APPLICATION_VARIABLE( PFEM_BASE_APPLICATION, WeakPointerVector< Node<3> >,    MASTER_NODES )

  //condition variables
  KRATOS_DEFINE_APPLICATION_VARIABLE( PFEM_BASE_APPLICATION, ConditionContainerType, CHILDREN_CONDITIONS)
    
  //modeler criteria
  KRATOS_DEFINE_APPLICATION_VARIABLE( PFEM_BASE_APPLICATION, double, MEAN_ERROR )


  ///@}

}

#endif	/* KRATOS_PFEM_BASE_APPLICATION_VARIABLES_H_INCLUDED */
