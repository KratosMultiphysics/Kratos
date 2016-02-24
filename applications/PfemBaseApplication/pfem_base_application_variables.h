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
#include "includes/kratos_application.h"
#include "includes/variables.h"


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

  //nodal dofs
  KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS( PFEM_BASE_APPLICATION, OFFSET )
  KRATOS_DEFINE_APPLICATION_VARIABLE( PFEM_BASE_APPLICATION, double, SHRINK_FACTOR )
  //KRATOS_DEFINE_APPLICATION_VARIABLE( PFEM_BASE_APPLICATION, Vector, NORMAL )
  //KRATOS_DEFINE_APPLICATION_VARIABLE( PFEM_BASE_APPLICATION, double, NODAL_H )


  //domain definition
  KRATOS_DEFINE_APPLICATION_VARIABLE( PFEM_BASE_APPLICATION, unsigned int, DOMAIN_LABEL )
  //KRATOS_DEFINE_APPLICATION_VARIABLE( PFEM_BASE_APPLICATION, int         , RIGID_WALL )
  //KRATOS_DEFINE_APPLICATION_VARIABLE( PFEM_BASE_APPLICATION, double      , WALL_TIP_RADIUS )
  //KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS( PFEM_BASE_APPLICATION, WALL_REFERENCE_POINT )
  //KRATOS_DEFINE_3D_APPLICATION_VARIABLE_WITH_COMPONENTS( PFEM_BASE_APPLICATION, WALL_VELOCITY )

  //boundary definition
  KRATOS_DEFINE_APPLICATION_VARIABLE( PFEM_BASE_APPLICATION, int,                               RIGID_WALL )
  KRATOS_DEFINE_APPLICATION_VARIABLE( PFEM_BASE_APPLICATION, Condition::Pointer,          MASTER_CONDITION )
  KRATOS_DEFINE_APPLICATION_VARIABLE( PFEM_BASE_APPLICATION, WeakPointerVector< Element >, MASTER_ELEMENTS )
  KRATOS_DEFINE_APPLICATION_VARIABLE( PFEM_BASE_APPLICATION, WeakPointerVector< Node<3> >,    MASTER_NODES )


  //modeler criteria
  KRATOS_DEFINE_APPLICATION_VARIABLE( PFEM_BASE_APPLICATION, double, MEAN_ERROR )


  ///@}

}

#endif	/* KRATOS_PFEM_BASE_APPLICATION_VARIABLES_H_INCLUDED */
