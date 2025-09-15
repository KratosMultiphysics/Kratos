//------------------------------------------------------------------
//           _   _            _                                    .
//   KRATOS | | | |_ __  __ _| |_                                  .
//          | |_| | '  \/ _` |  _|                                 .
//           \___/|_|_|_\__,_|\__| INTERFACE                       .                             
//			                                           .
//   License:(BSD)	  UmatApplication/license.txt              .
//   Main authors:        LlMonforte, JMCarbonell                  .
//                        ..                                       .
//------------------------------------------------------------------
//
//   Project Name:        KratosUmatApplication      $
//   Created by:          $Author:      JMCarbonell  $
//   Last modified by:    $Co-Author:                $
//   Date:                $Date:      September 2017 $
//   Revision:            $Revision:             0.0 $
//
//


#if !defined(KRATOS_UMAT_APPLICATION_VARIABLES_H_INCLUDED )
#define  KRATOS_UMAT_APPLICATION_VARIABLES_H_INCLUDED

// System includes

// External includes

// Project includes
#include "constitutive_models_application_variables.h"

namespace Kratos
{
   KRATOS_DEFINE_APPLICATION_VARIABLE( UMAT_APPLICATION, double, ALPHA )
   KRATOS_DEFINE_APPLICATION_VARIABLE( UMAT_APPLICATION, double, BETA )
   KRATOS_DEFINE_APPLICATION_VARIABLE( UMAT_APPLICATION, double, MF )
   KRATOS_DEFINE_APPLICATION_VARIABLE( UMAT_APPLICATION, double, CC )
   KRATOS_DEFINE_APPLICATION_VARIABLE( UMAT_APPLICATION, double, MM )
   KRATOS_DEFINE_APPLICATION_VARIABLE( UMAT_APPLICATION, double, KSIS )
   KRATOS_DEFINE_APPLICATION_VARIABLE( UMAT_APPLICATION, double, RHOM )
   KRATOS_DEFINE_APPLICATION_VARIABLE( UMAT_APPLICATION, double, PC0 )
   KRATOS_DEFINE_APPLICATION_VARIABLE( UMAT_APPLICATION, double, VOID_RATIO )
   KRATOS_DEFINE_APPLICATION_VARIABLE( UMAT_APPLICATION, double, PLASTIC_MULTIPLIER )

}

#endif	/* KRATOS_UMAT_APPLICATION_VARIABLES_H_INCLUDED */
