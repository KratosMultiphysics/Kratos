//--------------------------------------------------------
//          ___  __                                      .
//  KRATOS | _ \/ _|___ _ __                             .
//         |  _/  _/ -_) '  \                            .
//         |_| |_| \___|_|_|_| APPLICATION               .
//                                                       .
//  License:(BSD)         PfemApplication/license.txt    .
//  Main authors:         Josep Maria Carbonell          .
//                        ..                             .
//--------------------------------------------------------
//
//   Project Name:        KratosPfemApplication     $
//   Created by:          $Author:      JMCarbonell $
//   Last modified by:    $Co-Author:               $
//   Date:                $Date:           May 2018 $
//   Revision:            $Revision:            0.0 $
//
//

#if !defined(KRATOS_PFEM_APPLICATION_VARIABLES_H_INCLUDED )
#define  KRATOS_PFEM_APPLICATION_VARIABLES_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/cfd_variables.h"
#include "includes/kratos_flags.h"

#include "solid_mechanics_application_variables.h"

namespace Kratos
{
  ///@name Type Definitions
  ///@{
  ///@}

  typedef PointerVectorSet<Properties, IndexedObject> PropertiesContainerType;
  typedef typename PropertiesContainerType::Pointer   PropertiesContainerPointerType;

  ///@name Kratos Globals
  ///@{

  //Define Variables
  KRATOS_DEFINE_APPLICATION_VARIABLE( PFEM_APPLICATION, double, FLUID_PRESSURE )
  KRATOS_DEFINE_APPLICATION_VARIABLE( PFEM_APPLICATION, double, FLUID_PRESSURE_VELOCITY )
  KRATOS_DEFINE_APPLICATION_VARIABLE( PFEM_APPLICATION, double, FLUID_PRESSURE_ACCELERATION )
  KRATOS_DEFINE_APPLICATION_VARIABLE( PFEM_APPLICATION, double, FLUID_PRESSURE_REACTION )
  KRATOS_DEFINE_APPLICATION_VARIABLE( PFEM_APPLICATION, double, VOLUME_WEAR )

  KRATOS_DEFINE_APPLICATION_VARIABLE( PFEM_APPLICATION, PropertiesContainerPointerType, PROPERTIES_VECTOR )
  KRATOS_DEFINE_APPLICATION_VARIABLE( PFEM_APPLICATION, Vector, MATERIAL_PERCENTAGE )

  //Adaptive time step (review needed)
  KRATOS_DEFINE_APPLICATION_VARIABLE( PFEM_APPLICATION, double, INITIAL_DELTA_TIME )
  KRATOS_DEFINE_APPLICATION_VARIABLE( PFEM_APPLICATION, double, CURRENT_DELTA_TIME )
  KRATOS_DEFINE_APPLICATION_VARIABLE( PFEM_APPLICATION,   bool, TIME_INTERVAL_CHANGED )
  KRATOS_DEFINE_APPLICATION_VARIABLE( PFEM_APPLICATION,   bool, BAD_VELOCITY_CONVERGENCE )
  KRATOS_DEFINE_APPLICATION_VARIABLE( PFEM_APPLICATION,   bool, BAD_PRESSURE_CONVERGENCE )

  //Material variables
  KRATOS_DEFINE_APPLICATION_VARIABLE( PFEM_APPLICATION, double, WEAR_COEFFICIENT )
  KRATOS_DEFINE_APPLICATION_VARIABLE( PFEM_APPLICATION, double, INDENTATION_HARDNESS )



  ///@}

}

#endif	// KRATOS_PFEM_APPLICATION_VARIABLES_H_INCLUDED defined
