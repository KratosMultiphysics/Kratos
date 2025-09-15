//--------------------------------------------------------------------
//    |  /           |                                               .
//    ' /   __| _` | __|  _ \   __|                                  .
//    . \  |   (   | |   (   |\__ \                                  .
//   _|\_\_|  \__,_|\__|\___/ ____/                                  .
//                        __  __      _           _      _           .
//   KRATOS CONSTITUTIVE |  \/  |__ _| |_ ___ _ _(_)__ _| |          .
//                       | |\/| / _` |  _/ -_) '_| / _` | |          .
//                       |_|  |_\__,_|\__\___|_| |_\__,_|_| MODELS   .
//			                                             .
//   License:(BSD)	  ConstitutiveModelsApplication/license.txt  .
//   Main authors:        Josep Maria Carbonell                      .
//                        ..                                         .
//--------------------------------------------------------------------
//
//   Project Name:        KratosConstitutiveModelsApplication $
//   Created by:          $Author:                JMCarbonell $
//   Last modified by:    $Co-Author:                         $
//   Date:                $Date:                   April 2017 $
//   Revision:            $Revision:                      0.0 $
//
//

#if !defined(KRATOS_CONSTITUTIVE_MODELS_APPLICATION_VARIABLES_H_INCLUDED )
#define  KRATOS_CONSTITUTIVE_MODELS_APPLICATION_VARIABLES_H_INCLUDED

// System includes

// External includes

// Project includes
#include "includes/define.h"
#include "includes/variables.h"
#include "includes/mat_variables.h"
#include "includes/cfd_variables.h"
#include "includes/kratos_application.h"
#include "includes/checks.h"
#include "custom_utilities/properties_layout.hpp"

namespace Kratos
{
  //specific constitutive models variables must be DEFINED here
   KRATOS_DEFINE_APPLICATION_VARIABLE( CONSTITUTIVE_MODELS_APPLICATION, double, RHOS )
   KRATOS_DEFINE_APPLICATION_VARIABLE( CONSTITUTIVE_MODELS_APPLICATION, double, RHOT )
   KRATOS_DEFINE_APPLICATION_VARIABLE( CONSTITUTIVE_MODELS_APPLICATION, double, CHIS )
   KRATOS_DEFINE_APPLICATION_VARIABLE( CONSTITUTIVE_MODELS_APPLICATION, double, CHIT )
   KRATOS_DEFINE_APPLICATION_VARIABLE( CONSTITUTIVE_MODELS_APPLICATION, double, REFERENCE_PRESSURE )

   KRATOS_DEFINE_APPLICATION_VARIABLE( CONSTITUTIVE_MODELS_APPLICATION, double, PLASTIC_VOL_DEF )
   KRATOS_DEFINE_APPLICATION_VARIABLE( CONSTITUTIVE_MODELS_APPLICATION, double, NONLOCAL_PLASTIC_VOL_DEF )
   KRATOS_DEFINE_APPLICATION_VARIABLE( CONSTITUTIVE_MODELS_APPLICATION, double, PLASTIC_VOL_DEF_ABS )
   KRATOS_DEFINE_APPLICATION_VARIABLE( CONSTITUTIVE_MODELS_APPLICATION, double, NONLOCAL_PLASTIC_VOL_DEF_ABS )
   KRATOS_DEFINE_APPLICATION_VARIABLE( CONSTITUTIVE_MODELS_APPLICATION, double, PLASTIC_DEV_DEF )
   KRATOS_DEFINE_APPLICATION_VARIABLE( CONSTITUTIVE_MODELS_APPLICATION, double, NONLOCAL_PLASTIC_DEV_DEF )

   KRATOS_DEFINE_APPLICATION_VARIABLE( CONSTITUTIVE_MODELS_APPLICATION, double, KSIM )
   KRATOS_DEFINE_APPLICATION_VARIABLE( CONSTITUTIVE_MODELS_APPLICATION, double, PS )
   KRATOS_DEFINE_APPLICATION_VARIABLE( CONSTITUTIVE_MODELS_APPLICATION, double, PT )
   KRATOS_DEFINE_APPLICATION_VARIABLE( CONSTITUTIVE_MODELS_APPLICATION, double, PM )
  KRATOS_DEFINE_APPLICATION_VARIABLE( CONSTITUTIVE_MODELS_APPLICATION, std::string, TEMPERATURE_VARIABLE )
  KRATOS_DEFINE_APPLICATION_VARIABLE( CONSTITUTIVE_MODELS_APPLICATION, std::string, PRESSURE_VARIABLE )

  KRATOS_DEFINE_APPLICATION_VARIABLE( CONSTITUTIVE_MODELS_APPLICATION, PropertiesLayout, PROPERTIES_LAYOUT )
}

#endif	/* KRATOS_CONSTITUTIVE_MODELS_APPLICATION_VARIABLES_H_INCLUDED */
