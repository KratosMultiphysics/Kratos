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

  ///@name Kratos Globals
  ///@{

  //Create Variables
  KRATOS_CREATE_VARIABLE( Vector, MATERIAL_PERCENT_COMPOSITION )
  KRATOS_CREATE_VARIABLE( double, PRESSURE_VELOCITY )
  KRATOS_CREATE_VARIABLE( double, PRESSURE_ACCELERATION )

  ///@}

}
