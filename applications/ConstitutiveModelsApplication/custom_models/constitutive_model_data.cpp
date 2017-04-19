//
//   Project Name:        KratosConstitutiveModelsApplication $
//   Created by:          $Author:                JMCarbonell $
//   Last modified by:    $Co-Author:                         $
//   Date:                $Date:                   April 2017 $
//   Revision:            $Revision:                      0.0 $
//
//

// System includes

// External includes

// Project includes
#include "custom_models/constitutive_model_data.hpp"

namespace Kratos
{
  
  KRATOS_CREATE_LOCAL_FLAG( ConstitutiveModelData, IMPLEX_ACTIVE,                0 );

  KRATOS_CREATE_LOCAL_FLAG( ConstitutiveModelData, COMPUTED_STRAIN,              1 );
  KRATOS_CREATE_LOCAL_FLAG( ConstitutiveModelData, COMPUTED_STRESS,              2 );
  KRATOS_CREATE_LOCAL_FLAG( ConstitutiveModelData, COMPUTED_CONSTITUTIVE_MATRIX, 3 );

  KRATOS_CREATE_LOCAL_FLAG( ConstitutiveModelData, PLASTIC_REGION,               4 );
  KRATOS_CREATE_LOCAL_FLAG( ConstitutiveModelData, PLASTIC_RATE_REGION,          5 );
  KRATOS_CREATE_LOCAL_FLAG( ConstitutiveModelData, COMPUTED_RETURN_MAPPING,      6 );
  KRATOS_CREATE_LOCAL_FLAG( ConstitutiveModelData, UPDATE_INTERNAL_VARIABLES,    7 );
  
}  // namespace Kratos.




