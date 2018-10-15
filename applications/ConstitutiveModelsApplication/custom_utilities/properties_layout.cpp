//
//   Project Name:        KratosConstitutiveModelsApplication $
//   Created by:          $Author:                JMCarbonell $
//   Last modified by:    $Co-Author:                         $
//   Date:                $Date:                 October 2018 $
//   Revision:            $Revision:                      0.0 $
//
//

// System includes

// External includes

// Project includes
#include "custom_utilities/properties_layout.hpp"

#include "constitutive_models_application_variables.h"


namespace Kratos
{

void PropertiesLayout::Configure(const Properties& rProperties, const GeometryType& rGeometry, const Vector& rShapeFunctions)
{
  mpData = &(rProperties.Data());
  mpTables = &(rProperties.Tables());
   
  const TableKeyVariables<double,double>::Pointer& ScalarVariables = mpData->GetValue(TABLES_SCALAR_VARIABLES);
    
  for(auto it = mpTables->begin(); it != mpTables->end(); ++it)
  {
    double Variable = 0.0;
    for(std::size_t j=1; j<rShapeFunctions.size(); ++j)
      Variable = rShapeFunctions[j] * rGeometry[j].FastGetSolutionStepValue(ScalarVariables->GetXVariable(std::size_t(it->first)));
    mTableArguments.push_back(ScalarTableArgumentsType(it->first, VariableKeyArgumentsType(ScalarVariables->GetYVariable(std::size_t(it->first)).Key(), Variable)));
  }
}

} // Namespace Kratos
