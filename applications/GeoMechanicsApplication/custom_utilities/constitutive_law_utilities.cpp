// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Anne van de Graaf
//                   Wijtze Pieter Kikstra
//

#include "custom_utilities/constitutive_law_utilities.hpp"

namespace Kratos
{

int ConstitutiveLawUtilities::GetStateVariableIndex(const Variable<double>& rThisVariable)
{
   int index = -1;
   if (const std::string prefix{"STATE_VARIABLE_"};
       rThisVariable.Name().substr(0, prefix.length()) == prefix) {
       index = std::stoi(rThisVariable.Name().substr(prefix.length()));
   }

   return index - 1;
}

}
