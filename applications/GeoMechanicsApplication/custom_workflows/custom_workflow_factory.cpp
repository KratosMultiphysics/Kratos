// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Richard Faasse
//

#include "custom_workflow_factory.h"
#include "dgeosettlement.h"
#include "custom_utilities/input_utility.h"
#include <memory>

namespace Kratos
{

KratosGeoSettlement* CustomWorkflowFactory::CreateKratosGeoSettlement()
{
    auto inputUtility = std::make_unique<InputUtility>();
    return new KratosGeoSettlement(std::move(inputUtility));
}
}
