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
//                   Richard Faasse
//

#include "geo_mechanics_fast_suite.h"

namespace Kratos::Testing
{

KratosGeoMechanicsFastSuite::KratosGeoMechanicsFastSuite() : KratosCoreFastSuite()
{
    mpGeoApp = std::make_shared<KratosGeoMechanicsApplication>();
    this->ImportApplicationIntoKernel(mpGeoApp);
    mpLinearSolversApp = std::make_shared<KratosLinearSolversApplication>();
    this->ImportApplicationIntoKernel(mpLinearSolversApp);
}

KratosGeoMechanicsIntegrationSuite::KratosGeoMechanicsIntegrationSuite() : KratosCoreFastSuite()
{
    mpGeoApp = std::make_shared<KratosGeoMechanicsApplication>();
    this->ImportApplicationIntoKernel(mpGeoApp);
    mpLinearSolversApp = std::make_shared<KratosLinearSolversApplication>();
    this->ImportApplicationIntoKernel(mpLinearSolversApp);
}

} // namespace Kratos::Testing
