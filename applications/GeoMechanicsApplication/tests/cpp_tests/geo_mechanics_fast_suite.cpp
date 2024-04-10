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
    if (!this->mKernel.IsImported("GeoMechanicsApplication")) {
        mpGeoApp = std::make_shared<KratosGeoMechanicsApplication>();
        this->mKernel.ImportApplication(mpGeoApp);
    }
    if (!this->mKernel.IsImported("LinearSolversApplication")) {
        mpLinearSolversApp = std::make_shared<KratosLinearSolversApplication>();
        this->mKernel.ImportApplication(mpLinearSolversApp);
    }
}

KratosGeoMechanicsIntegrationSuite::KratosGeoMechanicsIntegrationSuite() : KratosCoreFastSuite()
{
    if (!this->mKernel.IsImported("GeoMechanicsApplication")) {
        mpGeoApp = std::make_shared<KratosGeoMechanicsApplication>();
        this->mKernel.ImportApplication(mpGeoApp);
    }
    if (!this->mKernel.IsImported("LinearSolversApplication")) {
        mpLinearSolversApp = std::make_shared<KratosLinearSolversApplication>();
        this->mKernel.ImportApplication(mpLinearSolversApp);
    }
}

} // namespace Kratos::Testing
