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

// namespace Kratos::Testing
// {

// KratosGeoMechanicsFastSuite::KratosGeoMechanicsFastSuite() : KratosCoreFastSuite()
// {
//     if (!this->mKernel.IsImported("GeoMechanicsApplication")) {
//         mpGeoApp = std::make_shared<KratosGeoMechanicsApplication>();
//         this->mKernel.ImportApplication(mpGeoApp);
//     }
//     if (!this->mKernel.IsImported("LinearSolversApplication")) {
//         mpLinearSolversApp = std::make_shared<KratosLinearSolversApplication>();
//         this->mKernel.ImportApplication(mpLinearSolversApp);
//     }
// }

// KratosGeoMechanicsIntegrationSuite::KratosGeoMechanicsIntegrationSuite() : KratosCoreFastSuite()
// {
//     if (!this->mKernel.IsImported("GeoMechanicsApplication")) {
//         mpGeoApp = std::make_shared<KratosGeoMechanicsApplication>();
//         this->mKernel.ImportApplication(mpGeoApp);
//     }
//     if (!this->mKernel.IsImported("LinearSolversApplication")) {
//         mpLinearSolversApp = std::make_shared<KratosLinearSolversApplication>();
//         this->mKernel.ImportApplication(mpLinearSolversApp);
//     }
// }

// } // namespace Kratos::Testing

// External includes
#include <gtest/gtest.h>
#include <gmock/gmock.h>

// Project includes
#include "testing/testing.h"
#include "geo_mechanics_fast_suite.h"

int main(int argc, char* argv[]) 
{
    ::testing::InitGoogleTest(&argc, argv);

    Kratos::Testing::mApplicationInitializerList.push_back([](std::vector<Kratos::KratosApplication::Pointer> & rRegisteredApplications, Kratos::Kernel & rKernel) {
      if (!rKernel.IsImported("GeoMechanicsApplication")) {
        auto pGeoApp = std::make_shared<Kratos::KratosGeoMechanicsApplication>();
        auto pLinearSolversApp = std::make_shared<Kratos::KratosLinearSolversApplication>();
        rKernel.ImportApplication(pGeoApp);
        rKernel.ImportApplication(pLinearSolversApp);
        rRegisteredApplications.push_back(std::move(pGeoApp));
        rRegisteredApplications.push_back(std::move(pLinearSolversApp));
      }
    });

    return RUN_ALL_TESTS();
}
