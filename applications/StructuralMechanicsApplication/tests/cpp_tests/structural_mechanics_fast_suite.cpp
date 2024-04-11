// KRATOS  ___|  |                   |                   |
//       \___ \  __|  __| |   |  __| __| |   |  __| _` | |
//             | |   |    |   | (    |   |   | |   (   | |
//       _____/ \__|_|   \__,_|\___|\__|\__,_|_|  \__,_|_| MECHANICS
//
//  License:         BSD License
//                   license: StructuralMechanicsApplication/license.txt
//
//  Main authors:    Richard Faasse
//

// External includes
#include <gtest/gtest.h>
#include <gmock/gmock.h>

// Project includes
#include "structural_mechanics_fast_suite.h"

int main(int argc, char* argv[]) 
{
    ::testing::InitGoogleTest(&argc, argv);

    Kratos::Testing::mApplicationInitializerList.push_back([](std::vector<Kratos::KratosApplication::Pointer> & rRegisteredApplications, Kratos::Kernel & rKernel) {
      if (!rKernel.IsImported("StructuralMechanicsApplication")) {
        auto pApplication = std::make_shared<Kratos::KratosStructuralMechanicsApplication>();
        rKernel.ImportApplication(pApplication);
        rRegisteredApplications.push_back(std::move(pApplication));
      }
    });

    return RUN_ALL_TESTS();
}