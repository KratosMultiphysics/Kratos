//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Richard Faasse
//

// External includes
#include <gtest/gtest.h>
#include <gmock/gmock.h>

// Project includes
#include "testing/testing.h"
#include "exaqute_sandbox_application.h"

int main(int argc, char* argv[]) 
{
    ::testing::InitGoogleTest(&argc, argv);

    Kratos::Testing::mApplicationInitializerList.push_back([](std::vector<Kratos::KratosApplication::Pointer> & rRegisteredApplications, Kratos::Kernel & rKernel) {
      if (!rKernel.IsImported("ExaquteSandboxApplication")) {
        auto pApplication = std::make_shared<Kratos::KratosExaquteSandboxApplication>();
        rKernel.ImportApplication(pApplication);
        rRegisteredApplications.push_back(std::move(pApplication));
      }
    });

    return RUN_ALL_TESTS();
}

} // namespace Kratos::Testing
