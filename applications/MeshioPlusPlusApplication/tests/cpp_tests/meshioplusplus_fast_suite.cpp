//  ██████   ██████ ██████████  █████████  █████   █████ █████    ███████
// ░░██████ ██████ ░░███░░░░░█ ███░░░░░███░░███   ░░███ ░░███   ███░░░░░███      ███         ███
//  ░███░█████░███  ░███  █ ░ ░███    ░░░  ░███    ░███  ░███  ███     ░░███    ░███        ░███
//  ░███░░███ ░███  ░██████   ░░█████████  ░███████████  ░███ ░███      ░███ ███████████ ███████████
//  ░███ ░░░  ░███  ░███░░█    ░░░░░░░░███ ░███░░░░░███  ░███ ░███      ░███░░░░░███░░░ ░░░░░███░░░
//  ░███      ░███  ░███ ░   █ ███    ░███ ░███    ░███  ░███ ░░███     ███     ░███        ░███
//  █████     █████ ██████████░░█████████  █████   █████ █████ ░░░███████░      ░░░         ░░░
// ░░░░░     ░░░░░ ░░░░░░░░░░  ░░░░░░░░░  ░░░░░   ░░░░░ ░░░░░    ░░░░░░░                        Application
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Vicente Mataix Ferrandiz
//

// External includes
#include <gtest/gtest.h>
#include <gmock/gmock.h>

// Project includes
#include "testing/testing.h"
#include "meshioplusplus_application.h"

int main(int argc, char* argv[])
{
    ::testing::InitGoogleTest(&argc, argv);

    Kratos::Testing::mApplicationInitializerList.push_back([](std::vector<Kratos::KratosApplication::Pointer>& rRegisteredApplications, Kratos::Kernel& rKernel) {
      if (!rKernel.IsImported("MeshioPlusPlusApplication")) {
        auto pApplication = std::make_shared<Kratos::KratosMeshioPlusPlusApplication>();
        rKernel.ImportApplication(pApplication);
        rRegisteredApplications.push_back(std::move(pApplication));
      }
    });

    return RUN_ALL_TESTS();
}
