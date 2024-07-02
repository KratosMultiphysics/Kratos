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
#include "constitutive_laws_application.h"
#include "structural_mechanics_application.h"

int main(int argc, char* argv[]) 
{
    ::testing::InitGoogleTest(&argc, argv);

    Kratos::Testing::mApplicationInitializerList.push_back([](std::vector<Kratos::KratosApplication::Pointer> & rRegisteredApplications, Kratos::Kernel & rKernel) {
      if (!rKernel.IsImported("StructuralMechanicsApplication")) {
        auto pApplication = std::make_shared<Kratos::KratosStructuralMechanicsApplication>();
        rKernel.ImportApplication(pApplication);
        rRegisteredApplications.push_back(std::move(pApplication));
      }

      if (!rKernel.IsImported("ConstitutiveLawsApplication")) {
        auto pApplication = std::make_shared<Kratos::KratosConstitutiveLawsApplication>();
        rKernel.ImportApplication(pApplication);
        rRegisteredApplications.push_back(std::move(pApplication));
      }
    });

    return RUN_ALL_TESTS();
}
