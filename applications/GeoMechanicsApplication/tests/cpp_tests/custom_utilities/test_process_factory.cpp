// KRATOS___
//     //   ) )
//    //         ___      ___
//   //  ____  //___) ) //   ) )
//  //    / / //       //   / /
// ((____/ / ((____   ((___/ /  MECHANICS
//
//  License:         geo_mechanics_application/license.txt
//
//  Main authors:    Wijtze Pieter Kikstra
//                   Anne van de Graaf
//

#include "custom_utilities/process_factory.hpp"
#include "includes/expect.h"
#include "tests/cpp_tests/geo_mechanics_fast_suite_without_kernel.h"

using namespace Kratos;

namespace
{

class FooProcess : public Process
{
};

class BarProcess : public Process
{
};

} // namespace

namespace Kratos::Testing
{

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, CreateNothingWhenNoCreatorWasAddedForRequestedProcess)
{
    ProcessFactory factory;

    const Parameters process_settings;
    const auto       process = factory.Create("UnknownProcess", process_settings);

    EXPECT_EQ(process.get(), nullptr);
}

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, CreateThrowsForUnknownProcess_WhenCallbackFunctionThrows)
{
    ProcessFactory factory;
    factory.SetCallBackWhenProcessIsUnknown([](const std::string& rProcessName) {
        KRATOS_ERROR << "Unexpected process (" << rProcessName << "), calculation is aborted";
    });

    const Parameters process_settings;
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(const auto process = factory.Create("UnknownProcess", process_settings),
                                      "Unexpected process (UnknownProcess), calculation is aborted")
}

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, CreateNothingWhenTheAddedCreatorIsEmpty)
{
    ProcessFactory factory;
    factory.AddCreator("TestProcess", {});

    const Parameters process_settings;
    const auto       process = factory.Create("TestProcess", process_settings);

    EXPECT_EQ(process.get(), nullptr);
}

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, CreateFooProcessAfterCreatorWasAdded)
{
    ProcessFactory factory;
    factory.AddCreator("FooProcess", [](const Parameters&) { return std::make_unique<FooProcess>(); });

    const Parameters process_settings;
    const auto       process = factory.Create("FooProcess", process_settings);

    EXPECT_NE(dynamic_cast<FooProcess*>(process.get()), nullptr);
    EXPECT_EQ(dynamic_cast<BarProcess*>(process.get()), nullptr);
}

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, CreateDifferentKindsOfProcessesAfterAddingCreators)
{
    ProcessFactory factory;
    factory.AddCreator("FooProcess", [](const Parameters&) { return std::make_unique<FooProcess>(); });
    factory.AddCreator("BarProcess", [](const Parameters&) { return std::make_unique<BarProcess>(); });

    const Parameters foo_process_settings;
    const auto       foo_process = factory.Create("FooProcess", foo_process_settings);
    const Parameters bar_process_settings;
    const auto       bar_process = factory.Create("BarProcess", bar_process_settings);

    EXPECT_NE(dynamic_cast<FooProcess*>(foo_process.get()), nullptr);
    EXPECT_NE(dynamic_cast<BarProcess*>(bar_process.get()), nullptr);
}

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, CreatorReceivesGivenSettingsWhenMakingNewProcess)
{
    ProcessFactory    factory;
    const Parameters* received_settings = nullptr;
    auto              creator           = [&received_settings](const Parameters& rProcessSettings) {
        received_settings = &rProcessSettings;
        return std::make_unique<FooProcess>();
    };
    factory.AddCreator("FooProcess", creator);

    const Parameters process_settings;
    const auto       foo_process = factory.Create("FooProcess", process_settings);

    EXPECT_EQ(&process_settings, received_settings);
}

TEST_F(KratosGeoMechanicsFastSuiteWithoutKernel, AddingCreatorWithSameProcessNameOverwritesPreviousOne)
{
    ProcessFactory factory;
    auto           creator_id = 0;
    auto           creator1   = [&creator_id](const Parameters&) {
        creator_id = 1;
        return std::make_unique<FooProcess>();
    };
    auto creator2 = [&creator_id](const Parameters&) {
        creator_id = 2;
        return std::make_unique<FooProcess>();
    };

    factory.AddCreator("FooProcess", creator1);
    auto foo_process = factory.Create("FooProcess", {});
    EXPECT_EQ(1, creator_id);

    factory.AddCreator("FooProcess", creator2);
    foo_process = factory.Create("FooProcess", {});
    EXPECT_EQ(2, creator_id);
}

} // namespace Kratos::Testing
