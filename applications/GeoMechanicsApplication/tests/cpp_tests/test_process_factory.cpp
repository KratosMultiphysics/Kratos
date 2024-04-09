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
#include "geo_mechanics_fast_suite.h"

using namespace Kratos;


namespace
{

class FooProcess : public Process
{
};

class BarProcess : public Process
{
};

}


namespace Kratos::Testing
{

KRATOS_TEST_CASE_IN_SUITE(CreateNothingWhenNoCreatorWasAddedForRequestedProcess, KratosGeoMechanicsFastSuite)
{
    ProcessFactory factory;

    const Parameters process_settings;
    const auto process = factory.Create("UnknownProcess", process_settings);

    KRATOS_EXPECT_EQ(process.get(), nullptr);
}

KRATOS_TEST_CASE_IN_SUITE(CreateThrowsForUnknownProcess_WhenCallbackFunctionThrows, KratosGeoMechanicsFastSuite)
{
    ProcessFactory factory;
    factory.SetCallBackWhenProcessIsUnknown([](const std::string& rProcessName){ KRATOS_ERROR << "Unexpected process (" << rProcessName << "), calculation is aborted";});

    const Parameters process_settings;
    KRATOS_EXPECT_EXCEPTION_IS_THROWN(const auto process = factory.Create("UnknownProcess", process_settings),
                                      "Unexpected process (UnknownProcess), calculation is aborted")
}

KRATOS_TEST_CASE_IN_SUITE(CreateNothingWhenTheAddedCreatorIsEmpty, KratosGeoMechanicsFastSuite)
{
    ProcessFactory factory;
    factory.AddCreator("TestProcess", {});

    const Parameters process_settings;
    const auto process = factory.Create("TestProcess", process_settings);

    KRATOS_EXPECT_EQ(process.get(), nullptr);
}

KRATOS_TEST_CASE_IN_SUITE(CreateFooProcessAfterCreatorWasAdded, KratosGeoMechanicsFastSuite)
{
    ProcessFactory factory;
    factory.AddCreator("FooProcess", [](const Parameters&){return std::make_unique<FooProcess>();});

    const Parameters process_settings;
    const auto process = factory.Create("FooProcess", process_settings);

    KRATOS_EXPECT_NE(dynamic_cast<FooProcess*>(process.get()), nullptr);
    KRATOS_EXPECT_EQ(dynamic_cast<BarProcess*>(process.get()), nullptr);
}

KRATOS_TEST_CASE_IN_SUITE(CreateDifferentKindsOfProcessesAfterAddingCreators, KratosGeoMechanicsFastSuite)
{
    ProcessFactory factory;
    factory.AddCreator("FooProcess", [](const Parameters&){return std::make_unique<FooProcess>();});
    factory.AddCreator("BarProcess", [](const Parameters&){return std::make_unique<BarProcess>();});

    const Parameters foo_process_settings;
    const auto foo_process = factory.Create("FooProcess", foo_process_settings);
    const Parameters bar_process_settings;
    const auto bar_process = factory.Create("BarProcess", bar_process_settings);

    KRATOS_EXPECT_NE(dynamic_cast<FooProcess*>(foo_process.get()), nullptr);
    KRATOS_EXPECT_NE(dynamic_cast<BarProcess*>(bar_process.get()), nullptr);
}

KRATOS_TEST_CASE_IN_SUITE(CreatorReceivesGivenSettingsWhenMakingNewProcess, KratosGeoMechanicsFastSuite)
{
    ProcessFactory factory;
    const Parameters* received_settings = nullptr;
    auto creator = [&received_settings](const Parameters& rProcessSettings){
        received_settings = &rProcessSettings;
        return std::make_unique<FooProcess>();
    };
    factory.AddCreator("FooProcess", creator);

    const Parameters process_settings;
    const auto foo_process = factory.Create("FooProcess", process_settings);

    KRATOS_EXPECT_EQ(&process_settings, received_settings);
}

KRATOS_TEST_CASE_IN_SUITE(AddingCreatorWithSameProcessNameOverwritesPreviousOne, KratosGeoMechanicsFastSuite)
{
    ProcessFactory factory;
    auto creator_id = 0;
    auto creator1 = [&creator_id](const Parameters&){
        creator_id = 1;
        return std::make_unique<FooProcess>();
    };
    auto creator2 = [&creator_id](const Parameters&){
        creator_id = 2;
        return std::make_unique<FooProcess>();
    };

    factory.AddCreator("FooProcess", creator1);
    auto foo_process = factory.Create("FooProcess", {});
    KRATOS_EXPECT_EQ(1, creator_id);

    factory.AddCreator("FooProcess", creator2);
    foo_process = factory.Create("FooProcess", {});
    KRATOS_EXPECT_EQ(2, creator_id);
}

}
