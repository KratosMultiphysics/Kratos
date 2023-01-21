//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   license: HDF5Application/license.txt
//
//  Main author:     Máté Kelemen
//

// Core includes
#include "includes/define.h"
#include "testing/testing.h"

// Internal includes
#include "custom_utilities/pipe.h"


namespace Kratos::Testing {


struct NoOpTestPipe : public Pipes::Traits<char,char>
{
    char operator()(char argument) const {return argument;}
}; // struct NoOpTestPipe


struct BoolIntTestPipe : public Pipes::Traits<bool,int>
{
    BoolIntTestPipe()
        : mFactor(1)
    {}

    BoolIntTestPipe(const Parameters& rParameters)
        : mFactor(rParameters["factor"].GetInt())
    {}

    BoolIntTestPipe(BoolIntTestPipe&& rOther) noexcept = default;

    BoolIntTestPipe(const BoolIntTestPipe& rOther) = default;

    int operator()(bool argument) const {return mFactor * argument;}

private:
    int mFactor;
}; // struct struct BoolIntTestPipe


struct IntBoolTestPipe : public Pipes::Traits<int,bool>
{
    IntBoolTestPipe()
        : mNegate(false)
    {}

    IntBoolTestPipe(const Parameters& rParameters)
        : mNegate(rParameters["negate"].GetBool())
    {}

    IntBoolTestPipe(IntBoolTestPipe&& rOther) noexcept = default;

    IntBoolTestPipe(const IntBoolTestPipe& rOther) = default;

    bool operator()(int argument) const {return mNegate ? !bool(argument) : argument;}

private:
    bool mNegate;
}; // struct IntBoolTestPipe


KRATOS_TEST_CASE_IN_SUITE(PipeTraits, KratosHDF5TestSuite)
{
    static_assert(Pipes::IsPipe<BoolIntTestPipe>::value);
    static_assert(Pipes::IsPipe<IntBoolTestPipe>::value);
}


KRATOS_TEST_CASE_IN_SUITE(PipeOperator, KratosHDF5TestSuite)
{
    BoolIntTestPipe bool_int;
    IntBoolTestPipe int_bool;

    {
        const bool output = false >> bool_int >> int_bool;
        KRATOS_CHECK_EQUAL(output, false);
    }

    {
        const bool output = true >> bool_int >> int_bool;
        KRATOS_CHECK_EQUAL(output, true);
    }
}


KRATOS_TEST_CASE_IN_SUITE(CompoundPipe, KratosHDF5TestSuite)
{
    BoolIntTestPipe bool_int;
    IntBoolTestPipe int_bool;

    constexpr const bool is_same = std::is_same_v<
        Pipes::Pipeline<BoolIntTestPipe,IntBoolTestPipe,BoolIntTestPipe>,
        decltype(bool_int | int_bool | bool_int)
    >;
    KRATOS_CHECK(is_same);

    // Instantiate compound pipes
    {
        const auto compound_test_pipe = NoOpTestPipe() | IntBoolTestPipe(); // NoOpTestPipe is not constructible from Parameters
        static_assert(Pipes::IsPipe<decltype(compound_test_pipe)>::value);
        KRATOS_CHECK_EQUAL(compound_test_pipe('a'), true);
        KRATOS_CHECK_EQUAL(compound_test_pipe(0), false);
    }

    {
        const auto compound_test_pipe = BoolIntTestPipe() | IntBoolTestPipe();
        static_assert(Pipes::IsPipe<decltype(compound_test_pipe)>::value);
        KRATOS_CHECK_EQUAL(compound_test_pipe(false), false);
        KRATOS_CHECK_EQUAL(compound_test_pipe(true), true);
    }

    {
        const auto compound_test_pipe = bool_int | IntBoolTestPipe();
        static_assert(Pipes::IsPipe<decltype(compound_test_pipe)>::value);
        KRATOS_CHECK_EQUAL(compound_test_pipe(false), false);
        KRATOS_CHECK_EQUAL(compound_test_pipe(true), true);
    }

    {
        const auto compound_test_pipe = BoolIntTestPipe() | int_bool;
        static_assert(Pipes::IsPipe<decltype(compound_test_pipe)>::value);
        KRATOS_CHECK_EQUAL(compound_test_pipe(false), false);
        KRATOS_CHECK_EQUAL(compound_test_pipe(true), true);
    }

    {
        const auto compound_test_pipe = bool_int | int_bool;
        static_assert(Pipes::IsPipe<decltype(compound_test_pipe)>::value);
        KRATOS_CHECK_EQUAL(compound_test_pipe(false), false);
        KRATOS_CHECK_EQUAL(compound_test_pipe(true), true);
    }

    {
        const auto compound_test_pipe = bool_int | int_bool | BoolIntTestPipe() | IntBoolTestPipe();
        static_assert(Pipes::IsPipe<decltype(compound_test_pipe)>::value);
        KRATOS_CHECK_EQUAL(compound_test_pipe(false), false);
        KRATOS_CHECK_EQUAL(compound_test_pipe(true), true);
    }
}


KRATOS_TEST_CASE_IN_SUITE(PipeFactory, KratosHDF5TestSuite)
{
    {
        using Pipeline = decltype(BoolIntTestPipe() | IntBoolTestPipe());
        Parameters parameters(R"([
            {"factor" : 2},
            {"negate" : true}
        ])");
        const auto pipe = Pipeline(parameters);
        KRATOS_CHECK_EQUAL(pipe(false), true);
        KRATOS_CHECK_EQUAL(pipe(true), false);
    }

    {
        using Pipeline = decltype(BoolIntTestPipe() | IntBoolTestPipe() | BoolIntTestPipe() | IntBoolTestPipe());
        Parameters parameters(R"([
            {"factor" : 2},
            {"negate" : true},
            {"factor" : 2},
            {"negate" : true}
        ])");
        const auto pipe = Pipeline(parameters);
        KRATOS_CHECK_EQUAL(pipe(false), false);
        KRATOS_CHECK_EQUAL(pipe(true), true);
    }
}


} // namespace Kratos::Testing
