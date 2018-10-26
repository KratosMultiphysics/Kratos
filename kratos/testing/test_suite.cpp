//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//
//

// Project includes
#include "testing/test_suite.h"
#include "includes/checks.h"

namespace Kratos {
namespace Testing {

TestSuite::TestSuite(std::string const& Name) : TestCase(Name), mTestCases() {}

TestSuite::~TestSuite() {
}  // The test cases are belong to tester class and deleted by it

void TestSuite::AddTestCase(TestCase* pTestCase) {
    mTestCases.push_back(pTestCase);
}

void TestSuite::Reset() {
    for (auto i_test = mTestCases.begin(); i_test != mTestCases.end(); i_test++)
        (*i_test)->Reset();
}

void TestSuite::ResetResult() {
    for (auto i_test = mTestCases.begin(); i_test != mTestCases.end(); i_test++)
        (*i_test)->ResetResult();
}

void TestSuite::Run() {
    for (auto i_test = mTestCases.begin(); i_test != mTestCases.end(); i_test++)
    {
        (*i_test)->Run();
    }
}

void TestSuite::Profile() {
    for (auto i_test = mTestCases.begin(); i_test != mTestCases.end(); i_test++)
    {
        (*i_test)->Profile();
    }
}

void TestSuite::Enable() {
    for (auto i_test = mTestCases.begin(); i_test != mTestCases.end(); i_test++)
        (*i_test)->Enable();
}

void TestSuite::Disable() {
    for (auto i_test = mTestCases.begin(); i_test != mTestCases.end(); i_test++)
        (*i_test)->Disable();
}

void TestSuite::Select() {
    for (auto i_test = mTestCases.begin(); i_test != mTestCases.end(); i_test++)
        if ((*i_test)->IsEnabled())
            (*i_test)->Select();
        else
            (*i_test)->UnSelect();
}

void TestSuite::UnSelect() {
    for (auto i_test = mTestCases.begin(); i_test != mTestCases.end(); i_test++)
        (*i_test)->UnSelect();
}

std::string TestSuite::Info() const { return "Test suite " + Name(); }

void TestSuite::PrintInfo(std::ostream& rOStream) const { rOStream << Info(); }

void TestSuite::PrintData(std::ostream& rOStream) const {
    rOStream << "This test suite includes following thest cases:" << std::endl;
    for (auto i_test = mTestCases.begin(); i_test != mTestCases.end(); i_test++)
        rOStream << "    " << (*i_test)->Info();
}

void TestSuite::TestFunction() {
    KRATOS_ERROR
        << "Attempting to call the TestFunction of the the test suite \""
        << Name() << "\"" << std::endl;
}

}  // manespace Testing.
}  // namespace Kratos.
