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


namespace Kratos
{
	namespace Testing
	{

		TestSuite::TestSuite(std::string const& Name) : TestCase(Name), mTestCases() {}

		TestSuite::~TestSuite() {} // The test cases are belong to tester class and deleted by it

		void TestSuite::AddTestCase(TestCase* pTestCase) {
			mTestCases.push_back(pTestCase);
		}


		void TestSuite::Reset()	{
			for (auto & mTestCase : mTestCases)
				mTestCase->Reset();
		} 

		void TestSuite::ResetResult() {
			for (auto & mTestCase : mTestCases)
				mTestCase->ResetResult();
		}

		void TestSuite::Run() {
			for (auto & mTestCase : mTestCases)
				mTestCase->Run();
		}

		void TestSuite::Profile() {
			for (auto & mTestCase : mTestCases)
				mTestCase->Profile();
		}

		void TestSuite::Enable() {
			for (auto & mTestCase : mTestCases)
				mTestCase->Enable();
		}

		void TestSuite::Disable() {
			for (auto & mTestCase : mTestCases)
				mTestCase->Disable();
		}

		void TestSuite::Select() {
			for (auto & mTestCase : mTestCases)
				mTestCase->Select();
		}

		void TestSuite::UnSelect() {
			for (auto & mTestCase : mTestCases)
				mTestCase->UnSelect();
		}


		std::string TestSuite::Info() const
		{
			return "Test suite " + Name();
		}

		void TestSuite::PrintInfo(std::ostream& rOStream) const
		{
			rOStream << Info();
		}

		void TestSuite::PrintData(std::ostream& rOStream) const
		{
			rOStream << "This test suite includes following thest cases:" << std::endl;
			for (auto mTestCase : mTestCases)
				rOStream << "    " << mTestCase->Info();
		}


		void TestSuite::TestFunction()
		{
			KRATOS_ERROR << R"(Attempting to call the TestFunction of the the test suite ")" << Name() << R"(")" << std::endl;
		}

	} // manespace Testing.
}  // namespace Kratos.
