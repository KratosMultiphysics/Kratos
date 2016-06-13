//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ \
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics 
//
//  License:		 BSD License 
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//                    
//
	           

// System includes


// External includes 


// Project includes
#include "testing/test_case.h"
#include "includes/exception.h"


namespace Kratos
{
	namespace Testing
	{
		TestCase::TestCase(std::string Name)
			:mName(Name), mIsEnabled(true) {}

		TestCase::~TestCase() {}

		void TestCase::Reset()
		{
			mResult.Reset();
			mIsEnabled = true;
		}

		void TestCase::ResetResult()
		{
			mResult.Reset();
		}

		void TestCase::Setup() {}

		void TestCase::Run()
		{
			try {
				Setup();
				TestFunction();
				TearDown();
				mResult.SetToSucceed();
			}
			catch (Exception& e) {
				std::stringstream buffer;
				buffer << e.what() << " in " << e.where() << std::endl;
				mResult.SetErrorMessage(buffer.str());
			}
			catch (std::exception& e) {
				mResult.SetErrorMessage(e.what());
			}
			catch (...) {
				mResult.SetErrorMessage("Unknown error");
			}
		}

		void TestCase::Profile()
		{
			try {
				Setup();
				TestFunction();
				TearDown();
			}
			catch (...) {
				std::cout << "Test " << Name() << " Failed!" << std::endl;
			}
		}

		void TestCase::TearDown() {}

		void TestCase::Enable()
		{
			mIsEnabled = true;
		}

		void TestCase::Disable() 
		{
			mIsEnabled = false;
		}

		const std::string& TestCase::Name()
		{
			return mName;
		}

		const TestCaseResult& TestCase::GetResult() const
		{
			return mResult;
		}

		void TestCase::SetResult(TestCaseResult const& TheResult)
		{
			mResult = TheResult;
		}

		bool TestCase::IsEnabled()
		{
			return mIsEnabled;
		}

		bool TestCase::IsDisabled()
		{
			return !mIsEnabled;
		}

		std::string TestCase::Info() const
		{
			return "Test case " + mName;
		}

		void TestCase::PrintInfo(std::ostream& rOStream) const
		{
			rOStream << Info();
		}

		void TestCase::PrintData(std::ostream& rOStream) const
		{
		}




	} // manespace Testing.
}  // namespace Kratos.


