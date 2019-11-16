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


// System includes
#include <chrono>
#include <sstream>

// External includes


// Project includes
#include "testing/test_case.h"
#include "testing/test_skipped_exception.h"
#include "includes/exception.h"


namespace Kratos
{
	namespace Testing
	{
		TestCase::TestCase(std::string const& Name)
			:mName(Name), mIsEnabled(true), mIsSelected(false) {}

		TestCase::~TestCase() {}

		void TestCase::Reset()
		{
			mResult.Reset();
			mIsEnabled = true;
			mIsSelected = false;
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
			catch (TestSkippedException& e) {
				mResult.SetToSkipped();
				mResult.SetErrorMessage(e.what());
			}
			catch (Exception& e) {
				mResult.SetToFailed();
				mResult.SetErrorMessage(e.what());
			}
			catch (std::exception& e) {
				mResult.SetToFailed();
				mResult.SetErrorMessage(e.what());
			}
			catch (...) {
				mResult.SetToFailed();
				mResult.SetErrorMessage("Unknown error");
			}
		}

		void TestCase::Profile()
		{
			try {
				auto start = std::chrono::steady_clock::now();
				Setup();
				auto end_setup = std::chrono::steady_clock::now();
				TestFunction();
				auto end_run = std::chrono::steady_clock::now();
				TearDown();
				auto end_tear_down = std::chrono::steady_clock::now();
				mResult.SetToSucceed();
				auto end = std::chrono::steady_clock::now();
				std::chrono::duration<double> setup_elapsed = end_setup - start;
				std::chrono::duration<double> run_elapsed = end_run - end_setup;
				std::chrono::duration<double> tear_down_elapsed = end_tear_down - end_run;
				std::chrono::duration<double> elapsed = end - start;
				mResult.SetSetupElapsedTime(setup_elapsed.count());
				mResult.SetRunElapsedTime(run_elapsed.count());
				mResult.SetTearDownElapsedTime(tear_down_elapsed.count());
				mResult.SetElapsedTime(elapsed.count());
			}
			catch (TestSkippedException& e) {
				mResult.SetToSkipped();
				mResult.SetErrorMessage(e.what());
			}
			catch (Exception& e) {
				mResult.SetToFailed();
				mResult.SetErrorMessage(e.what());
			}
			catch (std::exception& e) {
				mResult.SetToFailed();
				mResult.SetErrorMessage(e.what());
			}
			catch (...) {
				mResult.SetToFailed();
				mResult.SetErrorMessage("Unknown error");
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

		void TestCase::Select()
		{
			mIsSelected = true;
		}

		void TestCase::UnSelect()
		{
			mIsSelected = false;
		}

		const std::string& TestCase::Name() const
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


		void TestCase::SetResultOutput(std::string const& TheResultOutput) {
			mResult.SetOutput(TheResultOutput);
		}

		bool TestCase::IsEnabled() const
		{
			return mIsEnabled;
		}

		bool TestCase::IsDisabled() const
		{
			return !mIsEnabled;
		}

		bool TestCase::IsSelected() const
		{
			return mIsSelected;
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
