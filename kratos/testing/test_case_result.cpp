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


// External includes


// Project includes
#include "testing/test_case_result.h"


namespace Kratos
{
	namespace Testing
	{
		TestCaseResult::TestCaseResult() : mStatus(Result::DidNotRun), mOutput(""), mErrorMessage("")
			, mSetupElapsedTime(0.00), mRunElapsedTime(0.00), mTearDownElapsedTime(0.00), mElapsedTime(0.00) {}

		TestCaseResult::TestCaseResult(TestCaseResult const& rOther)
			: mStatus(rOther.mStatus), mOutput(rOther.mOutput), mErrorMessage(rOther.mErrorMessage)
			, mSetupElapsedTime(rOther.mSetupElapsedTime), mRunElapsedTime(rOther.mRunElapsedTime)
			, mTearDownElapsedTime(rOther.mTearDownElapsedTime), mElapsedTime(rOther.mElapsedTime) {}

		TestCaseResult::~TestCaseResult() {}

		TestCaseResult& TestCaseResult::operator=(TestCaseResult const& rOther)
		{
			mStatus = rOther.mStatus;
			mOutput = rOther.mOutput;
			mErrorMessage = rOther.mErrorMessage;
			mSetupElapsedTime = rOther.mSetupElapsedTime;
			mRunElapsedTime = rOther.mRunElapsedTime;
			mTearDownElapsedTime = rOther.mTearDownElapsedTime;
			mElapsedTime = rOther.mElapsedTime;

			return *this;
		}

		void TestCaseResult::Reset()
		{
			mStatus = Result::DidNotRun;
			mOutput = "";
			mErrorMessage = "";
			mSetupElapsedTime = 0.00;
			mRunElapsedTime = 0.00;
			mTearDownElapsedTime = 0.00;
			mElapsedTime = 0.00;
		}

		void TestCaseResult::SetToSucceed()	{
			mStatus = Result::Passed;
		}

		void TestCaseResult::SetToFailed() {
			mStatus = Result::Failed;
		}

		void TestCaseResult::SetToSkipped() {
			mStatus = Result::Skipped;
		}

		void TestCaseResult::SetOutput(const std::string& TheOutput) {
			mOutput = TheOutput;
		}

		const std::string& TestCaseResult::GetOutput() const {
			return mOutput;
		}

		void TestCaseResult::SetErrorMessage(const std::string& TheMessage)	{
			mErrorMessage = TheMessage;
		}

		const std::string& TestCaseResult::GetErrorMessage() const {
			return mErrorMessage;
		}

		void TestCaseResult::SetSetupElapsedTime(double ElapsedTime)
		{
			mSetupElapsedTime = ElapsedTime;
		}


		double TestCaseResult::GetSetupElapsedTime() const
		{
			return mSetupElapsedTime;
		}

		void TestCaseResult::SetRunElapsedTime(double ElapsedTime)
		{
			mRunElapsedTime = ElapsedTime;
		}

		double TestCaseResult::GetRunElapsedTime() const
		{
			return mRunElapsedTime;
		}

		void TestCaseResult::SetTearDownElapsedTime(double ElapsedTime)
		{
			mTearDownElapsedTime = ElapsedTime;
		}

		double TestCaseResult::GetTearDownElapsedTime() const
		{
			return mTearDownElapsedTime;
		}

		void TestCaseResult::SetElapsedTime(double ElapsedTime)
		{
			mElapsedTime = ElapsedTime;
		}

		double TestCaseResult::GetElapsedTime() const
		{
			return mElapsedTime;
		}

		bool TestCaseResult::IsSucceed() const
		{
			return mStatus == Result::Passed;
		}

		bool TestCaseResult::IsFailed() const
		{
			return mStatus == Result::Failed;
		}

		bool TestCaseResult::IsSkipped() const
		{
			return mStatus == Result::Skipped;
		}

		bool TestCaseResult::IsRun() const
		{
			return mStatus != Result::DidNotRun;
		}

		std::string TestCaseResult::Info() const
		{
			if (mStatus == Result::DidNotRun)
			{
				return "Test case not run";
			}
			else if (mStatus == Result::Passed)
			{
				return "Test case successful";
			}
			else if (mStatus == Result::Failed)
			{
				return "Test case failed";
			}
			else if (mStatus == Result::Skipped)
			{
				return "Test case skipped";
			}

			return "Unknown test case state";
		}

		/// Print information about this object.
		void TestCaseResult::PrintInfo(std::ostream& rOStream) const
		{
			rOStream << Info();
		}

		/// Print object's data.
		void TestCaseResult::PrintData(std::ostream& rOStream) const
		{
			rOStream << "Test output: " << mOutput << std::endl;
			if (IsFailed())
				rOStream << "Test Error message : " << mErrorMessage;
		}

	} // manespace Testing.
}  // namespace Kratos.
