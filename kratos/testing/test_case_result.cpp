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
		TestCaseResult::TestCaseResult() : mSucceed(true), mOutput(""), mErrorMessage("")
			, mSetupElapsedTime(0.00), mRunElapsedTime(0.00), mTearDownElapsedTime(0.00), mElapsedTime(0.00) {}

		TestCaseResult::TestCaseResult(TestCaseResult const& rOther)
			: mSucceed(rOther.mSucceed), mOutput(rOther.mOutput), mErrorMessage(rOther.mErrorMessage)
			, mSetupElapsedTime(rOther.mSetupElapsedTime), mRunElapsedTime(rOther.mRunElapsedTime)
			, mTearDownElapsedTime(rOther.mTearDownElapsedTime), mElapsedTime(rOther.mElapsedTime) {}

		TestCaseResult::~TestCaseResult() {}

		TestCaseResult& TestCaseResult::operator=(TestCaseResult const& rOther)
		{
			mSucceed = rOther.mSucceed;
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
			mSucceed = true;
			mOutput = "";
			mErrorMessage = "";
			mSetupElapsedTime = 0.00;
			mRunElapsedTime = 0.00;
			mTearDownElapsedTime = 0.00;
			mElapsedTime = 0.00;
		}

		void TestCaseResult::SetToSucceed()	{
			mSucceed = true;
		}

		void TestCaseResult::SetToFailed() {
			mSucceed = false;
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
			return mSucceed;
		}

		bool TestCaseResult::IsFailed() const
		{
			return !mSucceed;
		}

		std::string TestCaseResult::Info() const
		{
			if(mSucceed)
				return "Test case succeed";
			return "Test case failed";
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
