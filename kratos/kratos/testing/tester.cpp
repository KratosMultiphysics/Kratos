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

// External includes


// Project includes
#include "testing/tester.h"
#include "testing/test_case.h"
#include "includes/checks.h"
//#include "utilities/logger.h"


namespace Kratos
{
	namespace Testing
	{
		Tester::Tester() : mTestCases(), mVerbosity(Verbosity::QUITE)
			 {}

		Tester::~Tester()
		{
			for (auto i_test = GetInstance().mTestCases.begin();
			i_test != GetInstance().mTestCases.end(); i_test++)
				delete i_test->second;
		}

		void Tester::ResetAllTestCasesResults()
		{
			for (auto i_test = GetInstance().mTestCases.begin();
			i_test != GetInstance().mTestCases.end(); i_test++)
				i_test->second->ResetResult();
		}

		void Tester::RunAllTestCases()
		{
			auto start = std::chrono::steady_clock::now();
			ResetAllTestCasesResults();
			auto number_of_run_tests = NumberOfEnabledTestCases();

			std::size_t test_number = 0;
			for (auto i_test = GetInstance().mTestCases.begin();
			i_test != GetInstance().mTestCases.end(); i_test++)
			{
				i_test->second->Run();
				ShowProgress(++test_number, number_of_run_tests, i_test->second);
			}


			auto end = std::chrono::steady_clock::now();
			std::chrono::duration<double> elapsed = end - start;

			ReportResults(std::cout, number_of_run_tests, elapsed.count());

		}

		void Tester::ProfileAllTestCases()
		{
			auto start = std::chrono::steady_clock::now();
			ResetAllTestCasesResults();
			auto number_of_run_tests = NumberOfEnabledTestCases();

			for (auto i_test = GetInstance().mTestCases.begin();
			i_test != GetInstance().mTestCases.end(); i_test++)
				i_test->second->Profile();

			auto end = std::chrono::steady_clock::now();
			std::chrono::duration<double> elapsed = end - start;

			ReportResults(std::cout, number_of_run_tests, elapsed.count());
		}

		std::size_t Tester::NumberOfEnabledTestCases()
		{
			std::size_t result = 0;
			for (TestCasesContainerType::iterator i_test = GetInstance().mTestCases.begin();
			i_test != GetInstance().mTestCases.end(); i_test++)
				if (i_test->second->IsEnabled())
					result++;

			return result;
		}

		std::size_t Tester::NumberOfFailedTestCases()
		{
			std::size_t result = 0;
			for (auto i_test = GetInstance().mTestCases.begin();
			i_test != GetInstance().mTestCases.end(); i_test++)
			{
				TestCaseResult const& test_case_result = i_test->second->GetResult();
				if (test_case_result.IsFailed())
					result++;
			}

			return result;
		}

		void Tester::AddTestCase(TestCase* pHeapAllocatedTestCase)
		{
			KRATOS_CHECK(IsTestCaseNotAddedBefore(pHeapAllocatedTestCase)) << "A duplicated test case found! The test case \"" << pHeapAllocatedTestCase->Name() << "\" is already added." << std::endl;
			GetInstance().mTestCases[pHeapAllocatedTestCase->Name()] = pHeapAllocatedTestCase;
		}

		Tester& Tester::GetInstance()
		{
			static Tester instance;
			return instance;
		}

		void Tester::SetVerbosity(Verbosity TheVerbosity)
		{
			GetInstance().mVerbosity = TheVerbosity;
		}

		std::string Tester::Info() const
		{
			return "Tester" ;
		}

		void Tester::PrintInfo(std::ostream& rOStream) const
		{
			rOStream << Info();
		}

		void Tester::PrintData(std::ostream& rOStream) const
		{
			rOStream << "    Test cases:" << std::endl;
			for (TestCasesContainerType::iterator i_test = GetInstance().mTestCases.begin();
			i_test != GetInstance().mTestCases.end(); i_test++)
				rOStream << "        " << i_test->first << std::endl;
			rOStream << "    Total of " << GetInstance().mTestCases.size() << " Test cases";
		}

		bool Tester::IsTestCaseNotAddedBefore(TestCase* pHeapAllocatedTestCase)
		{
			return (GetInstance().mTestCases.find(pHeapAllocatedTestCase->Name()) == GetInstance().mTestCases.end());
		}

		void Tester::ShowProgress(std::size_t Current, std::size_t Total, const TestCase* const pTheTestCase)
		{
			if (GetInstance().mVerbosity == Verbosity::PROGRESS)
				if (pTheTestCase->GetResult().IsSucceed())
					std::cout << ".";
				else
					std::cout << "F";
			else if (GetInstance().mVerbosity >= Verbosity::TESTS_LIST)
			{
				if (pTheTestCase->GetResult().IsSucceed())
				{
					std::cout << "OK.     : "<< pTheTestCase->Name() << std::endl;
					if (GetInstance().mVerbosity == Verbosity::TESTS_OUTPUTS)
						std::cout << pTheTestCase->GetResult().GetOutput() << std::endl;
				}
				else
				{
					std::cout << "FAILED! : "<< pTheTestCase->Name()  << std::endl;
					if(GetInstance().mVerbosity >= Verbosity::FAILED_TESTS_OUTPUTS)
						std::cout << pTheTestCase->GetResult().GetOutput() << std::endl;
				}
			}
		}

		void Tester::ReportResults(std::ostream& rOStream, std::size_t NumberOfRunTests,double ElapsedTime)
		{
			if (GetInstance().mVerbosity == Verbosity::PROGRESS)
				rOStream << std::endl;
			std::size_t number_of_failed_tests = NumberOfFailedTestCases();
			std::string test_cases = " Test Case";
			if (NumberOfRunTests > 1)
				test_cases += "s";

			if (number_of_failed_tests == 0)
				rOStream << NumberOfRunTests << test_cases << " run in " << ElapsedTime << "s. OK." << std::endl;
			else
			{
				rOStream << NumberOfRunTests << test_cases << " run in " << ElapsedTime << "s. " << number_of_failed_tests << " failed:" << std::endl;
				ReportFailures(rOStream);
			}

		}


		void Tester::ReportFailures(std::ostream& rOStream)
		{
			for (auto i_test = GetInstance().mTestCases.begin();
			i_test != GetInstance().mTestCases.end(); i_test++)
			{
				TestCaseResult const& test_case_result = i_test->second->GetResult();
				if (test_case_result.IsFailed())
				{
					rOStream << "    " << i_test->first << " Failed";
					if (test_case_result.GetErrorMessage().size() == 0)
						rOStream << std::endl;
					else
					{
						rOStream << " with message: " << std::endl;
						rOStream << "        " << test_case_result.GetErrorMessage() << std::endl;
					}
				}
			}
		}


	} // manespace Testing.
}  // namespace Kratos.
