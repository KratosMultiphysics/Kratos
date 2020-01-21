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
#include <regex>

// External includes


// Project includes
#include "testing/tester.h"
#include "testing/test_suite.h" // it includes the test_case.h
#include "includes/exception.h"
#include "includes/parallel_environment.h"
#include "includes/data_communicator.h"


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

			for (auto i_test = GetInstance().mTestSuites.begin();
			i_test != GetInstance().mTestSuites.end(); i_test++)
				delete i_test->second;
		}

		void Tester::ResetAllTestCasesResults()
		{
			for (auto i_test = GetInstance().mTestCases.begin();
			i_test != GetInstance().mTestCases.end(); i_test++)
				i_test->second->ResetResult();
		}

		int Tester::RunAllTestCases()
		{
			// TODO: Including the initialization time in the timing.
			ResetAllTestCasesResults();
			SelectOnlyEnabledTestCases();
			return RunSelectedTestCases();
		}

		int Tester::ProfileAllTestCases()
		{
			ResetAllTestCasesResults();
			SelectOnlyEnabledTestCases();
			return ProfileSelectedTestCases();
		}

		int Tester::RunTestSuite(std::string const& TestSuiteName)
		{
			// TODO: Including the initialization time in the timing.
			ResetAllTestCasesResults();
			UnSelectAllTestCases();
			GetTestSuite(TestSuiteName).Select();
			return RunSelectedTestCases();
		}

		int Tester::RunTestCases(std::string const& TestCasesNamePattern)
		{
			ResetAllTestCasesResults();
			UnSelectAllTestCases();
			SelectTestCasesByPattern(TestCasesNamePattern);
			return RunSelectedTestCases();

			//KRATOS_CHECK(std::regex_match(s, std::regex(buffer.str())));

		}

		int Tester::ProfileTestSuite(std::string const& TestSuiteName)
		{
            KRATOS_ERROR << "Profile test suite is not implmented yet" << std::endl;
            return 0;
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

		std::size_t Tester::NumberOfSkippedTestCases()
		{
			std::size_t result = 0;
			for (auto i_test = GetInstance().mTestCases.begin();
			i_test != GetInstance().mTestCases.end(); i_test++)
			{
				TestCaseResult const& test_case_result = i_test->second->GetResult();
				if (test_case_result.IsSkipped())
					result++;
			}

			return result;
		}

		void Tester::AddTestCase(TestCase* pHeapAllocatedTestCase)
		{
			KRATOS_ERROR_IF(HasTestCase(pHeapAllocatedTestCase->Name())) << "A duplicated test case found! The test case \"" << pHeapAllocatedTestCase->Name() << "\" is already added." << std::endl;
			GetInstance().mTestCases[pHeapAllocatedTestCase->Name()] = pHeapAllocatedTestCase;
		}

		TestCase& Tester::GetTestCase(std::string const& TestCaseName) {
			return *pGetTestCase(TestCaseName);
		}

		TestCase* Tester::pGetTestCase(std::string const& TestCaseName) {
			KRATOS_ERROR_IF_NOT(HasTestCase(TestCaseName)) << "The test case \"" << TestCaseName << "\" is not registered in tester" << std::endl;
			return GetInstance().mTestCases[TestCaseName];
		}

		TestSuite* Tester::CreateTestSuite(std::string const& TestSuiteName)
		{
			if (HasTestSuite(TestSuiteName))
				return pGetTestSuite(TestSuiteName);
			TestSuite* p_new_test_suite = new TestSuite(TestSuiteName);
			AddTestSuite(p_new_test_suite);
			return p_new_test_suite;
		}

		TestSuite* Tester::CreateNewTestSuite(std::string const& TestSuiteName)
		{
			TestSuite* p_new_test_suite = new TestSuite(TestSuiteName);
			AddTestSuite(p_new_test_suite);
			return p_new_test_suite;
		}

		void Tester::AddTestSuite(TestSuite* pHeapAllocatedTestSuite)
		{
			KRATOS_ERROR_IF(HasTestSuite(pHeapAllocatedTestSuite->Name())) << "A duplicated test suite found! The test suite \"" << pHeapAllocatedTestSuite->Name() << "\" is already added." << std::endl;
			GetInstance().mTestSuites[pHeapAllocatedTestSuite->Name()] = pHeapAllocatedTestSuite;
		}

		TestSuite& Tester::GetTestSuite(std::string const& TestSuiteName) {
			return *pGetTestSuite(TestSuiteName);
		}

		TestSuite* Tester::pGetTestSuite(std::string const& TestSuiteName) {
			KRATOS_ERROR_IF_NOT(HasTestSuite(TestSuiteName)) << "The test suite \"" << TestSuiteName << "\" is not registered in tester" << std::endl;
			return GetInstance().mTestSuites[TestSuiteName];
		}

		void Tester::AddTestToTestSuite(std::string const& TestName, std::string const& TestSuiteName)
		{
			TestCase* p_test_case = nullptr;
			if (HasTestCase(TestName))
				p_test_case = pGetTestCase(TestName);
			else if (HasTestSuite(TestName))
				p_test_case = pGetTestSuite(TestName);
			else
				KRATOS_ERROR << "The test \"" << TestName << "\" is not registered in the tester" << std::endl;

			TestSuite* p_test_suite = CreateTestSuite(TestSuiteName);
			p_test_suite->AddTestCase(p_test_case);
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

		bool Tester::HasTestCase(std::string const& TestCaseName) {
			return (GetInstance().mTestCases.find(TestCaseName) != GetInstance().mTestCases.end());
		}

		bool Tester::HasTestSuite(std::string const& TestSuiteName) {
			return (GetInstance().mTestSuites.find(TestSuiteName) != GetInstance().mTestSuites.end());
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
			for (auto i_test = GetInstance().mTestCases.begin();
			i_test != GetInstance().mTestCases.end(); i_test++)
				rOStream << "        " << i_test->first << std::endl;
			rOStream << "    Total of " << GetInstance().mTestCases.size() << " Test cases" << std::endl;
			rOStream << "    Test suites:" << std::endl;
			for (auto i_test = GetInstance().mTestSuites.begin();
			i_test != GetInstance().mTestSuites.end(); i_test++)
				rOStream << "        " << i_test->first << std::endl;
			rOStream << "    Total of " << GetInstance().mTestSuites.size() << " Test Suites";
		}

		void Tester::UnSelectAllTestCases()
		{
			for (auto i_test = GetInstance().mTestCases.begin();
			i_test != GetInstance().mTestCases.end(); i_test++)
				i_test->second->UnSelect();
		}

		void Tester::SelectOnlyEnabledTestCases()
		{
			for (auto i_test = GetInstance().mTestCases.begin();
			i_test != GetInstance().mTestCases.end(); i_test++)
				if (i_test->second->IsEnabled())
					i_test->second->Select();
				else
					i_test->second->UnSelect();
		}

		void Tester::SelectTestCasesByPattern(std::string const& TestCasesNamePattern)
		{
			// creating the regex pattern replacing * with ".*"
#if defined(__GNUC__) && !defined(__clang__) && (__GNUC__ < 4 || (__GNUC__ == 4 && (__GNUC_MINOR__ < 9)))
			KRATOS_ERROR << "This method is not compiled well. You should use a GCC 4.9 or higher" << std::endl;
#else
			std::regex replace_star("\\*");
			std::stringstream regex_pattern_string;
			std::regex_replace(std::ostreambuf_iterator<char>(regex_pattern_string),
				TestCasesNamePattern.begin(), TestCasesNamePattern.end(), replace_star, ".*");
			for (auto i_test = GetInstance().mTestCases.begin();
				i_test != GetInstance().mTestCases.end(); i_test++)
				if (std::regex_match(i_test->second->Name(), std::regex(regex_pattern_string.str()))) {
					if (i_test->second->IsEnabled()) {
						i_test->second->Select();
                    } else {
                        i_test->second->UnSelect();
                    }
                }
#endif
		}

		int Tester::RunSelectedTestCases()
		{
			auto start = std::chrono::steady_clock::now();
			auto number_of_run_tests = NumberOfSelectedTestCases();

			std::stringstream secondary_stream;
			std::streambuf* original_buffer = nullptr;
			if (ParallelEnvironment::GetDefaultRank() != 0) {
				original_buffer = std::cout.rdbuf(secondary_stream.rdbuf());
			}

			std::size_t test_number = 0;
			for (auto i_test = GetInstance().mTestCases.begin();
			i_test != GetInstance().mTestCases.end(); i_test++)
			{
				if (i_test->second->IsSelected()) {
					StartShowProgress(test_number, number_of_run_tests, i_test->second);
					if (GetInstance().mVerbosity != Verbosity::TESTS_OUTPUTS) {
						std::stringstream output_stream;
						auto old_buffer(std::cout.rdbuf(output_stream.rdbuf()));
						i_test->second->Run();
						i_test->second->SetResultOutput(output_stream.str());
						std::cout.rdbuf(old_buffer);
					}
					else
						i_test->second->Run();
					EndShowProgress(++test_number, number_of_run_tests, i_test->second);
				}
			}

			auto end = std::chrono::steady_clock::now();
			std::chrono::duration<double> elapsed = end - start;

			auto tmp = ReportResults(std::cout, number_of_run_tests, elapsed.count());

			if (ParallelEnvironment::GetDefaultRank() != 0) {
				std::cout.rdbuf(original_buffer);
			}
			return tmp;
		}



		int Tester::ProfileSelectedTestCases()
		{
			KRATOS_ERROR << "Profile test cases is not implemented yet" << std::endl;
			auto start = std::chrono::steady_clock::now();
			auto number_of_run_tests = NumberOfSelectedTestCases();

			std::size_t test_number = 0;
			for (auto i_test = GetInstance().mTestCases.begin();
			i_test != GetInstance().mTestCases.end(); i_test++)
			{
				if (i_test->second->IsSelected()) {
					StartShowProgress(test_number, number_of_run_tests, i_test->second);
					i_test->second->Profile();
					EndShowProgress(++test_number, number_of_run_tests, i_test->second);
				}
			}

			auto end = std::chrono::steady_clock::now();
			std::chrono::duration<double> elapsed = end - start;

			return ReportResults(std::cout, number_of_run_tests, elapsed.count());
		}

		std::size_t Tester::NumberOfSelectedTestCases()
		{
			std::size_t result = 0;
			for (auto i_test = GetInstance().mTestCases.begin();
			i_test != GetInstance().mTestCases.end(); i_test++)
				if (i_test->second->IsEnabled())
					if (i_test->second->IsSelected())
						result++;

			return result;
		}

		void Tester::StartShowProgress(std::size_t Current, std::size_t Total, const TestCase* const pTheTestCase)
		{
			if (GetInstance().mVerbosity >= Verbosity::TESTS_LIST)
			{
				std::cout << pTheTestCase->Name();
			}
		}

		void Tester::EndShowProgress(std::size_t Current, std::size_t Total, const TestCase* const pTheTestCase)
		{
			constexpr std::size_t ok_culumn = 72;
			if (GetInstance().mVerbosity == Verbosity::PROGRESS) {
				if (pTheTestCase->GetResult().IsSucceed())
					std::cout << ".";
				else if (pTheTestCase->GetResult().IsFailed())
					std::cout << "F";
				else if (pTheTestCase->GetResult().IsSkipped())
					std::cout << "s";
			}
			else if (GetInstance().mVerbosity >= Verbosity::TESTS_LIST)
			{
				for (std::size_t i = pTheTestCase->Name().size(); i < ok_culumn; i++)
					std::cout << " ";

				if (pTheTestCase->GetResult().IsSucceed())
				{
					std::cout << "OK." << std::endl;
					if (GetInstance().mVerbosity == Verbosity::TESTS_OUTPUTS)
						std::cout << pTheTestCase->GetResult().GetOutput() << std::endl;
				}
				else if (pTheTestCase->GetResult().IsFailed())
				{
					std::cout << "FAILED!" << std::endl;
					if (GetInstance().mVerbosity >= Verbosity::FAILED_TESTS_OUTPUTS)
						std::cout << pTheTestCase->GetResult().GetOutput() << std::endl;
				}
				else if (pTheTestCase->GetResult().IsSkipped())
				{
					std::cout << "SKIPPED." << std::endl;
					if (GetInstance().mVerbosity == Verbosity::TESTS_OUTPUTS)
						std::cout << pTheTestCase->GetResult().GetErrorMessage() << std::endl;
				}
			}
		}

		int Tester::ReportResults(std::ostream& rOStream, std::size_t NumberOfRunTests,double ElapsedTime)
		{
            int exit_code = 0;

			if (GetInstance().mVerbosity == Verbosity::PROGRESS)
				rOStream << std::endl;

			auto number_of_failed_tests = NumberOfFailedTestCases();
			auto number_of_skipped_tests = NumberOfSkippedTestCases();

			std::string total_test_cases = " test case";
			auto total_test_cases_size = GetInstance().mTestCases.size();
			if (total_test_cases_size > 1)
				total_test_cases += "s";

			rOStream << "Ran " << NumberOfRunTests << " of " << total_test_cases_size << total_test_cases << " in " << ElapsedTime << "s";
			if (number_of_skipped_tests > 0) {
				rOStream << " (" << number_of_skipped_tests << " skipped)";
			}

			if (number_of_failed_tests == 0) {
				rOStream << ". OK" << std::endl;
			}
			else
			{
				rOStream << ". " << number_of_failed_tests << " failed:" << std::endl;
				ReportFailures(rOStream);
                exit_code = 1;
			}

            return exit_code;
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
					if (ParallelEnvironment::GetDefaultSize() == 1)
					{
						if (test_case_result.GetErrorMessage().size() == 0)
							rOStream << std::endl;
						else
						{
							rOStream << " with message: " << std::endl;
							rOStream << "        " << test_case_result.GetErrorMessage() << std::endl;
						}
					}
					else
					{
						Tester::ReportDistributedFailureDetails(rOStream, i_test->second);
					}
				}
			}
		}

		void Tester::ReportDistributedFailureDetails(std::ostream& rOStream, const TestCase* const pTheTestCase)
		{
			TestCaseResult const& r_test_case_result = pTheTestCase->GetResult();
			rOStream << " with messages: " << std::endl;
			rOStream << "From rank 0:" << std::endl << r_test_case_result.GetErrorMessage() << std::endl;
			const DataCommunicator& r_comm = ParallelEnvironment::GetDefaultDataCommunicator();
			const int parallel_rank = r_comm.Rank();
			const int parallel_size = r_comm.Size();
			if (parallel_rank == 0)
			{
				for (int i = 1; i < parallel_size; i++)
				{
					std::string remote_message;
					r_comm.Recv(remote_message, i,i);
					rOStream << "From rank " << i << ":" << std::endl << remote_message << std::endl;
				}
			}
			else
			{
				r_comm.Send(r_test_case_result.GetErrorMessage(), 0, parallel_rank);
			}
		}


	} // manespace Testing.
}  // namespace Kratos.
