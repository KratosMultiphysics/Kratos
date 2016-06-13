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
#include "testing/tester.h"
#include "testing/test_case.h"
//#include "includes/checks.h"
//#include "utilities/logger.h"


namespace Kratos
{
	namespace Testing
	{
		Tester::Tester() : mTestCases()
			 {}

		Tester::~Tester()
		{
			for (TestCasesContainerType::iterator i_test = GetInstance().mTestCases.begin();
			i_test != GetInstance().mTestCases.end(); i_test++)
				delete i_test->second;
		}

		void Tester::ResetAllTestCasesResults()
		{
			for (TestCasesContainerType::iterator i_test = GetInstance().mTestCases.begin();
			i_test != GetInstance().mTestCases.end(); i_test++)
				i_test->second->ResetResult();
		}

		void Tester::RunAllTestCases()
		{
			ResetAllTestCasesResults();
			std::size_t number_of_run_tests = NumberOfEnabledTestCases();

			for (TestCasesContainerType::iterator i_test = GetInstance().mTestCases.begin();
			i_test != GetInstance().mTestCases.end(); i_test++)
				i_test->second->Run();

			std::size_t number_of_failed_tests = NumberOfFailedTestCases();
			if (number_of_failed_tests == 0)
				std::cout << number_of_run_tests << " Test cases run. OK.";
			else
			{
				std::cout << number_of_run_tests << " Test cases run. " << number_of_failed_tests << " failed:" << std::endl;
				ReportFailures(std::cout);
			}

		}

		void Tester::ProfileAllTestCases()
		{
			ResetAllTestCasesResults();
			std::size_t number_of_run_tests = NumberOfEnabledTestCases();

			for (TestCasesContainerType::iterator i_test = GetInstance().mTestCases.begin();
			i_test != GetInstance().mTestCases.end(); i_test++)
				i_test->second->Profile();

			std::size_t number_of_failed_tests = NumberOfFailedTestCases();
			if (number_of_failed_tests == 0)
				std::cout << number_of_run_tests << " Test cases run. OK.";
			else
			{
				std::cout << number_of_run_tests << " Test cases run. " << number_of_failed_tests << " failed:" << std::endl;
				ReportFailures(std::cout);
			}
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
			for (TestCasesContainerType::iterator i_test = GetInstance().mTestCases.begin();
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
			//KRATOS_CHECK(IsTestCaseNotAddedBefore(pHeapAllocatedTestCase)) << "A duplicated test case found! The test case \"" << pHeapAllocatedTestCase->Name() << "\" is already added." << std::endl;
			GetInstance().mTestCases[pHeapAllocatedTestCase->Name()] = pHeapAllocatedTestCase;
		}

		Tester& Tester::GetInstance()
		{
			static Tester instance;
			return instance;
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

		void Tester::ReportFailures(std::ostream& rOStream)
		{
			for (TestCasesContainerType::iterator i_test = GetInstance().mTestCases.begin();
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
