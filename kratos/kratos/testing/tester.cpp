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
#include "testing/tester.h"
#include "testing/test_case.h"


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

		void Tester::RunAllTests()
		{
			for (TestCasesContainerType::iterator i_test = GetInstance().mTestCases.begin(); 
			i_test != GetInstance().mTestCases.end(); i_test++)
				i_test->second->Run();
		}

		void Tester::ProfileAllTests()
		{
			for (TestCasesContainerType::iterator i_test = GetInstance().mTestCases.begin(); 
			i_test != GetInstance().mTestCases.end(); i_test++)
				i_test->second->Profile();
		}

		void Tester::AddTestCase(TestCase* pHeapAllocatedTestCase)
		{
		// TODO: Check if test cases already added.
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




	} // manespace Testing.
}  // namespace Kratos.


