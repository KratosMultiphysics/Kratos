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


#if !defined(KRATOS_TESTER_H_INCLUDED )
#define  KRATOS_TESTER_H_INCLUDED



// System includes
#include <string>
#include <iostream>
#include <map>


// External includes


// Project includes
#include "includes/kratos_export_api.h"

namespace Kratos
{
	namespace Testing
	{

		///@addtogroup KratosCore
		///@{

		///@name Kratos Classes
		///@{

		class TestCase;

		/// Tester class manages all tests and provide interface to run them.
		/** Tester is a singletone class which registers the test cases and
		    test suits and gives interface to register test cases and run all
			or some of them.
		*/
		class KRATOS_API(KRATOS_CORE) Tester
		{
		public:
			///@name Type Definitions
			///@{

			typedef std::map<std::string, TestCase*> TestCasesContainerType;

			///@}
			///@name Life Cycle
			///@{

			/// The Tester cannot be copied to avoid duplications
			Tester(Tester const& rOther) = delete;

			/// Destructor.
			virtual ~Tester();


			///@}
			///@name Operators
			///@{

			/// Preventing the assignment of the tests
			Tester& operator=(Tester const& rOther) = delete;

			///@}
			///@name Operations
			///@{

			static void ResetAllTestCasesResults();

			static void RunAllTestCases();

			static void ProfileAllTestCases();

			static std::size_t NumberOfEnabledTestCases();

			static std::size_t NumberOfFailedTestCases();

			/// This method assumes that the given test case is allocated
			/// via new. So it will delete it at the end of the program
			static void AddTestCase(TestCase* pHeapAllocatedTestCase);

			///@}
			///@name Access
			///@{

			static Tester& GetInstance();




			///@}
			///@name Inquiry
			///@{


			///@}
			///@name Input and output
			///@{

			/// Turn back information as a string.
			virtual std::string Info() const;

			/// Print information about this object.
			virtual void PrintInfo(std::ostream& rOStream) const;

			/// Print object's data.
			virtual void PrintData(std::ostream& rOStream) const;


			///@}

		private:
			///@name Life Cycle
			///@{

			/// Tester cannot be created from outside. To ensure that the one created by instance is the only one.
			Tester();

			///@}
			///@name Static Member Variables
			///@{

			TestCasesContainerType mTestCases;

			///@}
			static bool IsTestCaseNotAddedBefore(TestCase* pHeapAllocatedTestCase);

			static void ReportFailures(std::ostream& rOStream);

		}; // Class Tester

	  ///@}

	  ///@name Input and output
	  ///@{

		/// output stream function
		inline std::ostream& operator << (std::ostream& rOStream,
			const Tester& rThis)
		{
			rThis.PrintInfo(rOStream);
			rOStream << std::endl;
			rThis.PrintData(rOStream);

			return rOStream;
		}
		///@}
		///@name macros
		///@{


		///@}

		///@} addtogroup block
	} // manespace Testing.
}  // namespace Kratos.

#endif // KRATOS_TESTER_H_INCLUDED  defined
