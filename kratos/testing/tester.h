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
		class TestSuite;

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

			typedef std::map<std::string, TestSuite*> TestSuitesContainerType;

			///@}
			///@name Enums
			///@{

			enum class Verbosity {QUITE, PROGRESS, TESTS_LIST, FAILED_TESTS_OUTPUTS, TESTS_OUTPUTS};

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

			static int RunAllTestCases();

			static int ProfileAllTestCases();

			static int RunTestSuite(std::string const& TestSuiteName);

			static int ProfileTestSuite(std::string const& TestSuiteName);

			/// The test case pattern can apply * as any number of any character
			static int RunTestCases(std::string const& TestCasesNamePattern);

			/// The test case pattern can apply * as any number of any character
			static int ProfileTestCases(std::string const& TestCasesNamePattern);

			static std::size_t NumberOfFailedTestCases();

			static std::size_t NumberOfSkippedTestCases();

			/// This method assumes that the given test case is allocated
			/// via new. So it will delete it at the end of the program
			static void AddTestCase(TestCase* pHeapAllocatedTestCase);

			static TestCase& GetTestCase(std::string const& TestCaseName);

			static TestCase* pGetTestCase(std::string const& TestCaseName);

			/// Creates a new test suite or return the already created one.
			static TestSuite* CreateTestSuite(std::string const& TestSuiteName);

			/// This method gives an error if test suite already exists.
			static TestSuite* CreateNewTestSuite(std::string const& TestSuiteName);

			/// This method assumes that the given test suite is allocated
			/// via new. So it will delete it at the end of the program
			static void AddTestSuite(TestSuite* pHeapAllocatedTestSuite);

			static TestSuite& GetTestSuite(std::string const& TestSuiteName);

			static TestSuite* pGetTestSuite(std::string const& TestSuiteName);

			/// This method creates the suite if is no exist
			static void AddTestToTestSuite(std::string const& TestName, std::string const& TestSuiteName);

			///@}
			///@name Access
			///@{

			static Tester& GetInstance();

			static void SetVerbosity(Verbosity TheVerbosity);



			///@}
			///@name Inquiry
			///@{

			static bool HasTestCase(std::string const& TestCaseName);

			static bool HasTestSuite(std::string const& TestSuiteName);

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
			///@name Member Variables
			///@{

			TestCasesContainerType mTestCases;

			TestSuitesContainerType mTestSuites;

			Verbosity mVerbosity;

			///@}

			static void UnSelectAllTestCases();

			static void SelectOnlyEnabledTestCases();

			static void SelectTestCasesByPattern(std::string const& TestCasesNamePattern);

			static int RunSelectedTestCases();

			static int ProfileSelectedTestCases();

			static std::size_t NumberOfSelectedTestCases();

			static void StartShowProgress(std::size_t Current, std::size_t Total, const TestCase* const pTheTestCase);

			static void EndShowProgress(std::size_t Current, std::size_t Total, const TestCase* const pTheTestCase);

			static int ReportResults(std::ostream& rOStream, std::size_t NumberOfRunTests, double ElapsedTime);

			static void ReportFailures(std::ostream& rOStream);

			static void ReportDistributedFailureDetails(std::ostream& rOStream, const TestCase* const pTheTestCase);


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
