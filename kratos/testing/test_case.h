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


#if !defined(KRATOS_TEST_CASE_H_INCLUDED )
#define  KRATOS_TEST_CASE_H_INCLUDED

#include <string>
#include <iostream>


// External includes


// Project includes
#include "testing/tester.h"
#include "testing/test_case_result.h"

namespace Kratos
{
	namespace Testing
	{
		namespace Internals
		{
			template <typename TestType> class RegisterThisTest {
			public:
				RegisterThisTest(bool IsDisabled=false)
				{
					TestType* p_test = new TestType;
					if (IsDisabled)
						p_test->Disable();
					Tester::AddTestCase(p_test);
				}
			};

		}

		///@addtogroup KratosCore
		///@{

		///@name Kratos Classes
		///@{

		/// The test case base class.
		/** Defines the interface for all the test cases and also fixtures
		*/
		class KRATOS_API(KRATOS_CORE) TestCase
		{
		public:
			///@name Type Definitions
			///@{

			///@}
			///@name Life Cycle
			///@{

			/// TestCase cannot be created without a name
			TestCase() = delete;

			/// The TestCase cannot be copied to avoid duplications
			TestCase(TestCase const& rOther) = delete;

			/// The constructor to be called
			TestCase(std::string const& Name);

			/// Destructor.
			virtual ~TestCase();


			///@}
			///@name Operators
			///@{

			/// Preventing the assignment of the tests
			TestCase& operator=(TestCase const& rOther) = delete;

			///@}
			///@name Operations
			///@{

			virtual void Reset();

			virtual void ResetResult();

			virtual void Setup();

			virtual void Run();

			virtual void Profile();

			virtual void TearDown();

			virtual void Enable();

			virtual void Disable();

			virtual void Select();

			virtual void UnSelect();

			///@}
			///@name Access
			///@{

			const std::string& Name() const;

			const TestCaseResult& GetResult() const;

			void SetResult(TestCaseResult const& TheResult);

			void SetResultOutput(std::string const& TheResultOutput);


			///@}
			///@name Inquiry
			///@{

			bool IsEnabled();
			bool IsDisabled();
			bool IsSelected();


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
			///@name Static Member Variables
			///@{


			///@}
			///@name Member Variables
			///@{

			const std::string mName;

			bool mIsEnabled;

			bool mIsSelected;

			TestCaseResult mResult;

			///@}
			///@name Operations
			///@{

			virtual void TestFunction() = 0;

			///@}


		}; // Class TestCase

	  ///@}

	  ///@name Input and output
	  ///@{

		/// output stream function
		inline std::ostream& operator << (std::ostream& rOStream,
			const TestCase& rThis)
		{
			rThis.PrintInfo(rOStream);
			rOStream << std::endl;
			rThis.PrintData(rOStream);

			return rOStream;
		}
		///@}
		///@name macros
		///@{


#define KRATOS_TESTING_CREATE_CLASS_NAME(TestCaseName) \
  Test##TestCaseName

#define KRATOS_TESTING_CONVERT_TO_STRING(Name) #Name

#define KRATOS_TESTING_TEST_CASE_CLASS_BODY(TestCaseName,ParentName) \
 public:\
  KRATOS_TESTING_CREATE_CLASS_NAME(TestCaseName)() : ParentName(KRATOS_TESTING_CONVERT_TO_STRING(Test##TestCaseName)) {}\
 private: \
  void TestFunction() override; \
  static const Internals::RegisterThisTest< KRATOS_TESTING_CREATE_CLASS_NAME(TestCaseName) > mDummy; 

#define KRATOS_TESTING_TEST_CASE_CLASS(TestCaseName,ParentName) \
class KRATOS_TESTING_CREATE_CLASS_NAME(TestCaseName) : public ParentName \
 {\
KRATOS_TESTING_TEST_CASE_CLASS_BODY(TestCaseName,ParentName) \
};


// This macro creates a sub class of TestCase for given TestCaseName prepending a Test to it
// For example if the TestCaseName is "ModelPartConstruction" then the result is:
// class TestModelPartConstruction : public TestCase
// {
//	public:
//		TestModelPartConstruction() : TestCase("TestModelPartConstruction") {}
//	private:
//		void TestFunction() override;
//		static Internals::RegisterThisTest<TestModelPartConstruction> mDummy;
//	};
//	Kratos::Testing::Internals::RegisterThisTest<TestModelPartConstruction>
//		TestModelPartConstruction::mDummy;
//	void TestModelPartConstruction::TestFunction()
//
#define KRATOS_TEST_CASE(TestCaseName) \
KRATOS_TESTING_TEST_CASE_CLASS(TestCaseName, TestCase) \
const Kratos::Testing::Internals::RegisterThisTest< KRATOS_TESTING_CREATE_CLASS_NAME(TestCaseName) > \
		KRATOS_TESTING_CREATE_CLASS_NAME(TestCaseName)::mDummy; \
\
void KRATOS_TESTING_CREATE_CLASS_NAME(TestCaseName)::TestFunction()

#define KRATOS_DISABLED_TEST_CASE(TestCaseName) \
KRATOS_TESTING_TEST_CASE_CLASS(TestCaseName, TestCase) \
const Kratos::Testing::Internals::RegisterThisTest< KRATOS_TESTING_CREATE_CLASS_NAME(TestCaseName) > \
		KRATOS_TESTING_CREATE_CLASS_NAME(TestCaseName)::mDummy(true); \
\
void KRATOS_TESTING_CREATE_CLASS_NAME(TestCaseName)::TestFunction()

// This macro creates a fixture sub class of TestCase for given FixtureName
// For example if the FixtureName is "ModelPartFixture" then the result is:
// class ModelPartFixture : public TestCase
// {
//	public:
//		ModelPartFixture(std::string const& Name) : TestCase(Name) {}
//	private:
//		void Setup() override;
//		void TearDown() override;
//	};
//
#define KRATOS_TEST_FIXTURE(TestFixtureName) \
class TestFixtureName : public TestCase \
 {\
 public:\
  TestFixtureName(std::string const& Name) : TestCase(Name) {}\
 private: \
  void Setup() override; \
  void TearDown() override; \
}; 

#define KRATOS_TEST_FIXTURE_SETUP(TestFixtureName) \
void TestFixtureName::Setup()

#define KRATOS_TEST_FIXTURE_TEAR_DOWN(TestFixtureName) \
void TestFixtureName::TearDown()

#define KRATOS_TEST_CASE_WITH_FIXTURE(TestCaseName,TestFixtureName) \
KRATOS_TESTING_TEST_CASE_CLASS(TestCaseName, TestFixtureName)  \
const Kratos::Testing::Internals::RegisterThisTest< KRATOS_TESTING_CREATE_CLASS_NAME(TestCaseName) > \
		KRATOS_TESTING_CREATE_CLASS_NAME(TestCaseName)::mDummy; \
\
void KRATOS_TESTING_CREATE_CLASS_NAME(TestCaseName)::TestFunction()

#define KRATOS_DISABLED_TEST_CASE_WITH_FIXTURE(TestCaseName,TestFixtureName) \
KRATOS_TESTING_TEST_CASE_CLASS(TestCaseName, TestFixtureName)  \
const Kratos::Testing::Internals::RegisterThisTest< KRATOS_TESTING_CREATE_CLASS_NAME(TestCaseName) > \
		KRATOS_TESTING_CREATE_CLASS_NAME(TestCaseName)::mDummy(true); \
\
void KRATOS_TESTING_CREATE_CLASS_NAME(TestCaseName)::TestFunction()


///@}

		///@} addtogroup block
	} // manespace Testing.
}  // namespace Kratos.

#endif // KRATOS_TEST_CASE_H_INCLUDED  defined
