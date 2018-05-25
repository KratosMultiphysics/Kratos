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

#if !defined(KRATOS_TEST_SUITE_H_INCLUDED)
#define KRATOS_TEST_SUITE_H_INCLUDED

// System includes
#include <vector>

// Project includes
#include "testing/test_case.h"

namespace Kratos {
namespace Testing {
namespace Internals {
class AddThisTestToTestSuite {
   public:
    AddThisTestToTestSuite(
        std::string const& TestName, std::string const& TestSuiteName) {
        Tester::AddTestToTestSuite(TestName, TestSuiteName);
    }
};
}

///@addtogroup KratosCore
///@{

///@name Kratos Classes
///@{

/// This class holds an array of test cases and run them one by one in its Run method
/** this class implements a composite pattern. Derived from TestCase and has an array pointers to the 
			TestCase to be run. The Run and Profile methods are overridden to call the corresponidng methods 
			of the TestCases.
		*/
class KRATOS_API(KRATOS_CORE) TestSuite : public TestCase {
   public:
    ///@name Type Definitions
    ///@{

    typedef std::vector<TestCase*> TestCasesContainerType;

    ///@}
    ///@name Life Cycle
    ///@{

    /// Default constructor is deleted becuase the Name should always passed to the constructor.
    TestSuite() = delete;

    /// The TestSuite cannot be copied to avoid duplications
    TestSuite(TestSuite const& rOther) = delete;

    /// The constructor to be called
    TestSuite(std::string const& Name);

    /// Destructor.
    virtual ~TestSuite();

    ///@}
    ///@name Operators
    ///@{

    /// Preventing the assignment of the test suites
    TestSuite& operator=(TestSuite const& rOther) = delete;

    ///@}
    ///@name Operations
    ///@{

    void AddTestCase(TestCase* pTestCase);

    void Reset() override;

    void ResetResult() override;

    void Run() override;

    void Profile() override;

    void Enable() override;

    void Disable() override;

    void Select() override;

    void UnSelect() override;

    ///@}
    ///@name Access
    ///@{

    ///@}
    ///@name Inquiry
    ///@{

    ///@}
    ///@name Input and output
    ///@{

    /// Turn back information as a string.
    virtual std::string Info() const override;

    /// Print information about this object.
    virtual void PrintInfo(std::ostream& rOStream) const override;

    /// Print object's data.
    virtual void PrintData(std::ostream& rOStream) const override;

    ///@}
   private:
    ///@name Member Variables
    ///@{

    TestCasesContainerType mTestCases;

    ///@}
    ///@name Private Operations
    ///@{

    void TestFunction() override;

    ///@}

};  // Class TestSuite

///@}

///@name Input and output
///@{

/// input stream function
inline std::istream& operator>>(std::istream& rIStream, TestSuite& rThis);

/// output stream function
inline std::ostream& operator<<(
    std::ostream& rOStream, const TestSuite& rThis) {
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);

    return rOStream;
}
///@}
///@name macros
///@{

#define KRATOS_TESTING_CONCATENATE(a, b) a##b

#define KRATOS_TESTING_CREATE_DUMMY_VARIABLE_NAME(prefix, UniqueNumber) \
    KRATOS_TESTING_CONCATENATE(prefix, UniqueNumber)

#define KRATOS_TESTING_ADD_TEST_TO_TEST_SUITE(TestName, TestSuiteName) \
    Kratos::Testing::Internals::AddThisTestToTestSuite                 \
        KRATOS_TESTING_CREATE_DUMMY_VARIABLE_NAME(dummy, __LINE__)(    \
            TestName, TestSuiteName)

#define KRATOS_TESTING_TEST_CASE_IN_SUITE_CLASS(TestCaseName, ParentName)      \
    class KRATOS_TESTING_CREATE_CLASS_NAME(TestCaseName) : public ParentName { \
        KRATOS_TESTING_TEST_CASE_CLASS_BODY(TestCaseName, ParentName)          \
        static const Kratos::Testing::Internals::AddThisTestToTestSuite        \
            mAnotherDummy;                                                     \
    };

// This is the macro to use
#define KRATOS_TEST_CASE_IN_SUITE(TestCaseName, TestSuiteName)          \
    KRATOS_TESTING_TEST_CASE_IN_SUITE_CLASS(TestCaseName, TestCase)     \
    const Kratos::Testing::Internals::RegisterThisTest<                 \
        KRATOS_TESTING_CREATE_CLASS_NAME(TestCaseName)>                 \
        KRATOS_TESTING_CREATE_CLASS_NAME(TestCaseName)::mDummy;         \
    const Kratos::Testing::Internals::AddThisTestToTestSuite            \
        KRATOS_TESTING_CREATE_CLASS_NAME(TestCaseName)::mAnotherDummy = \
            Kratos::Testing::Internals::AddThisTestToTestSuite(         \
                KRATOS_TESTING_CONVERT_TO_STRING(Test##TestCaseName),   \
                KRATOS_TESTING_CONVERT_TO_STRING(TestSuiteName));       \
                                                                        \
    void KRATOS_TESTING_CREATE_CLASS_NAME(TestCaseName)::TestFunction()

// Disabled version of the macro to used
#define KRATOS_DISABLED_TEST_CASE_IN_SUITE(TestCaseName, TestSuiteName) \
    KRATOS_TESTING_TEST_CASE_IN_SUITE_CLASS(TestCaseName, TestCase)     \
    const Kratos::Testing::Internals::RegisterThisTest<                 \
        KRATOS_TESTING_CREATE_CLASS_NAME(TestCaseName)>                 \
        KRATOS_TESTING_CREATE_CLASS_NAME(TestCaseName)::mDummy(true);   \
    const Kratos::Testing::Internals::AddThisTestToTestSuite            \
        KRATOS_TESTING_CREATE_CLASS_NAME(TestCaseName)::mAnotherDummy = \
            Kratos::Testing::Internals::AddThisTestToTestSuite(         \
                KRATOS_TESTING_CONVERT_TO_STRING(Test##TestCaseName),   \
                KRATOS_TESTING_CONVERT_TO_STRING(TestSuiteName));       \
                                                                        \
    void KRATOS_TESTING_CREATE_CLASS_NAME(TestCaseName)::TestFunction()

#define KRATOS_TEST_CASE_WITH_FIXTURE_IN_SUITE(                            \
    TestCaseName, TestFixtureName, TestSuiteName)                          \
    KRATOS_TESTING_TEST_CASE_IN_SUITE_CLASS(TestCaseName, TestFixtureName) \
    const Kratos::Testing::Internals::RegisterThisTest<                    \
        KRATOS_TESTING_CREATE_CLASS_NAME(TestCaseName)>                    \
        KRATOS_TESTING_CREATE_CLASS_NAME(TestCaseName)::mDummy;            \
    const Kratos::Testing::Internals::AddThisTestToTestSuite               \
        KRATOS_TESTING_CREATE_CLASS_NAME(TestCaseName)::mAnotherDummy =    \
            Kratos::Testing::Internals::AddThisTestToTestSuite(            \
                KRATOS_TESTING_CONVERT_TO_STRING(Test##TestCaseName),      \
                KRATOS_TESTING_CONVERT_TO_STRING(TestSuiteName));          \
                                                                           \
    void KRATOS_TESTING_CREATE_CLASS_NAME(TestCaseName)::TestFunction()

#define KRATOS_DISABLED_TEST_CASE_WITH_FIXTURE_IN_SUITE(                   \
    TestCaseName, TestFixtureName, TestSuiteName)                          \
    KRATOS_TESTING_TEST_CASE_IN_SUITE_CLASS(TestCaseName, TestFixtureName) \
    const Kratos::Testing::Internals::RegisterThisTest<                    \
        KRATOS_TESTING_CREATE_CLASS_NAME(TestCaseName)>                    \
        KRATOS_TESTING_CREATE_CLASS_NAME(TestCaseName)::mDummy(true);      \
    const Kratos::Testing::Internals::AddThisTestToTestSuite               \
        KRATOS_TESTING_CREATE_CLASS_NAME(TestCaseName)::mAnotherDummy =    \
            Kratos::Testing::Internals::AddThisTestToTestSuite(            \
                KRATOS_TESTING_CONVERT_TO_STRING(Test##TestCaseName),      \
                KRATOS_TESTING_CONVERT_TO_STRING(TestSuiteName));          \
                                                                           \
    void KRATOS_TESTING_CREATE_CLASS_NAME(TestCaseName)::TestFunction()

///@}

///@} addtogroup block
}  // namespace Testing
}  // namespace Kratos.

#endif  // KRATOS_TEST_SUITE_H_INCLUDED  defined
