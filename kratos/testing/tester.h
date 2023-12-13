//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand
//
//

#pragma once

// System includes
#include <string>
#include <iostream>
#include <map>

// External includes

// Project includes
#include "includes/kratos_export_api.h"

namespace Kratos::Testing
{

///@addtogroup KratosCore
///@{

///@name Kratos Classes
///@{

class TestCase;  // Forward declaration of TestCase class
class TestSuite; // Forward declaration of TestSuite class

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

    using TestCasesContainerType = std::map<std::string, TestCase*>;

    using TestSuitesContainerType = std::map<std::string, TestSuite*>;

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

    static int RunAllDistributedTestCases();

    static int ProfileAllTestCases();

    static int ProfileAllDistributedTestCases();

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
    ///@name Private Operations
    ///@{

    /**
    * @brief Unselects all test cases.
    * @details This function iterates over all test cases in the 'mTestCases' map and calls the 'UnSelect()' method on each of them, effectively unselecting them.
    */
    static void UnSelectAllTestCases();

    /**
    * @brief Selects only the enabled test cases.
    * @details This function iterates over all the test cases in the `mTestCases` map of the `Tester` singleton instance, and selects only the ones that are enabled by calling their `Select()` method. The test cases that are disabled are unselected by calling their `UnSelect()` method.
    */
    static void SelectOnlyEnabledTestCases();

    /**
    * @brief Selects only the distributed test cases.
    * @details This function iterates over all the test cases in the `mTestCases` map of the `Tester` singleton instance, and selects only the ones that are enabled by calling their `Select()` method. The test cases that are disabled are unselected by calling their `UnSelect()` method.
    */
    static void SelectOnlyDistributedTestCases();

    /**
    * @brief Select test cases whose names match a given pattern.
    * @param TestCasesNamePattern a string representing a regex pattern, where * is replaced with ".*"
    */
    static void SelectTestCasesByPattern(std::string const& rTestCasesNamePattern);

    /**
    * @brief Runs all selected test cases and reports the results
    * @return The number of selected test cases
    */
    static int RunSelectedTestCases();

    /**
    * @brief Profiles all selected test cases and reports the results
    * @return The number of selected test cases
    * @todo May be unified with benchmarking
    */
    static int ProfileSelectedTestCases();

    /**
    * @brief Calculates the number of selected test cases.
    * @details This method:
    * - Loops through all the test cases in the mTestCases list of the Tester instance.
    * - If a test case is enabled and selected, increments the result counter.
    * - Returns the final count of selected test cases.
    * @return An integer representing the number of selected test cases.
    */
    static std::size_t NumberOfSelectedTestCases();

    /**
    * @brief Starts the progress show of test cases
    * @param Current: a std::size_t representing the current progress
    * @param Total: a std::size_t representing the total progress
    * @param pTheTestCase: a pointer to a TestCase object representing the current test case
    */
    static void StartShowProgress(
        const std::size_t Current,
        const std::size_t Total,
        const TestCase* pTheTestCase
        );

    /**
    * @brief Ends the progress show of test cases
    * @param Current: a std::size_t representing the current progress
    * @param Total: a std::size_t representing the total progress
    * @param pTheTestCase: a pointer to a TestCase object representing the current test case
    */
    static void EndShowProgress(
        const std::size_t Current,
        const std::size_t Total,
        const TestCase* pTheTestCase
        );

    /**
    * @brief Reports test results to a given output stream.
    * @param rOStream Output stream to report results to.
    * @param NumberOfRunTests Number of test cases run.
    * @param ElapsedTime Time elapsed during test execution.
    * @return An integer representing the exit code of the function.
    * @return 0 if all tests passed, 1 if any tests failed.
    */
    static int ReportResults(
        std::ostream& rOStream,
        const std::size_t NumberOfRunTests,
        const double ElapsedTime
        );

    /**
    * @brief A description of the entire function, its parameters, and its return types.
    * @param rOStream The output stream to write test failure results to.
    */
    static void ReportFailures(std::ostream& rOStream);

    /**
    * Reports detailed failure information in a distributed setting.
    * @param rOStream The output stream to write the failure details to.
    * @param pTheTestCase A pointer to the test case object being reported.
    */
    static void ReportDistributedFailureDetails(
        std::ostream& rOStream,
        const TestCase* pTheTestCase
        );

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
} // namespace Kratos::Testing.
