//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand, Carlos A. Roig
//
//

// External includes

// Project includes
#include "includes/define_python.h"
#include "testing/testing.h"
#include "add_testing_to_python.h"


namespace Kratos
{

namespace Python
{

void ListOfAllTestCases() {
	std::cout << Testing::Tester::GetInstance() << std::endl;
}

void  AddTestingToPython(pybind11::module& m) {
	namespace py = pybind11;

    py::class_<Testing::Tester, Kratos::shared_ptr<Testing::Tester> > TesterPyBind(m, "Tester");

    // Properties
    TesterPyBind
        .def_static("SetVerbosity",&Testing::Tester::SetVerbosity)
        // Run methods
        .def_static("RunAllTestCases", &Testing::Tester::RunAllTestCases)
        .def_static("RunTestSuite", &Testing::Tester::RunTestSuite)
        .def_static("RunTestCases", &Testing::Tester::RunTestCases)

        // Profile tests
        .def_static("ProfileAllTestCases", &Testing::Tester::ProfileAllTestCases)
        .def_static("ProfileTestSuite", &Testing::Tester::ProfileTestSuite)

        // Utils
        .def_static("NumberOfFailedTestCases", &Testing::Tester::NumberOfFailedTestCases)
        .def_static("ResetAllTestCasesResults", &Testing::Tester::ResetAllTestCasesResults)

        // Info
        .def_static("ListOfAllTestCases", ListOfAllTestCases)
    ;

    py::enum_<Testing::Tester::Verbosity>(TesterPyBind, "Verbosity")
        .value("QUITE", Testing::Tester::Verbosity::QUITE)
        .value("PROGRESS", Testing::Tester::Verbosity::PROGRESS)
        .value("TESTS_LIST", Testing::Tester::Verbosity::TESTS_LIST)
        .value("FAILED_TESTS_OUTPUTS", Testing::Tester::Verbosity::FAILED_TESTS_OUTPUTS)
        .value("TESTS_OUTPUTS", Testing::Tester::Verbosity::TESTS_OUTPUTS)
    ;
}

}  // namespace Python.
} // Namespace Kratos
