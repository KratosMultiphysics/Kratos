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
#include <boost/python.hpp>

// Project includes
#include "testing/testing.h"


namespace Kratos
{

namespace Python
{

void ListOfAllTestCases() {
	std::cout << Testing::Tester::GetInstance() << std::endl;
}

void  AddTestingToPython() {
	using namespace boost::python;

  scope tester_scope = class_<Testing::Tester, boost::shared_ptr<Testing::Tester>, boost::noncopyable>("Tester", no_init)

  // Properties
  .def("SetVerbosity",&Testing::Tester::SetVerbosity)
  .staticmethod("SetVerbosity")

  // Run
  .def("RunAllTestCases", &Testing::Tester::RunAllTestCases)
  .staticmethod("RunAllTestCases")
  .def("RunTestSuite", &Testing::Tester::RunTestSuite)
  .staticmethod("RunTestSuite")
  .def("RunTestCases", &Testing::Tester::RunTestCases)
  .staticmethod("RunTestCases")

  // Profile tests
  .def("ProfileAllTestCases", &Testing::Tester::ProfileAllTestCases)
  .staticmethod("ProfileAllTestCases")
  .def("ProfileTestSuite", &Testing::Tester::ProfileTestSuite)
  .staticmethod("ProfileTestSuite")

  // Utils
  .def("NumberOfFailedTestCases", &Testing::Tester::NumberOfFailedTestCases)
  .staticmethod("NumberOfFailedTestCases")
  .def("ResetAllTestCasesResults", &Testing::Tester::ResetAllTestCasesResults)
  .staticmethod("ResetAllTestCasesResults")

  // Info
  .def("ListOfAllTestCases", ListOfAllTestCases)
  .staticmethod("ListOfAllTestCases")

  ;

  enum_<Testing::Tester::Verbosity>("Verbosity")
  // Enums
  .value("QUITE", Testing::Tester::Verbosity::QUITE)
  .value("PROGRESS", Testing::Tester::Verbosity::PROGRESS)
  .value("TESTS_LIST", Testing::Tester::Verbosity::TESTS_LIST)
  .value("FAILED_TESTS_OUTPUTS", Testing::Tester::Verbosity::FAILED_TESTS_OUTPUTS)
  .value("TESTS_OUTPUTS", Testing::Tester::Verbosity::TESTS_OUTPUTS);

}
}  // namespace Python.

} // Namespace Kratos
