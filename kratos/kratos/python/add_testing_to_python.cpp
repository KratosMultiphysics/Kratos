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

// External includes
#include <boost/python.hpp>

// Project includes
#include "testing/tester.h"


namespace Kratos
{

namespace Python
{

void ListOfAllTestCases()
{
	std::cout << Testing::Tester::GetInstance() << std::endl;
}

void  AddTestingToPython()
{
	using namespace boost::python;


	def("RunAllTestCases", &Testing::Tester::RunAllTests);
	def("ListOfAllTestCases", ListOfAllTestCases);
}
}  // namespace Python.

} // Namespace Kratos

