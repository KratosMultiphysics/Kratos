//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Pooyan Dadvand, Carlos A. Roig
//
//

// External includes

// Project includes
#include "includes/define_python.h"
#include "includes/parallel_environment.h"
// #include "testing/testing.h"
#include "add_testing_to_python.h"

namespace Kratos::Python
{

// void ListOfAllTestCases() {
//     std::cout << Testing::Tester::GetInstance() << std::endl;
// }

void  AddTestingToPython(pybind11::module& m) {
	namespace py = pybind11;

    auto m_testing = m.def_submodule("Testing");
    m_testing.def("GetDefaultDataCommunicator", []() -> DataCommunicator& {
        return ParallelEnvironment::GetDefaultDataCommunicator();
    }, py::return_value_policy::reference);
}

}  // namespace Kratos::Python.
