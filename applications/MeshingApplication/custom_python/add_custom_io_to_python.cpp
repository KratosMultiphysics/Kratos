// KRATOS  __  __ _____ ____  _   _ ___ _   _  ____
//        |  \/  | ____/ ___|| | | |_ _| \ | |/ ___|
//        | |\/| |  _| \___ \| |_| || ||  \| | |  _
//        | |  | | |___ ___) |  _  || || |\  | |_| |
//        |_|  |_|_____|____/|_| |_|___|_| \_|\____| APPLICATION
//
//  License:         BSD License
//                   license: MeshingApplication/license.txt
//
//  Main authors:    Pooyan Dadvand
//                   Josep Maria Carbonell Puigbo
//                   Vicente Matix Ferrandiz
//

// System includes

// External includes

// Project includes
#include "includes/define_python.h"
#include "custom_python/add_custom_io_to_python.h"

#ifdef INCLUDE_MMG
    #include "custom_io/mmg/mmg_io.h"
#endif

namespace Kratos::Python
{
namespace py = pybind11;

void  AddCustomIOToPython(pybind11::module& m)
{
#ifdef INCLUDE_MMG
    py::class_<MmgIO<MMGLibrary::MMG2D>, MmgIO<MMGLibrary::MMG2D>::Pointer, IO>(m, "MmgIO2D")
    .def(py::init<std::string const&>())
    .def(py::init<std::string const&, Parameters>())
    .def(py::init<std::string const&, Parameters, const Flags>())
    .def("GetMmgVersion", &MmgIO<MMGLibrary::MMG2D>::GetMmgVersion)
    ;
    py::class_<MmgIO<MMGLibrary::MMG3D>, MmgIO<MMGLibrary::MMG3D>::Pointer, IO>(m, "MmgIO3D")
    .def(py::init<std::string const&>())
    .def(py::init<std::string const&, Parameters>())
    .def(py::init<std::string const&, Parameters, const Flags>())
    .def("GetMmgVersion", &MmgIO<MMGLibrary::MMG3D>::GetMmgVersion)
    ;
    py::class_<MmgIO<MMGLibrary::MMGS>, MmgIO<MMGLibrary::MMGS>::Pointer, IO>(m, "MmgIOS")
    .def(py::init<std::string const&>())
    .def(py::init<std::string const&, Parameters>())
    .def(py::init<std::string const&, Parameters, const Flags>())
    .def("GetMmgVersion", &MmgIO<MMGLibrary::MMGS>::GetMmgVersion)
    ;
#endif

}

} // Namespace Kratos::Python.

