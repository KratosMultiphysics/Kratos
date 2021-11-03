//     ______     _____ _           ________
//    / ____/___ / ___/(_)___ ___  /  _/ __ |
//   / /   / __ \\__ \/ / __ `__ \ / // / / /
//  / /___/ /_/ /__/ / / / / / / // // /_/ /
//  \____/\____/____/_/_/ /_/ /_/___/\____/
//  Kratos CoSimulationApplication
//
//  License:         BSD License, see license.txt
//
//  Main authors:    Philipp Bucher (https://github.com/philbucher)
//

#ifndef CO_SIM_IO_INFO_TO_PYHON_INCLUDED
#define CO_SIM_IO_INFO_TO_PYHON_INCLUDED

// Exposure of the CoSimIO to Python

// System includes
#include <string>
#include <sstream>

// pybind includes
#include <pybind11/pybind11.h>

// CoSimIO include
#include "co_sim_io.hpp"

namespace {

template<typename TDataType>
void AddGetSetInterface(pybind11::class_<CoSimIO::Info>& PythonInfo, const std::string& Name)
{
    PythonInfo.def(std::string("Get"+Name).c_str(), [](const CoSimIO::Info& I_Info, const std::string& I_Key){return I_Info.Get<TDataType>(I_Key);});
    PythonInfo.def(std::string("Get"+Name).c_str(), [](const CoSimIO::Info& I_Info, const std::string& I_Key, const TDataType& I_Default){return I_Info.Get<TDataType>(I_Key, I_Default);});
    PythonInfo.def(std::string("Set"+Name).c_str(), &CoSimIO::Info::Set<TDataType>);
}

}

namespace CoSimIO {

void AddCoSimIOInfoToPython(pybind11::module& m)
{
    namespace py = pybind11;

    auto py_info = py::class_<CoSimIO::Info>(m,"Info")
        .def(py::init<>())
        .def(py::init<const CoSimIO::Info&>())
        .def("Has",       &CoSimIO::Info::Has)
        .def("Erase",     &CoSimIO::Info::Erase)
        .def("Clear",     &CoSimIO::Info::Clear)
        .def("Size",      &CoSimIO::Info::Size)
        .def("__len__",   [](const CoSimIO::Info& I_Info){return I_Info.Size();})
        .def("__str__",   [](const CoSimIO::Info& I_Info)
            { std::stringstream ss; ss << I_Info; return ss.str(); } )
        ;

    AddGetSetInterface<int>(py_info, "Int");
    AddGetSetInterface<double>(py_info, "Double");
    AddGetSetInterface<bool>(py_info, "Bool");
    AddGetSetInterface<std::string>(py_info, "String");
    AddGetSetInterface<Info>(py_info, "Info");
}

} // namespace CoSimIO

#endif // CO_SIM_IO_INFO_TO_PYHON_INCLUDED
