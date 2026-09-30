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

// Exposure of the CoSimIO to Python

// System includes
#include <functional>
#include <vector>
#include <string>
#include <tuple>

// External includes
#include "intrusive_ptr/intrusive_ptr.hpp"

// pybind includes
#include <pybind11/pybind11.h>
#include <pybind11/functional.h>
#include <pybind11/numpy.h>

// must be done before CoSimIO includes
PYBIND11_DECLARE_HOLDER_TYPE(T, CoSimIO::intrusive_ptr<T>)

// CoSimIO includes
#include "co_sim_io.hpp"
#include "info_to_python.hpp"
#include "model_part_to_python.hpp"
#include "vector_to_python.hpp"
#include "connection_status_to_python.hpp"
#include "version_to_python.hpp"


PYBIND11_MODULE(PyCoSimIO, m)
{
    namespace py = pybind11;

    m.def("Hello", & CoSimIO::Hello);

    m.def("Connect",    &CoSimIO::Connect);
    m.def("Disconnect", &CoSimIO::Disconnect);

    m.def("ImportMesh", &CoSimIO::ImportMesh);
    m.def("ExportMesh", &CoSimIO::ExportMesh);

    m.def("ImportData", [](const CoSimIO::Info& I_Info, CoSimIO::VectorWrapper<double>& rValues){
        return CoSimIO::ImportData(
        I_Info,
        rValues.Vector());
    });
    m.def("ExportData", [](const CoSimIO::Info& I_Info, const CoSimIO::VectorWrapper<double>& rValues){
        return CoSimIO::ExportData(
        I_Info,
        rValues.Vector());
    });

    m.def("ImportInfo", &CoSimIO::ImportInfo);
    m.def("ExportInfo", &CoSimIO::ExportInfo);

    // functions for CoSim-orchestrated CoSimulation
    m.def("Run", &CoSimIO::Run);

    m.def("Register", [](
        const CoSimIO::Info& I_Info,
        std::function<CoSimIO::Info(const CoSimIO::Info&)> FunctionPointer)
        { return CoSimIO::Register(I_Info, FunctionPointer); } );

    CoSimIO::AddCoSimIOInfoToPython(m);
    CoSimIO::AddCoSimIOModelPartToPython(m);
    CoSimIO::AddCoSimIOVectorToPython(m);
    CoSimIO::AddCoSimIOConnectionStatusToPython(m);
    CoSimIO::AddCoSimIOVersionToPython(m);
}
