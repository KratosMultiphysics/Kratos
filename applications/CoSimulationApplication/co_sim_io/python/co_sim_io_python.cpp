// KRATOS  / ___|___/ ___|(_)_ __ ___  _   _| | __ _| |_(_) ___  _ ___
//        | |   / _ \___ \| | '_ ` _ \| | | | |/ _` | __| |/ _ \| '_  |
//        | |__| (_) |__) | | | | | | | |_| | | (_| | |_| | (_) | | | |
//         \____\___/____/|_|_| |_| |_|\__,_|_|\__,_|\__|_|\___/|_| |_|
//
//  License:		 BSD License
//                   license: CoSimulationApplication/license.txt
//
//  Main authors:    Philipp Bucher (https://github.com/philbucher)
//

// Exposure of the CoSimIO to Python

#include <functional>

#include <pybind11/pybind11.h>
#include <pybind11/functional.h>

#include "../co_sim_io.h"


PYBIND11_MODULE(CoSimIO, m)
{
    void (*ConnectWithSettingsFileName)(const std::string&, const std::string&) = &CoSimIO::Connect;
    void (*ConnectWithSettings)(const std::string&, CoSimIO::SettingsType) = &CoSimIO::Connect;

    m.def("Connect", ConnectWithSettingsFileName);
    m.def("Connect", ConnectWithSettings);

    m.def("Disconnect", &CoSimIO::Disconnect);

    m.def("IsConverged", &CoSimIO::IsConverged);

    // mCoSimIO.def("ImportMesh", CoSimIO_Wrappers::ImportMesh);
    // mCoSimIO.def("ExportMesh", CoSimIO_Wrappers::ExportMesh);

    // mCoSimIO.def("ImportData", CoSimIO_Wrappers::ImportData);
    // mCoSimIO.def("ExportData", CoSimIO_Wrappers::ExportData);

    // functions for CoSim-orchestrated CoSimulation
    m.def("Run", &CoSimIO::Run);

    m.def("Register_AdvanceInTime",          [](const std::string& rConnectionName, std::function<double(double)> FunctionPointer){
        CoSimIO::Register(rConnectionName, "AdvanceInTime", FunctionPointer);});
    m.def("Register_InitializeSolutionStep", [](const std::string& rConnectionName, std::function<void()> FunctionPointer){
        CoSimIO::Register(rConnectionName, "InitializeSolutionStep", FunctionPointer);});
    m.def("Register_SolveSolutionStep",      [](const std::string& rConnectionName, std::function<void()> FunctionPointer){
        CoSimIO::Register(rConnectionName, "SolveSolutionStep", FunctionPointer);});
    m.def("Register_FinalizeSolutionStep",   [](const std::string& rConnectionName, std::function<void()> FunctionPointer){
        CoSimIO::Register(rConnectionName, "FinalizeSolutionStep", FunctionPointer);});

    m.def("Register_ImportData",   [](const std::string& rConnectionName, std::function<void(const std::string&, const std::string&)> FunctionPointer){
        CoSimIO::Register(rConnectionName, "ImportData", FunctionPointer);});;
    m.def("Register_ExportData",   [](const std::string& rConnectionName, std::function<void(const std::string&, const std::string&)> FunctionPointer){
        CoSimIO::Register(rConnectionName, "ExportData", FunctionPointer);});

    m.def("Register_ImportMesh",   [](const std::string& rConnectionName, std::function<void(const std::string&, const std::string&)> FunctionPointer){
        CoSimIO::Register(rConnectionName, "ImportMesh", FunctionPointer);});;
    m.def("Register_ExportMesh",   [](const std::string& rConnectionName, std::function<void(const std::string&, const std::string&)> FunctionPointer){
        CoSimIO::Register(rConnectionName, "ExportMesh", FunctionPointer);});

}
