// KRATOS  / ___|___/ ___|(_)_ __ ___  _   _| | __ _| |_(_) ___  _ ___
//        | |   / _ \___ \| | '_ ` _ \| | | | |/ _` | __| |/ _ \| '_  |
//        | |__| (_) |__) | | | | | | | |_| | | (_| | |_| | (_) | | | |
//         \____\___/____/|_|_| |_| |_|\__,_|_|\__,_|\__|_|\___/|_| |_|
//
//  License:		 BSD License
//                   license: CoSimulationApplication/license.txt
//
//  Main authors:    Philipp Bucher
//

#ifndef KRATOS_CO_SIM_IO_H_INCLUDED
#define KRATOS_CO_SIM_IO_H_INCLUDED

/*
This file defines the IO of Kratos-CoSimulation for the exchange of data
with external solvers.
By default the communication is done through files,
support for sockets and MPI can optionally be enabled
*/

// #define KRATOS_CO_SIM_IO_ENABLE_SOCKETS // uncomment for Sockets support
// #define KRATOS_CO_SIM_IO_ENABLE_MPI // uncomment for MPI support

// System includes
#include <string>
#include <memory>

// Project includes
#include "co_sim_io_define.h"
#include "co_sim_comm.h"
#include "co_sim_data_containers.h"

namespace CoSim {

class CoSimIO
{

public:
    typedef CoSimComm::SettingsType SettingsType;

    explicit CoSimIO(const std::string& rName, const std::string& rSettingsFileName, const bool IsConnectionMaster=false);
    explicit CoSimIO(const std::string& rName, SettingsType rSettings, const bool IsConnectionMaster=false);

    bool Connect();
    bool Disconnect();

    void SendControlSignal(const Internals::ControlSignal Signal, const std::string& rIdentifier);
    bool IsConverged();

    void RegisterAdvanceInTime(double (*pFuncPtr)(double));
    void RegisterInitializeSolutionStep(void (*pFuncPtr)());
    void RegisterSolveSolutionStep(void (*pFuncPtr)());
    void RegisterFinalizeSolutionStep(void (*pFuncPtr)());

    void RegisterImportGeometry(void (*pFuncPtr)(const std::string&));
    void RegisterExportGeometry(void (*pFuncPtr)(const std::string&));
    void RegisterImportMesh(void (*pFuncPtr)(const std::string&));
    void RegisterExportMesh(void (*pFuncPtr)(const std::string&));
    void RegisterImportData(void (*pFuncPtr)(const std::string&));
    void RegisterExportData(void (*pFuncPtr)(const std::string&));

    void Run();

    template<class DataContainer>
    bool Import(DataContainer& rDataContainer, const std::string& rIdentifier);

    template<class DataContainer>
    bool Export(const DataContainer& rDataContainer, const std::string& rIdentifier);

private:
    std::unique_ptr<CoSimComm> mpComm; // handles communication (File, Sockets, MPI, ...)

    double (*mpAdvInTime)(double) = nullptr;
    void (*mpInitSolStep)() = nullptr;
    void (*mpSolSolStep)() = nullptr;
    void (*mpFinSolStep)() = nullptr;

    void (*mpImportGeom)(const std::string&) = nullptr;
    void (*mpExportGeom)(const std::string&) = nullptr;

    void (*mpImportMesh)(const std::string&) = nullptr;
    void (*mpExportMesh)(const std::string&) = nullptr;

    void (*mpImportData)(const std::string&) = nullptr;
    void (*mpExportData)(const std::string&) = nullptr;

    void Initialize(const std::string& rName, SettingsType& rSettings, const bool IsConnectionMaster);
    Internals::ControlSignal RecvControlSignal(std::string& rIdentifier);

}; // class CoSimIO

} // namespace CoSim

// #include "co_sim_io_impl.h"

#endif /* KRATOS_CO_SIM_IO_H_INCLUDED */
