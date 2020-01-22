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

#ifndef KRATOS_CO_SIM_COMM_H_INCLUDED
#define KRATOS_CO_SIM_COMM_H_INCLUDED

// System includes

// Project includes
#include "co_sim_io_internals.h"

namespace CoSimIO {
namespace Internals {

class CoSimCommunication
{
public:
    explicit CoSimCommunication(const std::string& rName, SettingsType& rSettings, const bool IsConnectionMaster) : mrSettings(rSettings),mConnectionName(rName), mIsConnectionMaster(IsConnectionMaster)
    {
        const SettingsType default_settings {
            {"echo_level",   "1"},
            {"print_timing", "0"}
        };
        Internals::AddMissingSettings(default_settings, mrSettings);

        mEchoLevel = std::stoi(mrSettings.at("echo_level"));
        mPrintTiming = (mrSettings.at("print_timing") == "1");
    }

    virtual ~CoSimCommunication() = default; // impl of disconnect has to be in derived class due to order of class destruction

    bool Connect()
    {
        KRATOS_CO_SIM_INFO_IF("CoSimIO", GetEchoLevel()>0) << "Connecting \"" << mConnectionName << "\" as Connection-" << (mIsConnectionMaster ? "MASTER" : "SLAVE") << " ..." << std::endl;

        KRATOS_CO_SIM_ERROR_IF(mIsConnected) << "A connection was already established!" << std::endl;

        mIsConnected = ConnectDetail();

        KRATOS_CO_SIM_ERROR_IF_NOT(mIsConnected) << "Connection was not successful!" << std::endl;

        KRATOS_CO_SIM_INFO_IF("CoSimIO", GetEchoLevel()>0) << "Connection established" << std::endl;

        return mIsConnected;
    }

    bool Disconnect()
    {
        KRATOS_CO_SIM_INFO_IF("CoSimIO", GetEchoLevel()>0) << "Disconnecting \"" << mConnectionName << "\" ..." << std::endl;

        if (mIsConnected) {
            mIsConnected = !DisconnectDetail();
            if (mIsConnected) {
                KRATOS_CO_SIM_INFO("CoSimIO") << "Warning: Disconnect was not successful!" << std::endl;
                return false;
            }
        } else {
            KRATOS_CO_SIM_INFO("CoSimIO") << "Warning: Calling Disconnect but there was no active connection!" << std::endl;
            return false;
        }

        KRATOS_CO_SIM_INFO_IF("CoSimIO", GetEchoLevel()>0) << "Disconnecting successful" << std::endl;

        return true;
    }

    void SendControlSignal(const std::string& rIdentifier, const CoSimIO::ControlSignal Signal)
    {
        CheckConnection(); SendControlSignalDetail(rIdentifier, Signal);
    }
    CoSimIO::ControlSignal RecvControlSignal(std::string& rIdentifier)
    {
        CheckConnection(); return RecvControlSignalDetail(rIdentifier);
    }

    template<class... Args>
    void ImportData(Args&&... args)
    {
        CheckConnection(); ImportDataImpl(std::forward<Args>(args)...);
    }

    template<class... Args>
    void ExportData(Args&&... args)
    {
        CheckConnection(); ExportDataImpl(std::forward<Args>(args)...);
    }

    template<class... Args>
    void ImportMesh(Args&&... args)
    {
        CheckConnection(); ImportMeshImpl(std::forward<Args>(args)...);
    }

    template<class... Args>
    void ExportMesh(Args&&... args)
    {
        CheckConnection(); ExportMeshImpl(std::forward<Args>(args)...);
    }

    template<class... Args>
    void ImportGeometry(Args&&... args)
    {
        CheckConnection(); // ImportGeometryImpl(std::forward<Args>(args)...);
    }

    template<class... Args>
    void ExportGeometry(Args&&... args)
    {
        CheckConnection(); // ExportGeometryImpl(std::forward<Args>(args)...);
    }

protected:
    SettingsType& mrSettings;

    std::string GetConnectionName() const {return mConnectionName;}
    int GetEchoLevel() const              {return mEchoLevel;}
    bool GetIsConnectionMaster() const    {return mIsConnectionMaster;}
    bool GetPrintTiming() const           {return mPrintTiming;}
    bool GetIsConnected() const           {return mIsConnected;}

private:
    std::string mConnectionName;
    int mEchoLevel = 1;
    bool mIsConnectionMaster = false;
    bool mPrintTiming = false;
    bool mIsConnected = false;

    virtual bool ConnectDetail() = 0;
    virtual bool DisconnectDetail() = 0;

    virtual void SendControlSignalDetail(const std::string& rIdentifier, CoSimIO::ControlSignal Signal)
    {
        KRATOS_CO_SIM_ERROR << "SendControlSignalDetail not implemented for this comm-type" << std::endl;
    }
    virtual CoSimIO::ControlSignal RecvControlSignalDetail(std::string& rIdentifier)
    {
        KRATOS_CO_SIM_ERROR << "RecvControlSignalDetail not implemented for this comm-type" << std::endl;
        return CoSimIO::ControlSignal::Dummy;
    }

    virtual void ImportDataImpl(
        const std::string& rIdentifier,
        CoSimIO::Internals::DataContainer<double>& rData)
    {
        KRATOS_CO_SIM_ERROR << "ImportDataImpl not implemented for this comm-type!" << std::endl;
    }

    virtual void ExportDataImpl(
        const std::string& rIdentifier,
        const CoSimIO::Internals::DataContainer<double>& rData)
    {
        KRATOS_CO_SIM_ERROR << "ExportDataImpl not implemented for this comm-type!" << std::endl;
    }

    virtual void ImportMeshImpl(
        const std::string& rIdentifier,
        CoSimIO::Internals::DataContainer<double>& rNodalCoordinates,
        CoSimIO::Internals::DataContainer<int>& rElementConnectivities,
        CoSimIO::Internals::DataContainer<int>& rElementTypes)
    {
        KRATOS_CO_SIM_ERROR << "ImportDataImpl not implemented for this comm-type!" << std::endl;
    }

    virtual void ExportMeshImpl(
        const std::string& rIdentifier,
        CoSimIO::Internals::DataContainer<double>& rNodalCoordinates,
        CoSimIO::Internals::DataContainer<int>& rElementConnectivities,
        CoSimIO::Internals::DataContainer<int>& rElementTypes)
    {
        KRATOS_CO_SIM_ERROR << "ImportDataImpl not implemented for this comm-type!" << std::endl;
    }

    void CheckConnection()
    {
        KRATOS_CO_SIM_ERROR_IF_NOT(mIsConnected) << "No active connection exists!" << std::endl;;
    }
};

} // namespace Internals
} // namespace CoSimIO

#endif /* KRATOS_CO_SIM_COMM_H_INCLUDED */
