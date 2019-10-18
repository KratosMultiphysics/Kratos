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

#ifndef KRATOS_CO_SIM_COMM_H_INCLUDED
#define KRATOS_CO_SIM_COMM_H_INCLUDED

// System includes
#include <unordered_map>
#include <string>
#include <stdexcept>

// Project includes
#include "co_sim_api_internals.h"
#include "co_sim_data_containers.h"

namespace CoSim {

#define CO_SIM_COMM_REGISTER_DATA_CONTAINER_TYPE(DataContainerType)                                    \
    virtual bool ImportDetail(DataContainerType& rDataContainer, const std::string& rIdentifier)       \
        { KRATOS_CO_SIM_ERROR << "Type of data not yet supported" << std::endl; return false;}         \
    virtual bool ExportDetail(const DataContainerType& rDataContainer, const std::string& rIdentifier) \
        { KRATOS_CO_SIM_ERROR << "Type of data not yet supported" << std::endl; return false;}         \

class CoSimComm
{
public:
    typedef Internals::SettingsType SettingsType;

    explicit CoSimComm(const std::string& rName, SettingsType& rSettings, const bool IsConnectionMaster) : mrSettings(rSettings),mName(rName), mIsConnectionMaster(IsConnectionMaster)
    {
        const SettingsType default_settings = {
            {"echo_level",   "1"},
            {"print_timing", "0"}
        };
        Internals::AddMissingSettings(default_settings, mrSettings);

        mEchoLevel = std::stoi(mrSettings.at("echo_level"));
        mPrintTiming = (mrSettings.at("print_timing") == "1");
    }

    virtual ~CoSimComm() = default; // impl of disconnect has to be in derived class due to order of class destruction

    bool Connect()
    {
        CS_LOG << "Connecting \"" << mName << "\" as Connection-" << (mIsConnectionMaster ? "MASTER" : "SLAVE") << " ..." << std::endl;

        KRATOS_CO_SIM_ERROR_IF(mIsConnected) << "A connection was already established!" << std::endl;

        mIsConnected = ConnectDetail();

        KRATOS_CO_SIM_ERROR_IF_NOT(mIsConnected) << "Connection was not successful!" << std::endl;

        CS_LOG << "Connection established" << std::endl;

        return mIsConnected;
    }

    bool Disconnect()
    {
        CS_LOG << "Disconnecting \"" << mName << "\" ..." << std::endl;

        if (mIsConnected) {
            mIsConnected = !DisconnectDetail();
            if (mIsConnected) {
                CS_LOG << "Warning: Disconnect was not successful!" << std::endl;
                return false;
            }
        } else {
            CS_LOG << "Warning: Calling Disconnect but there was no active connection!" << std::endl;
            return false;
        }

        CS_LOG << "Disconnecting successful" << std::endl;

        return true;
    }

    void SendControlSignal(Internals::ControlSignal Signal, const std::string& rIdentifier)
    {
        CheckConnection(); return SendControlSignalDetail(Signal, rIdentifier);
    }
    Internals::ControlSignal RecvControlSignal(std::string& rIdentifier)
    {
        CheckConnection(); return RecvControlSignalDetail(rIdentifier);
    }

    template<class DataContainer>
    bool Import(DataContainer& rDataContainer, const std::string& rIdentifier)
        { CheckConnection(); return ImportDetail(rDataContainer, rIdentifier); }

    template<class DataContainer>
    bool Export(const DataContainer& rDataContainer, const std::string& rIdentifier)
        { CheckConnection(); return ExportDetail(rDataContainer, rIdentifier); }

protected:
    SettingsType& mrSettings;

    std::string GetName() const        {return mName;}
    int GetEchoLevel() const           {return mEchoLevel;}
    bool GetIsConnectionMaster() const {return mIsConnectionMaster;}
    bool GetPrintTiming() const        {return mPrintTiming;}
    bool GetIsConnected() const        {return mIsConnected;}

private:
    std::string mName;
    int mEchoLevel = 1;
    bool mIsConnectionMaster = false;
    bool mPrintTiming = false;
    bool mIsConnected = false;

    virtual bool ConnectDetail() = 0;
    virtual bool DisconnectDetail() = 0;

    virtual void SendControlSignalDetail(Internals::ControlSignal Signal, const std::string& rIdentifier)
    {
        KRATOS_CO_SIM_ERROR << "SendControlSignalDetail not implemented for this comm-type" << std::endl;
    }
    virtual Internals::ControlSignal RecvControlSignalDetail(std::string& rIdentifier)
    {
        KRATOS_CO_SIM_ERROR << "RecvControlSignalDetail not implemented for this comm-type" << std::endl;
        return Internals::ControlSignal::Dummy;
    }

    CO_SIM_COMM_REGISTER_DATA_CONTAINER_TYPE(DataContainers::Geometry);
    CO_SIM_COMM_REGISTER_DATA_CONTAINER_TYPE(DataContainers::Mesh);
    CO_SIM_COMM_REGISTER_DATA_CONTAINER_TYPE(DataContainers::Data);

    void CheckConnection()
    {
        KRATOS_CO_SIM_ERROR_IF_NOT(mIsConnected) << "No active connection exists!" << std::endl;;
    }

};

#undef CO_SIM_COMM_REGISTER_DATA_CONTAINER_TYPE // this macro should only be used for this class

} // namespace CoSim

#endif /* KRATOS_CO_SIM_COMM_H_INCLUDED */
