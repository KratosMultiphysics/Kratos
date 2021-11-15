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

#ifndef CO_SIM_IO_COMMUNICATION_INCLUDED
#define CO_SIM_IO_COMMUNICATION_INCLUDED

// System includes
#include <utility>

// Project includes
#include "includes/info.hpp"
#include "includes/data_container.hpp"
#include "includes/model_part.hpp"
#include "includes/data_communicator.hpp"
#include "includes/filesystem_inc.hpp"

namespace CoSimIO {
namespace Internals {


class CO_SIM_IO_API Communication
{
public:
    Communication(
        const Info& I_Settings,
        std::shared_ptr<DataCommunicator> I_DataComm);

    // might throw when trying to remove files!
    virtual ~Communication() noexcept(false) {}; // impl of disconnect has to be in derived class due to order of class destruction

    Info Connect(const Info& I_Info);

    Info Disconnect(const Info& I_Info);

    template<class... Args>
    Info ExportInfo(Args&&... args)
    {
        CheckConnection(); return ExportInfoImpl(std::forward<Args>(args)...);
    }

    template<class... Args>
    Info ImportInfo(Args&&... args)
    {
        CheckConnection(); return ImportInfoImpl(std::forward<Args>(args)...);
    }

    template<class... Args>
    Info ImportData(Args&&... args)
    {
        CheckConnection(); return ImportDataImpl(std::forward<Args>(args)...);
    }

    template<class... Args>
    Info ExportData(Args&&... args)
    {
        CheckConnection(); return ExportDataImpl(std::forward<Args>(args)...);
    }

    template<class... Args>
    Info ImportMesh(Args&&... args)
    {
        CheckConnection(); return ImportMeshImpl(std::forward<Args>(args)...);
    }

    template<class... Args>
    Info ExportMesh(Args&&... args)
    {
        CheckConnection(); return ExportMeshImpl(std::forward<Args>(args)...);
    }

protected:
    std::string GetConnectionName() const      {return mConnectionName;}
    fs::path GetWorkingDirectory() const       {return mWorkingDirectory;}
    fs::path GetCommunicationDirectory() const {return mCommFolder;}
    int GetEchoLevel() const                   {return mEchoLevel;}
    bool GetIsPrimaryConnection() const        {return mIsPrimaryConnection;}
    bool GetPrintTiming() const                {return mPrintTiming;}
    bool GetIsConnected() const                {return mIsConnected;}
    const DataCommunicator& GetDataCommunicator() const                {return *mpDataComm;}

    Info GetMyInfo() const;
    Info GetPartnerInfo() const {return mPartnerInfo;};

    fs::path GetTempFileName(const fs::path& rPath) const;

    fs::path GetFileName(const fs::path& rPath, const std::string& rExtension) const;

    fs::path GetFileName(const fs::path& rPath, const int Rank, const std::string& rExtension) const;

    void WaitForPath(const fs::path& rPath) const;

    void WaitUntilFileIsRemoved(const fs::path& rPath) const;

    void MakeFileVisible(const fs::path& rPath) const;

    void RemovePath(const fs::path& rPath) const;

    void SynchronizeAll() const;

private:
    std::shared_ptr<DataCommunicator> mpDataComm;

    std::string mConnectionName;
    std::string mMyName;
    std::string mConnectTo;

    Info mPartnerInfo;

    fs::path mCommFolder;
    bool mCommInFolder = true;
    bool mUseAuxFileForFileAvailability = false;

    fs::path mWorkingDirectory;
    int mEchoLevel = 1;
    bool mIsPrimaryConnection;
    bool mPrimaryWasExplicitlySpecified;
    bool mPrintTiming = false;
    bool mIsConnected = false;

    void CheckConnection()
    {
        CO_SIM_IO_ERROR_IF_NOT(mIsConnected) << "No active connection exists!" << std::endl;;
    }

    virtual std::string GetCommunicationName() const = 0;
    virtual Info GetCommunicationSettings() const {return Info();}

    virtual void BaseConnectDetail(const Info& I_Info);
    virtual void BaseDisconnectDetail(const Info& I_Info);
    virtual Info ConnectDetail(const Info& I_Info){return Info();}
    virtual Info DisconnectDetail(const Info& I_Info){return Info();}

    virtual Info ImportInfoImpl(const Info& I_Info)
    {
        CO_SIM_IO_ERROR << "ImportInfo not implemented for this comm-type" << std::endl;
        return Info();
    }

    virtual Info ExportInfoImpl(const Info& I_Info)
    {
        CO_SIM_IO_ERROR << "ExportInfo not implemented for this comm-type" << std::endl;
        return Info();
    }

    virtual Info ImportDataImpl(
        const Info& I_Info,
        Internals::DataContainer<double>& rData)
    {
        CO_SIM_IO_ERROR << "ImportDataImpl not implemented for this comm-type!" << std::endl;
        return Info();
    }

    virtual Info ExportDataImpl(
        const Info& I_Info,
        const Internals::DataContainer<double>& rData)
    {
        CO_SIM_IO_ERROR << "ExportDataImpl not implemented for this comm-type!" << std::endl;
        return Info();
    }

    virtual Info ImportMeshImpl(
        const Info& I_Info,
        ModelPart& O_ModelPart)
    {
        CO_SIM_IO_ERROR << "ImportMeshImpl not implemented for this comm-type!" << std::endl;
        return Info();
    }

    virtual Info ExportMeshImpl(
        const Info& I_Info,
        const ModelPart& I_ModelPart)
    {
        CO_SIM_IO_ERROR << "ExportMeshImpl not implemented for this comm-type!" << std::endl;
        return Info();
    }

    void HandShake(const Info& I_Info);

    virtual void DerivedHandShake() {};
};

} // namespace Internals
} // namespace CoSimIO

#endif // CO_SIM_IO_COMMUNICATION_INCLUDED
