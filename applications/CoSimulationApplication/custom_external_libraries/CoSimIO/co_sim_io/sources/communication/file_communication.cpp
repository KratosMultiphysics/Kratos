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

// System includes
#include <chrono>

// Project includes
#include "includes/communication/file_communication.hpp"
#include "includes/utilities.hpp"
#include "includes/file_serializer.hpp"

namespace CoSimIO {
namespace Internals {

namespace {

// Important: having this in a function also makes sure that the FileSerializer releases its resources (i.e. the file) at destruction
template<class TObject>
void SerializeToFile(const fs::path& rPath, const std::string& rTag, const TObject& rObject, const Serializer::TraceType SerializerTrace)
{
    CO_SIM_IO_TRY

    FileSerializer serializer(rPath.string(), SerializerTrace);
    serializer.save(rTag, rObject);

    CO_SIM_IO_CATCH
}

// important: having this in a function also makes sure that the FileSerializer releases its resources (i.e. the file) at destruction
template<class TObject>
void SerializeFromFile(const fs::path& rPath, const std::string& rTag, TObject& rObject, const Serializer::TraceType SerializerTrace)
{
    CO_SIM_IO_TRY

    FileSerializer serializer(rPath.string(), SerializerTrace);
    serializer.load(rTag, rObject);

    CO_SIM_IO_CATCH
}

}

FileCommunication::FileCommunication(
    const Info& I_Settings,
    std::shared_ptr<DataCommunicator> I_DataComm)
    : Communication(I_Settings, I_DataComm)
{
}

FileCommunication::~FileCommunication()
{
    CO_SIM_IO_TRY

    if (GetIsConnected()) {
        CO_SIM_IO_INFO("CoSimIO") << "Warning: Disconnect was not performed, attempting automatic disconnection!" << std::endl;
        Info tmp;
        Disconnect(tmp);
    }

    CO_SIM_IO_CATCH
}

Info FileCommunication::ImportInfoImpl(const Info& I_Info)
{
    CO_SIM_IO_TRY

    const std::string identifier = I_Info.Get<std::string>("identifier");
    Utilities::CheckEntry(identifier, "identifier");

    const fs::path file_name(GetFileName("CoSimIO_info_" + GetConnectionName() + "_" + identifier + "_" + std::to_string(GetDataCommunicator().Rank()), "dat"));

    CO_SIM_IO_INFO_IF("CoSimIO", GetEchoLevel()>1) << "Attempting to import Info in file " << file_name << " ..." << std::endl;

    WaitForPath(file_name);

    Info info;
    SerializeFromFile(file_name, "info", info, Serializer::TraceType::SERIALIZER_NO_TRACE);

    RemovePath(file_name);

    CO_SIM_IO_INFO_IF("CoSimIO", GetEchoLevel()>1) << "Finished importing Info" << std::endl;

    return info;

    CO_SIM_IO_CATCH
}

Info FileCommunication::ExportInfoImpl(const Info& I_Info)
{
    CO_SIM_IO_TRY

    const std::string identifier = I_Info.Get<std::string>("identifier");
    Utilities::CheckEntry(identifier, "identifier");

    const fs::path file_name(GetFileName("CoSimIO_info_" + GetConnectionName() + "_" + identifier + "_" + std::to_string(GetDataCommunicator().Rank()), "dat"));

    CO_SIM_IO_INFO_IF("CoSimIO", GetEchoLevel()>1) << "Attempting to export Info in file " << file_name << " ..." << std::endl;

    WaitUntilFileIsRemoved(file_name); // TODO maybe this can be queued somehow ... => then it would not block the sender

    SerializeToFile(GetTempFileName(file_name), "info", I_Info, Serializer::TraceType::SERIALIZER_NO_TRACE);

    MakeFileVisible(file_name);

    CO_SIM_IO_INFO_IF("CoSimIO", GetEchoLevel()>1) << "Finished exporting Info" << std::endl;

    return Info(); // TODO use

    CO_SIM_IO_CATCH
}

Info FileCommunication::ImportDataImpl(
    const Info& I_Info,
    Internals::DataContainer<double>& rData)
{
    CO_SIM_IO_TRY

    const std::string identifier = I_Info.Get<std::string>("identifier");
    Utilities::CheckEntry(identifier, "identifier");

    const fs::path file_name(GetFileName("CoSimIO_data_" + GetConnectionName() + "_" + identifier + "_" + std::to_string(GetDataCommunicator().Rank()), "dat"));

    CO_SIM_IO_INFO_IF("CoSimIO", GetEchoLevel()>1) << "Attempting to import array \"" << identifier << "\" in file " << file_name << " ..." << std::endl;

    WaitForPath(file_name);

    const auto start_time(std::chrono::steady_clock::now());

    SerializeFromFile(file_name, "data", rData, Serializer::TraceType::SERIALIZER_NO_TRACE);

    RemovePath(file_name);

    CO_SIM_IO_INFO_IF("CoSimIO", GetEchoLevel()>1) << "Finished importing array with size: " << rData.size() << std::endl;

    CO_SIM_IO_INFO_IF("CoSimIO", GetPrintTiming()) << "Importing Array \"" << identifier << "\" took: " << Utilities::ElapsedSeconds(start_time) << " [sec]" << std::endl;

    return Info(); // TODO use

    CO_SIM_IO_CATCH
}

Info FileCommunication::ExportDataImpl(
    const Info& I_Info,
    const Internals::DataContainer<double>& rData)
{
    CO_SIM_IO_TRY

    const std::string identifier = I_Info.Get<std::string>("identifier");
    Utilities::CheckEntry(identifier, "identifier");

    const fs::path file_name(GetFileName("CoSimIO_data_" + GetConnectionName() + "_" + identifier + "_" + std::to_string(GetDataCommunicator().Rank()), "dat"));

    WaitUntilFileIsRemoved(file_name); // TODO maybe this can be queued somehow ... => then it would not block the sender

    const std::size_t size = rData.size();
    CO_SIM_IO_INFO_IF("CoSimIO", GetEchoLevel()>1) << "Attempting to export array \"" << identifier << "\" with size: " << size << " in file " << file_name << " ..." << std::endl;

    const auto start_time(std::chrono::steady_clock::now());

    SerializeToFile(GetTempFileName(file_name), "data", rData, Serializer::TraceType::SERIALIZER_NO_TRACE);

    MakeFileVisible(file_name);

    CO_SIM_IO_INFO_IF("CoSimIO", GetEchoLevel()>1) << "Finished exporting array" << std::endl;

    CO_SIM_IO_INFO_IF("CoSimIO", GetPrintTiming()) << "Exporting Array \"" << identifier << "\" took: " << Utilities::ElapsedSeconds(start_time) << " [sec]" << std::endl;

    return Info(); // TODO use

    CO_SIM_IO_CATCH
}

Info FileCommunication::ImportMeshImpl(
    const Info& I_Info,
    ModelPart& O_ModelPart)
{
    CO_SIM_IO_TRY

    const std::string identifier = I_Info.Get<std::string>("identifier");
    Utilities::CheckEntry(identifier, "identifier");

    const fs::path file_name(GetFileName("CoSimIO_mesh_" + GetConnectionName() + "_" + identifier + "_" + std::to_string(GetDataCommunicator().Rank()), "vtk"));

    CO_SIM_IO_INFO_IF("CoSimIO", GetEchoLevel()>1) << "Attempting to import mesh \"" << identifier << "\" in file " << file_name << " ..." << std::endl;

    WaitForPath(file_name);

    const auto start_time(std::chrono::steady_clock::now());

    SerializeFromFile(file_name, "model_part", O_ModelPart, Serializer::TraceType::SERIALIZER_NO_TRACE);

    RemovePath(file_name);

    CO_SIM_IO_INFO_IF("CoSimIO", GetEchoLevel()>1) << "Finished importing mesh" << std::endl;

    CO_SIM_IO_INFO_IF("CoSimIO", GetPrintTiming()) << "Importing Mesh \"" << identifier << "\" took: " << Utilities::ElapsedSeconds(start_time) << " [sec]" << std::endl;

    return Info(); // TODO use

    CO_SIM_IO_CATCH
}

Info FileCommunication::ExportMeshImpl(
    const Info& I_Info,
    const ModelPart& I_ModelPart)
{
    CO_SIM_IO_TRY

    const std::string identifier = I_Info.Get<std::string>("identifier");
    Utilities::CheckEntry(identifier, "identifier");

    const fs::path file_name(GetFileName("CoSimIO_mesh_" + GetConnectionName() + "_" + identifier + "_" + std::to_string(GetDataCommunicator().Rank()), "vtk"));

    WaitUntilFileIsRemoved(file_name); // TODO maybe this can be queued somehow ... => then it would not block the sender

    CO_SIM_IO_INFO_IF("CoSimIO", GetEchoLevel()>1) << "Attempting to export mesh \"" << identifier << "\" with " << I_ModelPart.NumberOfNodes() << " Nodes | " << I_ModelPart.NumberOfElements() << " Elements in file " << file_name << " ..." << std::endl;

    const auto start_time(std::chrono::steady_clock::now());

    SerializeToFile(GetTempFileName(file_name), "model_part", I_ModelPart, Serializer::TraceType::SERIALIZER_NO_TRACE);

    MakeFileVisible(file_name);

    CO_SIM_IO_INFO_IF("CoSimIO", GetEchoLevel()>1) << "Finished exporting mesh" << std::endl;

    CO_SIM_IO_INFO_IF("CoSimIO", GetPrintTiming()) << "Exporting Mesh \"" << identifier << "\" took: " << Utilities::ElapsedSeconds(start_time) << " [sec]" << std::endl;

    return Info(); // TODO use

    CO_SIM_IO_CATCH
}

} // namespace Internals
} // namespace CoSimIO
