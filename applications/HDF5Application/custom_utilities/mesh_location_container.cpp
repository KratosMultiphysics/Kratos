//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:        BSD License
//                  license: HDF5Application/license.txt
//
//  Main author:    Suneth Warnakulasuriya
//

// System includes

// Project includes

// Application includes
#include "hdf5_application_variables.h"

// Include base h
#include "mesh_location_container.h"

namespace Kratos
{

namespace HDF5
{

void MeshLocationContainer::Set(
    const int HDF5RankId,
    const int HDF5ProcessId,
    const std::string& rMeshLocation)
{
    for (IndexType i = 0; i < mProcessIds.size(); ++i) {
        const auto& current_id = mProcessIds[i];
        if (std::get<0>(current_id) == HDF5RankId && std::get<1>(current_id) == HDF5ProcessId) {
            mMeshLocations[i] = rMeshLocation;
            return;
        }
    }

    mProcessIds.push_back(std::make_pair(HDF5RankId, HDF5ProcessId));
    mMeshLocations.push_back(rMeshLocation);
}

bool MeshLocationContainer::Has(
    const int HDF5RankId,
    const int HDF5ProcessId) const
{
    for (IndexType i = 0; i < mProcessIds.size(); ++i) {
        const auto& current_id = mProcessIds[i];
        if (std::get<0>(current_id) == HDF5RankId && std::get<1>(current_id) == HDF5ProcessId) {
            return true;
        }
    }
    return false;
}

std::string MeshLocationContainer::Get(
    const int HDF5RankId,
    const int HDF5ProcessId) const
{
    KRATOS_TRY

    for (IndexType i = 0; i < mProcessIds.size(); ++i) {
        const auto& current_id = mProcessIds[i];
        if (std::get<0>(current_id) == HDF5RankId && std::get<1>(current_id) == HDF5ProcessId) {
            return mMeshLocations[i];
        }
    }

    KRATOS_ERROR
        << "Mesh location for HDF5RankId = "
        << HDF5RankId << ", and HDF5ProcessId = "
        << HDF5ProcessId << " not found.";

    KRATOS_CATCH("");
}

std::string MeshLocationContainer::Info() const
{
    return "MeshLocationContainer";
}

void MeshLocationContainer::PrintInfo(std::ostream& rOStream) const
{
    rOStream << Info();
}

void MeshLocationContainer::PrintData(std::ostream& rOStream) const
{
    rOStream << "Available (rank_id, process_id): mesh_locations :" << std::endl;
    for (IndexType i = 0; i < mProcessIds.size(); ++i) {
        const auto& current_id = mProcessIds[i];
        rOStream << "\t(" << std::get<0>(current_id) << ", "
                 << std::get<1>(current_id) << "): "
                 << mMeshLocations[i] << std::endl;
    }
}

void MeshLocationContainer::save(Kratos::Serializer &rSerializer) const {
    rSerializer.save("ProcessIds", mProcessIds);
    rSerializer.save("MeshLocations", mMeshLocations);
}
void MeshLocationContainer::load(Kratos::Serializer &rSerializer) {
    rSerializer.load("ProcessIds", mProcessIds);
    rSerializer.load("MeshLocations", mMeshLocations);
}

void SetMeshLocationContainer(
    ModelPart& rModelPart,
    MeshLocationContainer::Pointer pMeshLocationContainer)
{
    rModelPart.SetValue(HDF5_MESH_LOCATION_CONTAINER, pMeshLocationContainer);
}

bool HasMeshLocationContainer(const ModelPart& rModelPart)
{
    return rModelPart.Has(HDF5_MESH_LOCATION_CONTAINER);
}

MeshLocationContainer& GetMeshLocationContainer(ModelPart& rModelPart)
{
    return *rModelPart.GetValue(HDF5_MESH_LOCATION_CONTAINER);
}

const MeshLocationContainer& GetMeshLocationContainer(const ModelPart& rModelPart)
{
    return *rModelPart.GetValue(HDF5_MESH_LOCATION_CONTAINER);
}

void AddProcessId(
    Parameters Settings,
    const int HDF5RankId,
    const int HDF5ProcessId)
{
    KRATOS_TRY

    KRATOS_ERROR_IF(Settings.Has("__hdf5_process_id"))
        << "The \"__hdf5_process_id\" is a reserved attribute. Please remove it from attributes. Settings = \n"
        << Settings << "\n.";

    Settings.AddEmptyArray("__hdf5_process_id");
    Settings["__hdf5_process_id"].Append(HDF5RankId);
    Settings["__hdf5_process_id"].Append(HDF5ProcessId);

    KRATOS_CATCH("");
}

bool HasProcessId(const Parameters Settings)
{
    return (Settings.Has("__hdf5_process_id") && Settings["__hdf5_process_id"].IsArray());
}

std::pair<int, int> GetProcessId(const Parameters Settings)
{
    const auto& process_params = Settings["__hdf5_process_id"];
    return std::make_pair(process_params.GetArrayItem(0).GetInt(), process_params.GetArrayItem(1).GetInt());
}

std::ostream& operator << (
    std::ostream& rOStream,
    const MeshLocationContainer& rThis)
{
    rThis.PrintInfo(rOStream);
    rOStream << std::endl;
    rThis.PrintData(rOStream);
    return rOStream;
}

} // namespace HDF5
} // namespace Kratos