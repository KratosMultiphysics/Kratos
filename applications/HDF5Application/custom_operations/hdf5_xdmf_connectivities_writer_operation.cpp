//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:        BSD License
//                  license: HDF5Application/license.txt
//
//  Main author:    Michael Andre, https://github.com/msandre
//                  Suneht Warnakulasuriya
//

// System includes

// Project includes
#include "includes/kratos_parameters.h"
#include "utilities/parallel_utilities.h"

// Application includes
#include "custom_io/hdf5_file.h"

// Include base h
#include "hdf5_xdmf_connectivities_writer_operation.h"

namespace Kratos
{
namespace HDF5
{

namespace Detail
{

void GetModelDataGroups(
    std::vector<std::string>& rModelDataPaths,
    File& rFile,
    const std::string& rCurrentPath)
{
    KRATOS_TRY

    if (rCurrentPath != "/" && rFile.HasAttribute(rCurrentPath, "__model_part_name") && rFile.HasAttributeType<char>(rCurrentPath, "__model_part_name")) {
        // found a model data path. add it. No model data within a model data group is allowed.
        // Hence, recursive check is not required for this branch.
        rModelDataPaths.push_back(rCurrentPath);
    } else {
        // this does not contain model data, hence recursively searching.
        const auto& group_names = rFile.GetGroupNames(rCurrentPath);
        for (const auto& group_name : group_names) {
            GetModelDataGroups(rModelDataPaths, rFile, rCurrentPath + "/" + group_name);
        }
    }

    KRATOS_CATCH("Path: " << rCurrentPath);
}

} // namespace Detail

XdmfConnectivitiesWriterOperation::XdmfConnectivitiesWriterOperation(
    const std::string& rFileName,
    const std::vector<std::string>& rModelPrefixes)
    : mListOfModelDataPaths(rModelPrefixes)
{
    KRATOS_TRY;

    Parameters file_params(R"({
        "file_name"       : "",
        "file_access_mode": "read_write"
    })");
    file_params["file_name"].SetString(rFileName);

    // since the current implementation can only create the xdmf data in serial, we will
    // use the serial data communicator
    mpFile = File::Pointer(new File(mSerialDataCommunicator, file_params));

    if (mListOfModelDataPaths.empty()) {
        // generate the possible model data paths in the file.
        Detail::GetModelDataGroups(mListOfModelDataPaths, *mpFile, "/");
    }

    KRATOS_CATCH("");
}

void XdmfConnectivitiesWriterOperation::Execute()
{
    KRATOS_TRY;

    for (const auto& model_data_path : mListOfModelDataPaths) {
        CreateXdmfModelData(model_data_path);
    }

    KRATOS_CATCH("");
}

void XdmfConnectivitiesWriterOperation::CreateXdmfModelData(const std::string& rModelDataPath) const
{
    KRATOS_TRY

    // generate the node indices map
    std::string node_ids_path = rModelDataPath + "/Nodes/Local/Ids";

    KRATOS_ERROR_IF_NOT(mpFile->HasPath(node_ids_path))
        << "Path \"" << node_ids_path << "\" was not found." << std::endl;

    Vector<int> node_ids;
    const int num_points = mpFile->GetDataDimensions(node_ids_path)[0];
    mpFile->ReadDataSet(node_ids_path, node_ids, 0, num_points);

    // Set the parametric coordinate ids.
    IdMapType kratos_to_xdmf_id_map;
    kratos_to_xdmf_id_map.reserve(num_points);
    for (int i = 0; i < num_points; ++i) {
        const auto p = kratos_to_xdmf_id_map.insert(IdMapType::value_type(node_ids[i], i));
        KRATOS_ERROR_IF_NOT(p.second) << "Node #" << node_ids[i] << " at position "
            << i << " already exists at position " << p.first->second << std::endl;
    }

    // create the root level model part element connectivities
    KRATOS_ERROR_IF(mpFile->HasPath(rModelDataPath + "/Xdmf/Elements")) << "Path \"" << rModelDataPath + "/Xdmf/Elements\" exists." << std::endl;
    for (const auto& name : mpFile->GetGroupNames(rModelDataPath + "/Elements")) {
        CreateXdmfConnectivities(rModelDataPath + "/Elements/" + name, rModelDataPath + "/Xdmf/Elements/" + name, kratos_to_xdmf_id_map);
    }

    // create the root level model part condition connectivities
    KRATOS_ERROR_IF(mpFile->HasPath(rModelDataPath + "/Xdmf/Conditions")) << "Path \"" << rModelDataPath + "/Xdmf/Conditions\" exists." << std::endl;
    for (const auto& name : mpFile->GetGroupNames(rModelDataPath + "/Conditions")) {
        CreateXdmfConnectivities(rModelDataPath + "/Conditions/" + name, rModelDataPath + "/Xdmf/Conditions/" + name, kratos_to_xdmf_id_map);
    }

    // now create the submodel part connectivities
    KRATOS_ERROR_IF(mpFile->HasPath(rModelDataPath + "/Xdmf/SubModelParts")) << "Path \"" << rModelDataPath + "/Xdmf/SubModelParts\" exists." << std::endl;
    CreateXdmfConnectivitiesForSubModelParts(rModelDataPath + "/SubModelParts", rModelDataPath + "/Xdmf/SubModelParts", kratos_to_xdmf_id_map);

    KRATOS_CATCH("");
}

void XdmfConnectivitiesWriterOperation::CreateXdmfConnectivitiesForSubModelParts(
    const std::string& rPath,
    const std::string& rDestinationPrefix,
    const IdMapType& rKratosToXdmfIdMap) const
{
    CreateXdmfPoints(rPath, rDestinationPrefix, rKratosToXdmfIdMap);
    for (const auto& r_name : mpFile->GetGroupNames(rPath)) {
        if (r_name == "Conditions" ||  r_name == "Elements") {
            for (const auto& r_item_name : mpFile->GetGroupNames(rPath + "/" + r_name)) {
                CreateXdmfConnectivities(rPath + "/" + r_name + "/" + r_item_name,  rDestinationPrefix + "/" + r_name + "/" + r_item_name, rKratosToXdmfIdMap);
            }
        } else {
            CreateXdmfConnectivitiesForSubModelParts(rPath + "/" + r_name, rDestinationPrefix + "/" + r_name, rKratosToXdmfIdMap);
        }
    }
}

void XdmfConnectivitiesWriterOperation::CreateXdmfConnectivities(
    const std::string& rKratosConnectivitiesPath,
    const std::string& rXdmfConnectivitiesPath,
    const IdMapType& rKratosToXdmfIdMap) const
{
    KRATOS_TRY;

    Matrix<int> kratos_connectivities, xdmf_connectivities;

    const int block_size = mpFile->GetDataDimensions(rKratosConnectivitiesPath + "/Connectivities")[0];
    mpFile->ReadDataSet(rKratosConnectivitiesPath + "/Connectivities", kratos_connectivities, 0, block_size);

    xdmf_connectivities.resize(kratos_connectivities.size1(), kratos_connectivities.size2(), false);

    IndexPartition<int>(block_size).for_each([&kratos_connectivities, &rKratosToXdmfIdMap, &xdmf_connectivities](const auto i) {
        for (unsigned int j = 0; j < kratos_connectivities.size2(); ++j) {
            const auto kratos_id = kratos_connectivities(i, j);
            xdmf_connectivities(i, j) = rKratosToXdmfIdMap.at(kratos_id);
        }
    });

    WriteInfo info;
    mpFile->WriteDataSet(rXdmfConnectivitiesPath + "/Connectivities", xdmf_connectivities, info);

    int tmp;
    mpFile->ReadAttribute(rKratosConnectivitiesPath, "WorkingSpaceDimension", tmp);
    mpFile->WriteAttribute(rXdmfConnectivitiesPath, "WorkingSpaceDimension", tmp);
    mpFile->WriteAttribute(rXdmfConnectivitiesPath, "Dimension", tmp);
    mpFile->ReadAttribute(rKratosConnectivitiesPath, "NumberOfNodes", tmp);
    mpFile->WriteAttribute(rXdmfConnectivitiesPath, "NumberOfNodes", tmp);

    KRATOS_CATCH("");
}

void XdmfConnectivitiesWriterOperation::CreateXdmfPoints(
    const std::string& rKratosNodeIdsPath,
    const std::string& rXdmfNodeIdsPath,
    const IdMapType& rKratosToXdmfIdMap) const
{
    KRATOS_TRY

    if (mpFile->IsDataSet(rKratosNodeIdsPath + "/Nodes/Ids")) {
        Vector<int> kratos_ids, xdmf_ids;
        const int block_size = mpFile->GetDataDimensions(rKratosNodeIdsPath + "/Nodes/Ids")[0];
        mpFile->ReadDataSet(rKratosNodeIdsPath + "/Nodes/Ids", kratos_ids, 0, block_size);

        xdmf_ids.resize(kratos_ids.size(), false);

        IndexPartition<int>(block_size).for_each([&xdmf_ids, &rKratosToXdmfIdMap, &kratos_ids](const auto Index) {
            xdmf_ids[Index] = rKratosToXdmfIdMap.at(kratos_ids[Index]);
        });

        WriteInfo info;
        mpFile->WriteDataSet(rXdmfNodeIdsPath + "/Points", xdmf_ids, info);
    }

    KRATOS_CATCH("");
}

} // namespace HDF5.
} // namespace Kratos.