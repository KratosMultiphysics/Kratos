#include "custom_processes/hdf5_xdmf_connectivities_writer_process.h"

#include "custom_io/hdf5_file_serial.h"
#include "utilities/openmp_utils.h"

namespace Kratos
{
namespace HDF5
{
XdmfConnectivitiesWriterProcess::XdmfConnectivitiesWriterProcess(const std::string& rFileName, const std::string& rPrefix)
{
    KRATOS_TRY;

    Parameters file_params(R"(
            {
                "file_name" : "",
                "file_access_mode": "read_write"
            })");
    file_params["file_name"].SetString(rFileName);
    mpFile = FileSerial::Pointer(new FileSerial(file_params));
    mPrefix = rPrefix;
    std::string node_ids_path = mPrefix + "/Nodes/Local/Ids";

    KRATOS_ERROR_IF(mpFile->HasPath(node_ids_path) == false)
        << "Path \"" << node_ids_path << "\" was not found." << std::endl;

    Vector<int> node_ids;
    const int num_points = mpFile->GetDataDimensions(node_ids_path)[0];
    mpFile->ReadDataSet(node_ids_path, node_ids, 0, num_points);

    // Set the parametric coordinate ids.
    mKratosToXdmfIdTable.resize(num_points);
#pragma omp parallel for
    for (int i = 0; i < num_points; ++i)
    {
        // Here we expect the ids are one-based indices.
        mKratosToXdmfIdTable(node_ids[i] - 1) = i;
    }

    KRATOS_CATCH("");
}

void XdmfConnectivitiesWriterProcess::Execute()
{
    KRATOS_TRY;

    std::vector<std::string> labels;
    KRATOS_ERROR_IF(mpFile->HasPath(mPrefix + "/Xdmf/Elements")) << "Path \"" << mPrefix + "/Xdmf/Elements\" exists." << std::endl;
    mpFile->GetGroupNames(mPrefix + "/Elements", labels);
    for (unsigned i = 0; i < labels.size(); ++i)
        CreateXdmfConnectivities(mPrefix + "/Elements/" + labels[i], mPrefix + "/Xdmf/Elements/" + labels[i]);

    KRATOS_ERROR_IF(mpFile->HasPath(mPrefix + "/Xdmf/Conditions")) << "Path \"" << mPrefix + "/Xdmf/Conditions\" exists." << std::endl;
    mpFile->GetGroupNames(mPrefix + "/Conditions", labels);
    for (unsigned i = 0; i < labels.size(); ++i)
        CreateXdmfConnectivities(mPrefix + "/Conditions/" + labels[i], mPrefix + "/Xdmf/Conditions/" + labels[i]);

    KRATOS_CATCH("");
}

void XdmfConnectivitiesWriterProcess::CreateXdmfConnectivities(const std::string& rKratosConnectivitiesPath, const std::string& rXdmfConnectivitiesPath) const
{
    KRATOS_TRY;

    Matrix<int> kratos_connectivities, xdmf_connectivities;

    const int block_size = mpFile->GetDataDimensions(rKratosConnectivitiesPath + "/Connectivities")[0];
    mpFile->ReadDataSet(rKratosConnectivitiesPath + "/Connectivities", kratos_connectivities, 0, block_size);

    xdmf_connectivities.resize(kratos_connectivities.size1(), kratos_connectivities.size2(), false);

#pragma omp parallel for
    for (int i = 0; i < block_size; ++i)
    {
        for (unsigned j = 0; j < kratos_connectivities.size2(); ++j)
        {
            const int kratos_id = kratos_connectivities(i, j);
            // Here we expect the ids are one-based indices.
            xdmf_connectivities(i, j) = mKratosToXdmfIdTable(kratos_id - 1);
        }
    }

    WriteInfo info;
    mpFile->WriteDataSet(rXdmfConnectivitiesPath + "/Connectivities", xdmf_connectivities, info);

    int tmp;
    mpFile->ReadAttribute(rKratosConnectivitiesPath, "WorkingSpaceDimension", tmp);
    mpFile->WriteAttribute(rXdmfConnectivitiesPath, "WorkingSpaceDimension", tmp);
    mpFile->ReadAttribute(rKratosConnectivitiesPath, "Dimension", tmp);
    mpFile->WriteAttribute(rXdmfConnectivitiesPath, "Dimension", tmp);
    mpFile->ReadAttribute(rKratosConnectivitiesPath, "NumberOfNodes", tmp);
    mpFile->WriteAttribute(rXdmfConnectivitiesPath, "NumberOfNodes", tmp);

    KRATOS_CATCH("");
}

} // namespace HDF5.
} // namespace Kratos.
