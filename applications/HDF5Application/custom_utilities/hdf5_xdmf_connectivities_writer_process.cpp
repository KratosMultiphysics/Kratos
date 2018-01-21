#include "custom_utilities/hdf5_xdmf_connectivities_writer_process.h"

#include "custom_io/hdf5_file_serial.h"
#include "utilities/openmp_utils.h"

namespace Kratos
{
namespace HDF5
{
XdmfConnectivitiesWriterProcess::XdmfConnectivitiesWriterProcess(std::string FileName, std::string Prefix)
{
    KRATOS_TRY;

    Parameters file_params(R"(
            {
                "file_name" : "",
                "file_access_mode": "read_write"
            })");
    file_params["file_name"].SetString(FileName);
    mpFile = FileSerial::Pointer(new FileSerial(file_params));
    mPrefix = Prefix;
    std::string node_ids_path = mPrefix + "/Nodes/Local/Ids";

    KRATOS_ERROR_IF(mpFile->HasPath(node_ids_path) == false)
        << "Path \"" << node_ids_path << "\" was not found." << std::endl;

    Vector<int> node_ids;
    const unsigned num_points = mpFile->GetDataDimensions(node_ids_path)[0];
    mpFile->ReadDataSet(node_ids_path, node_ids, 0, num_points);

    // Set the parametric coordinate ids.
    mKratosToXdmfIdTable.resize(num_points);
    const int num_threads = OpenMPUtils::GetNumThreads();
    OpenMPUtils::PartitionVector partition;
    OpenMPUtils::DivideInPartitions(num_points, num_threads, partition);
#pragma omp parallel
    {
        const int thread_id = OpenMPUtils::ThisThread();
        for (auto i = partition[thread_id]; i < partition[thread_id + 1]; ++i)
        {
            // Here we expect the ids are one-based indices.
            mKratosToXdmfIdTable(node_ids[i] - 1) = i;
        }
    }

    KRATOS_CATCH("");
}

void XdmfConnectivitiesWriterProcess::Execute()
{
    KRATOS_TRY;

    std::vector<std::string> labels;
    Matrix<int> connectivities;
    KRATOS_ERROR_IF(mpFile->HasPath(mPrefix + "/Xdmf/Elements")) << "Path \"" << mPrefix + "/Xdmf/Elements\" exists." << std::endl;
    mpFile->GetLinkNames(mPrefix + "/Elements", labels);
    for (unsigned i = 0; i < labels.size(); ++i)
        CreateXdmfConnectivities(mPrefix + "/Elements/" + labels[i], mPrefix + "/Xdmf/Elements/" + labels[i]);

    KRATOS_ERROR_IF(mpFile->HasPath(mPrefix + "/Xdmf/Conditions")) << "Path \"" << mPrefix + "/Xdmf/Conditions\" exists." << std::endl;
    mpFile->GetLinkNames(mPrefix + "/Conditions", labels);
    for (unsigned i = 0; i < labels.size(); ++i)
        CreateXdmfConnectivities(mPrefix + "/Conditions/" + labels[i], mPrefix + "/Xdmf/Conditions/" + labels[i]);

    KRATOS_CATCH("");
}

void XdmfConnectivitiesWriterProcess::CreateXdmfConnectivities(std::string KratosConnectivitiesPath, std::string XdmfConnectivitiesPath) const
{
    KRATOS_TRY;

    Matrix<int> kratos_connectivities, xdmf_connectivities;

    const unsigned block_size = mpFile->GetDataDimensions(KratosConnectivitiesPath + "/Connectivities")[0];
    mpFile->ReadDataSet(KratosConnectivitiesPath + "/Connectivities", kratos_connectivities, 0, block_size);

    xdmf_connectivities.resize(kratos_connectivities.size1(), kratos_connectivities.size2(), false);

    const int num_threads = OpenMPUtils::GetNumThreads();
    OpenMPUtils::PartitionVector partition;
    OpenMPUtils::DivideInPartitions(block_size, num_threads, partition);
#pragma omp parallel
    {
        const int thread_id = OpenMPUtils::ThisThread();
        for (auto i = partition[thread_id]; i < partition[thread_id + 1]; ++i)
        {
            for (unsigned j = 0; j < kratos_connectivities.size2(); ++j)
            {
                int kratos_id = kratos_connectivities(i, j);
                // Here we expect the ids are one-based indices.
                xdmf_connectivities(i, j) = mKratosToXdmfIdTable(kratos_id - 1);
            }
        }
    }

    mpFile->WriteDataSet(XdmfConnectivitiesPath + "/Connectivities", xdmf_connectivities);

    int tmp;
    mpFile->ReadAttribute(KratosConnectivitiesPath, "WorkingSpaceDimension", tmp);
    mpFile->WriteAttribute(XdmfConnectivitiesPath, "WorkingSpaceDimension", tmp);
    mpFile->ReadAttribute(KratosConnectivitiesPath, "Dimension", tmp);
    mpFile->WriteAttribute(XdmfConnectivitiesPath, "Dimension", tmp);
    mpFile->ReadAttribute(KratosConnectivitiesPath, "NumberOfNodes", tmp);
    mpFile->WriteAttribute(XdmfConnectivitiesPath, "NumberOfNodes", tmp);

    KRATOS_CATCH("");
}

} // namespace HDF5.
} // namespace Kratos.
