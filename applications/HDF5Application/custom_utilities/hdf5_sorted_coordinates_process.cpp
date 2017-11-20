#include "custom_utilities/hdf5_sorted_coordinates_process.h"

#include "utilities/openmp_utils.h"
#include "containers/array_1d.h"
#include "hdf5_application_define.h"
#include "custom_io/hdf5_file_serial.h"


namespace Kratos
{
namespace HDF5
{
SortedCoordinatesProcess::SortedCoordinatesProcess(std::string FileName, std::string Prefix)
    : Process(), mFileName(FileName), mPrefix(Prefix)
{
}

void SortedCoordinatesProcess::Execute()
{
    KRATOS_TRY;

    Parameters file_params(R"(
            {
                "file_name" : "",
                "file_access_mode": "read_write"
            })");
    file_params["file_name"].SetString(mFileName);
    FileSerial hdf5_file(file_params);

    KRATOS_ERROR_IF(hdf5_file.HasPath(mPrefix) == false)
        << "Path \"" << mPrefix << "\" was not found." << std::endl;

    if (hdf5_file.HasPath(mPrefix + "/SortedCoordinates"))
        return;

    Vector<int> ids, pcs;
    Vector<array_1d<double, 3>> unsorted_coords, sorted_coords;

    const unsigned num_points = hdf5_file.GetDataDimensions(mPrefix + "/Ids")[0];
    hdf5_file.ReadDataSet(mPrefix + "/Ids", ids, 0, num_points);
    hdf5_file.ReadDataSet(mPrefix + "/Coordinates", unsorted_coords, 0, num_points);

    // Set the parametric coordinate ids.
    pcs.resize(num_points);
    sorted_coords.resize(num_points);
    const int num_threads = OpenMPUtils::GetNumThreads();
    OpenMPUtils::PartitionVector partition;
    OpenMPUtils::DivideInPartitions(num_points, num_threads, partition);
#pragma omp parallel
    {
        const int thread_id = OpenMPUtils::ThisThread();
        for (auto i = partition[thread_id]; i < partition[thread_id + 1]; ++i)
        {
            pcs[ids[i] - 1] = i; // Here we expect the ids are one-based indices.
        }

        for (auto i = partition[thread_id]; i < partition[thread_id + 1]; ++i)
        {
            sorted_coords[i] = unsorted_coords[pcs[i]];
        }
    }

    hdf5_file.WriteDataSet(mPrefix + "/SortedCoordinates", sorted_coords);

    KRATOS_CATCH("");
}

} // namespace HDF5.
} // namespace Kratos.
