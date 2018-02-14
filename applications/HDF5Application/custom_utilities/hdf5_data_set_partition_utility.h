//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 license: HDF5Application/license.txt
//
//  Main author:    Michael Andre, https://github.com/msandre
//

#if !defined(KRATOS_HDF5_DATA_SET_PARTITION_UTILITY_H_INCLUDED)
#define KRATOS_HDF5_DATA_SET_PARTITION_UTILITY_H_INCLUDED

// System includes
#include <string>
#include <tuple>

// External includes

// Project includes
#include "includes/define.h"

// Application includes
#include "hdf5_application_define.h"
#include "custom_io/hdf5_file.h"

namespace Kratos
{
namespace HDF5
{
///@addtogroup HDF5Application
///@{

/// Write the start and end indices of data blocks (by process rank).
/**
 * Performs collective write.
 */
void WritePartitionTable(File& rFile, std::string const& rPath, WriteInfo const& rInfo);

/// Write a user-defined partition table of start and end indices (by process rank).
/**
 * Performs independent write.
 */
void WritePartitionTableIndependent(File& rFile, std::string const& rPath, Vector<int> const& rPartition);

// Check if a path has a data set partition.
bool HasPartitionTable(File& rFile, std::string const& rPath);

// Get the start index and block size from an existing partition for this PID.
std::tuple<unsigned, unsigned> StartIndexAndBlockSize(File& rFile, std::string const& rPath);

///@} addtogroup
} // namespace HDF5.
} // namespace Kratos.

#endif // KRATOS_HDF5_DATA_SET_PARTITION_UTILITY_H_INCLUDED defined
