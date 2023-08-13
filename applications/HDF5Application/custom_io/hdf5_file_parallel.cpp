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
//                  Suneth Warnakulasuriya
//

// System includes
#include <numeric>

// Project includes
#include "includes/kratos_parameters.h"
#include "utilities/builtin_timer.h"
#include "input_output/logger.h"
#include "utilities/data_type_traits.h"

// Application includes
#include "custom_utilities/h5_data_type_traits.h"

// Include base h
#include "hdf5_file_parallel.h"

namespace Kratos
{
namespace HDF5
{

FileParallel::FileParallel(Parameters& rSettings) : File(rSettings)
{
}

FileParallel::FileParallel(
    const DataCommunicator& rDataCommunicator,
    Parameters Settings)
    : File(rDataCommunicator, Settings)
{
}

} // namespace HDF5.
} // namespace Kratos.