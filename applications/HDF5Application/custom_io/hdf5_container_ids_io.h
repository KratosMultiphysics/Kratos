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

#if !defined(KRATOS_HDF5_CONTAINER_IDS_IO_H_INCLUDED)
#define KRATOS_HDF5_CONTAINER_IDS_IO_H_INCLUDED

// System includes
#include <string>
#include <vector>

// External includes

// Project includes
#include "includes/define.h"
#include "utilities/openmp_utils.h"

// Application includes
#include "hdf5_application_define.h"
#include "custom_io/hdf5_file.h"

namespace Kratos
{
namespace HDF5
{
namespace Internals
{
///@addtogroup HDF5Application
///@{

template <typename TContainer>
void WriteContainerIds(File& rFile, std::string const& rPath, TContainerType const& rContainer, WriteInfo& rInfo)
{
    KRATOS_TRY;
    Vector<int> ids;
    ids.resize(rContainer.size());
    #pragma omp parallel for
    for (std::size_t i = 0; i < rContainer.size(); ++i)
    {
        const auto it = rContainer.begin() + i;
        ids[i] = it->Id();
    }
    rFile.WriteDataSet(rPath, ids, rInfo);
    KRATOS_CATCH("");
}

///@} addtogroup
} // namespace Internals.
} // namespace HDF5.
} // namespace Kratos.

#endif // KRATOS_HDF5_CONTAINER_IDS_IO_H_INCLUDED defined