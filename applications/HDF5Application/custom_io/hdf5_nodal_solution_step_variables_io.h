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

/** @file hdf5_nodal_solution_step_variables_io.h
 *  @brief Methods for storing and retrieving the storage layout for nodal variables in an HDF5 file.
 */

#if !defined(KRATOS_HDF5_NODAL_SOLUTION_STEP_VARIABLES_IO_H_INCLUDED)
#define KRATOS_HDF5_NODAL_SOLUTION_STEP_VARIABLES_IO_H_INCLUDED

// System includes
#include <string>

// External includes

// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "containers/variables_list.h"

// Application includes
#include "hdf5_application_define.h"

namespace Kratos
{
namespace HDF5
{

class File;

namespace Internals
{
///@addtogroup HDF5Application
///@{

void WriteVariablesList(File& rFile, std::string const& rPrefix, ModelPart const& rModelPart);

void ReadAndAssignVariablesList(File& rFile, std::string const& rPrefix, ModelPart& rModelPart);

void WriteBufferSize(File& rFile, std::string const& rPrefix, int BufferSize);

void ReadAndAssignBufferSize(File& rFile, std::string const& rPrefix, ModelPart& rModelPart);

///@} addtogroup
} // namespace Internals.
} // namespace HDF5.
} // namespace Kratos.

#endif // KRATOS_HDF5_NODAL_SOLUTION_STEP_VARIABLES_IO_H_INCLUDED defined
