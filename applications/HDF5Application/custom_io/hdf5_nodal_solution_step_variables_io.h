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
//

/** @file hdf5_nodal_solution_step_variables_io.h
 *  @brief Methods for storing and retrieving the storage layout for nodal variables in an HDF5 file.
 */

#pragma once

// System includes
#include <string>

// External includes

// Project includes
#include "includes/model_part.h"

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

void WriteVariablesList(
    File& rFile,
    const std::string& rPrefix,
    const ModelPart& rModelPart);

void ReadAndAssignVariablesList(
    File& rFile,
    const std::string& rPrefix,
    ModelPart& rModelPart);

void WriteBufferSize(
    File& rFile,
    const std::string& rPrefix,
    const int BufferSize);

void ReadAndAssignBufferSize(
    File& rFile,
    const std::string& rPrefix,
    ModelPart& rModelPart);

///@} addtogroup
} // namespace Internals.
} // namespace HDF5.
} // namespace Kratos.