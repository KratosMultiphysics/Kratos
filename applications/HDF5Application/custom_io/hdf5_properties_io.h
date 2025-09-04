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

/** @file hdf5_properties_io.h
 *  @brief Methods for storing and retrieving a properties in an HDF5 file.
 */

#pragma once

// System includes
#include <string>

// External includes

// Project includes
#include "includes/define.h"

// Application includes
#include "custom_io/hdf5_file.h"
#include "hdf5_application_define.h"

namespace Kratos::HDF5::Internals
{
///@addtogroup HDF5Application
///@{

void ReadProperties(
    File& rFile,
    const std::string& rPrefix,
    PropertiesContainerType& rProperties);

void WriteProperties(
    File& rFile,
    const std::string& rPrefix,
    const PropertiesContainerType& rProperties);

///@}
} // namespace Kratos::HDF5::Internals