// KRATOS  __  __          _    _                _ _           _   _
//        |  \/  | ___  __| |  / \   _ __  _ __ | (_) ___ __ _| |_(_) ___  _ ___
//        | |\/| |/ _ \/ _` | / _ \ | '_ \| '_ \| | |/ __/ _` | __| |/ _ \| '_  |
//        | |  | |  __/ (_| |/ ___ \| |_) | |_) | | | (_| (_| | |_| | (_) | | | |
//        |_|  |_|\___|\__,_/_/   \_\ .__/| .__/|_|_|\___\__,_|\__|_|\___/|_| |_|
//                                  |_|   |_|
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Philipp Bucher (https://github.com/philbucher)
//

#pragma once

// System includes
#include <cstdint>
#include <string>
#include <vector>

// External includes
#include "hdf5.h"

// Project includes

namespace Kratos {

/// Minimal, dependency-free helpers around the raw HDF5 C API.
/**
 * These wrap the handful of HDF5 calls (H5G/H5D/H5A/H5L) needed by MedModelPartIO to
 * read and write the MED-on-HDF5 group/dataset layout by hand, without linking the
 * external "med" library.
 *
 * HDF5Application's own "File" class was deliberately not reused here: its path
 * validation only accepts alphanumeric/underscore/parenthesis characters, which
 * rejects MED's own mandatory computation-step group name (e.g.
 * "-0000000000000000001-0000000000000000001"), and it has no notion of the
 * fixed-width string-array dataset MED uses for SubModelPart/group names.
 */
namespace MedHdf5Utilities {

/// Returns whether every segment of "rPath" exists, without erroring on missing
/// intermediate groups (unlike a bare H5Lexists on a multi-component path).
bool LinkExists(
    hid_t FileId,
    const std::string& rPath);

/// Returns the names of the direct children (links) of the group at "rGroupPath".
std::vector<std::string> ListGroupChildren(
    hid_t FileId,
    const std::string& rGroupPath);

/// Creates the group at "rPath", creating missing intermediate groups as needed.
/// The caller is responsible for closing the returned handle.
hid_t CreateGroup(
    hid_t FileId,
    const std::string& rPath);

/// Same as "CreateGroup", but closes the group immediately.
void CreateAndCloseGroup(
    hid_t FileId,
    const std::string& rPath);

int64_t ReadInt64Attribute(
    hid_t FileId,
    const std::string& rObjectPath,
    const std::string& rName);

void WriteInt64Attribute(
    hid_t FileId,
    const std::string& rObjectPath,
    const std::string& rName,
    int64_t Value);

void WriteFloat64Attribute(
    hid_t FileId,
    const std::string& rObjectPath,
    const std::string& rName,
    double Value);

/// Writes a scalar string attribute, sized exactly to the content (plus null
/// terminator), matching how the libmed C library sizes its string attributes.
void WriteStringAttribute(
    hid_t FileId,
    const std::string& rObjectPath,
    const std::string& rName,
    const std::string& rValue);

std::vector<int64_t> ReadInt64Dataset(
    hid_t FileId,
    const std::string& rPath);

std::vector<double> ReadFloat64Dataset(
    hid_t FileId,
    const std::string& rPath);

/// Writing zero entries is a no-op (no dataset is created), matching the observed
/// behavior of libmed's family-number/global-numbering write functions.
void WriteInt64Dataset(
    hid_t FileId,
    const std::string& rPath,
    const std::vector<int64_t>& rData);

void WriteFloat64Dataset(
    hid_t FileId,
    const std::string& rPath,
    const std::vector<double>& rData);

/// Reads a dataset whose element type is a fixed-width, null-terminated (or
/// null-padded) array of "Width" chars - the layout MED uses for its
/// SubModelPart/group name lists.
std::vector<std::string> ReadFixedWidthStringDataset(
    hid_t FileId,
    const std::string& rPath,
    std::size_t Width);

/// Writes a dataset whose element type is a fixed-width array of "Width" chars, one
/// entry per name in "rNames". Throws if any name is at least "Width" characters long.
void WriteFixedWidthStringDataset(
    hid_t FileId,
    const std::string& rPath,
    const std::vector<std::string>& rNames,
    std::size_t Width);

} // namespace MedHdf5Utilities
} // namespace Kratos
