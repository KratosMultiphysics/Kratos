//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Máté Kelemen
//


#pragma once

// STL includes
#include <filesystem>
#include <fstream>

#include "includes/kratos_export_api.h"

namespace Kratos::Testing {


/**
 *  A base class for filesystem entities that get created upon construction
 *  and get deleted when the object gets destroyed.
 */
class KRATOS_API(KRATOS_TEST_UTILS) ScopedEntry
{
public:
    ScopedEntry(const std::filesystem::path& rPath);

    ScopedEntry(ScopedEntry&& rOther) = default;

    ScopedEntry(const ScopedEntry& rOther) = delete;

    ScopedEntry& operator=(ScopedEntry&& rOther) = delete;

    ScopedEntry& operator=(const ScopedEntry& rOther) = delete;

    virtual ~ScopedEntry();

    operator const std::filesystem::path& () const;

private:
    const std::filesystem::path mPath;
}; // class ScopedEntry


/// Class representing a directory that follows RAII.
struct KRATOS_API(KRATOS_TEST_UTILS) ScopedDirectory final : public ScopedEntry
{
    ScopedDirectory(const std::filesystem::path& rPath);
}; // struct ScopedDirectory


/// Class representing a file that follows RAII.
class KRATOS_API(KRATOS_TEST_UTILS) ScopedFile final : public ScopedEntry
{
public:
    ScopedFile(const std::filesystem::path& rPath);

    ~ScopedFile() override;

    template <class T>
    friend ScopedFile& operator<<(ScopedFile& rFile, const T& rContent)
    {
        rFile.mStream << rContent;
        rFile.mStream.flush();
        return rFile;
    }

private:
    std::ofstream mStream;
}; // class ScopedFile


/// Class representing a symlink that follows RAII.
struct KRATOS_API(KRATOS_TEST_UTILS) ScopedSymlink final : public ScopedEntry
{
    ScopedSymlink(const std::filesystem::path& rSymlinkPath, const std::filesystem::path& rTargetPath);
}; // struct ScopedSymlink


} // namespace Kratos::Testing
