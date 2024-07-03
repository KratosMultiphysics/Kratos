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


// Project includes
#include "testing/scoped_file.h"


namespace Kratos::Testing {


ScopedEntry::ScopedEntry(const std::filesystem::path& rPath)
    : mPath(rPath)
{}


ScopedEntry::~ScopedEntry()
{
    std::filesystem::is_directory(mPath) ? std::filesystem::remove_all(mPath) : std::filesystem::remove(mPath);
}


ScopedEntry::operator const std::filesystem::path& () const
{
    return mPath;
}


ScopedDirectory::ScopedDirectory(const std::filesystem::path& rPath)
    : ScopedEntry(rPath)
{
    std::filesystem::create_directories(rPath);
}


ScopedFile::ScopedFile(const std::filesystem::path& rPath)
    : ScopedEntry(rPath),
      mStream(rPath)
{
}


ScopedFile::~ScopedFile()
{
    mStream.flush();
    mStream.close();
}


ScopedSymlink::ScopedSymlink(const std::filesystem::path& rSymlinkPath, const std::filesystem::path& rTargetPath)
    : ScopedEntry(rSymlinkPath)
{
    std::filesystem::exists(rTargetPath) && std::filesystem::is_directory(rTargetPath) ? std::filesystem::create_directory_symlink(rTargetPath, rSymlinkPath) : std::filesystem::create_symlink(rTargetPath, rSymlinkPath);
}


} // namespace Kratos::Testing
