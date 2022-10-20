//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   license: HDF5Application/license.txt
//
//  Main author:     Máté Kelemen
//

// HDF includes
#include "custom_utilities/registry_file.h"
#include "includes/define.h"

// STL includes
#include <filesystem>
#include <cstdio>


namespace Kratos
{


RegistryFile::iterator::iterator(std::shared_ptr<RegistryFile::RegistryAccess>&& rRegistryAccess,
                                 std::size_t EntryIndex)
    : mpRegistryAccess(std::move(rRegistryAccess)),
      mFilePointer(0)
{
    for (std::size_t i=0; i<EntryIndex; ++i) {
        KRATOS_ERROR_IF_NOT(this->SeekNext('\n'))
        << "Index " << EntryIndex << " is out of range";
    }
    this->StoreState();
}


RegistryFile::iterator& RegistryFile::iterator::operator++()
{
    this->RestoreState();
    this->SeekNext('\n');
    this->StoreState();
    return *this;
}


RegistryFile::iterator RegistryFile::iterator::operator++(int)
{
    this->RestoreState();
    auto copy = *this;
    ++(*this);
    return copy;
}


Parameters RegistryFile::iterator::operator*() const
{
    this->RestoreState();
    auto& r_stream = mpRegistryAccess->value().first;
    std::string string;
    std::getline(r_stream, string, '\n');
    return Parameters(string);
}


bool operator==(const RegistryFile::iterator& rLeft, const RegistryFile::iterator& rRight)
{
    return rLeft.mpRegistryAccess.get() == rRight.mpRegistryAccess.get()
           && rLeft.mpRegistryAccess.get() != nullptr
           && rLeft.mFilePointer == rRight.mFilePointer;
}


bool operator!=(const RegistryFile::iterator& rLeft, const RegistryFile::iterator& rRight)
{
    return !(rLeft == rRight);
}


bool operator<(const RegistryFile::iterator& rLeft, const RegistryFile::iterator& rRight)
{
    KRATOS_ERROR_IF((rLeft.mpRegistryAccess.get() != rRight.mpRegistryAccess.get()) || rLeft.mpRegistryAccess.get() == nullptr)
    << "Invalid iterator comparison";
    return rLeft.mFilePointer < rRight.mFilePointer;
}


void RegistryFile::iterator::StoreState()
{
    KRATOS_TRY;
    this->mFilePointer = this->mpRegistryAccess->value().first.tellg();
    KRATOS_CATCH("");
}


void RegistryFile::iterator::RestoreState() const
{
    KRATOS_TRY;
    this->mpRegistryAccess->value().first.seekg(mFilePointer);
    KRATOS_CATCH("");
}


bool RegistryFile::iterator::SeekNext(char Target)
{
    this->RestoreState();
    auto& r_stream = this->mpRegistryAccess->value().first;

    if (r_stream.eof()) {
        return false;
    }

    char c;
    while (r_stream.get(c) && c != Target) {}
    return c == Target;
}


void RegistryFile::iterator::SeekEnd()
{
    this->mpRegistryAccess->value().first.seekg(0, std::ios::end);
    this->StoreState();
}


RegistryFile::RegistryFile()
    : RegistryFile("")
{
}


RegistryFile::RegistryFile(const std::filesystem::path& rRegistryFilePath)
    : RegistryFile(rRegistryFilePath, [](const Model& rModel){return Parameters();})
{
}


RegistryFile::RegistryFile(const std::filesystem::path& rRegistryFilePath,
                           const Extractor& rExtractor)
    : mRegistryFilePath(rRegistryFilePath),
      mExtractor(rExtractor)
{
}


RegistryFile::RegistryFile(const RegistryFile& rOther)
    : RegistryFile()
{
    *this = rOther;
}


RegistryFile& RegistryFile::operator=(const RegistryFile& rOther)
{
    KRATOS_ERROR_IF(this->IsOpen())
    << "Attempt to overwrite RegistryFile at '" << mRegistryFilePath << "' while it's open";

    // Delete the stored file if it is not the same as the incoming one.
    if (std::filesystem::is_regular_file(this->mRegistryFilePath) && this->mRegistryFilePath != rOther.mRegistryFilePath) {
        std::remove(this->mRegistryFilePath.string().c_str());
        this->mRegistryFilePath = rOther.mRegistryFilePath;
    }

    this->mExtractor = rOther.mExtractor;
    this->mpRegistryAccess = rOther.mpRegistryAccess;

    return *this;
}


const std::filesystem::path& RegistryFile::GetFilePath() const noexcept
{
    return this->mRegistryFilePath;
}


void RegistryFile::SetExtractor(Extractor&& rExtractor)
{
    this->mExtractor = std::move(rExtractor);
}


void RegistryFile::SetExtractor(const Extractor& rExtractor)
{
    this->mExtractor = rExtractor;
}


void RegistryFile::Push(const Model& rModel)
{
    KRATOS_TRY;

    auto p_access = this->Open(std::ios::app | std::ios::out);

    const auto output = this->mExtractor(rModel);
    KRATOS_ERROR_IF_NOT(this->IsValidEntry(output))
    << "Extractor returned invalid output: " << output;

    p_access->value().first << output.WriteJsonString() << std::endl;

    KRATOS_CATCH("");
}


void RegistryFile::Clear()
{
    KRATOS_ERROR_IF(this->IsOpen())
    << "Attempt to clear RegistryFile while it's open";

    KRATOS_TRY;
    if (std::filesystem::is_regular_file(this->mRegistryFilePath)) {
        std::remove(this->mRegistryFilePath.string().c_str());
    }
    KRATOS_CATCH("");
}


bool RegistryFile::IsOpen() const noexcept
{
    if (auto p_registry_access = mpRegistryAccess.lock()) {
        if (*p_registry_access) {
            return true;
        }
    }
    return false;
}


RegistryFile::const_iterator RegistryFile::begin() const
{
    std::shared_ptr<RegistryAccess> p_access;
    if (this->IsOpen()) {
        p_access = this->mpRegistryAccess.lock();
    } else {
        p_access = this->Open();
    }

    return const_iterator(std::move(p_access));
}


RegistryFile::const_iterator RegistryFile::end() const
{
    std::shared_ptr<RegistryAccess> p_access;
    if (this->IsOpen()) {
        p_access = this->mpRegistryAccess.lock();
    } else {
        p_access = this->Open();
    }

    auto it = const_iterator(std::move(p_access));
    it.SeekEnd();
    return it;
}


RegistryFile::size_type RegistryFile::size() const
{
    std::shared_ptr<RegistryAccess> p_access;
    if (this->IsOpen()) {
        p_access = this->mpRegistryAccess.lock();
    } else {
        p_access = this->Open();
    }

    auto& r_stream = p_access->value().first;
    char c;
    size_type counter = 0;
    while (r_stream.get(c)) {
        if (c == '\n') ++counter;
    }

    return counter;
}


bool RegistryFile::IsValidEntry(const Parameters& rEntry)
{
    if (!rEntry.IsArray()) {
        return false;
    } else {
        for (const auto& r_item : rEntry) {
            if (!r_item.IsString()) {
                return false;
            } else {
                if (r_item.GetString().find('\n') != std::string::npos) return false;
            }
        }
    }

    return true;
}


std::shared_ptr<RegistryFile::RegistryAccess> RegistryFile::Open(std::ios::openmode OpenMode) const
{
    KRATOS_ERROR_IF(this->IsOpen())
    << "Attempt to repoen RegistryFile at '" << this->mRegistryFilePath  << "'";

    KRATOS_TRY;

    auto p_access = std::make_shared<RegistryAccess>();
    this->mpRegistryAccess = p_access;

    auto pair = std::make_pair(
        std::fstream(this->mRegistryFilePath, OpenMode),
        std::unique_lock<LockObject>(this->mMutex)
    );
    p_access->emplace(std::move(pair));

    return p_access;

    KRATOS_CATCH("");
}


} // namespace Kratos
