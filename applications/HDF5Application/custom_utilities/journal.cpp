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
#include "custom_utilities/journal.h"
#include "includes/define.h"

// STL includes
#include <filesystem>
#include <cstdio>


namespace Kratos
{


Journal::iterator::iterator(std::shared_ptr<Journal::FileAccess>&& rFileAccess,
                                 std::size_t EntryIndex)
    : mpFileAccess(std::move(rFileAccess)),
      mFilePointer(0)
{
    for (std::size_t i=0; i<EntryIndex; ++i) {
        KRATOS_ERROR_IF_NOT(this->SeekNext('\n'))
        << "Index " << EntryIndex << " is out of range";
    }
    this->StoreState();
}


Journal::iterator& Journal::iterator::operator++()
{
    this->RestoreState();
    this->SeekNext('\n');
    this->StoreState();
    return *this;
}


Journal::iterator Journal::iterator::operator++(int)
{
    this->RestoreState();
    auto copy = *this;
    ++(*this);
    return copy;
}


Parameters Journal::iterator::operator*() const
{
    this->RestoreState();
    auto& r_stream = mpFileAccess->value().first;
    std::string string;
    std::getline(r_stream, string, '\n');
    return Parameters(string);
}


bool operator==(const Journal::iterator& rLeft, const Journal::iterator& rRight)
{
    return rLeft.mpFileAccess.get() == rRight.mpFileAccess.get()
           && rLeft.mpFileAccess.get() != nullptr
           && rLeft.mFilePointer == rRight.mFilePointer;
}


bool operator!=(const Journal::iterator& rLeft, const Journal::iterator& rRight)
{
    return !(rLeft == rRight);
}


bool operator<(const Journal::iterator& rLeft, const Journal::iterator& rRight)
{
    KRATOS_ERROR_IF((rLeft.mpFileAccess.get() != rRight.mpFileAccess.get()) || rLeft.mpFileAccess.get() == nullptr)
    << "Invalid iterator comparison";
    return rLeft.mFilePointer < rRight.mFilePointer;
}


void Journal::iterator::StoreState()
{
    KRATOS_TRY;
    this->mFilePointer = this->mpFileAccess->value().first.tellg();
    KRATOS_CATCH("");
}


void Journal::iterator::RestoreState() const
{
    KRATOS_TRY;
    this->mpFileAccess->value().first.seekg(mFilePointer);
    KRATOS_CATCH("");
}


bool Journal::iterator::SeekNext(char Target)
{
    this->RestoreState();
    auto& r_stream = this->mpFileAccess->value().first;

    if (r_stream.eof()) {
        return false;
    }

    char c;
    while (r_stream.get(c) && c != Target) {}
    return c == Target;
}


void Journal::iterator::SeekEnd()
{
    this->mpFileAccess->value().first.seekg(0, std::ios::end);
    this->StoreState();
}


Journal::Journal()
    : Journal("")
{
}


Journal::Journal(const std::filesystem::path& rJournalPath)
    : Journal(rJournalPath, [](const Model& rModel){return Parameters();})
{
}


Journal::Journal(const std::filesystem::path& rJournalPath,
                           const Extractor& rExtractor)
    : mJournalPath(rJournalPath),
      mExtractor(rExtractor)
{
}


Journal::Journal(const Journal& rOther)
    : Journal()
{
    *this = rOther;
}


Journal& Journal::operator=(const Journal& rOther)
{
    KRATOS_ERROR_IF(this->IsOpen())
    << "Attempt to overwrite Journal at '" << mJournalPath << "' while it's open";

    // Delete the stored file if it is not the same as the incoming one.
    if (std::filesystem::is_regular_file(this->mJournalPath) && this->mJournalPath != rOther.mJournalPath) {
        std::remove(this->mJournalPath.string().c_str());
        this->mJournalPath = rOther.mJournalPath;
    }

    this->mExtractor = rOther.mExtractor;
    this->mpFileAccess = rOther.mpFileAccess;

    return *this;
}


const std::filesystem::path& Journal::GetFilePath() const noexcept
{
    return this->mJournalPath;
}


void Journal::SetExtractor(Extractor&& rExtractor)
{
    this->mExtractor = std::move(rExtractor);
}


void Journal::SetExtractor(const Extractor& rExtractor)
{
    this->mExtractor = rExtractor;
}


void Journal::Push(const Model& rModel)
{
    KRATOS_TRY;

    auto p_access = this->Open(std::ios::app | std::ios::out);

    const auto output = this->mExtractor(rModel);
    KRATOS_ERROR_IF_NOT(this->IsValidEntry(output))
    << "Extractor returned invalid output: " << output;

    p_access->value().first << output.WriteJsonString() << std::endl;

    KRATOS_CATCH("");
}


void Journal::Clear()
{
    KRATOS_ERROR_IF(this->IsOpen())
    << "Attempt to clear Journal while it's open";

    KRATOS_TRY;
    if (std::filesystem::is_regular_file(this->mJournalPath)) {
        std::remove(this->mJournalPath.string().c_str());
    }
    KRATOS_CATCH("");
}


bool Journal::IsOpen() const noexcept
{
    if (auto p_registry_access = mpFileAccess.lock()) {
        if (*p_registry_access) {
            return true;
        }
    }
    return false;
}


Journal::const_iterator Journal::begin() const
{
    std::shared_ptr<FileAccess> p_access;
    if (this->IsOpen()) {
        p_access = this->mpFileAccess.lock();
    } else {
        p_access = this->Open();
    }

    return const_iterator(std::move(p_access));
}


Journal::const_iterator Journal::end() const
{
    std::shared_ptr<FileAccess> p_access;
    if (this->IsOpen()) {
        p_access = this->mpFileAccess.lock();
    } else {
        p_access = this->Open();
    }

    auto it = const_iterator(std::move(p_access));
    it.SeekEnd();
    return it;
}


Journal::size_type Journal::size() const
{
    std::shared_ptr<FileAccess> p_access;
    if (this->IsOpen()) {
        p_access = this->mpFileAccess.lock();
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


bool Journal::IsValidEntry(const Parameters& rEntry)
{
    if (!rEntry.IsArray()) {
        return false;
    } else {
        for (const auto& r_item : rEntry) {
            if (r_item.IsString()) {
                if (r_item.GetString().find('\n') != std::string::npos) return false;
            } else {
                if (!(r_item.IsNumber() || r_item.IsBool())) {
                    return false;
                }
            }
        }
    }

    return true;
}


std::shared_ptr<Journal::FileAccess> Journal::Open(std::ios::openmode OpenMode) const
{
    KRATOS_ERROR_IF(this->IsOpen())
    << "Attempt to repoen Journal at '" << this->mJournalPath  << "'";

    KRATOS_TRY;

    auto p_access = std::make_shared<FileAccess>();
    this->mpFileAccess = p_access;

    auto pair = std::make_pair(
        std::fstream(this->mJournalPath, OpenMode),
        std::unique_lock<LockObject>(this->mMutex)
    );
    p_access->emplace(std::move(pair));

    return p_access;

    KRATOS_CATCH("");
}


} // namespace Kratos
