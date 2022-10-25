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

// External includes
#ifdef KRATOS_SMP_OPENMP
#include <omp.h>
#endif

// HDF includes
#include "custom_utilities/journal.h"
#include "includes/define.h"

// Core includes
#include "input_output/logger.h"

// STL includes
#include <filesystem>
#include <sstream>


namespace Kratos
{


JournalBase::JournalBase()
    : JournalBase("")
{
}


JournalBase::JournalBase(const std::filesystem::path& rJournalPath)
    : JournalBase(rJournalPath, [](const Model& rModel){return "";})
{
}


JournalBase::JournalBase(const std::filesystem::path& rJournalPath,
                         const Extractor& rExtractor)
    : mJournalPath(rJournalPath),
      mExtractor(rExtractor, ThreadID())
{
}


JournalBase::JournalBase(const JournalBase& rOther)
    : JournalBase()
{
    *this = rOther;
}


JournalBase& JournalBase::operator=(const JournalBase& rOther)
{
    KRATOS_ERROR_IF(this->IsOpen())
    << "Attempt to overwrite JournalBase at '" << mJournalPath << "' while it's open";

    // Delete the stored file if it is not the same as the incoming one.
    if (std::filesystem::is_regular_file(this->mJournalPath) && this->mJournalPath != rOther.mJournalPath) {
        std::remove(this->mJournalPath.string().c_str());
        this->mJournalPath = rOther.mJournalPath;
    }

    this->mExtractor = rOther.mExtractor;
    this->mpFileAccess = rOther.mpFileAccess;

    return *this;
}


const std::filesystem::path& JournalBase::GetFilePath() const noexcept
{
    return this->mJournalPath;
}


void JournalBase::SetExtractor(Extractor&& rExtractor)
{
    this->mExtractor = std::make_pair(std::move(rExtractor), ThreadID());
}


void JournalBase::SetExtractor(const Extractor& rExtractor)
{
    this->mExtractor = std::make_pair(rExtractor, ThreadID());
}


void JournalBase::Push(const Model& rModel)
{
    KRATOS_TRY;

    // The extractor can only be invoked by the thread that set it.
    if (mExtractor.second != ThreadID()) {
        std::stringstream message;
        message << "Journal::Push can only be invoked by the thread that set the extractor. "
        << "The extractor was set by thread " << mExtractor.second << " but invoked by " << ThreadID() << '\n';
        KRATOS_WARNING(message.str());
        return;
    }

    auto p_access = this->Open(std::ios::app | std::ios::out);

    const auto output = this->mExtractor.first(rModel);
    KRATOS_ERROR_IF_NOT(this->IsValidEntry(output))
    << "Extractor returned invalid output: " << output;

    p_access->value().first << output << std::endl;

    KRATOS_CATCH("");
}


void JournalBase::Clear()
{
    KRATOS_ERROR_IF(this->IsOpen())
    << "Attempt to clear JournalBase while it's open";

    KRATOS_TRY;
    if (std::filesystem::is_regular_file(this->mJournalPath)) {
        std::remove(this->mJournalPath.string().c_str());
    }
    KRATOS_CATCH("");
}


bool JournalBase::IsOpen() const noexcept
{
    if (auto p_registry_access = mpFileAccess.lock()) {
        if (*p_registry_access) {
            return true;
        }
    }
    return false;
}


JournalBase::const_iterator JournalBase::begin() const
{
    std::shared_ptr<iterator::FileAccess> p_access;
    if (this->IsOpen()) {
        p_access = this->mpFileAccess.lock();
    } else {
        p_access = this->Open();
    }

    return const_iterator(std::move(p_access), '\n');
}


JournalBase::const_iterator JournalBase::end() const
{
    std::shared_ptr<iterator::FileAccess> p_access;
    if (this->IsOpen()) {
        p_access = this->mpFileAccess.lock();
    } else {
        p_access = this->Open();
    }

    auto it = const_iterator(std::move(p_access), '\n');
    it.SeekEOF();
    return it;
}


JournalBase::size_type JournalBase::size() const
{
    std::shared_ptr<iterator::FileAccess> p_access;
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


JournalBase::JournalBase::ThreadID::ThreadID()
    :
    #ifdef KRATOS_SMP_OPENMP
        mID(omp_get_thread_num())
    #elif defined(KRATOS_CMP_CXX11)
        mID(std::this_thread::get_id())
    #endif
{
}


bool JournalBase::IsValidEntry(const value_type& rEntry)
{
    return rEntry.find('\n') == rEntry.npos;
}


std::shared_ptr<JournalBase::iterator::FileAccess> JournalBase::Open(std::ios::openmode OpenMode) const
{
    KRATOS_ERROR_IF(this->IsOpen())
    << "Attempt to repoen JournalBase at '" << this->mJournalPath  << "'";

    KRATOS_TRY;

    auto p_access = std::make_shared<iterator::FileAccess>();
    this->mpFileAccess = p_access;

    p_access->emplace();
    p_access->value().first = std::fstream(this->mJournalPath, OpenMode);

    return p_access;

    KRATOS_CATCH("");
}


Journal::iterator::iterator(JournalBase::iterator&& rWrapped)
    : mWrapped(std::move(rWrapped))
{
}


Journal::iterator& Journal::iterator::operator++()
{
    ++mWrapped;
    return *this;
}


Journal::iterator Journal::iterator::operator++(int)
{
    iterator copy = *this;
    ++(*this);
    return copy;
}


Journal::iterator::value_type Journal::iterator::operator*() const
{
    return value_type(*mWrapped);
}


bool operator==(const Journal::iterator& rLeft, const Journal::iterator& rRight)
{
    return rLeft.mWrapped == rRight.mWrapped;
}


bool operator!=(const Journal::iterator& rLeft, const Journal::iterator& rRight)
{
    return rLeft.mWrapped != rRight.mWrapped;
}


Journal::Journal()
    : Journal("")
{
}


Journal::Journal(const std::filesystem::path& rJournalPath)
    : Journal(rJournalPath, [](const Model& rModel){return Parameters();})
{
}


Journal::Journal(const std::filesystem::path& rJournalPath, const Extractor& rExtractor)
    : mBase(rJournalPath)
{
    this->SetExtractor(rExtractor);
}


const std::filesystem::path& Journal::GetFilePath() const noexcept
{
    return mBase.GetFilePath();
}


void Journal::SetExtractor(Extractor&& rExtractor)
{
    mBase.SetExtractor([rExtractor = std::move(rExtractor)](const Model& rModel){
        return rExtractor(rModel).WriteJsonString();
    });
}


void Journal::SetExtractor(const Extractor& rExtractor)
{
    this->SetExtractor(Extractor(rExtractor));
}


void Journal::Push(const Model& rModel)
{
    this->mBase.Push(rModel);
}


void Journal::Clear()
{
    mBase.Clear();
}


bool Journal::IsOpen() const noexcept
{
    return mBase.IsOpen();
}


Journal::const_iterator Journal::begin() const
{
    return const_iterator(mBase.begin());
}


Journal::const_iterator Journal::end() const
{
    return const_iterator(mBase.end());
}


Journal::size_type Journal::size() const
{
    return mBase.size();
}


} // namespace Kratos
