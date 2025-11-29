/// @author Máté Kelemen

#pragma once

// --- WRApplication Includes ---
#include "custom_utilities/file_content_iterator.hpp"

// --- Core Includes ---
#include "includes/kratos_parameters.h"
#include "containers/model.h"
#include "includes/smart_pointers.h"

// --- STL Includes ---
#include <filesystem>
#include <functional>

#ifdef KRATOS_SMP_CXX11
#include <thread>
#endif


namespace Kratos {


/// @addtogroup HDF5Application
/// @{


/**
 *  @brief A class for keeping administrational text data throughout an analysis.
 *  @details This class is primarily meant to keep track of output files generated during a simulation,
 *           though extra flexibility is provided to extend its purpose as necessary. An associated
 *           text file can be written to via @ref JournalBase::Push. On each push, a @ref Model is
 *           taken as input, parsed, and a string consisting of a single line is appended to the file.
 *           Generating the output string from the input @ref Model is the job of a functor that can be
 *           set in the constructor, or in @ref JournalBase::SetExtractor.
 *  @warning The associated file must be empty or end with an empty line.
 */
class JournalBase
{
public:
    using iterator = FileStringIterator;

    using const_iterator = iterator;

    using value_type = iterator::value_type;

    using size_type = std::size_t;

    using Extractor = std::function<value_type(const Model&)>;

    KRATOS_CLASS_POINTER_DEFINITION(JournalBase);

public:
    /// @brief Construct the object with an empty file name (invalid).
    JournalBase();

    /// @brief Construct a registry given an associated file path and a no-op extractor.
    /// @note Should the file already exist, it will be appended to instead of overwritten.
    JournalBase(const std::filesystem::path& rJournalPath);

    /// @brief Construct a registry given an associated file path and extractor.
    /// @note Should the file already exist, it will be appended to instead of overwritten.
    JournalBase(const std::filesystem::path& rJournalPath,
                const Extractor& rExtractor);

    /// @brief The move constructor is deleted because it's ambiguous what should happen to the associated file.
    JournalBase(JournalBase&& rOther) = delete;

    /// @brief Copy constructor.
    /// @warning Deletes the associated file if the incoming instance has another one.
    /// @warning The associated file must not be open.
    explicit JournalBase(const JournalBase& rOther);

    /// @brief The move assignment operator is deleted because it's ambiguous what should happen to the associated file.
    JournalBase& operator=(JournalBase&& rOther) = delete;

    /// @brief Copy assignment operator.
    /// @warning Deletes the associated file if the incoming instance has another one.
    /// @warning The associated file must not be open.
    JournalBase& operator=(const JournalBase& rOther);

    /// @brief Get the path to the associated file.
    const std::filesystem::path& GetFilePath() const noexcept;

    /**
     *  @brief Set the extractor handling the conversion from @ref Model to a string.
     *  @details The extractor must generate a string without any line breaks.
     *  @warning Once @a SetExtractor is called, the ID of the caller thread is
     *           stored and only that thread is allowed to invoke @ref JournalBase::Push.
     *           This is necessary to support python function objects (avoid hanging
     *           at the global interpreter lock).
     */
    void SetExtractor(Extractor&& rExtractor);

    /// @copydoc JournalBase::SetExtractor(const Extractor&)
    void SetExtractor(const Extractor& rExtractor);

    /**
     *  @brief Call the extractor and append the results to the associated file.
     *  @warning The associated file must not be open.
     *  @warning Only the thread that assigned the extractor can successfully
     *           invoke this function. Otherwise a warning is issued but no
     *           further operation is taken (associated file is not written to).
     */
    void Push(const Model& rModel);

    /// @copydoc JournalBase::Push(const Model&)
    void Push(const value_type& rEntry);

    /**
     *  @brief Erase an entry from the associated file.
     *  @warning The associated file must not be open.
     *  @details This operation is extremely inefficient, as a new
     *           new file has to be written that contains every line
     *           of the currently associated file except the erased one.
     *  @todo Make this function thread safe.
    */
    void Erase(const_iterator itEntry);

    /**
     *  @brief Erase a range of entries from the associated file.
     *  @warning The associated file must not be open.
     *  @details This operation is extremely inefficient, as a new
     *           new file has to be written that contains every line
     *           of the currently associated file except those to be
     *           erased.
     *  @todo Make this function thread safe.
    */
    void Erase(const_iterator Begin, const_iterator End);

    /**
     *  @brief Erase all lines from the associated file matching the input predicate.
     *  @arg rPredicate Functor returning @a true for lines to be deleted.
     *  @warning The associated file must not be open.
     *  @details This function iterates through all lines in the associated file one
     *           by one and calls JournalBase::Erase on each that matches the predicate.
    */
    void EraseIf(const std::function<bool(const value_type&)>& rPredicate);

    /// @brief Delete the associated file.
    /// @warning The associated file must not be open.
    void Clear();

    /// @brief Check whether the associated file is opened by this object.
    bool IsOpen() const noexcept;

    /// @brief Create an iterator to the first line of the associated file.
    /// @warning The file remains open for the iterator's lifetime.
    [[nodiscard]] const_iterator begin() const;

    /// @brief Create an iterator to the last (empty) line of the associated file.
    /// @warning The file remains open for the iterator's lifetime.
    [[nodiscard]] const_iterator end() const;

    /// @brief Count the number of lines in the file (including the last empty one).
    size_type size() const;

protected:
    /// @brief A wrapper class for identifying threads.
    class ThreadID
    {
    private:
        #ifdef KRATOS_SMP_OPENMP
            using ID = int;
            ID mID;
        #elif defined(KRATOS_SMP_CXX11)
            using ID = std::thread::thread_id;
            ID mID;
        #endif

    public:
        ThreadID();

        friend bool operator==(ThreadID Left, ThreadID Right)
        {return Left.mID == Right.mID;}

        friend bool operator!=(ThreadID Left, ThreadID Right)
        {return Left.mID != Right.mID;}

        friend std::ostream& operator<<(std::ostream& rStream, ThreadID threadID)
        {return rStream << threadID.mID;}
    }; // class ThreadID

private:
    /// @brief Erase all lines in the provided set.
    void Erase(const std::set<const_iterator>& rLines);

    /// @brief Check whether the result of an extractor invocation is valid (has no line breaks).
    static bool IsValidEntry(const value_type& rEntry);

    /// @brief Open the associated file and lock its mutex.
    [[nodiscard]] std::shared_ptr<iterator::FileAccess> Open(std::ios::openmode OpenMode = std::ios::in) const;

    std::filesystem::path mJournalPath;

    std::pair<Extractor,ThreadID> mExtractor;

    mutable std::weak_ptr<iterator::FileAccess> mpFileAccess;
}; // class JournalBase


/**
 *  @brief A class for keeping administrational text data throughout an analysis.
 *  @details This is a wrapper class around @ref JournalBase, but writes entries
 *           via @ref Parameters objects instead of raw strings, offering a bit
 *           more convenience.
 */
class Journal
{
public:
    /// @brief An iterator wrapping @ref JournalBase::iterator but with @ref Parameters as value_type.
    class iterator
    {
    public:
        friend class Journal;

        using iterator_category = std::forward_iterator_tag;

        using value_type = Parameters;

        using difference_type = std::ptrdiff_t;

        using pointer = void;

        using reference = void;

    public:
        iterator() = delete;

        iterator(iterator&& rOther) = default;

        iterator(const iterator& rOther) = default;

        iterator& operator=(iterator&& rOther) = default;

        iterator& operator=(const iterator& rOther) = default;

        iterator& operator++();

        iterator operator++(int);

        value_type operator*() const;

        JournalBase::iterator GetBase() const
        {return mWrapped;}

        friend bool operator==(const iterator& rLeft, const iterator& rRight);

        friend bool operator!=(const iterator& rLeft, const iterator& rRight);

    private:
        iterator(JournalBase::iterator&& rWrapped);

        JournalBase::iterator mWrapped;
    }; // class iterator

    using const_iterator = iterator;

    using value_type = iterator::value_type;

    using size_type = std::size_t;

    using Extractor = std::function<Parameters(const Model&)>;

    KRATOS_CLASS_POINTER_DEFINITION(Journal);

public:
    /// @copydoc JournalBase::JournalBase()
    Journal();

    /// @copydoc JournalBase::JournalBase(const std::filesystem::path&)
    Journal(const std::filesystem::path& rJournalPath);

    /// @copydoc JournalBase::JournalBase(const std::filesystem::path&, const Extractor&)
    Journal(const std::filesystem::path& rJournalPath,
            const Extractor& rExtractor);

    /// @copydoc JournalBase::JournalBase(JournalBase&&)
    Journal(Journal&& rOther) = delete;

    /// @copydoc JournalBase::JournalBase(const JournalBase&)
    explicit Journal(const Journal& rOther) = default;

    /// @copydoc JournalBase::operator=(JournalBase&&)
    Journal& operator=(Journal&& rOther) = delete;

    /// @copydoc JournalBase::operator=(const JournalBase&)
    Journal& operator=(const Journal& rOther) = default;

    /// @copydoc JournalBase::GetFilePath()
    const std::filesystem::path& GetFilePath() const noexcept;

    /// @copydoc JournalBase::SetExtractor(Extractor&&)
    void SetExtractor(Extractor&& rExtractor);

    /// @copydoc JournalBase::SetExtractor(const Extractor&)
    void SetExtractor(const Extractor& rExtractor);

    /// @copydoc JournalBase::Push(const Model&)
    void Push(const Model& rModel);

    /// @copydoc JournalBase::Erase(const_iterator)
    void Erase(const_iterator itEntry);

    /// @copydoc JournalBase::Erase(const_iterator, const_iterator)
    void Erase(const_iterator Begin, const_iterator End);

    /// @copydoc JournalBase::EraseIf
    void EraseIf(const std::function<bool(const value_type&)>& rPredicate);

    /// JournalBase::Clear()
    void Clear();

    /// @copydoc JournalBase::IsOpen()
    bool IsOpen() const noexcept;

    /// @copydoc JournalBase::begin()
    [[nodiscard]] const_iterator begin() const;

    /// @copydoc JournalBase::end()
    [[nodiscard]] const_iterator end() const;

    /// @copydoc JournalBase::size()
    size_type size() const;

private:
    JournalBase mBase;
}; // class Journal


/// @}


} // namespace Kratos
