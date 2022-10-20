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

#pragma once

// Core includes
#include "includes/lock_object.h"
#include "includes/kratos_parameters.h"
#include "containers/model.h"
#include "includes/smart_pointers.h"
#include "utilities/parallel_utilities.h"

// STL includes
#include <filesystem>
#include <fstream>
#include <iterator>
#include <optional>
#include <functional>


namespace Kratos
{


/**
 *  @brief A class for keeping administrational text data throughout an analysis.
 *  @details This class is primarily meant to keep track of output files generated during a simulation,
 *           though extra flexibility is provided to extend its purpose as necessary. An associated
 *           text file can be written to via @ref RegistryFile::Push. On each push, a @ref Model is
 *           taken as input, parsed, and a @ref Parameters object representing a scalar / string array
 *           is appended to the file. Generating the output @ref Parameters from the input @ref Model
 *           is the job of a functor that can be set in the constructor, or in
 *           @ref RegistryFile::SetExtractor.
 *  @warning The associated file must be empty or end with an empty line.
 */
class RegistryFile
{
private:
    using RegistryAccess = std::optional<std::pair<
        std::fstream,
        std::unique_lock<LockObject>
    >>;

    friend class iterator;

public:
    /**
     *  @brief An iterator providing access to all stored data from pushes.
     *  @details Each push is allowed to write exclusively to a single line,
     *           so the iterator essentially loops through the lines of the
     *           output file.
     */
    class iterator
    {
    public:
        friend class RegistryFile;

        /// Minimal forward iterator interface.
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

        Parameters operator*() const;

        friend bool operator==(const iterator& rLeft, const iterator& rRight);

        friend bool operator!=(const iterator& rLeft, const iterator& rRight);

        friend bool operator<(const iterator& rLeft, const iterator& rRight);

    private:
        iterator(std::shared_ptr<RegistryFile::RegistryAccess>&& rRegistryAccess,
                 std::size_t EntryIndex = 0);

        /// Store the current read position of the stream.
        void StoreState();

        /// Restore the saved read position of the stream.
        void RestoreState() const;

        /// Find the next position in the stream with a new line character.
        bool SeekNext(char Target);

        /// Set the stream position to EOF.
        void SeekEnd();

        mutable std::shared_ptr<RegistryFile::RegistryAccess> mpRegistryAccess;

        std::ifstream::pos_type mFilePointer;
    }; // class iterator

    using const_iterator = iterator;

    using value_type = iterator::value_type;

    using size_type = std::size_t;

    using Extractor = std::function<Parameters(const Model&)>;

    KRATOS_CLASS_POINTER_DEFINITION(RegistryFile);

public:
    /// @brief Construct the object with an empty file name (invalid).
    RegistryFile();

    /// @brief Construct a registry given an associated file path and a no-op extractor.
    /// @note Should the file already exist, it will be appended to instead of overwritten.
    RegistryFile(const std::filesystem::path& rRegistryFilePath);

    /// @brief Construct a registry given an associated file path and extractor.
    /// @note Should the file already exist, it will be appended to instead of overwritten.
    RegistryFile(const std::filesystem::path& rRegistryFilePath,
                 const Extractor& rExtractor);

    /// @brief The move constructor is deleted because it's ambiguous what should happen to the associated file.
    RegistryFile(RegistryFile&& rOther) = delete;

    /// @brief Copy constructor.
    /// @warning Deletes the associated file if the incoming instance has another one.
    /// @warning The associated file must not be open.
    explicit RegistryFile(const RegistryFile& rOther);

    /// @brief The move assignment operator is deleted because it's ambiguous what should happen to the associated file.
    RegistryFile& operator=(RegistryFile&& rOther) = delete;

    /// @brief Copy assignment operator.
    /// @warning Deletes the associated file if the incoming instance has another one.
    /// @warning The associated file must not be open.
    RegistryFile& operator=(const RegistryFile& rOther);

    /// @brief Get the path to the associated file.
    const std::filesystem::path& GetFilePath() const noexcept;

    /**
     *  @brief Set the extractor handling the conversion from @ref Model to @ref Parameters.
     *  @details The extractor must generate a @ref Parameters object that consists of a
     *           (possibly heterogeneous) array holding bools, numeric values, or strings
     *           without new line characters ('\n').
     */
    void SetExtractor(Extractor&& rExtractor);

    /// @brief Set the extractor handling the conversion from @ref Model to @ref Parameters.
    void SetExtractor(const Extractor& rExtractor);

    /// @brief Call the extractor and append the results to the associated file.
    /// @warning The associated file must not be open.
    void Push(const Model& rModel);

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

private:
    /// @brief Check whether the result of an extractor invocation is valid.
    static bool IsValidEntry(const Parameters& rEntry);

    /// @brief Open the associated file and lock the mutex.
    [[nodiscard]] std::shared_ptr<RegistryAccess> Open(std::ios::openmode OpenMode = std::ios::in) const;

    std::filesystem::path mRegistryFilePath;

    Extractor mExtractor;

    mutable std::weak_ptr<RegistryAccess> mpRegistryAccess;

    mutable LockObject mMutex;
}; // class RegistryFile


} // namespace Kratos
