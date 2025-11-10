/// @author Máté Kelemen

#pragma once

// --- Core Includes ---
#include "includes/lock_object.h" // Kratos::LockObject

// --- STL Includes ---
#include <iterator> // std::forward_iterator_tag, std::ptrdiff_t
#include <optional> // std::optional
#include <fstream> // std::fstream
#include <utility> // std::pair


namespace Kratos {


/// @addtogroup HDF5Application
/// @{


/**
 *  @brief Iterator interface around a file pointer.
 *  @details This class represents a read-only iterator over the items
 *           of a delimited file. Dereferencing this iterator yields a
 *           a range instead of reading the file to a string; if you
 *           need the value_type to be a string, use @ref FileStringIterator.
 *           Alternatively, you can use @ref FileContentIterator::value to
 *           read the pointed-to range into a string.
 *
 *  @details Access to the file is provided by an std::fstream and an associated
 *           lock that ensures that no other thread has access to the file while
 *           any iterators to it are alive. The iterators position within the file
 *           is stored separately, and is restored every time some operation is
 *           performed on the iterator. This makes sure multiple iterators can
 *           use the file.
 *
 *  @warning Even though the file stream is tied to the lifetime of a lock, the
 *           user can still create a separate stream+lock pair to the same file
 *           on a different thread. Please don't do that ... .
 */
class FileContentIterator
{
public:
    using FileAccess = std::optional<std::pair<
        std::fstream,
        LockObject
    >>;

    using iterator_category = std::forward_iterator_tag;

    /// @brief An iterator pair {begin,end} representing a range within the file.
    /// @todo Change this to the Range class if #10114 ever gets merged (@matekelemen).
    using value_type = std::pair<FileContentIterator,FileContentIterator>;

    using difference_type = std::ptrdiff_t;

    using pointer = void;

    using reference = void;

public:
    FileContentIterator(std::shared_ptr<FileAccess>&& rpFileAccess,
                        char delimiter = '\n');

    FileContentIterator(const std::shared_ptr<FileAccess>& rpFileAccess,
                        char delimiter = '\n');

    FileContentIterator() = delete;

    FileContentIterator(FileContentIterator&& rOther) = default;

    FileContentIterator(const FileContentIterator& rOther) = default;

    FileContentIterator& operator=(FileContentIterator&& rOther) = default;

    FileContentIterator& operator=(const FileContentIterator& rOther) = default;

    FileContentIterator& operator++();

    FileContentIterator operator++(int);

    value_type operator*() const;

    /// @brief Read the associated range from the file into a string.
    std::string value() const;

    friend bool operator==(const FileContentIterator& rLeft, const FileContentIterator& rRight);

    friend bool operator!=(const FileContentIterator& rLeft, const FileContentIterator& rRight);

    friend bool operator<(const FileContentIterator& rLeft, const FileContentIterator& rRight);

    /// @brief Set the iterator to the end of the file.
    void SeekEOF();

private:
    using Position = std::fstream::pos_type;

    /// @brief Get the current position of the iterator in the file.
    Position GetPosition() const noexcept;

    /// @brief Find the next position in the stream with the given character.
    bool SeekNext(char Target);

    /// @brief Store the current read position of the stream.
    void StoreState();

    /// @brief Restore the saved read position of the stream.
    void RestoreState() const;

    mutable std::shared_ptr<FileAccess> mpFileAccess;

    char mDelimiter;

    Position mPosition;
}; // class FileContentIterator


/**
 *  @brief A read-only iterator over the items of a delimited file.
 *  @details This is a wrapper class around @ref FileContentIterator but with
 *           std::string as value_type;
 */
class FileStringIterator
{
public:
    using iterator_category = std::forward_iterator_tag;

    using value_type = std::string;

    using difference_type = std::ptrdiff_t;

    using pointer = void;

    using reference = void;

    using FileAccess = FileContentIterator::FileAccess;

public:
    FileStringIterator(std::shared_ptr<FileAccess>&& rpFileAccess,
                       char delimiter = '\n');

    FileStringIterator(const std::shared_ptr<FileAccess>& rpFileAccess,
                       char delimiter = '\n');

    FileStringIterator(FileStringIterator&& rOther) = default;

    FileStringIterator(const FileStringIterator& rOther) = default;

    FileStringIterator& operator=(FileStringIterator&& rOther) = default;

    FileStringIterator& operator=(const FileStringIterator& rOther) = default;

    FileStringIterator& operator++();

    FileStringIterator operator++(int);

    value_type operator*() const;

    /// @copydoc FileContentIterator::SeekEOF
    void SeekEOF();

    friend bool operator==(const FileStringIterator& rLeft, const FileStringIterator& rRight);

    friend bool operator!=(const FileStringIterator& rLeft, const FileStringIterator& rRight);

    friend bool operator<(const FileStringIterator& rLeft, const FileStringIterator& rRight);

private:
    FileContentIterator mWrapped;
}; // class FileStringIterator


/// @}


} // namespace Kratos
