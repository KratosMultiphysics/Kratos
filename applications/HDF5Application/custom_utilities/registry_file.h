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


class RegistryFile
{
private:
    using RegistryAccess = std::optional<std::pair<
        std::fstream,
        std::unique_lock<LockObject>
    >>;

    friend class iterator;

public:
    class iterator
    {
    public:
        friend class RegistryFile;

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

        void StoreState();

        void RestoreState() const;

        bool SeekNext(char Target);

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
    RegistryFile();

    RegistryFile(const std::filesystem::path& rRegistryFilePath);

    RegistryFile(const std::filesystem::path& rRegistryFilePath,
                 const Extractor& rExtractor);

    RegistryFile(RegistryFile&& rOther) = delete;

    RegistryFile(const RegistryFile& rOther);

    RegistryFile& operator=(RegistryFile&& rOther) = delete;

    RegistryFile& operator=(const RegistryFile& rOther);

    const std::filesystem::path& GetFilePath() const noexcept;

    void SetExtractor(Extractor&& rExtractor);

    void SetExtractor(const Extractor& rExtractor);

    void Push(const Model& rModel);

    void Clear();

    bool IsOpen() const noexcept;

    [[nodiscard]] const_iterator begin() const;

    [[nodiscard]] const_iterator end() const;

    size_type size() const;

private:
    static bool IsValidEntry(const Parameters& rEntry);

    [[nodiscard]] std::shared_ptr<RegistryAccess> Open(std::ios::openmode OpenMode = std::ios::in) const;

    std::filesystem::path mRegistryFilePath;

    Extractor mExtractor;

    mutable std::weak_ptr<RegistryAccess> mpRegistryAccess;

    mutable LockObject mMutex;
}; // class RegistryFile


} // namespace Kratos
