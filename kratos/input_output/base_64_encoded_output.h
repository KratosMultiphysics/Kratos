//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:         BSD License
//                   Kratos default license: kratos/license.txt
//
//  Main authors:    Suneth Warnakulasuriya
//

#pragma once

// System includes
#include <algorithm>
#include <array>
#include <iterator>
#include <string>

// External includes

// Project includes

namespace Kratos {

///@name Kratos Classes
///@{

/**
 * @brief Encodes given iterator data to base 64 string representation
 *
 *  @brief Encode incoming data to base 64 string representation.
 *  @details Input data is read from chunks of consecutive iterators. Non-consecutive iterators
 *           can be processed by repeatedly calling @ref Base64EncodedOutput::WriteOutputData.
 *  @note The destructor of the instance will write the remaining bytes of the data container with padding.
 *        Hence, scoping of the encoder should be done to properly finish writing of the data stream
 */
class Base64EncodedOutput
{
public:
    ///@name Type definitions
    //@{

    using IndexType = std::size_t;

    ///@}
    ///@name Public classes
    ///@{

    template<class TIteratorType>
    class ByteIterator
    {
    public:
        ///@name Type definitions
        ///@{

        using value_type = char;

        using pointer = char*;

        using reference = char&;

        using iterator_category = std::forward_iterator_tag;

        using difference_type = std::ptrdiff_t;

        ///@}
        ///@name Life cycle
        ///@{

        ByteIterator(TIteratorType it)
            : mIt(it),
            mValue(*it),
            mByteIndex(0u)
        {}

        ///@}
        ///@name Public operations
        ///@{

        value_type operator*() const
        {
            return reinterpret_cast<const value_type*>(&mValue)[mByteIndex];
        }

        ByteIterator& operator++()
        {
            using value_type = typename std::iterator_traits<TIteratorType>::value_type;
            constexpr IndexType value_size = sizeof(value_type);
            if ((++mByteIndex) == value_size) {
                mValue = *++mIt;
                mByteIndex = 0;
            }
            return *this;
        }

        ByteIterator operator++(int)
        {
            ByteIterator copy = *this;
            ++(*this);
            return copy;
        }

        ///@}

    private:
        ///@name Private member variables
        ///@{

        TIteratorType mIt;

        typename std::iterator_traits<TIteratorType>::value_type mValue;

        IndexType mByteIndex = 0u;

        ///@}
    };

    ///}
    ///@name Life cycle
    ///@{

    Base64EncodedOutput(std::ostream& rOStream)
        : mrOStream(rOStream)
    {
    }

    ~Base64EncodedOutput()
    {
        const IndexType padding = 3 - mByteTripletIndex;
        if (padding != 0 && mByteTripletIndex != 0) {
            for (; mByteTripletIndex < 3; ++mByteTripletIndex) {
                mByteTriplet[mByteTripletIndex] = '\0';
            }

            EncodeTriplet(mrOStream, mByteTriplet, padding);
        }
        mByteTripletIndex = 0;
    }

    ///@}
    ///@name Public operations
    ///@{

    /**
     * @brief writes the iterator data to the given output stream
     *
     * This method continues writing the data given by the @ref Begin iterator
     * to the speified @ref rOutput stream in the base64 format.
     *
     * @tparam TIteratorType            Type of the iterator. Should satisfy requirements for std::input_iterator.
     * @param Begin                     Begining of the iterator.
     * @param N                         Number of data items in the data collection represented by the Begin iterator.
     */
    template <typename TIteratorType,
              std::enable_if_t<std::is_base_of_v<std::input_iterator_tag, typename std::iterator_traits<TIteratorType>::iterator_category>, bool> = true>
    void WriteData(
        TIteratorType Begin,
        const IndexType N)
    {
        // skips writing if an empty container is passed.
        if (N == 0) {
            return;
        }

        using value_type = typename std::iterator_traits<TIteratorType>::value_type;

        ByteIterator itr_byte_triplet(Begin);

        const IndexType raw_bytes = N * sizeof(value_type);

        // first fill the existing mByteTriplet
        IndexType number_of_written_bytes = 0;
        for (; mByteTripletIndex < std::min(raw_bytes, IndexType{3}); ++mByteTripletIndex) {
            mByteTriplet[mByteTripletIndex] = *itr_byte_triplet++;
            ++number_of_written_bytes;
        }

        // check if the triplets are filled otherwise don't write, wait for the
        // closeOutputData call to write remaining bytes in mByteTriplet
        if (mByteTripletIndex == 3) {
            mByteTripletIndex = 0;
            // write the last byte triplet
            EncodeTriplet(mrOStream, mByteTriplet, 0);

            // in steps of 3
            const IndexType number_of_triplets = (raw_bytes - number_of_written_bytes) / 3;

            for (IndexType i = 0; i < number_of_triplets; ++i) {
                EncodeTriplet(mrOStream, {*itr_byte_triplet++, *itr_byte_triplet++, *itr_byte_triplet++}, 0);
            }

            number_of_written_bytes += number_of_triplets * 3;

            // now fill the mByteTriplet with the remaining bytes of the data
            const IndexType remaining_bytes = raw_bytes - number_of_written_bytes;
            for (; mByteTripletIndex < remaining_bytes; ++mByteTripletIndex) {
                mByteTriplet[mByteTripletIndex] = *itr_byte_triplet++;
            }
        }
    }

    ///@}

private:
    ///@name Private member variables
    ///@{

    std::ostream& mrOStream;

    IndexType mByteTripletIndex = 0;

    std::array<char, 3> mByteTriplet;

    constexpr static char base64Map[] = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";

    ///@}
    ///@name Private static operations
    ///@{

    static void EncodeTriplet(
        std::ostream& rOutput,
        const std::array<char, 3>& rBytes,
        const IndexType Padding)
    {
        char tmp[5] = {
            base64Map[(rBytes[0] & 0xfc) >> 2],
            base64Map[((rBytes[0] & 0x03) << 4) + ((rBytes[1] & 0xf0) >> 4)],
            base64Map[((rBytes[1] & 0x0f) << 2) + ((rBytes[2] & 0xc0) >> 6)],
            base64Map[rBytes[2] & 0x3f], '\0'};

        std::fill(tmp + 4 - Padding, tmp + 4, '=');

        rOutput << tmp;
    }

    ///@}
};

///@}

} // namespace Kratos
