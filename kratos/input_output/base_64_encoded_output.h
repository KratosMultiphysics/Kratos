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
#include <ostream>

// External includes

// Project includes
#include "includes/define.h"

namespace Kratos
{
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
class KRATOS_API(KRATOS_CORE) Base64EncodedOutput
{
public:
    ///@name Type definitions
    //@{

    /// The index type
    using IndexType = std::size_t;

    ///@name Public classes
    ///@{

    /**
    * @class ByteIterator
    * @brief A forward iterator that iterates over bytes in a sequence.
    * @details This class provides an iterator interface for iterating over bytes in a sequence specified by the underlying iterator type.
    * It allows accessing individual bytes of the sequence and advancing the iterator by one byte at a time.
    * @tparam TIteratorType The underlying iterator type that provides the sequence to iterate over.
    */
    template<class TIteratorType>
    class ByteIterator
    {
    public:
        ///@name Type definitions
        ///@{

        using value_type = char;                               //// The value type representing a byte.
        using pointer = char*;                                 //// Pointer to a byte.
        using reference = char&;                               //// Reference to a byte.
        using iterator_category = std::forward_iterator_tag;   //// Iterator category - forward iterator.
        using difference_type = std::ptrdiff_t;                //// Difference type between iterators.

        ///@}
        ///@name Life cycle
        ///@{

        /**
        * @brief Constructor.
        * @param it The underlying iterator representing the current position in the sequence.
        * @details This constructor initializes a ByteIterator object with the specified underlying iterator.
        * The current byte is determined from the value pointed to by the underlying iterator.
        * The byte index is set to the first byte (index 0).
        */
        ByteIterator(TIteratorType it)
            : mIt(it),
            mValue(*it),
            mByteIndex(0u)
        {}

        ///@}
        ///@name Public operations
        ///@{

        /**
        * @brief Dereference operator.
        * @details This operator returns the value of the current byte being pointed to by the iterator.
        * The byte value is obtained by interpreting the value pointed to by the underlying iterator as a sequence of bytes.
        * The byte at the current byte index is returned.
        * @return The current byte value.
        */
        value_type operator*() const
        {
            return reinterpret_cast<const value_type*>(&mValue)[mByteIndex];
        }

        /**
        * @brief Pre-increment operator.
        * @details This operator advances the iterator to the next byte in the sequence.
        * If the current byte index reaches the size of the value type, the underlying iterator is incremented and the byte index is reset to zero.
        * The updated iterator is returned by reference.
        * @return Reference to the updated iterator.
        */
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

        /**
        * @brief Post-increment operator.
        * @details This operator advances the iterator to the next byte in the sequence and returns a copy of the iterator before the increment.
        * The pre-increment operator is invoked to perform the actual incrementation.
        * @return A copy of the iterator before incrementing.
        */
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

        TIteratorType mIt;                                               /// The underlying iterator representing the current position in the sequence.
        typename std::iterator_traits<TIteratorType>::value_type mValue; /// The value pointed to by the underlying iterator.
        IndexType mByteIndex = 0u;                                       /// The index of the current byte being processed.

        ///@}
    };

    ///}
    ///@name Life cycle
    ///@{

    /**
     * @brief Constructor.
     * @details This constructor initializes a Base64EncodedOutput object with the specified output stream.
     * The encoded data will be written to the provided output stream during encoding operations.
     * @param rOStream The output stream where the encoded data will be written.
     */
    Base64EncodedOutput(std::ostream& rOStream)
        : mrOStream(rOStream)
    {
    }

    /**
     * @brief Destructor.
     * @brief This destructor finalizes the encoding process by handling any remaining bytes and padding at the end of the data.
     * It ensures that the encoded data is correctly written to the output stream before the object is destroyed.
     */
    ~Base64EncodedOutput();

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

        // Following offset forces to write bytes to mByteTripliet array even
        // if raw_bytes count is less than the mByteTripletIndex. This is important
        // in the case if there are less number of raw_bytes (in case of UInt8 data type),
        // then raw bytes may be less than mByteTripletIndex, hence avoiding
        // reading data from the data stream.
        const IndexType initial_byte_triplet_offset = raw_bytes + mByteTripletIndex;

        // first fill the existing mByteTriplet
        IndexType number_of_written_bytes = 0;
        for (; mByteTripletIndex < std::min(initial_byte_triplet_offset, IndexType{3}); ++mByteTripletIndex) {
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

    /**
    * @brief The output stream for encoding operations.
    * @details This member variable represents the output stream where the encoded data will be written.
    * It is used by the encoding functions in this class to write the encoded data.
    */
    std::ostream& mrOStream;

    /**
    * @brief The index of the current byte triplet being processed.
    * @details This member variable keeps track of the current byte triplet being processed during encoding.
    * It is used to determine the position of the bytes in the triplet and to handle padding correctly.
    */
    IndexType mByteTripletIndex = 0;

    /**
    * @brief The byte triplet buffer for encoding operations.
    * @details This member variable is an array of three characters representing the current byte triplet being encoded.
    * It is used to store the bytes before they are processed and encoded into base64 format.
    */
    std::array<char, 3> mByteTriplet;

    /**
    * @brief The base64 encoding character mapping.
    * @details This static constexpr member variable is a character array that maps the values used for base64 encoding.
    * It contains the characters 'A' to 'Z', 'a' to 'z', '0' to '9', and the '+' and '/' characters.
    * It is used by the encoding functions to convert the binary data into base64-encoded characters.
    */
    constexpr static char base64Map[] = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";

    ///@}
    ///@name Private static operations
    ///@{

    /**
    * @brief Encodes a triplet of bytes into base64 format and writes it to the specified output stream.
    * @details This function takes three bytes and encodes them using the base64 encoding scheme.
    * The resulting encoded triplet is written to the provided output stream.
    * @param rOutput The output stream where the encoded triplet will be written.
    * @param rBytes An array of three bytes to encode.
    * @param Padding The number of padding characters to append to the encoded triplet.
    *                It should be a value between 0 and 2 (inclusive).
    * @note The output stream should be in a valid state and open for writing.
    */
    static void EncodeTriplet(
        std::ostream& rOutput,
        const std::array<char, 3>& rBytes,
        const IndexType Padding
        )
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
