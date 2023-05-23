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
 * This class allows encoding data given by the iterator (can be non-contiguous) to base 64
 * string representation. It allows having list of iterators as well. In the case if a list of
 * iterators needs to be processed for one base64 string representation, then they needs
 * to be passed in the order to the WriteOutputData method. Finally the output should be
 * closed with CloseOutputData call.
 *
 * The TIteratorType should satisfy requirements for a std::input_iterator.
 *
 */
class Base64EncodedOutput
{
public:
    ///@name Type definitions
    //@{

    using IndexType = std::size_t;

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
     * @param rOutput                   Output stream.
     * @param Begin                     Begining of the iterator.
     * @param N                         Number of data items in the data collection represented by the Begin iterator.
     */
    template <typename TIteratorType,
              std::enable_if_t<std::is_base_of_v<std::input_iterator_tag, typename std::iterator_traits<TIteratorType>::iterator_category>, bool> = true>
    void WriteOutputData(
        std::ostream& rOutput,
        TIteratorType Begin,
        const IndexType N)
    {
        using value_type = typename std::iterator_traits<TIteratorType>::value_type;

        constexpr IndexType size = sizeof(value_type);
        const IndexType raw_bytes = N * size;

        auto it = Begin;
        auto value = *it;
        IndexType byte_index = 0;

        auto next = [&value, &byte_index, &it]() {
            char byte = *(reinterpret_cast<const char*>(&value) + byte_index++);

            if (byte_index == size) {
                ++it;
                value = *it;
                byte_index = 0;
            }

            return byte;
        };

        // first fill the existing mByteTriplet
        IndexType number_of_written_bytes = 0;
        for (; mByteTripletIndex < std::min(raw_bytes, IndexType{3}); ++mByteTripletIndex) {
            mByteTriplet[mByteTripletIndex] = next();
            ++number_of_written_bytes;
        }

        // check if the triplets are filled otherwise don't write, wait for the
        // closeOutputData call to write remaining bytes in mByteTriplet
        if (mByteTripletIndex == 3) {
            mByteTripletIndex = 0;
            // write the last byte triplet
            EncodeTriplet(rOutput, mByteTriplet, 0);

            // in steps of 3
            const IndexType number_of_triplets = (raw_bytes - number_of_written_bytes) / 3;

            for (IndexType i = 0; i < number_of_triplets; ++i) {
                EncodeTriplet(rOutput, {next(), next(), next()}, 0);
            }

            number_of_written_bytes += number_of_triplets * 3;

            // now fill the mByteTriplet with the remaining bytes of the data
            const IndexType remaining_bytes = raw_bytes - number_of_written_bytes;
            for (; mByteTripletIndex < remaining_bytes; ++mByteTripletIndex) {
                mByteTriplet[mByteTripletIndex] = next();
            }
        }
    }

    /**
     * @brief Closes the base64 output.
     *
     * This method should be called when the data stream is finished to write
     * if required the padding. This will reset the writer as well so
     * a new data collection can be written next time.
     *
     * @param rOutput           Output stream.
     */
    void CloseOutputData(std::ostream& rOutput)
    {
        const IndexType padding = 3 - mByteTripletIndex;
        if (padding != 0 && mByteTripletIndex != 0) {
            for (; mByteTripletIndex < 3; ++mByteTripletIndex) {
                mByteTriplet[mByteTripletIndex] = '\0';
            }

            EncodeTriplet(rOutput, mByteTriplet, padding);
        }
        mByteTripletIndex = 0;
    }

    ///@}

private:
    ///@name Private member variables
    ///@{

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
