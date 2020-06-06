//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Vincent La (This an adaptation of https://github.com/vincentlaucsb/csv-parser)
//                   Vicente Mataix Ferrandiz
//
//

#ifdef KRATOS_BUILD_CSV_TESTING

// System includes

// External includes

// Project includes
#include "testing/testing.h"
#include "csv-parser/include/csv.hpp"

namespace Kratos {

namespace Testing {

using namespace csv::internals;

std::size_t split_at(BufferPtr buffer, ColumnPositions pos, int n) {
    return buffer->split_buffer[pos.start + n];
}

std::string get_string(
    const std::string& buffer,
    const std::pair<std::size_t, std::size_t>& row_str) {
    return std::string(
        buffer.c_str() + row_str.first, // Beginning
        row_str.second // Count
    );
}

KRATOS_TEST_CASE_IN_SUITE(GiantStringBufferTest, KratosExternalLibrariesFastSuite)
{
    BufferPtr buffer = BufferPtr(new RawRowBuffer());

    buffer->buffer.append("1234");
    std::string first_row = get_string(
        buffer->buffer,
        buffer->get_row().row_str
    );

    buffer->buffer.append("5678");
    std::string second_row = get_string(
        buffer->buffer,
        buffer->get_row().row_str
    );

    buffer = buffer->reset();
    buffer->buffer.append("abcd");
    std::string third_row = get_string(
        buffer->buffer,
        buffer->get_row().row_str
    );

    KRATOS_CHECK_STRING_EQUAL(first_row, "1234");
    KRATOS_CHECK_STRING_EQUAL(second_row, "5678");
    KRATOS_CHECK_STRING_EQUAL(third_row, "abcd");
}

KRATOS_TEST_CASE_IN_SUITE(GiantSplitBufferTest, KratosExternalLibrariesFastSuite)
{
    BufferPtr buffer_1 = BufferPtr(new RawRowBuffer());
    auto* splits_1 = &(buffer_1->split_buffer);

    splits_1->push_back(1);
    splits_1->push_back(2);
    splits_1->push_back(3);

    auto pos = buffer_1->get_row().col_pos;
    KRATOS_CHECK_EQUAL(split_at(buffer_1, pos, 0), 1);
    KRATOS_CHECK_EQUAL(split_at(buffer_1, pos, 1), 2);
    KRATOS_CHECK_EQUAL(split_at(buffer_1, pos, 2), 3);
    KRATOS_CHECK_EQUAL(pos.n_cols, 4);

    // Two Splits Test
    {
        splits_1->push_back(4);
        splits_1->push_back(5);

        pos = buffer_1->get_row().col_pos;
        KRATOS_CHECK_EQUAL(split_at(buffer_1, pos, 0), 4);
        KRATOS_CHECK_EQUAL(split_at(buffer_1, pos, 1), 5);
        KRATOS_CHECK_EQUAL(pos.n_cols, 3);
    }

    BufferPtr buffer_2 = BufferPtr(new RawRowBuffer());
    auto* splits_2 = &(buffer_2->split_buffer);

    splits_2->push_back(1);
    splits_2->push_back(2);
    splits_2->push_back(3);

    pos = buffer_2->get_row().col_pos;

    // Reset In Middle Test
    {
        splits_2->push_back(1);
        buffer_2 = buffer_2->reset();
        splits_2 = &(buffer_2->split_buffer);
        splits_2->push_back(2);
        splits_2->push_back(3);

        pos = buffer_2->get_row().col_pos;
        KRATOS_CHECK_EQUAL(split_at(buffer_2, pos, 0), 1);
        KRATOS_CHECK_EQUAL(split_at(buffer_2, pos, 1), 2);
        KRATOS_CHECK_EQUAL(split_at(buffer_2, pos, 2), 3);
        KRATOS_CHECK_EQUAL(pos.n_cols, 4);
    }
}

} // namespace Testing.
} // namespace Kratos.

#endif
