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

// System includes

// External includes

// Project includes
#include "testing/testing.h"
#include "csv-parser/include/csv.hpp"

namespace Kratos {

namespace Testing {

using namespace csv::internals;

size_t split_at(BufferPtr buffer, ColumnPositions pos, int n) {
    return buffer->split_buffer[pos.start + n];
}

std::string get_string(
    const std::string& buffer,
    const std::pair<size_t, size_t>& row_str) {
    return std::string(
        buffer.c_str() + row_str.first, // Beginning
        row_str.second // Count
    );
}

// KRATOS_TEST_CASE_IN_SUITE(GiantStringBufferTest, KratosExternalLibrariesFastSuite) {
//     BufferPtr buffer = BufferPtr(new RawRowBuffer());
//
//     buffer->buffer.append("1234");
//     std::string first_row = get_string(
//         buffer->buffer,
//         buffer->get_row().row_str
//     );
//
//     buffer->buffer.append("5678");
//     std::string second_row = get_string(
//         buffer->buffer,
//         buffer->get_row().row_str
//     );
//
//     buffer = buffer->reset();
//     buffer->buffer.append("abcd");
//     std::string third_row = get_string(
//         buffer->buffer,
//         buffer->get_row().row_str
//     );
//
//     KRATOS_CHECK(first_row == "1234");
//     KRATOS_CHECK(second_row == "5678");
//     KRATOS_CHECK(third_row == "abcd");
// }
//
// KRATOS_TEST_CASE_IN_SUITE(GiantSplitBufferTest, KratosExternalLibrariesFastSuite) {
//     BufferPtr buffer = BufferPtr(new RawRowBuffer());
//     auto * splits = &(buffer->split_buffer);
//
//     splits->push_back(1);
//     splits->push_back(2);
//     splits->push_back(3);
//
//     auto pos = buffer->get_row().col_pos;
//     KRATOS_CHECK(split_at(buffer, pos, 0) == 1);
//     KRATOS_CHECK(split_at(buffer, pos, 1) == 2);
//     KRATOS_CHECK(split_at(buffer, pos, 2) == 3);
//     KRATOS_CHECK(pos.n_cols == 4);
//
//     SECTION("Two Splits Test") {
//         splits->push_back(4);
//         splits->push_back(5);
//
//         pos = buffer->get_row().col_pos;
//         KRATOS_CHECK(split_at(buffer, pos, 0) == 4);
//         KRATOS_CHECK(split_at(buffer, pos, 1) == 5);
//         KRATOS_CHECK(pos.n_cols == 3);
//     }
//
//     SECTION("Reset In Middle Test") {
//         splits->push_back(1);
//         buffer = buffer->reset();
//         splits = &(buffer->split_buffer);
//         splits->push_back(2);
//         splits->push_back(3);
//
//         pos = buffer->get_row().col_pos;
//         KRATOS_CHECK(split_at(buffer, pos, 0) == 1);
//         KRATOS_CHECK(split_at(buffer, pos, 1) == 2);
//         KRATOS_CHECK(split_at(buffer, pos, 2) == 3);
//         KRATOS_CHECK(pos.n_cols == 4);
//     }
// }

} // namespace Testing.
} // namespace Kratos.
