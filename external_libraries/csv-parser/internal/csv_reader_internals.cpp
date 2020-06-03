#include "csv_reader_internals.hpp"
#include <iostream>

namespace csv {
    namespace internals {
        CSV_INLINE GuessScore calculate_score(csv::string_view head, CSVFormat format) {
            // Frequency counter of row length
            std::unordered_map<size_t, size_t> row_tally = { { 0, 0 } };

            // Map row lengths to row num where they first occurred
            std::unordered_map<size_t, size_t> row_when = { { 0, 0 } };

            // Parse the CSV
            auto buffer_ptr = internals::BufferPtr(new internals::RawRowBuffer());
            std::deque<CSVRow> rows;
            bool quote_escape = false;

            internals::parse({
                head,
                internals::make_parse_flags(format.get_delim(), '"'),
                internals::make_ws_flags({}, 0),
                buffer_ptr,
                quote_escape,
                rows
            });

            for (size_t i = 0; i < rows.size(); i++) {
                auto& row = rows[i];

                // Ignore zero-length rows
                if (row.size() > 0) {
                    if (row_tally.find(row.size()) != row_tally.end()) {
                        row_tally[row.size()]++;
                    }
                    else {
                        row_tally[row.size()] = 1;
                        row_when[row.size()] = i;
                    }
                }
            }

            // Most common numbers of columns
            auto max = std::max_element(row_tally.begin(), row_tally.end(),
                [](const std::pair<size_t, size_t>& x,
                    const std::pair<size_t, size_t>& y) {
                        return x.second < y.second; });

            return {
                (double)(max->first * max->second),
                row_when[max->first]
            };
        }

        /** Guess the delimiter used by a delimiter-separated values file */
        CSV_INLINE CSVGuessResult _guess_format(csv::string_view head, const std::vector<char>& delims) {
            /** For each delimiter, find out which row length was most common.
             *  The delimiter with the longest mode row length wins.
             *  Then, the line number of the header row is the first row with
             *  the mode row length.
             */

            CSVFormat format;
            size_t max_score = 0,
                   header = 0;
            char current_delim = delims[0];

            for (char cand_delim : delims) {
                auto result = calculate_score(head, format.delimiter(cand_delim));

                if (result.score > max_score) {
                    max_score = result.score;
                    current_delim = cand_delim;
                    header = result.header;
                }
            }

            return { current_delim, (int)header };
        }

        CSV_INLINE BufferPtr parse(const ParseData& data) {
            using internals::ParseFlags;

            // Optimizations
            auto * HEDLEY_RESTRICT parse_flags = data.parse_flags.data();
            auto * HEDLEY_RESTRICT ws_flags = data.ws_flags.data();
            auto& in = data.in;
            auto& row_buffer = *(data.row_buffer.get());
            auto& text_buffer = row_buffer.buffer;
            auto& split_buffer = row_buffer.split_buffer;
            text_buffer.reserve(data.in.size());
            split_buffer.reserve(data.in.size() / 10);

            for (size_t i = 0; i < in.size(); i++) {
                switch (parse_flags[data.in[i] + 128]) {
                case ParseFlags::DELIMITER:
                    if (!data.quote_escape) {
                        split_buffer.push_back((internals::StrBufferPos)row_buffer.size());
                        break;
                    }

                    HEDLEY_FALL_THROUGH;
                case ParseFlags::NEWLINE:
                    if (!data.quote_escape) {
                        // End of record -> Write record
                        if (i + 1 < in.size() && in[i + 1] == '\n') // Catches CRLF (or LFLF)
                            ++i;

                        data.records.push_back(CSVRow(data.row_buffer));
                        break;
                    }

                    // Treat as regular character
                    text_buffer += in[i];
                    break;
                case ParseFlags::NOT_SPECIAL: {
                    size_t start, end;

                    if (!parse_not_special(
                        in,
                        parse_flags,
                        ws_flags,
                        i,
                        start,
                        end
                    )) {
                        break;
                    }

                    // Finally append text
#ifdef CSV_HAS_CXX17
                    text_buffer += in.substr(start, end - start + 1);
#else
                    for (; start < end + 1; start++) {
                        text_buffer += in[start];
                    }
#endif

                    break;
                }
                default: // Quote
                    if (!data.quote_escape) {
                        // Don't deref past beginning
                        if (i && parse_flags[in[i - 1] + 128] >= ParseFlags::DELIMITER) {
                            // Case: Previous character was delimiter or newline
                            data.quote_escape = true;
                        }
                    }
                    else if (i + 1 < in.size()) {
                        auto next_ch = parse_flags[in[i + 1] + 128];
                        if (next_ch >= ParseFlags::DELIMITER) {
                            // Case: Delim or newline => end of field
                            data.quote_escape = false;
                            break;
                        }

                        // Case: Escaped quote
                        text_buffer += in[i];

                        // Note: Unescaped single quotes can be handled by the parser
                        if (next_ch == ParseFlags::QUOTE)
                            ++i;  // Case: Two consecutive quotes
                    }

                    break;
                }
            }

            return row_buffer.reset();
        }

        CSV_INLINE std::string get_csv_head(csv::string_view filename) {
            const size_t bytes = 500000;
            std::ifstream infile(filename.data());
            if (!infile.is_open()) {
                throw std::runtime_error("Cannot open file " + std::string(filename));
            }

            std::unique_ptr<char[]> buffer(new char[bytes + 1]);
            char * head_buffer = buffer.get();

            for (size_t i = 0; i < bytes + 1; i++) {
                head_buffer[i] = '\0';
            }

            infile.read(head_buffer, bytes);
            return std::string(head_buffer);
        }
    }
}