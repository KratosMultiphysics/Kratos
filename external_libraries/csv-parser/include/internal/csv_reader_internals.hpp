#pragma once
#include <array>
#include <deque>
#include <functional>
#include <fstream>
#include <memory>
#include <string>

#include "compatibility.hpp"
#include "csv_row.hpp"
#include "row_buffer.hpp"

namespace csv {
    namespace internals {
        /**  @typedef ParseFlags
         *   An enum used for describing the significance of each character
         *   with respect to CSV parsing
         */
        enum ParseFlags {
            NOT_SPECIAL, /**< Characters with no special meaning */
            QUOTE,       /**< Characters which may signify a quote escape */
            DELIMITER,   /**< Characters which may signify a new field */
            NEWLINE      /**< Characters which may signify a new row */
        };

        using ParseFlagMap = std::array<ParseFlags, 256>;
        using WhitespaceMap = std::array<bool, 256>;

        /** Create a vector v where each index i corresponds to the
         *  ASCII number for a character and, v[i + 128] labels it according to
         *  the CSVReader::ParseFlags enum
         */
        HEDLEY_CONST CONSTEXPR ParseFlagMap make_parse_flags(char delimiter) {
            std::array<ParseFlags, 256> ret = {};
            for (int i = -128; i < 128; i++) {
                const int arr_idx = i + 128;
                char ch = char(i);

                if (ch == delimiter)
                    ret[arr_idx] = DELIMITER;
                else if (ch == '\r' || ch == '\n')
                    ret[arr_idx] = NEWLINE;
                else
                    ret[arr_idx] = NOT_SPECIAL;
            }

            return ret;
        }

        /** Create a vector v where each index i corresponds to the
         *  ASCII number for a character and, v[i + 128] labels it according to
         *  the CSVReader::ParseFlags enum
         */
        HEDLEY_CONST CONSTEXPR ParseFlagMap make_parse_flags(char delimiter, char quote_char) {
            std::array<ParseFlags, 256> ret = make_parse_flags(delimiter);
            ret[(size_t)quote_char + 128] = QUOTE;
            return ret;
        }

        /** Create a vector v where each index i corresponds to the
         *  ASCII number for a character c and, v[i + 128] is true if
         *  c is a whitespace character
         */
        HEDLEY_CONST CONSTEXPR WhitespaceMap make_ws_flags(const char * ws_chars, size_t n_chars) {
            std::array<bool, 256> ret = {};
            for (int i = -128; i < 128; i++) {
                const int arr_idx = i + 128;
                char ch = char(i);
                ret[arr_idx] = false;

                for (size_t j = 0; j < n_chars; j++) {
                    if (ws_chars[j] == ch) {
                        ret[arr_idx] = true;
                    }
                }
            }

            return ret;
        }

        struct GuessScore {
            double score;
            size_t header;
        };

        CSV_INLINE GuessScore calculate_score(csv::string_view head, CSVFormat format);

        CSVGuessResult _guess_format(csv::string_view head, const std::vector<char>& delims = { ',', '|', '\t', ';', '^', '~' });

        /** Parse a CSV field until a delimiter is hit
         *  @return A value indicating whether or not text to be
         *          saved to the text buffer
         */
        CONSTEXPR bool parse_not_special(
            csv::string_view in,
            const csv::internals::ParseFlags* const parse_flags,
            const bool* const ws_flags,
            size_t& i,
            size_t& start,
            size_t& end) {
            // Trim off leading whitespace
            while (i < in.size() && ws_flags[in[i] + 128]) {
                i++;
            }

            start = i;

            // Case: This field is entirely whitespace
            if (parse_flags[in[start] + 128] >= ParseFlags::DELIMITER) {
                // Back the parser up one character so switch statement
                // can process the delimiter or newline
                i--;
                return false;
            }

            // Optimization: Since NOT_SPECIAL characters tend to occur in contiguous
            // sequences, use the loop below to avoid having to go through the outer
            // switch statement as much as possible
            while (i + 1 < in.size()
                && parse_flags[in[i + 1] + 128] == ParseFlags::NOT_SPECIAL) {
                i++;
            }

            // Trim off trailing whitespace
            end = i;
            while (ws_flags[in[end] + 128]) {
                end--;
            }

            return true;
        }

        struct ParseData {
            csv::string_view in;
            ParseFlagMap parse_flags;
            WhitespaceMap ws_flags;
            BufferPtr row_buffer;
            bool& quote_escape; /*< Whether or not we are currently in a quote escaped field */
            std::deque<CSVRow>& records;
        };

        CSV_INLINE BufferPtr parse(const ParseData& data);
        CSV_INLINE void write_record(const ParseData& data);

        /** Read the first 500KB of a CSV file */
        CSV_INLINE std::string get_csv_head(csv::string_view filename);
    }
}