/** @file
 *  Defines an object used to store CSV format settings
 */

#include <algorithm>
#include <set>

#include "csv_format.hpp"

namespace csv {
    CSV_INLINE CSVFormat& CSVFormat::delimiter(char delim) {
        this->possible_delimiters = { delim };
        this->assert_no_char_overlap();
        return *this;
    }

    CSV_INLINE CSVFormat& CSVFormat::delimiter(const std::vector<char> & delim) {
        this->possible_delimiters = delim;
        this->assert_no_char_overlap();
        return *this;
    }

    CSV_INLINE CSVFormat& CSVFormat::quote(char quote) {
        this->no_quote = false;
        this->quote_char = quote;
        this->assert_no_char_overlap();
        return *this;
    }

    CSV_INLINE CSVFormat& CSVFormat::trim(const std::vector<char> & chars) {
        this->trim_chars = chars;
        this->assert_no_char_overlap();
        return *this;
    }

    CSV_INLINE CSVFormat& CSVFormat::column_names(const std::vector<std::string>& names) {
        this->col_names = names;
        this->header = -1;
        return *this;
    }

    CSV_INLINE CSVFormat& CSVFormat::header_row(int row) {
        this->header = row;
        this->col_names = {};
        return *this;
    }

    CSV_INLINE void CSVFormat::assert_no_char_overlap()
    {
        auto delims = std::set<char>(
            this->possible_delimiters.begin(), this->possible_delimiters.end()),
            trims = std::set<char>(
                this->trim_chars.begin(), this->trim_chars.end());

        // Stores intersection of possible delimiters and trim characters
        std::vector<char> intersection = {};

        // Find which characters overlap, if any
        std::set_intersection(
            delims.begin(), delims.end(),
            trims.begin(), trims.end(),
            std::back_inserter(intersection));

        // Make sure quote character is not contained in possible delimiters
        // or whitespace characters
        if (delims.find(this->quote_char) != delims.end() ||
            trims.find(this->quote_char) != trims.end()) {
            intersection.push_back(this->quote_char);
        }

        if (!intersection.empty()) {
            std::string err_msg = "There should be no overlap between the quote character, "
                "the set of possible delimiters "
                "and the set of whitespace characters. Offending characters: ";

            // Create a pretty error message with the list of overlapping
            // characters
            for (size_t i = 0; i < intersection.size(); i++) {
                err_msg += "'";
                err_msg += intersection[i];
                err_msg += "'";

                if (i + 1 < intersection.size())
                    err_msg += ", ";
            }

            throw std::runtime_error(err_msg + '.');
        }
    }
}