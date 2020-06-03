/** @file
 *  @brief Defines functionality needed for basic CSV parsing
 */

#include <algorithm>
#include <cstdio>   // For read_csv()
#include <cstring>  // For read_csv()
#include <fstream>
#include <sstream>

#include "constants.hpp"
#include "csv_reader.hpp"

namespace csv {
    namespace internals {
        CSV_INLINE std::string format_row(const std::vector<std::string>& row, csv::string_view delim) {
            /** Print a CSV row */
            std::stringstream ret;
            for (size_t i = 0; i < row.size(); i++) {
                ret << row[i];
                if (i + 1 < row.size()) ret << delim;
                else ret << std::endl;
            }

            return ret.str();
        }

        /** Return a CSV's column names
         *
         *  @param[in] filename  Path to CSV file
         *  @param[in] format    Format of the CSV file
         *
         */
        CSV_INLINE std::vector<std::string> _get_col_names(csv::string_view head, CSVFormat format) {
            auto parse_flags = internals::make_parse_flags(format.get_delim());
            if (format.is_quoting_enabled()) {
                parse_flags = internals::make_parse_flags(format.get_delim(), format.get_quote_char());
            }

            // Parse the CSV
            auto buffer_ptr = internals::BufferPtr(new internals::RawRowBuffer());
            auto trim_chars = format.get_trim_chars();

            std::deque<CSVRow> rows;
            bool quote_escape = false;

            internals::parse({
                head,
                parse_flags,
                internals::make_ws_flags(trim_chars.data(), trim_chars.size()),
                buffer_ptr,
                quote_escape,
                rows
            });

            return rows[format.get_header()];
        }
    }

    /** Return a CSV's column names
     *
     *  @param[in] filename  Path to CSV file
     *  @param[in] format    Format of the CSV file
     *
     */
    CSV_INLINE std::vector<std::string> get_col_names(csv::string_view filename, CSVFormat format) {
        auto head = internals::get_csv_head(filename);

        /** Guess delimiter and header row */
        if (format.guess_delim()) {
            auto guess_result = guess_format(filename, format.get_possible_delims());
            format.delimiter(guess_result.delim).header_row(guess_result.header_row);
        }

        return internals::_get_col_names(head, format);
    }

    /** Guess the delimiter used by a delimiter-separated values file */
    CSV_INLINE CSVGuessResult guess_format(csv::string_view filename, const std::vector<char>& delims) {
        auto head = internals::get_csv_head(filename);
        return internals::_guess_format(head, delims);
    }

    /** Allows parsing in-memory sources (by calling feed() and end_feed()). */
    CSV_INLINE CSVReader::CSVReader(CSVFormat format) : 
        unicode_bom_scan(!format.unicode_detect), feed_state(new ThreadedReadingState)  {
        if (!format.col_names.empty()) {
            this->set_col_names(format.col_names);
        }
        
        this->set_parse_flags(format);
    }

    /** Allows reading a CSV file in chunks, using overlapped
     *  threads for simulatenously reading from disk and parsing.
     *  Rows should be retrieved with read_row() or by using
     *  CSVReader::iterator.
     *
     *  **Details:** Reads the first 500kB of a CSV file to infer file information
     *              such as column names and delimiting character.
     *
     *  @param[in] filename  Path to CSV file
     *  @param[in] format    Format of the CSV file
     *
     *  \snippet tests/test_read_csv.cpp CSVField Example
     *
     */
    CSV_INLINE CSVReader::CSVReader(csv::string_view filename, CSVFormat format) : feed_state(new ThreadedReadingState) {
        auto head = internals::get_csv_head(filename);

        /** Guess delimiter and header row */
        if (format.guess_delim()) {
            auto guess_result = internals::_guess_format(head, format.possible_delimiters);
            format.delimiter(guess_result.delim);
            format.header = guess_result.header_row;
        }

        if (format.col_names.empty()) {
            this->set_col_names(internals::_get_col_names(head, format));
        }
        else {
            this->set_col_names(format.col_names);
        }

        this->set_parse_flags(format);
        this->fopen(filename);
    }

    /** Return the format of the original raw CSV */
    CSV_INLINE CSVFormat CSVReader::get_format() const {
        CSVFormat new_format = this->format;

        // Since users are normally not allowed to set 
        // column names and header row simulatenously,
        // we will set the backing variables directly here
        new_format.col_names = this->col_names->get_col_names();
        new_format.header = this->format.header;

        return new_format;
    }

    /** Return the CSV's column names as a vector of strings. */
    CSV_INLINE std::vector<std::string> CSVReader::get_col_names() const {
        if (this->col_names) {
            return this->col_names->get_col_names();
        }

        return std::vector<std::string>();
    }

    /** Return the index of the column name if found or
     *         csv::CSV_NOT_FOUND otherwise.
     */
    CSV_INLINE int CSVReader::index_of(csv::string_view col_name) const {
        auto _col_names = this->get_col_names();
        for (size_t i = 0; i < _col_names.size(); i++)
            if (_col_names[i] == col_name) return (int)i;

        return CSV_NOT_FOUND;
    }

    CSV_INLINE void CSVReader::feed(WorkItem&& buff) {
        this->feed( csv::string_view(buff.first.get(), buff.second) );
    }

    /** Parse a CSV-formatted string.
     *
     *  @par Usage
     *  Incomplete CSV fragments can be joined together by calling feed() on them sequentially.
     *
     *  @note
     *  `end_feed()` should be called after the last string.
     */
    CSV_INLINE void CSVReader::feed(csv::string_view in) {
        /** Handle possible Unicode byte order mark */
        if (!this->unicode_bom_scan) {
            if (in[0] == '\xEF' && in[1] == '\xBB' && in[2] == '\xBF') {
                in.remove_prefix(3); // Remove BOM from input string
                this->utf8_bom = true;
            }

            this->unicode_bom_scan = true;
        }

        this->record_buffer = internals::parse({
            in,
            this->parse_flags,
            this->ws_flags,
            this->record_buffer,
            this->quote_escape,
            this->records
        });

        if (!this->header_trimmed) {
            for (int i = 0; i <= this->format.header && !this->records.empty(); i++) {
                if (i == this->format.header && this->col_names->empty()) {
                    this->set_col_names(this->records.front());
                }

                this->records.pop_front();
            }

            this->header_trimmed = true;
        }
    }

    CSV_INLINE void CSVReader::end_feed() {
        /** Indicate that there is no more data to receive,
         *  and handle the last row
         */
        if (this->record_buffer->size() > 0) {
            this->records.push_back(CSVRow(this->record_buffer));
        }
    }

    /** Worker thread for read_csv() which parses CSV rows (while the main
     *  thread pulls data from disk)
     */
    CSV_INLINE void CSVReader::read_csv_worker() {
        while (true) {
            std::unique_lock<std::mutex> lock{ this->feed_state->feed_lock }; // Get lock
            this->feed_state->feed_cond.wait(lock,                            // Wait
                [this] { return !(this->feed_state->feed_buffer.empty()); });

            // Wake-up
            auto in = std::move(this->feed_state->feed_buffer.front());
            this->feed_state->feed_buffer.pop_front();

            // Nullptr --> Die
            if (!in.first) break;

            lock.unlock();      // Release lock
            this->feed(std::move(in));
        }
    }

    CSV_INLINE void CSVReader::set_parse_flags(const CSVFormat& format)
    {
        this->format = format;
        if (format.no_quote) {
            this->parse_flags = internals::make_parse_flags(format.get_delim());
        }
        else {
            this->parse_flags = internals::make_parse_flags(format.get_delim(), format.quote_char);
        }

        this->ws_flags = internals::make_ws_flags(format.trim_chars.data(), format.trim_chars.size());
    }

    CSV_INLINE void CSVReader::fopen(csv::string_view filename) {
        if (!this->infile) {
#ifdef _MSC_BUILD
            // Silence compiler warnings in Microsoft Visual C++
            size_t err = fopen_s(&(this->infile), filename.data(), "rb");
            if (err)
                throw std::runtime_error("Cannot open file " + std::string(filename));
#else
            this->infile = std::fopen(filename.data(), "rb");
            if (!this->infile)
                throw std::runtime_error("Cannot open file " + std::string(filename));
#endif
        }
    }

    /**
     *  @param[in] names Column names
     */
    CSV_INLINE void CSVReader::set_col_names(const std::vector<std::string>& names)
    {
        this->col_names->set_col_names(names);
        this->n_cols = names.size();
    }

    /**
     * Parse a CSV file using multiple threads
     *
     * @pre CSVReader::infile points to a valid file handle, i.e. CSVReader::fopen was called
     *
     * @param[in] bytes Number of bytes to read.
     * @see CSVReader::read_row()
     */
    CSV_INLINE void CSVReader::read_csv(const size_t& bytes) {
        const size_t BUFFER_UPPER_LIMIT = std::min(bytes, (size_t)1000000);
        std::unique_ptr<char[]> buffer(new char[BUFFER_UPPER_LIMIT]);
        auto * HEDLEY_RESTRICT line_buffer = buffer.get();
        line_buffer[0] = '\0';

        std::thread worker(&CSVReader::read_csv_worker, this);

        for (size_t processed = 0; processed < bytes; ) {
            char * HEDLEY_RESTRICT result = std::fgets(line_buffer, internals::PAGE_SIZE, this->infile);
            if (result == NULL) break;
            line_buffer += std::strlen(line_buffer);
            size_t current_strlen = line_buffer - buffer.get();

            if (current_strlen >= 0.9 * BUFFER_UPPER_LIMIT) {
                processed += (line_buffer - buffer.get());
                std::unique_lock<std::mutex> lock{ this->feed_state->feed_lock };

                this->feed_state->feed_buffer.push_back(std::make_pair<>(std::move(buffer), current_strlen));

                buffer = std::unique_ptr<char[]>(new char[BUFFER_UPPER_LIMIT]); // New pointer
                line_buffer = buffer.get();
                line_buffer[0] = '\0';

                this->feed_state->feed_cond.notify_one();
            }
        }

        // Feed remaining bits
        std::unique_lock<std::mutex> lock{ this->feed_state->feed_lock };
        this->feed_state->feed_buffer.push_back(std::make_pair<>(std::move(buffer), line_buffer - buffer.get()));
        this->feed_state->feed_buffer.push_back(std::make_pair<>(nullptr, 0)); // Termination signal
        this->feed_state->feed_cond.notify_one();
        lock.unlock();
        worker.join();

        if (std::feof(this->infile)) {
            this->end_feed();
            this->close();
        }
    }

    /** Close the open file handle.
     *
     *  @note Automatically called by ~CSVReader().
     */
    CSV_INLINE void CSVReader::close() {
        if (this->infile) {
            std::fclose(this->infile);
            this->infile = nullptr;
        }
    }

    /**
     * Retrieve rows as CSVRow objects, returning true if more rows are available.
     *
     * **Performance Notes**:
     *  - The number of rows read in at a time is determined by csv::ITERATION_CHUNK_SIZE
     *  - For performance details, read the documentation for CSVRow and CSVField.
     *
     * @param[out] row The variable where the parsed row will be stored
     * @see CSVRow, CSVField
     *
     * **Example:**
     * \snippet tests/test_read_csv.cpp CSVField Example
     *
     */
    CSV_INLINE bool CSVReader::read_row(CSVRow &row) {
        if (this->records.empty()) {
            if (!this->eof()) {
                // TODO/Suggestion: Make this call non-blocking, 
                // i.e. move to it another thread
                this->read_csv(internals::ITERATION_CHUNK_SIZE);
            }
            else return false; // Stop reading
        }

        while (!this->records.empty()) {
            if (this->records.front().size() != this->n_cols &&
                this->format.variable_column_policy != VariableColumnPolicy::KEEP) {
                if (this->format.variable_column_policy == VariableColumnPolicy::THROW) {
                    if (this->records.front().size() < this->n_cols) {
                        throw std::runtime_error("Line too short " + internals::format_row(this->records.front()));
                    }

                    throw std::runtime_error("Line too long " + internals::format_row(this->records.front()));
                }

                // Silently drop row (default)
                this->records.pop_front();
            }
            else {
                row = std::move(this->records.front());
                this->num_rows++;
                this->records.pop_front();
                return true;
            }
        }
    
        return false;
    }
}