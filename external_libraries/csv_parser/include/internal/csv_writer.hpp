/** @file
  *  A standalone header file for writing delimiter-separated files
  */

#pragma once
#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include "compatibility.hpp"

namespace csv {
    /** @name CSV Writing */
    ///@{
    #ifndef DOXYGEN_SHOULD_SKIP_THIS
    template<char Delim = ',', char Quote = '"'>
    inline std::string csv_escape(csv::string_view in, const bool quote_minimal = true) {
        /** Format a string to be RFC 4180-compliant
         *  @param[in]  in              String to be CSV-formatted
         *  @param[out] quote_minimal   Only quote fields if necessary.
         *                              If False, everything is quoted.
         */

        // Sequence used for escaping quote characters that appear in text
        constexpr char double_quote[3] = { Quote, Quote };

        std::string new_string;
        bool quote_escape = false;     // Do we need a quote escape
        new_string += Quote;           // Start initial quote escape sequence

        for (size_t i = 0; i < in.size(); i++) {
            switch (in[i]) {
            case Quote:
                new_string += double_quote;
                quote_escape = true;
                break;
            case Delim:
                quote_escape = true;
                HEDLEY_FALL_THROUGH;
            default:
                new_string += in[i];
            }
        }

        if (quote_escape || !quote_minimal) {
            new_string += Quote; // Finish off quote escape
            return new_string;
        }

        return std::string(in);
    }
    #endif

    /** 
     *  Class for writing delimiter separated values files
     *
     *  To write formatted strings, one should
     *   -# Initialize a DelimWriter with respect to some output stream 
     *   -# Call write_row() on std::vector<std::string>s of unformatted text
     *
     *  @tparam OutputStream The output stream, e.g. `std::ofstream`, `std::stringstream`
     *  @tparam Delim        The delimiter character
     *  @tparam Quote        The quote character
     *
     *  @par Hint
     *  Use the aliases csv::CSVWriter<OutputStream> to write CSV
     *  formatted strings and csv::TSVWriter<OutputStream>
     *  to write tab separated strings
     *
     *  @par Example
     *  @snippet test_write_csv.cpp CSV Writer Example
     */
    template<class OutputStream, char Delim, char Quote>
    class DelimWriter {
    public:
        /** Construct a DelimWriter over the specified output stream */
        DelimWriter(OutputStream& _out) : out(_out) {};

        /** Construct a DelimWriter over the file
         *
         *  @param[out] filename  File to write to
         */
        DelimWriter(const std::string& filename) : DelimWriter(std::ifstream(filename)) {};

        /** Format a sequence of strings and write to CSV according to RFC 4180
         *
         *  @warning This does not check to make sure row lengths are consistent
         *
         *  @param[in]  record          Sequence of strings to be formatted
         *  @param      quote_minimal   Only quote fields if necessary
         */
        template<typename T, typename Alloc, template <typename, typename> class Container>
        void write_row(const Container<T, Alloc>& record, bool quote_minimal = true) {
            const size_t ilen = record.size();
            size_t i = 0;
            for (auto& field: record) {
                out << csv_escape<Delim, Quote>(field, quote_minimal);
                if (i + 1 != ilen) out << Delim;
                i++;
            }

            out << std::endl;
        }

        /** @copydoc write_row
         *  @return  The current DelimWriter instance (allowing for operator chaining)
         */
        template<typename T, typename Alloc, template <typename, typename> class Container>
        DelimWriter& operator<<(const Container<T, Alloc>& record) {
            this->write_row(record);
            return *this;
        }

    private:
        OutputStream & out;
    };

    /* Uncomment when C++17 support is better
    template<class OutputStream>
    DelimWriter(OutputStream&) -> DelimWriter<OutputStream>;
    */

    /** Class for writing CSV files
     *
     *  @sa csv::DelimWriter::write_row()
     *  @sa csv::DelimWriter::operator<<()
     *
     *  @note Use `csv::make_csv_writer()` to in instatiate this class over
     *        an actual output stream.
     */
    template<class OutputStream>
    using CSVWriter = DelimWriter<OutputStream, ',', '"'>;

    /** Class for writing tab-separated values files
*
     *  @sa csv::DelimWriter::write_row()
     *  @sa csv::DelimWriter::operator<<()
     *
     *  @note Use `csv::make_tsv_writer()` to in instatiate this class over
     *        an actual output stream.
     */
    template<class OutputStream>
    using TSVWriter = DelimWriter<OutputStream, '\t', '"'>;

    //
    // Temporary: Until more C++17 compilers support template deduction guides
    //
    template<class OutputStream>
    inline CSVWriter<OutputStream> make_csv_writer(OutputStream& out) {
        /** Return a CSVWriter over the output stream */
        return CSVWriter<OutputStream>(out);
    }

    template<class OutputStream>
    inline TSVWriter<OutputStream> make_tsv_writer(OutputStream& out) {
        /** Return a TSVWriter over the output stream */
        return TSVWriter<OutputStream>(out);
    }

    ///@}
}