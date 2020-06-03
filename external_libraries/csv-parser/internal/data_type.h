/** @file
 *  @brief Implements data type parsing functionality
 */

#pragma once
#include <math.h>
#include <cctype>
#include <string>
#include <cassert>

#include "compatibility.hpp"

namespace csv {
    /** Enumerates the different CSV field types that are
     *  recognized by this library
     *
     *  @note Overflowing integers will be stored and classified as doubles.
     *  @note Unlike previous releases, integer enums here are platform agnostic.
     */
    enum DataType {
        UNKNOWN = -1,
        CSV_NULL,   /**< Empty string */
        CSV_STRING, /**< Non-numeric string */
        CSV_INT8,   /**< 8-bit integer */
        CSV_INT16,  /**< 16-bit integer (short on MSVC/GCC) */
        CSV_INT32,  /**< 32-bit integer (int on MSVC/GCC) */
        CSV_INT64,  /**< 64-bit integer (long long on MSVC/GCC) */
        CSV_DOUBLE  /**< Floating point value */
    };

    static_assert(CSV_STRING < CSV_INT8, "String type should come before numeric types.");
    static_assert(CSV_INT8 < CSV_INT64, "Smaller integer types should come before larger integer types.");
    static_assert(CSV_INT64 < CSV_DOUBLE, "Integer types should come before floating point value types.");

    namespace internals {
        /** Compute 10 to the power of n */
        template<typename T>
        HEDLEY_CONST CONSTEXPR
        long double pow10(const T& n) noexcept {
            long double multiplicand = n > 0 ? 10 : 0.1,
                ret = 1;

            // Make all numbers positive
            T iterations = n > 0 ? n : -n;
            
            for (T i = 0; i < iterations; i++) {
                ret *= multiplicand;
            }

            return ret;
        }

        /** Compute 10 to the power of n */
        template<>
        HEDLEY_CONST CONSTEXPR
        long double pow10(const unsigned& n) noexcept {
            long double multiplicand = n > 0 ? 10 : 0.1,
                ret = 1;

            for (unsigned i = 0; i < n; i++) {
                ret *= multiplicand;
            }

            return ret;
        }

#ifndef DOXYGEN_SHOULD_SKIP_THIS
        /** Private site-indexed array mapping byte sizes to an integer size enum */
        constexpr DataType int_type_arr[8] = {
            CSV_INT8,  // 1
            CSV_INT16, // 2
            UNKNOWN,
            CSV_INT32, // 4
            UNKNOWN,
            UNKNOWN,
            UNKNOWN,
            CSV_INT64  // 8
        };

        template<typename T>
        inline DataType type_num() {
            static_assert(std::is_integral<T>::value, "T should be an integral type.");
            static_assert(sizeof(T) <= 8, "Byte size must be no greater than 8.");
            return int_type_arr[sizeof(T) - 1];
        }

        template<> inline DataType type_num<float>() { return CSV_DOUBLE; }
        template<> inline DataType type_num<double>() { return CSV_DOUBLE; }
        template<> inline DataType type_num<long double>() { return CSV_DOUBLE; }
        template<> inline DataType type_num<std::nullptr_t>() { return CSV_NULL; }
        template<> inline DataType type_num<std::string>() { return CSV_STRING; }

        CONSTEXPR DataType data_type(csv::string_view in, long double* const out = nullptr);
#endif

        /** Given a byte size, return the largest number than can be stored in
         *  an integer of that size
         *
         *  Note: Provides a platform-agnostic way of mapping names like "long int" to
         *  byte sizes
         */
        template<size_t Bytes>
        CONSTEXPR long double get_int_max() {
            static_assert(Bytes == 1 || Bytes == 2 || Bytes == 4 || Bytes == 8,
                "Bytes must be a power of 2 below 8.");

            IF_CONSTEXPR (sizeof(signed char) == Bytes) {
                return (long double)std::numeric_limits<signed char>::max();
            }

            IF_CONSTEXPR (sizeof(short) == Bytes) {
                return (long double)std::numeric_limits<short>::max();
            }

            IF_CONSTEXPR (sizeof(int) == Bytes) {
                return (long double)std::numeric_limits<int>::max();
            }

            IF_CONSTEXPR (sizeof(long int) == Bytes) {
                return (long double)std::numeric_limits<long int>::max();
            }

            IF_CONSTEXPR (sizeof(long long int) == Bytes) {
                return (long double)std::numeric_limits<long long int>::max();
            }

            HEDLEY_UNREACHABLE();
        }

        /** Given a byte size, return the largest number than can be stored in
         *  an unsigned integer of that size
         */
        template<size_t Bytes>
        CONSTEXPR long double get_uint_max() {
            static_assert(Bytes == 1 || Bytes == 2 || Bytes == 4 || Bytes == 8,
                "Bytes must be a power of 2 below 8.");

            IF_CONSTEXPR(sizeof(unsigned char) == Bytes) {
                return (long double)std::numeric_limits<unsigned char>::max();
            }

            IF_CONSTEXPR(sizeof(unsigned short) == Bytes) {
                return (long double)std::numeric_limits<unsigned short>::max();
            }

            IF_CONSTEXPR(sizeof(unsigned int) == Bytes) {
                return (long double)std::numeric_limits<unsigned int>::max();
            }

            IF_CONSTEXPR(sizeof(unsigned long int) == Bytes) {
                return (long double)std::numeric_limits<unsigned long int>::max();
            }

            IF_CONSTEXPR(sizeof(unsigned long long int) == Bytes) {
                return (long double)std::numeric_limits<unsigned long long int>::max();
            }

            HEDLEY_UNREACHABLE();
        }

        /** Largest number that can be stored in a 8-bit integer */
        CONSTEXPR_VALUE long double CSV_INT8_MAX = get_int_max<1>();

        /** Largest number that can be stored in a 16-bit integer */
        CONSTEXPR_VALUE long double CSV_INT16_MAX = get_int_max<2>();

        /** Largest number that can be stored in a 32-bit integer */
        CONSTEXPR_VALUE long double CSV_INT32_MAX = get_int_max<4>();

        /** Largest number that can be stored in a 64-bit integer */
        CONSTEXPR_VALUE long double CSV_INT64_MAX = get_int_max<8>();

        /** Largest number that can be stored in a 8-bit ungisned integer */
        CONSTEXPR_VALUE long double CSV_UINT8_MAX = get_uint_max<1>();

        /** Largest number that can be stored in a 16-bit unsigned integer */
        CONSTEXPR_VALUE long double CSV_UINT16_MAX = get_uint_max<2>();

        /** Largest number that can be stored in a 32-bit unsigned integer */
        CONSTEXPR_VALUE long double CSV_UINT32_MAX = get_uint_max<4>();

        /** Largest number that can be stored in a 64-bit unsigned integer */
        CONSTEXPR_VALUE long double CSV_UINT64_MAX = get_uint_max<8>();

        /** Given a pointer to the start of what is start of
         *  the exponential part of a number written (possibly) in scientific notation
         *  parse the exponent
         */
        HEDLEY_PRIVATE CONSTEXPR
        DataType _process_potential_exponential(
            csv::string_view exponential_part,
            const long double& coeff,
            long double * const out) {
            long double exponent = 0;
            auto result = data_type(exponential_part, &exponent);

            // Exponents in scientific notation should not be decimal numbers
            if (result >= CSV_INT8 && result < CSV_DOUBLE) {
                if (out) *out = coeff * pow10(exponent);
                return CSV_DOUBLE;
            }

            return CSV_STRING;
        }

        /** Given the absolute value of an integer, determine what numeric type
         *  it fits in
         */
        HEDLEY_PRIVATE HEDLEY_PURE CONSTEXPR
        DataType _determine_integral_type(const long double& number) noexcept {
            // We can assume number is always non-negative
            assert(number >= 0);

            if (number <= internals::CSV_INT8_MAX)
                return CSV_INT8;
            else if (number <= internals::CSV_INT16_MAX)
                return CSV_INT16;
            else if (number <= internals::CSV_INT32_MAX)
                return CSV_INT32;
            else if (number <= internals::CSV_INT64_MAX)
                return CSV_INT64;
            else // Conversion to long long will cause an overflow
                return CSV_DOUBLE;
        }

        /** Distinguishes numeric from other text values. Used by various
         *  type casting functions, like csv_parser::CSVReader::read_row()
         *
         *  #### Rules
         *   - Leading and trailing whitespace ("padding") ignored
         *   - A string of just whitespace is NULL
         *
         *  @param[in]  in  String value to be examined
         *  @param[out] out Pointer to long double where results of numeric parsing
         *                  get stored
         */
        CONSTEXPR
        DataType data_type(csv::string_view in, long double* const out) {
            // Empty string --> NULL
            if (in.size() == 0)
                return CSV_NULL;

            bool ws_allowed = true,
                neg_allowed = true,
                dot_allowed = true,
                digit_allowed = true,
                has_digit = false,
                prob_float = false;

            unsigned places_after_decimal = 0;
            long double integral_part = 0,
                decimal_part = 0;

            for (size_t i = 0, ilen = in.size(); i < ilen; i++) {
                const char& current = in[i];

                switch (current) {
                case ' ':
                    if (!ws_allowed) {
                        if (isdigit(in[i - 1])) {
                            digit_allowed = false;
                            ws_allowed = true;
                        }
                        else {
                            // Ex: '510 123 4567'
                            return CSV_STRING;
                        }
                    }
                    break;
                case '-':
                    if (!neg_allowed) {
                        // Ex: '510-123-4567'
                        return CSV_STRING;
                    }

                    neg_allowed = false;
                    break;
                case '.':
                    if (!dot_allowed) {
                        return CSV_STRING;
                    }

                    dot_allowed = false;
                    prob_float = true;
                    break;
                case 'e':
                case 'E':
                    // Process scientific notation
                    if (prob_float || (i && i + 1 < ilen && isdigit(in[i - 1]))) {
                        size_t exponent_start_idx = i + 1;
                        prob_float = true;

                        // Strip out plus sign
                        if (in[i + 1] == '+') {
                            exponent_start_idx++;
                        }

                        return _process_potential_exponential(
                            in.substr(exponent_start_idx),
                            neg_allowed ? integral_part + decimal_part : -(integral_part + decimal_part),
                            out
                        );
                    }

                    return CSV_STRING;
                    break;
                default:
                    short digit = current - '0';
                    if (digit >= 0 && digit <= 9) {
                        // Process digit
                        has_digit = true;

                        if (!digit_allowed)
                            return CSV_STRING;
                        else if (ws_allowed) // Ex: '510 456'
                            ws_allowed = false;

                        // Build current number
                        if (prob_float)
                            decimal_part += digit / pow10(++places_after_decimal);
                        else
                            integral_part = (integral_part * 10) + digit;
                    }
                    else {
                        return CSV_STRING;
                    }
                }
            }

            // No non-numeric/non-whitespace characters found
            if (has_digit) {
                long double number = integral_part + decimal_part;
                if (out) {
                    *out = neg_allowed ? number : -number;
                }

                return prob_float ? CSV_DOUBLE : _determine_integral_type(number);
            }

            // Just whitespace
            return CSV_NULL;
        }
    }
}