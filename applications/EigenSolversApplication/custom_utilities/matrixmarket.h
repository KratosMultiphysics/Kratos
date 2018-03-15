/*
//  KRATOS _______
//        / ____(_)___ ____  ____
//       / __/ / / __ `/ _ \/ __ \
//      / /___/ / /_/ /  __/ / / /
//     /_____/_/\__, /\___/_/ /_/ SolversApplication
//             /____/
//
//  Author: Thomas Oberbichler
*/

#pragma once

#include <cctype>
#include <complex>
#include <exception>
#include <fstream>
#include <string>
#include <tuple>
#include <vector>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>

class MatrixMarket
{
    static std::vector<std::string> split(
        const std::string& str,
        const std::string& delimiters = "\n\r\t "
    )
    {
        std::vector<std::string> tokens;

        // skip delimiters at beginning
        std::string::size_type lastPos = str.find_first_not_of(delimiters, 0);

        // find first "non-delimiter"
        std::string::size_type pos     = str.find_first_of(delimiters, lastPos);

        while (std::string::npos != pos || std::string::npos != lastPos)
        {
            // Found a token, add it to the vector.
            tokens.push_back(str.substr(lastPos, pos - lastPos));

            // Skip delimiters.  Note the "not_of"
            lastPos = str.find_first_not_of(delimiters, pos);

            // Find next "non-delimiter"
            pos = str.find_first_of(delimiters, lastPos);
        }

        return tokens;
    }

    template <typename TMatrix>
    struct matrix_builder { };

    template <typename Scalar>
    struct matrix_builder<boost::numeric::ublas::compressed_matrix<Scalar>>
    {
        using matrix_t = boost::numeric::ublas::compressed_matrix<Scalar>;

        using scalar_t = Scalar;

        using index_t = typename matrix_t::size_type;

        matrix_t& m_matrix;

        matrix_builder(
            matrix_t& matrix
        ) : m_matrix(matrix)
        { }

        void begin_coordinate(
            const index_t& rows,
            const index_t& cols,
            const index_t& nnz
        )
        {
            m_matrix.resize(rows, cols);
            m_matrix.reserve(nnz);
        }

        void begin_array(
            const index_t& rows,
            const index_t& cols
        )
        {
            m_matrix.resize(rows, cols);
            m_matrix.reserve(rows * cols);
        }

        void add_value(
            const index_t& row,
            const index_t& col,
            const scalar_t& value
        )
        {
            m_matrix.insert_element(row, col, value);
        }
    };

    template <typename Scalar>
    struct matrix_builder<boost::numeric::ublas::matrix<Scalar>>
    {
        using matrix_t = boost::numeric::ublas::matrix<Scalar>;

        using scalar_t = Scalar;

        using index_t = typename matrix_t::size_type;

        matrix_t& m_matrix;

        matrix_builder(
            matrix_t& matrix
        ) : m_matrix(matrix)
        { }

        void begin_coordinate(
            const index_t& rows,
            const index_t& cols,
            const index_t& nnz
        )
        {
            m_matrix.resize(rows, cols);
            m_matrix.clear();
        }

        void begin_array(
            const index_t& rows,
            const index_t& cols
        )
        {
            m_matrix.resize(rows, cols);
        }

        void add_value(
            const index_t& row,
            const index_t& col,
            const scalar_t& value
        )
        {
            m_matrix(row, col) = value;
        }
    };

    template <typename T>
    struct read;

    template <typename T>
    struct read<typename std::complex<T>>
    {
        template <typename TStream>
        static std::complex<T> from(TStream& input)
        {
            using complex_t = std::complex<T>;

            T re(1);
            T im(0);

            input >> re;
            input >> im;

            return complex_t(re, im);
        }
    };

    template <typename T>
    struct read
    {
        template <typename TStream>
        static T from(TStream& input)
        {
            T value(1);

            input >> value;

            return value;
        }
    };

    template <typename TIndex = int>
    struct read_coordinates
    {
        template <typename TStream>
        static std::pair<TIndex, TIndex> from(TStream& input)
        {
            TIndex row;
            TIndex col;

            input >> row;
            input >> col;

            return std::make_pair(row - 1, col - 1);
        }
    };

    template <typename TStream>
    static void get_data_line(TStream& input, std::string& line)
    {
        do {
            std::getline(input, line);
        } while (line[0] == '%');
    }

  public:
    template <typename TMatrix, typename TStream>
    static void read_stream(TMatrix& matrix, TStream& input)
    {
        using scalar_t = typename matrix_builder<TMatrix>::scalar_t;
        using index_t = typename matrix_builder<TMatrix>::index_t;

        // --- read banner

        std::string line;
        std::vector<std::string> tokens;

        std::getline(input, line);
        tokens = split(line);

        if (tokens.size() != 5 || tokens[0] != "%%MatrixMarket" || tokens[1] != "matrix") {
            throw std::runtime_error("MatrixMarket banner invalid");
        }

        const auto storage  = tokens[2];
        const auto type = tokens[3];
        const auto symmetry = tokens[4];

        if (storage != "array" && storage != "coordinate") {
            throw std::runtime_error("MatrixMarket storage format '" + storage + "' invalid");
        }

        if (type != "pattern" && type != "integer" && type != "real" && type != "complex") {
            throw std::runtime_error("MatrixMarket data type '" + type + "' invalid");
        }

        if (symmetry != "general" && symmetry != "symmetric" && symmetry != "hermitian" && symmetry != "skew-symmetric") {
            throw std::runtime_error("MatrixMarket symmetry '" + symmetry + "' invalid");
        }


        // --- read matrix size

        get_data_line(input, line);



        // --- read entries

        matrix_builder<TMatrix> builder(matrix);

        if (storage == "coordinate")
        {
            tokens = split(line);

            if (tokens.size() != 3) {
                throw std::runtime_error("MatrixMarket format invalid");
            }

            index_t rows, cols, entries;

            std::istringstream(tokens[0]) >> rows;
            std::istringstream(tokens[1]) >> cols;
            std::istringstream(tokens[2]) >> entries;

            builder.begin_coordinate(rows, cols, entries);

            size_t entries_read = 0;

            while (entries_read < entries && !input.eof()) {
                get_data_line(input, line);
                std::istringstream is(line);

                int row;
                int col;

                std::tie(row, col) = read_coordinates<index_t>::from(is);

                scalar_t value = read<scalar_t>::from(is);

                builder.add_value(row, col, value);

                entries_read += 1;
            }
        } else { // storage == "array"
            tokens = split(line);

            if (tokens.size() != 2) {
                throw std::runtime_error("MatrixMarket format invalid");
            }

            index_t rows, cols;

            std::istringstream(tokens[0]) >> rows;
            std::istringstream(tokens[1]) >> cols;

            builder.begin_array(rows, cols);

            for (size_t i = 0; i != cols && !input.eof(); ++i) {
                for (size_t j = 0; j != rows && !input.eof(); ++j) {
                    get_data_line(input, line);
                    std::istringstream is(line);

                    scalar_t value = read<scalar_t>::from(is);

                    builder.add_value(i, j, value);
                }
            }
        }
    }

    template <typename TMatrix>
    static void read_file(TMatrix& matrix, const std::string& filename)
    {
        std::ifstream file(filename.c_str());

        if (!file) {
            throw std::exception();
        }

        read_stream(matrix, file);
    }
};