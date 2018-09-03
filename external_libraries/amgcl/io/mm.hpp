#ifndef AMGCL_IO_MM_HPP
#define AMGCL_IO_MM_HPP

/*
The MIT License

Copyright (c) 2012-2018 Denis Demidov <dennis.demidov@gmail.com>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
*/

/**
 * \file   amgcl/io/mm.hpp
 * \author Denis Demidov <dennis.demidov@gmail.com>
 * \brief  Readers for Matrix Market sparse matrices and dense vectors.
 */

#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <numeric>

#include <type_traits>
#include <tuple>

#include <amgcl/util.hpp>
#include <amgcl/backend/interface.hpp>
#include <amgcl/value_type/interface.hpp>
#include <amgcl/detail/sort_row.hpp>

namespace amgcl {
namespace io {

/// Matrix market reader.
class mm_reader {
    public:
        /// Open the file by name
        mm_reader(const std::string &fname) : f(fname.c_str()) {
            precondition(f, "Failed to open file \"" + fname + "\"");

            // Read banner.
            std::string line;
            precondition(std::getline(f, line), format_error());

            std::istringstream is(line);
            std::string banner, mtx, coord, dtype, storage;

            precondition(
                    is >> banner >> mtx >> coord >> dtype >> storage,
                    format_error());

            precondition(banner  == "%%MatrixMarket", format_error("no banner"));
            precondition(mtx     == "matrix",         format_error("not a matrix"));

            if (storage == "general") {
                _symmetric = false;
            } else if (storage == "symmetric") {
                _symmetric = true;
            } else {
                precondition(false, "unsupported storage type");
            }

            if (coord == "coordinate") {
                _sparse = true;
            } else if (coord == "array") {
                _sparse = false;
            } else {
                precondition(false, format_error("unsupported coordinate type"));
            }

            if (dtype == "real") {
                _complex = false;
                _integer = false;
            } else if (dtype == "complex") {
                _complex = true;
                _integer = false;
            } else if (dtype == "integer") {
                _complex = false;
                _integer = true;
            } else {
                precondition(false, format_error("unsupported data type"));
            }

            // Skip comments.
            std::streampos pos;
            do {
                pos = f.tellg();
                precondition(std::getline(f, line), format_error("unexpected eof"));
            } while (line[0] == '%');

            // Get back to the first non-comment line.
            f.seekg(pos);

            // Read matrix size
            is.clear(); is.str(line);
            precondition(is >> nrows >> ncols, format_error());
        }

        /// Matrix in the file is symmetric.
        bool is_symmetric()  const { return _symmetric; }

        /// Matrix in the file is sparse.
        bool is_sparse()  const { return _sparse; }

        /// Matrix in the file is complex-valued.
        bool is_complex() const { return _complex; }

        /// Matrix in the file is integer-valued.
        bool is_integer() const { return _integer; }

        /// Number of rows.
        size_t rows() const { return nrows; }

        /// Number of rows.
        size_t cols() const { return ncols; }

        /// Read sparse matrix from the file.
        template <typename Idx, typename Val>
        std::tuple<size_t, size_t> operator()(
                std::vector<Idx> &ptr,
                std::vector<Idx> &col,
                std::vector<Val> &val,
                ptrdiff_t row_beg = -1,
                ptrdiff_t row_end = -1
                )
        {
            precondition(_sparse, format_error("not a sparse matrix"));
            precondition(amgcl::is_complex<Val>::value == _complex,
                    _complex ?
                        "attempt to read complex values into real vector" :
                        "attempt to read real values into complex vector"
                        );
            precondition(std::is_integral<Val>::value == _integer,
                    _integer ?
                        "attempt to read integer values into real vector" :
                        "attempt to read real values into integer vector"
                        );

            // Read sizes
            ptrdiff_t n, m;
            size_t nnz;
            std::string line;
            std::istringstream is;
            {
                precondition(std::getline(f, line), format_error("unexpected eof"));
                is.clear(); is.str(line);
                precondition(is >> n >> m >> nnz, format_error());
            }

            if (row_beg < 0) row_beg = 0;
            if (row_end < 0) row_end = n;

            precondition(row_beg >= 0 && row_end <= n,
                    "Wrong subset of rows is requested");

            ptrdiff_t _nnz = _symmetric ? 2 * nnz : nnz;

            if (row_beg != 0 || row_end != n)
                _nnz *= 1.2 * (row_end - row_beg) / n;

            std::vector<Idx> _row; _row.reserve(_nnz);
            std::vector<Idx> _col; _col.reserve(_nnz);
            std::vector<Val> _val; _val.reserve(_nnz);

            ptrdiff_t chunk = row_end - row_beg;

            ptr.resize(chunk + 1); std::fill(ptr.begin(), ptr.end(), 0);

            for(size_t k = 0; k < nnz; ++k) {
                precondition(std::getline(f, line), format_error("unexpected eof"));
                is.clear(); is.str(line);

                Idx i, j;
                Val v;

                precondition(is >> i >> j, format_error());

                i -= 1;
                j -= 1;

                v = read_value<Val>(is);

                if (row_beg <= i && i < row_end) {
                    ++ptr[i - row_beg + 1];

                    _row.push_back(i - row_beg);
                    _col.push_back(j);
                    _val.push_back(v);
                }

                if (_symmetric && i != j && row_beg <= j && j < row_end) {
                    ++ptr[j - row_beg + 1];

                    _row.push_back(j - row_beg);
                    _col.push_back(i);
                    _val.push_back(v);
                }
            }

            std::partial_sum(ptr.begin(), ptr.end(), ptr.begin());

            col.resize(ptr.back());
            val.resize(ptr.back());

            for(size_t k = 0, e = val.size(); k < e; ++k) {
                Idx i = _row[k];
                Idx j = _col[k];
                Val v = _val[k];

                Idx head = ptr[i]++;
                col[head] = j;
                val[head] = v;
            }

            std::rotate(ptr.begin(), ptr.end() - 1, ptr.end());
            ptr.front() = 0;

#pragma omp parallel for
            for(ptrdiff_t i = 0; i < chunk; ++i) {
                Idx beg = ptr[i];
                Idx end = ptr[i+1];

                amgcl::detail::sort_row(&col[0] + beg, &val[0] + beg, end - beg);
            }

            return std::make_tuple(chunk, m);
        }

        /// Read dense array from the file.
        template <typename Val>
        std::tuple<size_t, size_t> operator()(
                std::vector<Val> &val,
                ptrdiff_t row_beg = -1,
                ptrdiff_t row_end = -1
                )
        {
            precondition(!_sparse, format_error("not a dense array"));
            precondition(amgcl::is_complex<Val>::value == _complex,
                    _complex ?
                        "attempt to read complex values into real vector" :
                        "attempt to read real values into complex vector"
                        );
            precondition(std::is_integral<Val>::value == _integer,
                    _integer ?
                        "attempt to read integer values into real vector" :
                        "attempt to read real values into integer vector"
                        );

            // Read sizes
            ptrdiff_t n, m;
            std::string line;
            std::istringstream is;
            {
                precondition(std::getline(f, line), format_error("unexpected eof"));
                is.clear(); is.str(line);
                precondition(is >> n >> m, format_error());
            }

            if (row_beg < 0) row_beg = 0;
            if (row_end < 0) row_end = n;

            precondition(row_beg >= 0 && row_end <= n,
                    "Wrong subset of rows is requested");

            val.resize((row_end - row_beg) * m);

            for(ptrdiff_t j = 0; j < m; ++j) {
                for(ptrdiff_t i = 0; i < n; ++i) {
                    precondition(std::getline(f, line), format_error("unexpected eof"));
                    if (row_beg <= i && i < row_end) {
                        is.clear(); is.str(line);
                        val[(i - row_beg) * m + j] = read_value<Val>(is);
                    }
                }
            }

            return std::make_tuple(row_end - row_beg, m);
        }
    private:
        std::ifstream f;

        bool _sparse;
        bool _symmetric;
        bool _complex;
        bool _integer;

        size_t nrows, ncols;

        std::string format_error(const std::string &msg = "") const {
            std::string err_string = "MatrixMarket format error";
            if (!msg.empty())
                err_string += " (" + msg + ")";
            return err_string;
        }

        template <typename T>
        typename std::enable_if<amgcl::is_complex<T>::value, T>::type
        read_value(std::istream &s) {
            typename math::scalar_of<T>::type x,y;
            precondition(s >> x >> y, format_error());
            return T(x,y);
        }

        template <typename T>
        typename std::enable_if<!amgcl::is_complex<T>::value, T>::type
        read_value(std::istream &s) {
            T x;
            if (std::is_same<T, char>::value) {
                // Special case:
                // We want to read 8bit integers from MatrixMarket, not chars.
                int i;
                precondition(s >> i, format_error());
                x = static_cast<char>(i);
            } else {
                precondition(s >> x, format_error());
            }
            return x;
        }

};

namespace detail {
template <typename Val>
typename std::enable_if<is_complex<Val>::value, std::ostream&>::type
write_value(std::ostream &s, Val v) {
    return s << std::real(v) << " " << std::imag(v);
}

template <typename Val>
typename std::enable_if<!is_complex<Val>::value, std::ostream&>::type
write_value(std::ostream &s, Val v) {
    return s << v;
}

} // namespace detail

/// Write dense array in Matrix Market format.
template <typename Val>
void mm_write(
        const std::string &fname,
        const Val *data,
        size_t rows,
        size_t cols = 1
        )
{
    std::ofstream f(fname.c_str());
    precondition(f, "Failed to open file \"" + fname + "\" for writing");

    // Banner
    f << "%%MatrixMarket matrix array ";
    if (is_complex<Val>::value) {
        f << "complex ";
    } else if(std::is_integral<Val>::value) {
        f << "integer ";
    } else {
        f << "real ";
    }
    f << "general\n";

    // Sizes
    f << rows << " " << cols << "\n";

    // Data
    for(size_t j = 0; j < cols; ++j) {
        for(size_t i = 0; i < rows; ++i) {
            detail::write_value(f, data[i * cols + j]) << "\n";
        }
    }
}

/// Write sparse matrix in Matrix Market format.
template <class Matrix>
void mm_write(const std::string &fname, const Matrix &A) {
    typedef typename backend::value_type<Matrix>::type Val;

    const size_t rows = backend::rows(A);
    const size_t cols = backend::cols(A);
    const size_t nnz  = backend::nonzeros(A);

    std::ofstream f(fname.c_str());
    precondition(f, "Failed to open file \"" + fname + "\" for writing");

    // Banner
    f << "%%MatrixMarket matrix coordinate ";
    if (is_complex<Val>::value) {
        f << "complex ";
    } else if(std::is_integral<Val>::value) {
        f << "integer ";
    } else {
        f << "real ";
    }
    f << "general\n";

    // Sizes
    f << rows << " " << cols << " " << nnz << "\n";

    // Data
    for(size_t i = 0; i < rows; ++i) {
        for(auto a = backend::row_begin(A, i); a; ++a) {
            f << i + 1 << " " << a.col() + 1 << " ";
            detail::write_value(f, a.value()) << "\n";
        }
    }
}

} // namespace io
} // namespace amgcl


#endif
