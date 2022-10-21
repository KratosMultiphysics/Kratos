#ifndef AMGCL_ADAPTER_BLOCK_MATRIX_HPP
#define AMGCL_ADAPTER_BLOCK_MATRIX_HPP

/*
The MIT License

Copyright (c) 2012-2022 Denis Demidov <dennis.demidov@gmail.com>

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
\file    amgcl/adapter/block_matrix.hpp
\author  Denis Demidov <dennis.demidov@gmail.com>
\brief   On-the-fly conversion to block valued matrix.
\ingroup adapters
*/

#include <amgcl/util.hpp>
#include <amgcl/backend/detail/matrix_ops.hpp>
#include <amgcl/backend/builtin.hpp>
#include <amgcl/value_type/interface.hpp>

namespace amgcl {
namespace adapter {

template <class Matrix, class BlockType>
struct block_matrix_adapter {
    typedef BlockType value_type;
    static const int BlockSize = math::static_rows<BlockType>::value;

    const Matrix &A;

    block_matrix_adapter(const Matrix &A) : A(A) {
        precondition(
                backend::rows(A) % BlockSize == 0 &&
                backend::cols(A) % BlockSize == 0,
                "Matrix size is not divisible by block size!"
                );
    }

    size_t rows() const {
        return backend::rows(A) / BlockSize;
    }

    size_t cols() const {
        return backend::cols(A) / BlockSize;
    }

    size_t nonzeros() const {
        // Just an estimate:
        return backend::nonzeros(A) / (BlockSize * BlockSize);
    }

    struct row_iterator {
        typedef typename backend::row_iterator<Matrix>::type Base;
        typedef ptrdiff_t col_type;
        typedef BlockType val_type;

        std::array<char, sizeof(Base) * BlockSize> buf {};
        Base * base;

        bool done;
        col_type cur_col;
        val_type cur_val;

        row_iterator(const Matrix &A, col_type row)
            : base(reinterpret_cast<Base*>(buf.data())), done(true)
        {
            for(int i = 0; i < BlockSize; ++i) {
                new (base + i) Base(backend::row_begin(A, row * BlockSize + i));

                if (base[i]) {
                    col_type col = base[i].col() / BlockSize;
                    if (done) {
                        cur_col = col;
                        done = false;
                    } else {
                        cur_col = std::min<col_type>(cur_col, col);
                    }
                }
            }

            if (done) return;

            // While we are gathering the current value,
            // base iteratirs are advanced to the next block-column.
            cur_val = math::zero<val_type>();
            col_type end = (cur_col + 1) * BlockSize;
            for(int i = 0; i < BlockSize; ++i) {
                for(; base[i] && static_cast<ptrdiff_t>(base[i].col()) < end; ++base[i]) {
                    cur_val(i, base[i].col() % BlockSize) = base[i].value();
                }
            }
        }

        ~row_iterator() {
            for(int i = 0; i < BlockSize; ++i) base[i].~Base();
        }

        operator bool() const {
            return !done;
        }

        row_iterator& operator++() {
            // Base iterators are already at the next block-column.
            // We just need to gather the current column and value.
            done = true;

            col_type end = (cur_col + 1) * BlockSize;
            for(int i = 0; i < BlockSize; ++i) {
                if (base[i]) {
                    col_type col = base[i].col() / BlockSize;
                    if (done) {
                        cur_col = col;
                        done = false;
                    } else {
                        cur_col = std::min<col_type>(cur_col, col);
                    }
                }
            }

            if (done) return *this;

            cur_val = math::zero<val_type>();
            end = (cur_col + 1) * BlockSize;
            for(int i = 0; i < BlockSize; ++i) {
                for(; base[i] && static_cast<ptrdiff_t>(base[i].col()) < end; ++base[i]) {
                    cur_val(i, base[i].col() % BlockSize) = base[i].value();
                }
            }

            return *this;
        }

        col_type col() const {
            return cur_col;
        }

        val_type value() const {
            return cur_val;
        }
    };

    row_iterator row_begin(size_t i) const {
        return row_iterator(A, i);
    }
};

/// Convert scalar-valued matrix to a block-valued one.
template <class BlockType, class Matrix>
block_matrix_adapter<Matrix, BlockType> block_matrix(const Matrix &A) {
    return block_matrix_adapter<Matrix, BlockType>(A);
}

template <class Matrix>
std::shared_ptr<
    backend::crs<
        typename math::element_of<
            typename backend::value_type<Matrix>::type
            >::type,
        typename backend::col_type<Matrix>::type,
        typename backend::ptr_type<Matrix>::type
        >
    >
unblock_matrix(const Matrix &B) {
    typedef typename backend::value_type<Matrix>::type Block;
    typedef typename math::element_of<Block>::type Scalar;
    typedef typename backend::col_type<Matrix>::type Col;
    typedef typename backend::ptr_type<Matrix>::type Ptr;

    const int brows = math::static_rows<Block>::value;
    const int bcols = math::static_cols<Block>::value;

    static_assert(brows > 1 || bcols > 1, "Can not unblock scalar matrix!");

    auto A = std::make_shared<backend::crs<Scalar, Col, Ptr>>();

    A->set_size(backend::rows(B) * brows, backend::cols(B) * bcols);
    A->ptr[0] = 0;

    const ptrdiff_t nb = backend::rows(B);

#pragma omp for
    for (ptrdiff_t ib = 0; ib < nb; ++ib) {
        auto w = backend::row_nonzeros(B, ib);
        for (ptrdiff_t i = 0, ia = ib * brows; i < brows; ++i, ++ia) {
            A->ptr[ia + 1] = w * bcols;
        }
    }

    A->scan_row_sizes();
    A->set_nonzeros();

#pragma omp for
    for (ptrdiff_t ib = 0; ib < nb; ++ib) {
        for(auto b = backend::row_begin(B, ib); b; ++b) {
            auto c = b.col();
            auto v = b.value();

            for (ptrdiff_t i = 0, ia = ib * brows; i < brows; ++i, ++ia) {
                auto row_head = A->ptr[ia];
                for(int j = 0; j < bcols; ++j) {
                    A->col[row_head] = c * bcols + j;
                    A->val[row_head] = v(i,j);
                    ++row_head;
                }
                A->ptr[ia] = row_head;
            }
        }
    }

    std::rotate(A->ptr, A->ptr + A->nrows, A->ptr + A->nrows + 1);
    A->ptr[0] = 0;

    return A;
}

} // namespace adapter

namespace backend {
namespace detail {

template <class Matrix, class BlockType>
struct use_builtin_matrix_ops< adapter::block_matrix_adapter<Matrix, BlockType> >
    : std::true_type
{};

} // namespace detail
} // namespace backend
} // namespace amgcl

#endif
