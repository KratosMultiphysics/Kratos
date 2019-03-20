#ifndef VEXCL_BACKEND_CUDA_CSR_INL
#define VEXCL_BACKEND_CUDA_CSR_INL

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
 * \file   vexcl/backend/cuda/csr.inl
 * \author Denis Demidov <dennis.demidov@gmail.com>
 * \brief  CUSPARSE matrix in CSR format.
 */

struct SpMatCSR : public sparse_matrix {
    const backend::command_queue &queue;
    std::unique_ptr< backend::cuda::spmat_crs<val_t> > loc, rem;

    SpMatCSR(
            const backend::command_queue &queue,
            const idx_t *row_begin, const idx_t *row_end,
            const col_t *col, const val_t *val,
            col_t col_begin, col_t col_end,
            std::set<col_t> ghost_cols
            ) : queue(queue)
    {
        auto is_local = [col_begin, col_end](col_t c) {
            return c >= col_begin && c < col_end;
        };

        if (ghost_cols.empty()) {
            loc.reset(new backend::cuda::spmat_crs<val_t>(queue,
                        static_cast<int>(row_end - row_begin),
                        static_cast<int>(col_end - col_begin),
                        row_begin, col, val
                        ));
        } else {
            std::vector<idx_t> lrow;
            std::vector<col_t> lcol;
            std::vector<val_t> lval;

            std::vector<idx_t> rrow;
            std::vector<col_t> rcol;
            std::vector<val_t> rval;

            lrow.reserve(row_end - row_begin + 1);
            lrow.push_back(0);

            lcol.reserve(*row_end - *row_begin);
            lval.reserve(*row_end - *row_begin);

            if (!ghost_cols.empty()) {
                rrow.reserve(row_end - row_begin + 1);
                rrow.push_back(0);

                rcol.reserve(*row_end - *row_begin);
                rval.reserve(*row_end - *row_begin);
            }

            // Renumber columns.
            std::unordered_map<col_t, col_t> r2l(2 * ghost_cols.size());
            size_t nghost = 0;
            for(auto c = ghost_cols.begin(); c != ghost_cols.end(); ++c)
                r2l[*c] = static_cast<col_t>(nghost++);

            for(auto row = row_begin; row != row_end; ++row) {
                for(idx_t j = row[0]; j < row[1]; j++) {
                    if (is_local(col[j])) {
                        lcol.push_back(static_cast<col_t>(col[j] - col_begin));
                        lval.push_back(val[j]);
                    } else {
                        assert(r2l.count(col[j]));
                        rcol.push_back(r2l[col[j]]);
                        rval.push_back(val[j]);
                    }
                }

                lrow.push_back(static_cast<idx_t>(lcol.size()));
                rrow.push_back(static_cast<idx_t>(rcol.size()));
            }


            // Copy local part to the device.
            if (lrow.back()) {
                loc.reset(new backend::cuda::spmat_crs<val_t>(queue,
                        static_cast<int>(row_end - row_begin),
                        static_cast<int>(col_end - col_begin),
                        lrow.data(), lcol.data(), lval.data()
                        ));
            }

            // Copy remote part to the device.
            if (!ghost_cols.empty()) {
                rem.reset(new backend::cuda::spmat_crs<val_t>(queue,
                            static_cast<int>(row_end - row_begin),
                            static_cast<int>(nghost),
                            rrow.data(), rcol.data(), rval.data()
                            ));
            }
        }
    }

    template <class OP>
    void mul(const backend::cuda::spmat_crs<val_t> &A,
             const backend::device_vector<val_t> &in,
             backend::device_vector<val_t> &out,
             scalar_type scale
             ) const
    {
        A.mul(in, out, scale, std::is_same<OP, assign::ADD>::value);
    }

    void mul_local(
            const backend::device_vector<val_t> &in,
            backend::device_vector<val_t> &out,
            scalar_type scale, bool append) const
    {
        if (append) {
            if (loc) mul<assign::ADD>(*loc, in, out, scale);
        } else {
            if (loc)
                mul<assign::SET>(*loc, in, out, scale);
            else
                vector<val_t>(queue, out) = 0;
        }
    }

    void mul_remote(
            const backend::device_vector<val_t> &in,
            backend::device_vector<val_t> &out,
            scalar_type scale) const
    {
        if (rem) mul<assign::ADD>(*rem, in, out, scale);
    }
};

#endif
