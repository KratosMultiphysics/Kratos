#ifndef VEXCL_SPMAT_HYBRID_ELL_INL
#define VEXCL_SPMAT_HYBRID_ELL_INL

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
 * \file   vexcl/spmat/hybrid_ell.inl
 * \author Denis Demidov <dennis.demidov@gmail.com>
 * \brief  OpenCL sparse matrix in Hybrid ELL-CSR format.
 */

struct SpMatHELL : public sparse_matrix {
    const backend::command_queue &queue;
    size_t n, pitch;

    struct matrix_part {
        struct {
            size_t     width;
            backend::device_vector<col_t> col;
            backend::device_vector<val_t> val;
        } ell;

        struct {
            size_t     nnz;
            backend::device_vector<idx_t> row;
            backend::device_vector<col_t> col;
            backend::device_vector<val_t> val;
        } csr;
    } loc, rem;

    SpMatHELL(
            const backend::command_queue &queue,
            const idx_t *row_begin, const idx_t *row_end,
            const col_t *col, const val_t *val,
            size_t col_begin, size_t col_end,
            std::set<col_t> ghost_cols
            )
        : queue(queue), n(row_end - row_begin), pitch( alignup(n, 16U) )
    {
        auto is_local = [col_begin, col_end](size_t c) {
            return c >= col_begin && c < col_end;
        };

        /* 1. Get optimal ELL widths for local and remote parts. */
        {
            // Speed of ELL relative to CSR (e.g. 2.0 -> ELL is twice as fast):
            const double ell_vs_csr = 3.0;

            // Find maximum widths for local and remote parts:
            loc.ell.width = rem.ell.width = 0;
            for(auto row = row_begin; row != row_end; ++row) {
                size_t wl = 0, wr = 0;
                for(idx_t j = row[0]; j < row[1]; ++j) {
                    if (is_local(col[j]))
                        ++wl;
                    else
                        ++wr;
                }

                loc.ell.width = std::max(loc.ell.width, wl);
                rem.ell.width = std::max(rem.ell.width, wr);
            }

            // Build histograms for width distribution.
            std::vector<size_t> loc_hist(loc.ell.width + 1, 0);
            std::vector<size_t> rem_hist(rem.ell.width + 1, 0);

            for(auto row = row_begin; row != row_end; ++row) {
                size_t wl = 0, wr = 0;
                for(idx_t j = row[0]; j < row[1]; ++j) {
                    if (is_local(col[j]))
                        ++wl;
                    else
                        ++wr;
                }

                ++loc_hist[wl];
                ++rem_hist[wr];
            }

            auto optimal_width = [&](size_t max_width, const std::vector<size_t> &hist) -> size_t {
                for(size_t i = 0, rows = n; i < max_width; ++i) {
                    rows -= hist[i]; // Number of rows wider than i.
                    if (ell_vs_csr * rows < n) return i;
                }

                return max_width;
            };

            loc.ell.width = optimal_width(loc.ell.width, loc_hist);
            rem.ell.width = optimal_width(rem.ell.width, rem_hist);
        }

        /* 2. Count nonzeros in CSR parts of the matrix. */
        loc.csr.nnz = 0, rem.csr.nnz = 0;

        for(auto row = row_begin; row != row_end; ++row) {
            size_t wl = 0, wr = 0;
            for(idx_t j = row[0]; j < row[1]; ++j) {
                if (is_local(col[j]))
                    ++wl;
                else
                    ++wr;
            }

            if (wl > loc.ell.width) loc.csr.nnz += wl - loc.ell.width;
            if (wr > rem.ell.width) rem.csr.nnz += wr - rem.ell.width;
        }

        /* 3. Renumber columns. */
        std::unordered_map<col_t,col_t> r2l(2 * ghost_cols.size());
        size_t nghost = 0;
        for(auto c = ghost_cols.begin(); c != ghost_cols.end(); c++)
            r2l[*c] = static_cast<col_t>(nghost++);

        // Prepare ELL and COO formats for transfer to devices.
        const col_t not_a_column = static_cast<col_t>(-1);

        std::vector<col_t> lell_col(pitch * loc.ell.width, not_a_column);
        std::vector<val_t> lell_val(pitch * loc.ell.width, val_t());
        std::vector<col_t> rell_col(pitch * rem.ell.width, not_a_column);
        std::vector<val_t> rell_val(pitch * rem.ell.width, val_t());

        std::vector<idx_t> lcsr_row;
        std::vector<col_t> lcsr_col;
        std::vector<val_t> lcsr_val;

        std::vector<idx_t> rcsr_row;
        std::vector<col_t> rcsr_col;
        std::vector<val_t> rcsr_val;

        lcsr_row.reserve(n + 1);
        lcsr_col.reserve(loc.csr.nnz);
        lcsr_val.reserve(loc.csr.nnz);

        rcsr_row.reserve(n + 1);
        rcsr_col.reserve(rem.csr.nnz);
        rcsr_val.reserve(rem.csr.nnz);

        lcsr_row.push_back(0);
        rcsr_row.push_back(0);

        size_t k = 0;
        for(auto row = row_begin; row != row_end; ++row, ++k) {
            size_t lcnt = 0, rcnt = 0;
            for(idx_t j = row[0]; j < row[1]; ++j) {
                if (is_local(col[j])) {
                    if (lcnt < loc.ell.width) {
                        lell_col[k + pitch * lcnt] = static_cast<col_t>(col[j] - col_begin);
                        lell_val[k + pitch * lcnt] = val[j];
                        ++lcnt;
                    } else {
                        lcsr_col.push_back(static_cast<col_t>(col[j] - col_begin));
                        lcsr_val.push_back(val[j]);
                    }
                } else {
                    assert(r2l.count(col[j]));
                    if (rcnt < rem.ell.width) {
                        rell_col[k + pitch * rcnt] = r2l[col[j]];
                        rell_val[k + pitch * rcnt] = val[j];
                        ++rcnt;
                    } else {
                        rcsr_col.push_back(r2l[col[j]]);
                        rcsr_val.push_back(val[j]);
                    }
                }
            }

            lcsr_row.push_back(static_cast<idx_t>(lcsr_col.size()));
            rcsr_row.push_back(static_cast<idx_t>(rcsr_col.size()));
        }

        /* Copy data to device */
        if (loc.ell.width) {
            loc.ell.col = backend::device_vector<col_t>(queue, lell_col.size(), lell_col.data());
            loc.ell.val = backend::device_vector<val_t>(queue, lell_val.size(), lell_val.data());
        }

        if (loc.csr.nnz) {
            loc.csr.row = backend::device_vector<idx_t>(queue, lcsr_row.size(), lcsr_row.data());
            loc.csr.col = backend::device_vector<col_t>(queue, lcsr_col.size(), lcsr_col.data());
            loc.csr.val = backend::device_vector<val_t>(queue, lcsr_val.size(), lcsr_val.data());
        }

        if (rem.ell.width) {
            rem.ell.col = backend::device_vector<col_t>(queue, rell_col.size(), rell_col.data());
            rem.ell.val = backend::device_vector<val_t>(queue, rell_val.size(), rell_val.data());
        }

        if (rem.csr.nnz) {
            rem.csr.row = backend::device_vector<idx_t>(queue, rcsr_row.size(), rcsr_row.data());
            rem.csr.col = backend::device_vector<col_t>(queue, rcsr_col.size(), rcsr_col.data());
            rem.csr.val = backend::device_vector<val_t>(queue, rcsr_val.size(), rcsr_val.data());
        }
    }

    template <class OP>
    void mul(
            const matrix_part &part,
            const backend::device_vector<val_t> &in,
            backend::device_vector<val_t> &out,
            scalar_type scale
            ) const
    {
        using namespace detail;

        static kernel_cache cache;

        auto kernel = cache.find(queue);

        backend::select_context(queue);

        if (kernel == cache.end()) {
            backend::source_generator source(queue);

            source.begin_kernel("hybrid_ell_spmv");
            source.begin_kernel_parameters();
            source.template parameter<size_t>("n");
            source.template parameter<scalar_type>("scale");
            source.template parameter<size_t>("ell_w");
            source.template parameter<size_t>("ell_pitch");
            source.template parameter< global_ptr<const col_t> >("ell_col");
            source.template parameter< global_ptr<const val_t> >("ell_val");
            source.template parameter< global_ptr<const idx_t> >("csr_row");
            source.template parameter< global_ptr<const col_t> >("csr_col");
            source.template parameter< global_ptr<const val_t> >("csr_val");
            source.template parameter< global_ptr<const val_t> >("in");
            source.template parameter< global_ptr<val_t> >("out");
            source.end_kernel_parameters();
            source.grid_stride_loop("i").open("{");

            source.new_line() << type_name<val_t>() << " sum = 0;";
            source.new_line() << "for(size_t j = 0; j < ell_w; ++j)";
            source.open("{");
            source.new_line() << type_name<col_t>() << " c = ell_col[i + j * ell_pitch];";
            source.new_line() << "if (c != ("<< type_name<col_t>() << ")(-1))";
            source.open("{").new_line() << "sum += ell_val[i + j * ell_pitch] * in[c];";
            source.close("}").close("}");
            source.new_line() << "if (csr_row)";
            source.open("{");
            source.new_line() << "for(size_t j = csr_row[i], e = csr_row[i + 1]; j < e; ++j)";
            source.open("{");
            source.new_line() << "sum += csr_val[j] * in[csr_col[j]];";
            source.close("}").close("}");
            source.new_line() << "out[i] " << OP::string() << " scale * sum;";
            source.close("}");
            source.end_kernel();

            kernel = cache.insert(queue, backend::kernel(
                        queue, source.str(), "hybrid_ell_spmv"));
        }

        kernel->second.push_arg(n);
        kernel->second.push_arg(scale);
        kernel->second.push_arg(part.ell.width);
        kernel->second.push_arg(pitch);

        if (part.ell.width) {
            kernel->second.push_arg(part.ell.col);
            kernel->second.push_arg(part.ell.val);
        } else {
            kernel->second.push_arg(static_cast<size_t>(0));
            kernel->second.push_arg(static_cast<size_t>(0));
        }

        if (part.csr.nnz) {
            kernel->second.push_arg(part.csr.row);
            kernel->second.push_arg(part.csr.col);
            kernel->second.push_arg(part.csr.val);
        } else {
            kernel->second.push_arg(static_cast<size_t>(0));
            kernel->second.push_arg(static_cast<size_t>(0));
            kernel->second.push_arg(static_cast<size_t>(0));
        }
        kernel->second.push_arg(in);
        kernel->second.push_arg(out);

        kernel->second(queue);
    }

    void mul_local(
            const backend::device_vector<val_t> &in,
            backend::device_vector<val_t> &out,
            scalar_type scale, bool append) const
    {
        if (append)
            mul<assign::ADD>(loc, in, out, scale);
        else
            mul<assign::SET>(loc, in, out, scale);
    }

    void mul_remote(
            const backend::device_vector<val_t> &in,
            backend::device_vector<val_t> &out,
            scalar_type scale) const
    {
        mul<assign::ADD>(rem, in, out, scale);
    }

    static void inline_preamble(backend::source_generator &src,
        const std::string &prm_name)
    {
        src.begin_function<val_t>(prm_name + "_hell_spmv");
        src.begin_function_parameters();
        src.template parameter<size_t>("ell_w");
        src.template parameter<size_t>("ell_pitch");
        src.template parameter< global_ptr<const col_t> >("ell_col");
        src.template parameter< global_ptr<const val_t> >("ell_val");
        src.template parameter< global_ptr<const idx_t> >("csr_row");
        src.template parameter< global_ptr<const col_t> >("csr_col");
        src.template parameter< global_ptr<const val_t> >("csr_val");
        src.template parameter< global_ptr<const val_t> >("in");
        src.template parameter< size_t >("i");
        src.end_function_parameters();
        src.new_line() << type_name<val_t>() << " sum = 0;";
        src.new_line() << "for(size_t j = 0; j < ell_w; ++j)";
        src.open("{");
        src.new_line() << type_name<col_t>() << " c = ell_col[i + j * ell_pitch];";
        src.new_line() << "if (c != ("<< type_name<col_t>() << ")(-1))";
        src.open("{").new_line() << "sum += ell_val[i + j * ell_pitch] * in[c];";
        src.close("}").close("}");
        src.new_line() << "if (csr_row)";
        src.open("{");
        src.new_line() << "for(size_t j = csr_row[i], e = csr_row[i + 1]; j < e; ++j)";
        src.open("{").new_line() << "sum += csr_val[j] * in[csr_col[j]];";
        src.close("}").close("}");
        src.new_line() << "return sum;";
        src.end_function();
    }

    static void inline_expression(backend::source_generator &src,
            const std::string &prm_name)
    {
        src << prm_name << "_hell_spmv" << "("
            << prm_name << "_ell_w, "
            << prm_name << "_ell_pitch, "
            << prm_name << "_ell_col, "
            << prm_name << "_ell_val, "
            << prm_name << "_csr_row, "
            << prm_name << "_csr_col, "
            << prm_name << "_csr_val, "
            << prm_name << "_vec, idx)";
    }

    static void inline_parameters(backend::source_generator &src,
            const std::string &prm_name)
    {
        src.template parameter<size_t>(prm_name) << "_ell_w";
        src.template parameter<size_t>(prm_name) << "_ell_pitch";
        src.template parameter< global_ptr<const col_t> >(prm_name + "_ell_col");
        src.template parameter< global_ptr<const val_t> >(prm_name + "_ell_val");
        src.template parameter< global_ptr<const idx_t> >(prm_name + "_csr_row");
        src.template parameter< global_ptr<const col_t> >(prm_name + "_csr_col");
        src.template parameter< global_ptr<const val_t> >(prm_name + "_csr_val");
        src.template parameter< global_ptr<const val_t> >(prm_name + "_vec");
    }

    void setArgs(backend::kernel &krn, unsigned device, const vector<val_t> &x) const {
        krn.push_arg(loc.ell.width);
        krn.push_arg(pitch);
        if (loc.ell.width) {
            krn.push_arg(loc.ell.col);
            krn.push_arg(loc.ell.val);
        } else {
            krn.push_arg(static_cast<size_t>(0));
            krn.push_arg(static_cast<size_t>(0));
        }
        if (loc.csr.nnz) {
            krn.push_arg(loc.csr.row);
            krn.push_arg(loc.csr.col);
            krn.push_arg(loc.csr.val);
        } else {
            krn.push_arg(static_cast<size_t>(0));
            krn.push_arg(static_cast<size_t>(0));
            krn.push_arg(static_cast<size_t>(0));
        }
        krn.push_arg(x(device));
    }
};

#endif
