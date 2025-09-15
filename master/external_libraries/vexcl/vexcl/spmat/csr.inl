#ifndef VEXCL_SPMAT_CSR_INL
#define VEXCL_SPMAT_CSR_INL

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
 * \file   vexcl/spmat/csr.inl
 * \author Denis Demidov <dennis.demidov@gmail.com>
 * \brief  OpenCL sparse matrix in CSR format.
 */

struct SpMatCSR : public sparse_matrix {
    const backend::command_queue &queue;
    size_t n;

    struct matrix_part {
        size_t nnz;
        backend::device_vector<idx_t> row;
        backend::device_vector<col_t> col;
        backend::device_vector<val_t> val;
    } loc, rem;

    SpMatCSR(
            const backend::command_queue &queue,
            const idx_t *row_begin, const idx_t *row_end,
            const col_t *col, const val_t *val,
            col_t col_begin, col_t col_end,
            std::set<col_t> ghost_cols
            )
        : queue(queue), n(row_end - row_begin)
    {
        auto is_local = [col_begin, col_end](col_t c) {
            return c >= col_begin && c < col_end;
        };

        if (ghost_cols.empty()) {
            loc.nnz = *row_end - *row_begin;
            rem.nnz = 0;

            if (loc.nnz) {
                loc.row = backend::device_vector<idx_t>(queue, (n + 1), row_begin,        backend::MEM_READ_ONLY);
                loc.col = backend::device_vector<col_t>(queue, loc.nnz, col + *row_begin, backend::MEM_READ_ONLY);
                loc.val = backend::device_vector<val_t>(queue, loc.nnz, val + *row_begin, backend::MEM_READ_ONLY);

                if (*row_begin > 0) vector<idx_t>(queue, loc.row) -= *row_begin;
            }
        } else {
            std::vector<idx_t> lrow;
            std::vector<col_t> lcol;
            std::vector<val_t> lval;

            std::vector<idx_t> rrow;
            std::vector<col_t> rcol;
            std::vector<val_t> rval;

            lrow.reserve(n + 1);
            lrow.push_back(0);

            lcol.reserve(*row_end - *row_begin);
            lval.reserve(*row_end - *row_begin);

            if (!ghost_cols.empty()) {
                rrow.reserve(n + 1);
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
                loc.nnz = lrow.back();

                loc.row = backend::device_vector<idx_t>(queue, lrow.size(), lrow.data(), backend::MEM_READ_ONLY);
                loc.col = backend::device_vector<col_t>(queue, lcol.size(), lcol.data(), backend::MEM_READ_ONLY);
                loc.val = backend::device_vector<val_t>(queue, lval.size(), lval.data(), backend::MEM_READ_ONLY);
            }

            // Copy remote part to the device.
            if (!ghost_cols.empty()) {
                rem.nnz = rrow.back();

                rem.row = backend::device_vector<idx_t>(queue, rrow.size(), rrow.data(), backend::MEM_READ_ONLY);
                rem.col = backend::device_vector<col_t>(queue, rcol.size(), rcol.data(), backend::MEM_READ_ONLY);
                rem.val = backend::device_vector<val_t>(queue, rval.size(), rval.data(), backend::MEM_READ_ONLY);
            }
        }
    }

    template <class OP>
    void mul(const matrix_part &part,
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

            source.begin_kernel("csr_spmv");
            source.begin_kernel_parameters();
            source.template parameter<size_t>("n");
            source.template parameter<scalar_type>("scale");
            source.template parameter< global_ptr< const idx_t > >("row");
            source.template parameter< global_ptr< const col_t > >("col");
            source.template parameter< global_ptr< const val_t > >("val");
            source.template parameter< global_ptr< const val_t > >("in");
            source.template parameter< global_ptr< val_t > >("out");
            source.end_kernel_parameters();
            source.grid_stride_loop("i").open("{");
            source.new_line() << type_name<val_t>() << " sum = 0;";
            source.new_line() << "for(size_t j = row[i], e = row[i + 1]; j < e; ++j)";
            source.open("{");
            source.new_line() << "sum += val[j] * in[col[j]];";
            source.close("}");
            source.new_line() << "out[i] " << OP::string() << " scale * sum;";
            source.close("}");
            source.end_kernel();

            kernel = cache.insert(queue, backend::kernel(
                        queue, source.str(), "csr_spmv"));
        }

        kernel->second.push_arg(n);
        kernel->second.push_arg(scale);
        kernel->second.push_arg(part.row);
        kernel->second.push_arg(part.col);
        kernel->second.push_arg(part.val);
        kernel->second.push_arg(in);
        kernel->second.push_arg(out);

        kernel->second(queue);
    }

    void mul_local(
            const backend::device_vector<val_t> &in,
            backend::device_vector<val_t> &out,
            scalar_type scale, bool append) const
    {
        if (append) {
            if (loc.nnz) mul<assign::ADD>(loc, in, out, scale);
        } else {
            if (loc.nnz)
                mul<assign::SET>(loc, in, out, scale);
            else
                vector<val_t>(queue, out) = 0;
        }
    }

    void mul_remote(
            const backend::device_vector<val_t> &in,
            backend::device_vector<val_t> &out,
            scalar_type scale) const
    {
        if (rem.nnz) mul<assign::ADD>(rem, in, out, scale);
    }

    static void inline_preamble(backend::source_generator &src,
            const std::string &prm_name)
    {
        src.begin_function<val_t>(prm_name + "_csr_spmv");
        src.begin_function_parameters();
        src.template parameter< global_ptr<const idx_t> >("row");
        src.template parameter< global_ptr<const col_t> >("col");
        src.template parameter< global_ptr<const val_t> >("val");
        src.template parameter< global_ptr<const val_t> >("in");
        src.template parameter< size_t >("i");
        src.end_function_parameters();
        src.new_line() << type_name<val_t>() << " sum = 0;";
        src.new_line() << "for(size_t j = row[i], e = row[i + 1]; j < e; ++j)";
        src.open("{");
        src.new_line() << "sum += val[j] * in[col[j]];";
        src.close("}");
        src.new_line() << "return sum;";
        src.end_function();
    }

    static void inline_expression(backend::source_generator &src,
            const std::string &prm_name)
    {
        src << prm_name << "_csr_spmv" << "("
            << prm_name << "_row, "
            << prm_name << "_col, "
            << prm_name << "_val, "
            << prm_name << "_vec, idx)";
    }

    static void inline_parameters(backend::source_generator &src,
            const std::string &prm_name)
    {
        src.template parameter< global_ptr<const idx_t> >(prm_name + "_row");
        src.template parameter< global_ptr<const col_t> >(prm_name + "_col");
        src.template parameter< global_ptr<const val_t> >(prm_name + "_val");
        src.template parameter< global_ptr<const val_t> >(prm_name + "_vec");
    }

    void setArgs(backend::kernel &krn, unsigned device, const vector<val_t> &x) const {
        krn.push_arg(loc.row);
        krn.push_arg(loc.col);
        krn.push_arg(loc.val);
        krn.push_arg(x(device));
    }
};

#endif
