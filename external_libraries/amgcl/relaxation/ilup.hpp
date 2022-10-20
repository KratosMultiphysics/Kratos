#ifndef AMGCL_RELAXATION_ILUP_HPP
#define AMGCL_RELAXATION_ILUP_HPP

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
 * \file   amgcl/relaxation/ilup.hpp
 * \author Denis Demidov <dennis.demidov@gmail.com>
 * \brief  Incomplete LU with fill-in level.
 *
 * As opposed to the iluk, the fill-in is determined by taking a symbolic
 * power of the matrix.
 */

#include <vector>
#include <deque>
#include <queue>
#include <cmath>


#include <amgcl/backend/builtin.hpp>
#include <amgcl/util.hpp>
#include <amgcl/relaxation/ilu0.hpp>

namespace amgcl {
namespace relaxation {
namespace detail {

template <class Matrix>
std::shared_ptr<Matrix> symb_product(const Matrix &A, const Matrix &B) {
    auto C = std::make_shared<Matrix>();

    C->set_size(A.nrows, B.ncols);

    auto A_ptr = A.ptr;
    auto A_col = A.col;
    auto B_ptr = B.ptr;
    auto B_col = B.col;
    auto C_ptr = C->ptr;
    C_ptr[0] = 0;

#pragma omp parallel
    {
        std::vector<ptrdiff_t> marker(B.ncols, -1);

#pragma omp for
        for(ptrdiff_t ia = 0; ia < static_cast<ptrdiff_t>(A.nrows); ++ia) {
            ptrdiff_t C_cols = 0;
            for(ptrdiff_t ja = A_ptr[ia], ea = A_ptr[ia+1]; ja < ea; ++ja) {
                ptrdiff_t ca = A_col[ja];

                for(ptrdiff_t jb = B_ptr[ca], eb = B_ptr[ca+1]; jb < eb; ++jb) {
                    ptrdiff_t cb = B_col[jb];
                    if (marker[cb] != ia) {
                        marker[cb]  = ia;
                        ++C_cols;
                    }
                }
            }
            C_ptr[ia + 1] = C_cols;
        }
    }

    C->set_nonzeros(C->scan_row_sizes(), /*need_values = */false);
    auto C_col = C->col;

#pragma omp parallel
    {
        std::vector<ptrdiff_t> marker(B.ncols, -1);

#pragma omp for
        for(ptrdiff_t ia = 0; ia < static_cast<ptrdiff_t>(A.nrows); ++ia) {
            ptrdiff_t row_beg = C_ptr[ia];
            ptrdiff_t row_end = row_beg;

            for(ptrdiff_t ja = A_ptr[ia], ea = A_ptr[ia+1]; ja < ea; ++ja) {
                ptrdiff_t ca = A_col[ja];

                for(ptrdiff_t jb = B_ptr[ca], eb = B_ptr[ca+1]; jb < eb; ++jb) {
                    ptrdiff_t cb = B_col[jb];

                    if (marker[cb] < row_beg) {
                        marker[cb] = row_end;
                        C_col[row_end] = cb;
                        ++row_end;
                    }
                }
            }

            std::sort(C_col + row_beg, C_col + row_end);
        }
    }

    return C;
}

} // namespace detail

/// ILU(k) smoother.
template <class Backend>
struct ilup {
    typedef typename Backend::value_type      value_type;

    typedef ilu0<Backend> Base;

    /// Relaxation parameters.
    struct params : Base::params {
        typedef typename Base::params BasePrm;

        /// Level of fill-in.
        int k;

        params() : k(1) {}

#ifndef AMGCL_NO_BOOST
        params(const boost::property_tree::ptree &p)
            : BasePrm(p), AMGCL_PARAMS_IMPORT_VALUE(p, k)
        {
            check_params(p, {"k", "damping", "solve"});
        }

        void get(boost::property_tree::ptree &p, const std::string &path) const {
            BasePrm::get(p, path);
            AMGCL_PARAMS_EXPORT_VALUE(p, path, k);
        }
#endif
    } prm;

    /// \copydoc amgcl::relaxation::damped_jacobi::damped_jacobi
    template <class Matrix>
    ilup( const Matrix &A, const params &prm, const typename Backend::params &bprm)
      : prm(prm)
    {
        if (prm.k == 0) {
            base = std::make_shared<Base>(A, prm, bprm);
        } else {
            auto P = detail::symb_product(A, A);
            for(int k = 1; k < prm.k; ++k) {
                P = detail::symb_product(*P, A);
            }

            ptrdiff_t n = backend::rows(A);
            P->val = new value_type[P->nnz];

#pragma omp parallel for
            for(ptrdiff_t i = 0; i < n; ++i) {
                ptrdiff_t p_beg = P->ptr[i];
                ptrdiff_t p_end = P->ptr[i+1];
                ptrdiff_t a_beg = A.ptr[i];
                ptrdiff_t a_end = A.ptr[i+1];

                std::fill(P->val + p_beg, P->val + p_end, math::zero<value_type>());

                for(ptrdiff_t ja = a_beg, ea = a_end, jp = p_beg, ep = p_end; ja < ea; ++ja) {
                    ptrdiff_t ca = A.col[ja];
                    while(jp < ep && P->col[jp] < ca) ++jp;
                    if (P->col[jp] == ca) P->val[jp] = A.val[ja];
                }
            }

            base = std::make_shared<Base>(*P, prm, bprm);
        }
    }

    /// \copydoc amgcl::relaxation::damped_jacobi::apply_pre
    template <class Matrix, class VectorRHS, class VectorX, class VectorTMP>
    void apply_pre(
            const Matrix &A, const VectorRHS &rhs, VectorX &x, VectorTMP &tmp
            ) const
    {
        base->apply_pre(A, rhs, x, tmp);
    }

    /// \copydoc amgcl::relaxation::damped_jacobi::apply_post
    template <class Matrix, class VectorRHS, class VectorX, class VectorTMP>
    void apply_post(
            const Matrix &A, const VectorRHS &rhs, VectorX &x, VectorTMP &tmp
            ) const
    {
        base->apply_post(A, rhs, x, tmp);
    }

    template <class Matrix, class VectorRHS, class VectorX>
    void apply(const Matrix &A, const VectorRHS &rhs, VectorX &x) const
    {
        base->apply(A, rhs, x);
    }

    size_t bytes() const {
        return base->bytes();
    }

    private:
        std::shared_ptr<Base> base;
};

} // namespace relaxation
} // namespace amgcl

#endif
