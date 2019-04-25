#ifndef AMGCL_DETAIL_QR_HPP
#define AMGCL_DETAIL_QR_HPP

/*
The MIT License

Copyright (c) 2012-2019 Denis Demidov <dennis.demidov@gmail.com>

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
 * \file   amgcl/detail/qr.hpp
 * \author Denis Demidov <dennis.demidov@gmail.com>
 * \brief  QR decomposition.
 *
 * This is a port of ZGEQR2 procedure from LAPACK and its dependencies.
 * The original code included the following copyright notice:
 * \verbatim
   Copyright (c) 1992-2013 The University of Tennessee and The University
                           of Tennessee Research Foundation.  All rights
                           reserved.
   Copyright (c) 2000-2013 The University of California Berkeley. All
                           rights reserved.
   Copyright (c) 2006-2013 The University of Colorado Denver.  All rights
                           reserved.

   $COPYRIGHT$

   Additional copyrights may follow

   $HEADER$

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions are
   met:

   - Redistributions of source code must retain the above copyright
     notice, this list of conditions and the following disclaimer.

   - Redistributions in binary form must reproduce the above copyright
     notice, this list of conditions and the following disclaimer listed
     in this license in the documentation and/or other materials
     provided with the distribution.

   - Neither the name of the copyright holders nor the names of its
     contributors may be used to endorse or promote products derived from
     this software without specific prior written permission.

   The copyright holders provide no reassurances that the source code
   provided does not infringe any patent, copyright, or any other
   intellectual property rights of third parties.  The copyright holders
   disclaim any liability to any recipient for claims brought against
   recipient by any third party for infringement of that parties
   intellectual property rights.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
   A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
   OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
   SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
   LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
   DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
   THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
   (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
   OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 * \endverbatim
 */

#include <vector>
#include <complex>
#include <cmath>

#include <amgcl/util.hpp>
#include <amgcl/value_type/interface.hpp>

namespace amgcl {
namespace detail {

enum storage_order {
    row_major,
    col_major
};

template <class T>
inline T real(T a) {
    return a;
}

template <class T>
inline T real(std::complex<T> a) {
    return std::real(a);
}

/// In-place QR factorization.
template <typename value_type, class Enable = void>
class QR {
    public:
        QR() : m(0), n(0), row_stride(0), col_stride(0), r(NULL) {}

        void compute(int rows, int cols, int row_stride, int col_stride, value_type *A) {
            /*
             *  Ported from ZGEQR2
             *  ==================
             *
             *  Computes a QR factorization of an matrix A:
             *  A = Q * R.
             *
             *  Arguments
             *  =========
             *
             *  rows    The number of rows of the matrix A.
             *  cols    The number of columns of the matrix A.
             *
             *  A       On entry, the rows by cols matrix A.
             *          On exit, the elements on and above the diagonal of the
             *          array contain the min(m,n) by n upper trapezoidal
             *          matrix R (R is upper triangular if m >= n); the
             *          elements below the diagonal, with the array TAU,
             *          represent the unitary matrix Q as a product of
             *          elementary reflectors (see Further Details).
             *
             *  Further Details
             *  ===============
             *
             *  The matrix Q is represented as a product of elementary reflectors
             *
             *     Q = H(1) H(2) . . . H(k), where k = min(m,n).
             *
             *  Each H(i) has the form
             *
             *     H(i) = I - tau * v * v'
             *
             *  where tau is a value_type scalar, and v is a value_type vector
             *  with v[0:i) = 0 and v[i] = 1; v[i:m) is stored on exit in
             *  A[i+1:m)[i], and tau in tau[i].
             *  ==============================================================
             */
            const int m = rows;
            const int n = cols;
            const int k = std::min(m, n);

            if (k <= 0) return;

            r = A;

            tau.resize(k);

            for(int i = 0, ii = 0; i < k; ++i, ii += row_stride + col_stride) {
                // Generate elementary reflector H(i) to annihilate A[i+1:m)[i]
                tau[i] = gen_reflector(m-i, A[ii], A + ii + row_stride, row_stride);

                if (i+1 < n) {
                    // Apply H(i)' to A[i:m)[i+1:n) from the left
                    apply_reflector(m-i, n-i-1, A + ii, row_stride, math::adjoint(tau[i]),
                            A + ii + col_stride, row_stride, col_stride);
                }
            }
        }

        void compute(int rows, int cols, value_type *A, storage_order order = row_major) {
            int row_stride = (order == row_major ? cols : 1);
            int col_stride = (order == row_major ? 1 : rows);
            compute(rows, cols, row_stride, col_stride, A);
        }

        // Computes Q explicitly.
        void factorize(int rows, int cols, int row_stride, int col_stride, value_type *A) {
            /*
             *  Ported from ZUNG2R
             *  ==================
             *
             *  Generates an m by n matrix Q with orthonormal columns, which is
             *  defined as the first n columns of a product of k elementary
             *  reflectors of order m
             *
             *        Q  =  H(1) H(2) . . . H(k)
             *
             *  as returned by compute() [ZGEQR2].
             *
             *  ==============================================================
             */
            compute(rows, cols, row_stride, col_stride, A);

            m = rows;
            n = cols;

            int k = std::min(m, n);

            this->row_stride = row_stride;
            this->col_stride = col_stride;

            q.resize(m * n);

            // Initialise columns k+1:n to zero.
            // [In the original code these were initialized to the columns of
            // the unit matrix, but since k = min(n,m), the main diagonal is
            // never seen here].
            for(int i = 0, ia = 0; i < m; ++i, ia += row_stride)
                for(int j = k, ja = k * col_stride; j < n; ++j, ja += col_stride)
                    q[ia + ja] = (i == j ? math::identity<value_type>() : math::zero<value_type>());

            for(int i = k-1, ic = i * col_stride, ii = i*(row_stride + col_stride);
                    i >= 0; --i, ic -= col_stride, ii -= row_stride + col_stride)
            {
                // Apply H(i) to A[i:m)[i+1:n) from the left
                if (i < n-1)
                    apply_reflector(m-i, n-i-1, r+ii, row_stride, tau[i], &q[ii+col_stride], row_stride, col_stride);

                // Copy i-th reflector (including zeros and unit diagonal)
                // to the column of Q to be processed next
                for(int j = 0, jr = 0; j < i; ++j, jr += row_stride)
                    q[jr+ic] = math::zero<value_type>();

                q[ii] = math::identity<value_type>() - tau[i];

                for(int j = i + 1, jr=j*row_stride; j < m; ++j, jr += row_stride)
                    q[jr + ic] = -tau[i] * r[jr + ic];
            }
        }

        void factorize(int rows, int cols, value_type *A, storage_order order = row_major) {
            int row_stride = (order == row_major ? cols : 1);
            int col_stride = (order == row_major ? 1 : rows);
            factorize(rows, cols, row_stride, col_stride, A);
        }

        // Returns element of the matrix R.
        value_type R(int i, int j) const {
            if (j < i) return math::zero<value_type>();
            return r[i*row_stride + j*col_stride];
        }

        // Returns element of the matrix Q.
        value_type Q(int i, int j) const {
            return q[i*row_stride + j*col_stride];
        }

        // Solves the system Q R x = f
        void solve(
                int rows, int cols, int row_stride, int col_stride, value_type *A,
                const value_type *b, value_type *x, bool computed = false)
        {
            f.resize(rows);
            std::copy(b, b + rows, f.begin());

            if (rows >= cols) {
                // We are solving overdetermined (tall) system Ax = f by
                // writing the matrix A as A = QR and solving for x as
                // x = R^-1 Q^-1 f = R^-1 Q^T f.
                if (!computed) compute(rows, cols, row_stride, col_stride, A);

                for(int i = 0, ii = 0; i < cols; ++i, ii += row_stride + col_stride)
                    apply_reflector(rows-i, 1, r+ii, row_stride, math::adjoint(tau[i]), &f[i], 1, 1);

                std::copy(f.begin(), f.begin()+cols, x);

                for(int i = cols, ia = (cols-1) * col_stride; i --> 0; ia -= col_stride) {
                    value_type rii = r[i*(row_stride+col_stride)];
                    if (math::is_zero(rii)) continue;
                    x[i] = math::inverse(rii) * x[i];

                    for(int j = 0, ja = 0; j < i; ++j, ja += row_stride)
                        x[j] -= r[ia + ja] * x[i];
                }
            } else {
                // We are solving underdetermined (wide) system Ax = f by
                // writing the matrix A^T as A^T = QR and solving for x as
                // x = Q^-T R^-T f = Q R^-T f.
                if (!computed) {
                    for(int i = 0, n = cols * rows; i < n; ++i)
                        A[i] = math::adjoint(A[i]);
                    compute(cols, rows, col_stride, row_stride, A);
                }

                for(int i = 0, ia = 0; i < rows; ++i, ia += col_stride) {
                    value_type rii = math::adjoint(r[i*(row_stride+col_stride)]);
                    if (math::is_zero(rii)) continue;
                    f[i] = math::inverse(rii) * f[i];

                    for(int j = i+1, ja = j * row_stride; j < rows; ++j, ja += row_stride)
                        f[j] -= math::adjoint(r[ia + ja]) * f[i];
                }

                std::copy(f.begin(), f.end(), x);
                std::fill(x+rows, x+cols, math::zero<value_type>());

                for(int i = rows; i --> 0; ) {
                    int ii = i * (col_stride + row_stride);
                    apply_reflector(cols-i, 1, r+ii, col_stride, tau[i], x+i, 1, 1);
                }
            }
        }

        void solve(
                int rows, int cols, value_type *A, const value_type *b, value_type *x,
                storage_order order = row_major, bool computed = false
                )
        {
            int row_stride = (order == row_major ? cols : 1);
            int col_stride = (order == row_major ? 1 : rows);
            solve(rows, cols, row_stride, col_stride, A, b, x, computed);
        }

        size_t bytes() {
            return sizeof(value_type) * (tau.size() + f.size() + q.size());
        }

    private:
        typedef typename math::scalar_of<value_type>::type scalar_type;

        static scalar_type sqr(scalar_type x) { return x * x; }

        int m, n, row_stride, col_stride;

        value_type *r;
        std::vector<value_type> tau, f;
        std::vector<value_type> q;

        static value_type gen_reflector(int order, value_type &alpha, value_type *x, int stride) {
            /*
             *  Ported from ZLARFG
             *  ==================
             *
             *  Generates a value_type elementary reflector H of order n, such
             *  that
             *
             *        H' * ( alpha ) = ( beta ),   H' * H = I.
             *             (   x   )   (   0  )
             *
             *  where alpha and beta are scalars, with beta real, and x is an
             *  (n-1)-element value_type vector. H is represented in the form
             *
             *        H = I - tau * ( 1 ) * ( 1 v' ) ,
             *                      ( v )
             *
             *  where tau is a value_type scalar and v is a value_type
             *  (n-1)-element vector. Note that H is not hermitian.
             *
             *  If the elements of x are all zero and alpha is real,
             *  then tau = 0 and H is taken to be the unit matrix.
             *
             *  Otherwise  1 <= real(tau) <= 2  and  abs(tau-1) <= 1 .
             *
             *  Arguments
             *  =========
             *
             *  order   The order of the elementary reflector.
             *
             *  alpha   On entry, the value alpha.
             *          On exit, it is overwritten with the value beta.
             *
             *  x       dimension (1+(order-2)*abs(stride))
             *          On entry, the vector x.
             *          On exit, it is overwritten with the vector v.
             *
             *  stride  The increment between elements of x.
             *
             *  Returns the value tau.
             *
             *  ==============================================================
             */
            value_type tau = math::zero<value_type>();
            if (order <= 1) return tau;
            int n = order - 1;

            scalar_type xnorm2 = 0;
            for(int i = 0, ii = 0; i < n; ++i, ii += stride)
                xnorm2 += sqr(math::norm(x[ii]));

            if (math::is_zero(xnorm2)) return tau;

            scalar_type beta = -std::abs(sqrt(sqr(math::norm(alpha)) + xnorm2));
            if (amgcl::detail::real(alpha) < 0) beta = -beta;

            tau = math::identity<value_type>() - math::inverse(beta) * alpha;
            alpha = math::inverse(alpha - beta * math::identity<value_type>());

            for(int i = 0, ii = 0; i < n; ++i, ii += stride)
                x[ii] = alpha * x[ii];

            alpha = beta * math::identity<value_type>();
            return tau;
        }

        static void apply_reflector(
                int m, int n, const value_type *v, int v_stride, value_type tau,
                value_type *C, int row_stride, int col_stride
                )
        {
            /*
             *  Ported from ZLARF
             *  =================
             *
             *  Applies an elementary reflector H to an m-by-n matrix C from
             *  the left. H is represented in the form
             *
             *        H = I - v * tau * v'
             *
             *  where tau is a value_type scalar and v is a value_type vector.
             *
             *  If tau = 0, then H is taken to be the unit matrix.
             *
             *  To apply H' (the conjugate transpose of H), supply adjoint(tau)
             *  instead of tau.
             *
             *  Arguments
             *  =========
             *
             *  m          The number of rows of the matrix C.
             *
             *  n          The number of columns of the matrix C.
             *
             *  v          The vector v in the representation of H.
             *             v is not used if tau = 0.
             *             The value of v[0] is ignored and assumed to be 1.
             *
             *  v_stride   The increment between elements of v.
             *
             *  tau        The value tau in the representation of H.
             *
             *  C          On entry, the m-by-n matrix C.
             *             On exit, C is overwritten by the matrix H * C.
             *
             *  row_stride The increment between the rows of C.
             *  col_stride The increment between the columns of C.
             *
             *  ==============================================================
             */

            if (math::is_zero(tau)) return;

            // w = C` * v; C -= tau * v * w`
            for(int i = 0, ia=0; i < n; ++i, ia += col_stride) {
                value_type s = math::adjoint(C[ia]);
                for(int j = 1, jv = v_stride, ja=row_stride; j < m; ++j, jv += v_stride, ja += row_stride) {
                    s += math::adjoint(C[ja+ia]) * v[jv];
                }

                s = tau * math::adjoint(s);
                C[ia] -= s;
                for(int j = 1, jv = v_stride, ja=row_stride; j < m; ++j, jv += v_stride, ja += row_stride) {
                    C[ja+ia] -= v[jv] * s;
                }
            }
        }

};

template <class value_type>
class QR<value_type, typename std::enable_if<math::is_static_matrix<value_type>::value>::type>
{
    public:
        typedef typename amgcl::math::rhs_of<value_type>::type rhs_type;

        QR() {}

        void compute(int rows, int cols, int row_stride, int col_stride, value_type *A) {
            const int M = math::static_rows<value_type>::value;
            const int N = math::static_cols<value_type>::value;

            m = rows;
            n = cols;

            r = A;

            copy_to_scalar_buf(rows, cols, row_stride, col_stride, A);
            base.compute(rows * M, cols * N, 1, rows * M, buf.data());
        }

        void factorize(int rows, int cols, int row_stride, int col_stride, value_type *A) {
            const int M = math::static_rows<value_type>::value;
            const int N = math::static_cols<value_type>::value;

            m = rows * M;
            n = cols * N;

            r = A;

            copy_to_scalar_buf(rows, cols, row_stride, col_stride, A);
            base.factorize(m, n, 1, m, buf.data());
        }

        void factorize(int rows, int cols, value_type *A, storage_order order = row_major) {
            int row_stride = (order == row_major ? cols : 1);
            int col_stride = (order == row_major ? 1 : rows);
            factorize(rows, cols, row_stride, col_stride, A);
        }

        value_type R(int i, int j) const {
            const int N = math::static_rows<value_type>::value;
            const int M = math::static_cols<value_type>::value;

            value_type v;

            if (j < i) {
                v = math::zero<value_type>();
            } else {
                for(int ii = 0; ii < N; ++ii)
                    for(int jj = 0; jj < M; ++jj)
                        v(ii,jj) = base.R(i * N + ii, j * M + jj);
            }

            return v;
        }

        // Returns element of the matrix Q.
        value_type Q(int i, int j) const {
            const int N = math::static_rows<value_type>::value;
            const int M = math::static_cols<value_type>::value;

            value_type v;

            for(int ii = 0; ii < N; ++ii)
                for(int jj = 0; jj < M; ++jj)
                    v(ii,jj) = base.Q(i * N + ii, j * M + jj);

            return v;
        }

        // Solves the system Q R x = f
        void solve(
                int rows, int cols, int row_stride, int col_stride, value_type *A,
                const rhs_type *f, rhs_type *x, bool computed = false)
        {
            const int M = math::static_rows<value_type>::value;
            const int N = math::static_cols<value_type>::value;

            m = rows * M;
            n = cols * N;

            r = A;

            copy_to_scalar_buf(rows, cols, row_stride, col_stride, A);
            base.solve(m, n, 1, m, buf.data(),
                    reinterpret_cast<const scalar_type*>(f),
                    reinterpret_cast<scalar_type*>(x),
                    computed
                    );
        }

        void solve(
                int rows, int cols, value_type *A, const rhs_type *f, rhs_type *x,
                storage_order order = row_major, bool computed = false
                )
        {
            int row_stride = (order == row_major ? cols : 1);
            int col_stride = (order == row_major ? 1 : rows);
            solve(rows, cols, row_stride, col_stride, A, f, x, computed);
        }

        size_t bytes() const {
            return base.bytes() + sizeof(scalar_type) * buf.size();
        }

    private:
        typedef typename amgcl::math::scalar_of<value_type>::type scalar_type;

        int m, n;
        value_type *r;

        QR<scalar_type> base;
        std::vector<scalar_type> buf;

        void copy_to_scalar_buf(int rows, int cols, int row_stride, int col_stride, value_type *A) {
            const int M = math::static_rows<value_type>::value;
            const int N = math::static_cols<value_type>::value;

            buf.resize(M * rows * N * cols);

            const int scalar_rows = M * rows;

            for(int i = 0, ib = 0; i < rows; ++i)
                for(int ii = 0; ii < M; ++ii, ++ib)
                    for(int j = 0, jb = 0; j < cols; ++j)
                        for(int jj = 0; jj < N; ++jj, jb += scalar_rows)
                            buf[ib + jb] = A[i * row_stride + j * col_stride](ii, jj);
        }
};

} // namespace detail
} // namespace amgcl

#endif
