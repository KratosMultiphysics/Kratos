#ifndef AMGCL_DETAIL_QR_HPP
#define AMGCL_DETAIL_QR_HPP

///TAKEN FROM AMGCL

/*
The MIT License

Copyright (c) 2012-2016 Denis Demidov <dennis.demidov@gmail.com>

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

// #include <vector>
// #include <cmath>
// 
// #include <amgcl/util.hpp>
// #include <amgcl/value_type/interface.hpp>

namespace Kratos {

enum storage_order {
    row_major,
    col_major
};

/// In-place QR factorization.
/**
 * \tparam Order Storage order of the input matrix. Should be col_major for
 *               the best performance.
 */
template <typename value_type, storage_order Order>
class QR {
    public:
        QR() : m(0), n(0) {}

        void compute(int rows, int cols, value_type *A, int max_cols = -1) {
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
            m = rows;
            n = cols;
            k = std::min(m, n);

            nmax = (max_cols < 0 ? n : max_cols);

            r = A;

            tau.resize(std::min(m, nmax));

            if (k <= 0) return;

            const int row_stride = (Order == row_major ? n : 1);
            const int col_stride = (Order == row_major ? 1 : m);

            for(int i = 0, ii = 0; i < k; ++i, ii += row_stride + col_stride) {
                // Generate elementary reflector H(i) to annihilate A[i+1:m)[i]
                tau[i] = gen_reflector(m-i, A[ii], A + ii + row_stride, row_stride);

                if (i+1 < n) {
                    // Apply H(i)' to A[i:m)[i+1:n) from the left
                    apply_reflector(m-i, n-i-1, A + ii, row_stride, tau[i],
                            A + ii + col_stride, row_stride, col_stride);
                }
            }
        }
        
        //~ void append_cols(int cols) {
            //~ const int row_stride = (Order == row_major ? nmax : 1);
            //~ const int col_stride = (Order == row_major ? 1 : m);

            //~ int old_n = n;
            //~ n += cols;

            //~ precondition(n <= nmax, "Too many columns in QR::append_cols()");

            //~ int old_k = k;
            //~ k = std::min(m, n);

            //~ for(int i = 0, ii = 0; i < k; ++i, ii += row_stride + col_stride) {
                //~ if (i >= old_k) {
                    //~ // Generate elementary reflector H(i) to annihilate A[i+1:m)[i]
                    //~ tau[i] = gen_reflector(m-i, r[ii], r + ii + row_stride, row_stride);
                //~ }

                //~ if (i+1 < n) {
                    //~ // Apply H(i)' to A[i:m)[i+1:n) from the left
                    //~ int l = std::max(i, old_n-1);
                    //~ apply_reflector(m-i, n-l-1, r + ii, row_stride, tau[i],
                            //~ r + i * row_stride + (l + 1) * col_stride, row_stride, col_stride);
                //~ }
            //~ }
        //~ }

        // Returns element of the matrix R.
        value_type R(int i, int j) const {
            if (j < i) return 0.0;

            const int row_stride = (Order == row_major ? nmax : 1);
            const int col_stride = (Order == row_major ? 1 : m);

            return r[i*row_stride + j*col_stride];
        }

        // Returns the Frobenius norm of matrix R.
        value_type normR() const {
            value_type normR = 0;

            const int row_stride = (Order == row_major ? nmax : 1);
            const int col_stride = (Order == row_major ? 1 : m);

            for (int i = 0; i < n; ++i){
                for (int j = 0; j < n; ++j){
                    normR += std::pow(r[i*row_stride + j*col_stride], 2);
                }
            }

            return std::sqrt(normR);
        }

        // Returns element of the matrix Q.
        value_type Q(int i, int j) const {
            const int row_stride = (Order == row_major ? nmax : 1);
            const int col_stride = (Order == row_major ? 1 : m);

            return q[i*row_stride + j*col_stride];
        }

        // Solves the system Q R x = f
        void solve(value_type *f, value_type *x) const {
            const int row_stride = (Order == row_major ? nmax : 1);
            const int col_stride = (Order == row_major ? 1 : m);

            for(int i = 0, ii = 0; i < n; ++i, ii += row_stride + col_stride)
                apply_reflector(m-i, 1, r+ii, row_stride, tau[i], f+i, 1, 1);

            std::copy(f, f+n, x);

            for(int i = n; i --> 0; ) {
                value_type rii = r[i*(row_stride+col_stride)];
                if (rii==0.0) continue;
                x[i] = x[i]/rii;

                for(int j = 0, ja = 0; j < i; ++j, ja += row_stride)
                    x[j] -= r[ja+i*col_stride] * x[i];
            }
        }

        // Computes Q explicitly.
        void compute_q(int ncols = -1) {
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
            q.resize(m * nmax);

            ncols = (ncols < 0 ? n : ncols);

            const int row_stride = (Order == row_major ? nmax : 1);
            const int col_stride = (Order == row_major ? 1 : m);

            // Initialise columns k+1:n to zero.
            // [In the original code these were initialized to the columns of
            // the unit matrix, but since k = min(n,m), the main diagonal is
            // never seen here].
            for(int i = 0, ia = 0; i < m; ++i, ia += row_stride)
                for(int j = k, ja = k * col_stride; j < ncols; ++j, ja += col_stride)
                        q[ia + ja] = (i == j ? 1.0 : 0.0);

            for(int i = k-1, ic = i * col_stride, ii = i*(row_stride + col_stride);
                    i >= 0; --i, ic -= col_stride, ii -= row_stride + col_stride)
            {
                // Apply H(i) to A[i:m)[i+1:n) from the left
                if (i < ncols-1)
                    apply_reflector(m-i, ncols-i-1, r+ii, row_stride, tau[i], &q[ii+col_stride], row_stride, col_stride);

                // Copy i-th reflector (including zeros and unit diagonal)
                // to the column of Q to be processed next
                for(int j = 0, jr = 0; j < i; ++j, jr += row_stride)
                    q[jr+ic] = 0.0;

                q[ii] = 1.0 - tau[i];

                for(int j = i + 1, jr=j*row_stride; j < m; ++j, jr += row_stride)
                    q[jr + ic] = -tau[i] * r[jr + ic];
            }
        }
    private:
        typedef value_type scalar_type;

        static scalar_type sqr(scalar_type x) { return x * x; }

        int m, n, k, nmax;

        value_type *r;
        std::vector<value_type> tau;
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
            value_type tau = 0.0;
            if (order <= 1) return tau;
            int n = order - 1;

            scalar_type xnorm2 = 0;
            for(int i = 0, ii = 0; i < n; ++i, ii += stride)
                xnorm2 += sqr(x[ii]);

            if (xnorm2 == 0.0) return tau;

            scalar_type beta = sqrt(sqr(alpha) + xnorm2);

            tau = 1.0 - alpha/beta;
            alpha = 1.0/(alpha - beta * 1.0);

            for(int i = 0, ii = 0; i < n; ++i, ii += stride)
                x[ii] = alpha * x[ii];

            alpha = beta * 1.0;
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

            if (tau==0.0) return;

            // w = C` * v; C -= tau * v * w`
            for(int i = 0, ia=0; i < n; ++i, ia += col_stride) {
                value_type s = C[ia];
                for(int j = 1, jv = v_stride, ja=row_stride; j < m; ++j, jv += v_stride, ja += row_stride) {
                    s += C[ja+ia] * v[jv];
                }

                s = tau * s;
                C[ia] -= s;
                for(int j = 1, jv = v_stride, ja=row_stride; j < m; ++j, jv += v_stride, ja += row_stride) {
                    C[ja+ia] -= v[jv] * s;
                }
            }
        }

};

// template <class value_type, storage_order Order>
// class QR<value_type, Order >
// {
//     public:
//         typedef  value_type rhs_type;
// 
//         QR() {}
// 
//         void compute(int rows, int cols, value_type *A) {
//             const int N = rows; //math::static_rows<value_type>::value;
//             const int M = cols; //math::static_cols<value_type>::value;
// 
//             buf.resize(rows * cols * N * M);
// 
//             const int brows = M * rows;
//             const int row_stride = (Order == row_major ? cols : 1);
//             const int col_stride = (Order == row_major ? 1 : rows);
// 
//             for(int i = 0, ib = 0; i < rows; ++i)
//                 for(int ii = 0; ii < N; ++ii, ++ib)
//                     for(int j = 0, jb = 0; j < cols; ++j)
//                         for(int jj = 0; jj < M; ++jj, jb += brows)
//                             buf[ib + jb] = A[i * row_stride + j * col_stride](ii, jj);
// 
//             base.compute(rows * N, cols * M, buf.data());
//         }
// 
//         value_type R(int i, int j) const {
//             const int N = rows; //math::static_rows<value_type>::value;
//             const int M = cols; //math::static_cols<value_type>::value;
// 
//             value_type v;
// 
//             if (j < i) {
//                 v = 0.0;
//             } else {
//                 for(int ii = 0; ii < N; ++ii)
//                     for(int jj = 0; jj < M; ++jj)
//                         v(ii,jj) = base.R(i * N + ii, j * M + jj);
//             }
// 
//             return v;
//         }
// 
//         // Returns element of the matrix Q.
//         value_type Q(int i, int j) const {
//             const int N = rows; //math::static_rows<value_type>::value;
//             const int M = cols; //math::static_cols<value_type>::value;
// 
//             value_type v;
// 
//             for(int ii = 0; ii < N; ++ii)
//                 for(int jj = 0; jj < M; ++jj)
//                     v(ii,jj) = base.Q(i * N + ii, j * M + jj);
// 
//             return v;
//         }
// 
//         // Solves the system Q R x = f
//         void solve(rhs_type *f, rhs_type *x) const {
//             base.solve(
//                     reinterpret_cast<scalar_type*>(f),
//                     reinterpret_cast<scalar_type*>(x)
//                     );
//         }
// 
//         void compute_q() { base.compute_q(); }
// 
//     private:
//         typedef typename amgcl::math::scalar_of<value_type>::type scalar_type;
// 
//         QR<scalar_type, col_major> base;
//         std::vector<scalar_type> buf;
// };

} // namespace Kratos

#endif
