#ifndef VIENNACL_LINALG_QR_METHOD_HPP_
#define VIENNACL_LINALG_QR_METHOD_HPP_

/* =========================================================================
   Copyright (c) 2010-2014, Institute for Microelectronics,
                            Institute for Analysis and Scientific Computing,
                            TU Wien.
   Portions of this software are copyright by UChicago Argonne, LLC.

                            -----------------
                  ViennaCL - The Vienna Computing Library
                            -----------------

   Project Head:    Karl Rupp                   rupp@iue.tuwien.ac.at

   (A list of authors and contributors can be found in the PDF manual)

   License:         MIT (X11), see file LICENSE in the base directory
============================================================================= */

#include "viennacl/vector.hpp"
#include "viennacl/matrix.hpp"

#include "viennacl/linalg/qr-method-common.hpp"
#include "viennacl/linalg/prod.hpp"

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>

/** @file viennacl/linalg/qr-method.hpp
    @brief Implementation of the QR method for eigenvalue computations. Experimental.
*/

namespace viennacl
{
  namespace linalg
  {
    namespace detail
    {
        template<typename MatrixType, typename VectorType>
        void givens_next(MatrixType& matrix,
                        VectorType& tmp1,
                        VectorType& tmp2,
                        int l,
                        int m
                      )
        {
          viennacl::ocl::context & ctx = const_cast<viennacl::ocl::context &>(viennacl::traits::opencl_handle(matrix).context());

          typedef typename MatrixType::value_type                                   ScalarType;
          typedef typename viennacl::result_of::cpu_value_type<ScalarType>::type    CPU_ScalarType;

          viennacl::ocl::kernel& kernel = ctx.get_kernel(viennacl::linalg::opencl::kernels::svd<CPU_ScalarType>::program_name(), SVD_GIVENS_NEXT_KERNEL);

          kernel.global_work_size(0, viennacl::tools::align_to_multiple<cl_uint>(cl_uint(viennacl::traits::size1(matrix)), 256));
          kernel.local_work_size(0, 256);

          viennacl::ocl::enqueue(kernel(
                                        matrix,
                                        tmp1,
                                        tmp2,
                                        static_cast<cl_uint>(matrix.size1()),
                                        static_cast<cl_uint>(matrix.internal_size2()),
                                        static_cast<cl_uint>(l),
                                        static_cast<cl_uint>(m - 1)
                                ));
        }


        // Symmetric tridiagonal QL algorithm.
        // This is derived from the Algol procedures tql2, by Bowdler, Martin, Reinsch, and Wilkinson,
        // Handbook for Auto. Comp., Vol.ii-Linear Algebra, and the corresponding Fortran subroutine in EISPACK.
        template <typename SCALARTYPE, unsigned int ALIGNMENT>
        void tql2(viennacl::matrix<SCALARTYPE, row_major, ALIGNMENT> & Q,
                  boost::numeric::ublas::vector<SCALARTYPE> & d,
                  boost::numeric::ublas::vector<SCALARTYPE> & e)
        {
            int n = static_cast<int>(Q.size1());

            boost::numeric::ublas::vector<SCALARTYPE> cs(n), ss(n);
            viennacl::vector<SCALARTYPE> tmp1(n), tmp2(n);

            for (int i = 1; i < n; i++)
                e(i - 1) = e(i);

            e(n - 1) = 0;

            SCALARTYPE f = 0;
            SCALARTYPE tst1 = 0;
            SCALARTYPE eps = 2 * static_cast<SCALARTYPE>(EPS);

            for (int l = 0; l < n; l++)
            {
                // Find small subdiagonal element.
                tst1 = std::max<SCALARTYPE>(tst1, std::fabs(d(l)) + std::fabs(e(l)));
                int m = l;
                while (m < n)
                {
                    if (std::fabs(e(m)) <= eps * tst1)
                        break;
                    m++;
                }

                // If m == l, d(l) is an eigenvalue, otherwise, iterate.
                if (m > l)
                {
                    int iter = 0;
                    do
                    {
                        iter = iter + 1;  // (Could check iteration count here.)

                        // Compute implicit shift
                        SCALARTYPE g = d(l);
                        SCALARTYPE p = (d(l + 1) - g) / (2 * e(l));
                        SCALARTYPE r = pythag<SCALARTYPE>(p, 1);
                        if (p < 0)
                        {
                            r = -r;
                        }

                        d(l) = e(l) / (p + r);
                        d(l + 1) = e(l) * (p + r);
                        SCALARTYPE dl1 = d(l + 1);
                        SCALARTYPE h = g - d(l);
                        for (int i = l + 2; i < n; i++)
                        {
                            d(i) -= h;
                        }

                        f = f + h;

                        // Implicit QL transformation.
                        p = d(m);
                        SCALARTYPE c = 1;
                        SCALARTYPE c2 = c;
                        SCALARTYPE c3 = c;
                        SCALARTYPE el1 = e(l + 1);
                        SCALARTYPE s = 0;
                        SCALARTYPE s2 = 0;
                        for (int i = m - 1; i >= l; i--)
                        {
                            c3 = c2;
                            c2 = c;
                            s2 = s;
                            g = c * e(i);
                            h = c * p;
                            r = pythag(p, e(i));
                            e(i + 1) = s * r;
                            s = e(i) / r;
                            c = p / r;
                            p = c * d(i) - s * g;
                            d(i + 1) = h + s * (c * g + s * d(i));

                            cs[i] = c;
                            ss[i] = s;
                        }

                        p = -s * s2 * c3 * el1 * e(l) / dl1;
                        e(l) = s * p;
                        d(l) = c * p;

                        {
                            viennacl::copy(cs, tmp1);
                            viennacl::copy(ss, tmp2);

                            givens_next(Q, tmp1, tmp2, l, m);
                        }

                        // Check for convergence.
                    }
                    while (std::fabs(e(l)) > eps * tst1);
                }
                d(l) = d(l) + f;
                e(l) = 0;
            }
        }

        template <typename SCALARTYPE, typename MatrixT>
        void final_iter_update_gpu(MatrixT& A,
                                int n,
                                int last_n,
                                SCALARTYPE q,
                                SCALARTYPE p
                                )
        {
            viennacl::ocl::context & ctx = const_cast<viennacl::ocl::context &>(viennacl::traits::opencl_handle(A).context());

            viennacl::ocl::kernel& kernel = ctx.get_kernel(viennacl::linalg::opencl::kernels::svd<SCALARTYPE>::program_name(), SVD_FINAL_ITER_UPDATE_KERNEL);

            viennacl::ocl::enqueue(kernel(
                                          A,
                                          static_cast<cl_uint>(A.internal_size1()),
                                          static_cast<cl_uint>(n),
                                          static_cast<cl_uint>(last_n),
                                          q,
                                          p
                                  ));
        }

        template <typename SCALARTYPE, typename MatrixT>
        void update_float_QR_column_gpu(MatrixT& A,
                                const std::vector<SCALARTYPE>& buf,
                                viennacl::vector<SCALARTYPE>& buf_vcl,
                                int m,
                                int n,
                                int last_n,
                                bool //is_triangular
                                )
        {
            viennacl::ocl::context & ctx = const_cast<viennacl::ocl::context &>(viennacl::traits::opencl_handle(A).context());

            viennacl::fast_copy(buf, buf_vcl);

            viennacl::ocl::kernel& kernel = ctx.get_kernel(viennacl::linalg::opencl::kernels::svd<SCALARTYPE>::program_name(), SVD_UPDATE_QR_COLUMN_KERNEL);

            viennacl::ocl::enqueue(kernel(
                                          A,
                                          static_cast<cl_uint>(A.internal_size1()),
                                          buf_vcl,
                                          static_cast<cl_uint>(m),
                                          static_cast<cl_uint>(n),
                                          static_cast<cl_uint>(last_n)
                                  ));
        }

        template <typename SCALARTYPE, typename MatrixT>
        void final_iter_update(MatrixT& A,
                                int n,
                                int last_n,
                                SCALARTYPE q,
                                SCALARTYPE p
                                )
        {
            for (int i = 0; i < last_n; i++)
            {
                SCALARTYPE v_in = A(i, n);
                SCALARTYPE z = A(i, n - 1);
                A(i, n - 1) = q * z + p * v_in;
                A(i, n) = q * v_in - p * z;
            }
        }

        template <typename SCALARTYPE, typename MatrixT>
        void update_float_QR_column(MatrixT& A,
                                const std::vector<SCALARTYPE>& buf,
                                int m,
                                int n,
                                int last_i,
                                bool is_triangular
                                )
        {
            for (int i = 0; i < last_i; i++)
            {
                int start_k = is_triangular?std::max(i + 1, m):m;

                SCALARTYPE* a_row = A.row(i);

                SCALARTYPE a_ik   = a_row[start_k];
                SCALARTYPE a_ik_1 = 0;
                SCALARTYPE a_ik_2 = 0;

                if(start_k < n)
                    a_ik_1 = a_row[start_k + 1];

                for(int k = start_k; k < n; k++)
                {
                    bool notlast = (k != n - 1);

                    SCALARTYPE p = buf[5 * k] * a_ik + buf[5 * k + 1] * a_ik_1;

                    if (notlast)
                    {
                        a_ik_2 = a_row[k + 2];
                        p = p + buf[5 * k + 2] * a_ik_2;
                        a_ik_2 = a_ik_2 - p * buf[5 * k + 4];
                    }

                    a_row[k] = a_ik - p;
                    a_ik_1 = a_ik_1 - p * buf[5 * k + 3];

                    a_ik = a_ik_1;
                    a_ik_1 = a_ik_2;
                }

                if(start_k < n)
                    a_row[n] = a_ik;
            }
        }

        /** @brief Internal helper class representing a row-major dense matrix used for the QR method for the purpose of computing eigenvalues. */
        template <typename SCALARTYPE>
        class FastMatrix
        {
        public:
            FastMatrix()
            {
                size_ = 0;
            }

            FastMatrix(vcl_size_t sz, vcl_size_t internal_size) : size_(sz), internal_size_(internal_size)
            {
                data.resize(internal_size * internal_size);
            }

            SCALARTYPE& operator()(int i, int j)
            {
                return data[i * internal_size_ + j];
            }

            SCALARTYPE* row(int i)
            {
                return &data[i * internal_size_];
            }

            SCALARTYPE* begin()
            {
                return &data[0];
            }

            SCALARTYPE* end()
            {
                return &data[0] + data.size();
            }

            std::vector<SCALARTYPE> data;
        private:
            vcl_size_t size_;
            vcl_size_t internal_size_;
        };

        // Nonsymmetric reduction from Hessenberg to real Schur form.
        // This is derived from the Algol procedure hqr2, by Martin and Wilkinson, Handbook for Auto. Comp.,
        // Vol.ii-Linear Algebra, and the corresponding  Fortran subroutine in EISPACK.
        template <typename SCALARTYPE, unsigned int ALIGNMENT>
        void hqr2(viennacl::matrix<SCALARTYPE, row_major, ALIGNMENT>& vcl_H,
                    viennacl::matrix<SCALARTYPE, row_major, ALIGNMENT>& V,
                    boost::numeric::ublas::vector<SCALARTYPE>& d,
                    boost::numeric::ublas::vector<SCALARTYPE>& e)
        {
            transpose(V);

            int nn = static_cast<int>(vcl_H.size1());

            FastMatrix<SCALARTYPE> H(nn, vcl_H.internal_size2());//, V(nn);

            std::vector<SCALARTYPE> buf(5 * nn);
            viennacl::vector<SCALARTYPE> buf_vcl(5 * nn);

            viennacl::fast_copy(vcl_H, H.begin());


            int n = nn - 1;

            SCALARTYPE eps = 2 * static_cast<SCALARTYPE>(EPS);
            SCALARTYPE exshift = 0;
            SCALARTYPE p = 0;
            SCALARTYPE q = 0;
            SCALARTYPE r = 0;
            SCALARTYPE s = 0;
            SCALARTYPE z = 0;
            SCALARTYPE t;
            SCALARTYPE w;
            SCALARTYPE x;
            SCALARTYPE y;

            SCALARTYPE out1, out2;

            // compute matrix norm
            SCALARTYPE norm = 0;
            for (int i = 0; i < nn; i++)
            {
                for (int j = std::max(i - 1, 0); j < nn; j++)
                    norm = norm + std::fabs(H(i, j));
            }

            // Outer loop over eigenvalue index
            int iter = 0;
            while (n >= 0)
            {
                // Look for single small sub-diagonal element
                int l = n;
                while (l > 0)
                {
                    s = std::fabs(H(l - 1, l - 1)) + std::fabs(H(l, l));
                    if (s == 0) s = norm;
                    if (std::fabs(H(l, l - 1)) < eps * s)
                        break;

                    l--;
                }

                // Check for convergence
                if (l == n)
                {
                    // One root found
                    H(n, n) = H(n, n) + exshift;
                    d(n) = H(n, n);
                    e(n) = 0;
                    n--;
                    iter = 0;
                }
                else if (l == n - 1)
                {
                    // Two roots found
                    w = H(n, n - 1) * H(n - 1, n);
                    p = (H(n - 1, n - 1) - H(n, n)) / 2;
                    q = p * p + w;
                    z = static_cast<SCALARTYPE>(std::sqrt(std::fabs(q)));
                    H(n, n) = H(n, n) + exshift;
                    H(n - 1, n - 1) = H(n - 1, n - 1) + exshift;
                    x = H(n, n);

                    if (q >= 0)
                    {
                        // Real pair
                        z = (p >= 0) ? (p + z) : (p - z);
                        d(n - 1) = x + z;
                        d(n) = d(n - 1);
                        if (z != 0)
                            d(n) = x - w / z;
                        e(n - 1) = 0;
                        e(n) = 0;
                        x = H(n, n - 1);
                        s = std::fabs(x) + std::fabs(z);
                        p = x / s;
                        q = z / s;
                        r = static_cast<SCALARTYPE>(std::sqrt(p * p + q * q));
                        p = p / r;
                        q = q / r;

                        // Row modification
                        for (int j = n - 1; j < nn; j++)
                        {
                            SCALARTYPE h_nj = H(n, j);
                            z = H(n - 1, j);
                            H(n - 1, j) = q * z + p * h_nj;
                            H(n, j) = q * h_nj - p * z;
                        }

                        final_iter_update(H, n, n + 1, q, p);
                        final_iter_update_gpu(V, n, nn, q, p);
                    }
                    else
                    {
                        // Complex pair
                        d(n - 1) = x + p;
                        d(n) = x + p;
                        e(n - 1) = z;
                        e(n) = -z;
                    }

                    n = n - 2;
                    iter = 0;
                }
                else
                {
                    // No convergence yet

                    // Form shift
                    x = H(n, n);
                    y = 0;
                    w = 0;
                    if (l < n)
                    {
                        y = H(n - 1, n - 1);
                        w = H(n, n - 1) * H(n - 1, n);
                    }

                    // Wilkinson's original ad hoc shift
                    if (iter == 10)
                    {
                        exshift += x;
                        for (int i = 0; i <= n; i++)
                            H(i, i) -= x;

                        s = std::fabs(H(n, n - 1)) + std::fabs(H(n - 1, n - 2));
                        x = y = SCALARTYPE(0.75) * s;
                        w = SCALARTYPE(-0.4375) * s * s;
                    }

                    // MATLAB's new ad hoc shift
                    if (iter == 30)
                    {
                        s = (y - x) / 2;
                        s = s * s + w;
                        if (s > 0)
                        {
                            s = static_cast<SCALARTYPE>(std::sqrt(s));
                            if (y < x) s = -s;
                            s = x - w / ((y - x) / 2 + s);
                            for (int i = 0; i <= n; i++)
                                H(i, i) -= s;
                            exshift += s;
                            x = y = w = SCALARTYPE(0.964);
                        }
                    }

                    iter = iter + 1;

                    // Look for two consecutive small sub-diagonal elements
                    int m = n - 2;
                    while (m >= l)
                    {
                        SCALARTYPE h_m1_m1 = H(m + 1, m + 1);
                        z = H(m, m);
                        r = x - z;
                        s = y - z;
                        p = (r * s - w) / H(m + 1, m) + H(m, m + 1);
                        q = h_m1_m1 - z - r - s;
                        r = H(m + 2, m + 1);
                        s = std::fabs(p) + std::fabs(q) + std::fabs(r);
                        p = p / s;
                        q = q / s;
                        r = r / s;
                        if (m == l)
                            break;
                        if (std::fabs(H(m, m - 1)) * (std::fabs(q) + std::fabs(r)) < eps * (std::fabs(p) * (std::fabs(H(m - 1, m - 1)) + std::fabs(z) + std::fabs(h_m1_m1))))
                            break;
                        m--;
                    }

                    for (int i = m + 2; i <= n; i++)
                    {
                        H(i, i - 2) = 0;
                        if (i > m + 2)
                            H(i, i - 3) = 0;
                    }

                    // float QR step involving rows l:n and columns m:n
                    for (int k = m; k < n; k++)
                    {
                        bool notlast = (k != n - 1);
                        if (k != m)
                        {
                            p = H(k, k - 1);
                            q = H(k + 1, k - 1);
                            r = (notlast ? H(k + 2, k - 1) : 0);
                            x = std::fabs(p) + std::fabs(q) + std::fabs(r);
                            if (x != 0)
                            {
                                p = p / x;
                                q = q / x;
                                r = r / x;
                            }
                        }

                        if (x == 0) break;

                        s = static_cast<SCALARTYPE>(std::sqrt(p * p + q * q + r * r));
                        if (p < 0) s = -s;

                        if (s != 0)
                        {
                            if (k != m)
                                H(k, k - 1) = -s * x;
                            else
                                if (l != m)
                                    H(k, k - 1) = -H(k, k - 1);

                            p = p + s;
                            y = q / s;
                            z = r / s;
                            x = p / s;
                            q = q / p;
                            r = r / p;

                            buf[5 * k] = x;
                            buf[5 * k + 1] = y;
                            buf[5 * k + 2] = z;
                            buf[5 * k + 3] = q;
                            buf[5 * k + 4] = r;


                            SCALARTYPE* a_row_k = H.row(k);
                            SCALARTYPE* a_row_k_1 = H.row(k + 1);
                            SCALARTYPE* a_row_k_2 = H.row(k + 2);
                            // Row modification
                            for (int j = k; j < nn; j++)
                            {
                                SCALARTYPE h_kj = a_row_k[j];
                                SCALARTYPE h_k1_j = a_row_k_1[j];

                                p = h_kj + q * h_k1_j;
                                if (notlast)
                                {
                                    SCALARTYPE h_k2_j = a_row_k_2[j];
                                    p = p + r * h_k2_j;
                                    a_row_k_2[j] = h_k2_j - p * z;
                                }

                                a_row_k[j] = h_kj - p * x;
                                a_row_k_1[j] = h_k1_j - p * y;
                            }

                            //H(k + 1, nn - 1) = h_kj;


                            // Column modification
                            for (int i = k; i < std::min(nn, k + 4); i++)
                            {
                                p = x * H(i, k) + y * H(i, k + 1);
                                if (notlast)
                                {
                                    p = p + z * H(i, k + 2);
                                    H(i, k + 2) = H(i, k + 2) - p * r;
                                }

                                H(i, k) = H(i, k) - p;
                                H(i, k + 1) = H(i, k + 1) - p * q;
                            }
                        }
                        else
                        {
                            buf[5 * k] = 0;
                            buf[5 * k + 1] = 0;
                            buf[5 * k + 2] = 0;
                            buf[5 * k + 3] = 0;
                            buf[5 * k + 4] = 0;
                        }
                    }

                    // Timer timer;
                    // timer.start();

                    update_float_QR_column(H, buf, m, n, n, true);
                    update_float_QR_column_gpu(V, buf, buf_vcl, m, n, nn, false);

                    // std::cout << timer.get() << "\n";
                }
            }

            // Backsubstitute to find vectors of upper triangular form
            if (norm == 0)
            {
                return;
            }

            for (n = nn - 1; n >= 0; n--)
            {
                p = d(n);
                q = e(n);

                // Real vector
                if (q == 0)
                {
                    int l = n;
                    H(n, n) = 1;
                    for (int i = n - 1; i >= 0; i--)
                    {
                        w = H(i, i) - p;
                        r = 0;
                        for (int j = l; j <= n; j++)
                            r = r + H(i, j) * H(j, n);

                        if (e(i) < 0)
                        {
                            z = w;
                            s = r;
                        }
                        else
                        {
                            l = i;
                            if (e(i) == 0)
                            {
                                H(i, n) = (w != 0) ? (-r / w) : (-r / (eps * norm));
                            }
                            else
                            {
                                // Solve real equations
                                x = H(i, i + 1);
                                y = H(i + 1, i);
                                q = (d(i) - p) * (d(i) - p) + e(i) * e(i);
                                t = (x * s - z * r) / q;
                                H(i, n) = t;
                                H(i + 1, n) = (std::fabs(x) > std::fabs(z)) ? ((-r - w * t) / x) : ((-s - y * t) / z);
                            }

                            // Overflow control
                            t = std::fabs(H(i, n));
                            if ((eps * t) * t > 1)
                                for (int j = i; j <= n; j++)
                                    H(j, n) /= t;
                        }
                    }
                }
                else if (q < 0)
                {
                    // Complex vector
                    int l = n - 1;

                    // Last vector component imaginary so matrix is triangular
                    if (std::fabs(H(n, n - 1)) > std::fabs(H(n - 1, n)))
                    {
                        H(n - 1, n - 1) = q / H(n, n - 1);
                        H(n - 1, n) = -(H(n, n) - p) / H(n, n - 1);
                    }
                    else
                    {
                        cdiv<SCALARTYPE>(0, -H(n - 1, n), H(n - 1, n - 1) - p, q, out1, out2);

                        H(n - 1, n - 1) = out1;
                        H(n - 1, n) = out2;
                    }

                    H(n, n - 1) = 0;
                    H(n, n) = 1;
                    for (int i = n - 2; i >= 0; i--)
                    {
                        SCALARTYPE ra, sa, vr, vi;
                        ra = 0;
                        sa = 0;
                        for (int j = l; j <= n; j++)
                        {
                            SCALARTYPE h_ij = H(i, j);
                            ra = ra + h_ij * H(j, n - 1);
                            sa = sa + h_ij * H(j, n);
                        }

                        w = H(i, i) - p;

                        if (e(i) < 0)
                        {
                            z = w;
                            r = ra;
                            s = sa;
                        }
                        else
                        {
                            l = i;
                            if (e(i) == 0)
                            {
                                cdiv<SCALARTYPE>(-ra, -sa, w, q, out1, out2);
                                H(i, n - 1) = out1;
                                H(i, n) = out2;
                            }
                            else
                            {
                                // Solve complex equations
                                x = H(i, i + 1);
                                y = H(i + 1, i);
                                vr = (d(i) - p) * (d(i) - p) + e(i) * e(i) - q * q;
                                vi = (d(i) - p) * 2 * q;
                                if ( (vr == 0) && (vi == 0) )
                                    vr = eps * norm * (std::fabs(w) + std::fabs(q) + std::fabs(x) + std::fabs(y) + std::fabs(z));

                                cdiv<SCALARTYPE>(x * r - z * ra + q * sa, x * s - z * sa - q * ra, vr, vi, out1, out2);

                                H(i, n - 1) = out1;
                                H(i, n) = out2;


                                if (std::fabs(x) > (std::fabs(z) + std::fabs(q)))
                                {
                                    H(i + 1, n - 1) = (-ra - w * H(i, n - 1) + q * H(i, n)) / x;
                                    H(i + 1, n) = (-sa - w * H(i, n) - q * H(i, n - 1)) / x;
                                }
                                else
                                {
                                    cdiv<SCALARTYPE>(-r - y * H(i, n - 1), -s - y * H(i, n), z, q, out1, out2);

                                    H(i + 1, n - 1) = out1;
                                    H(i + 1, n) = out2;
                                }
                            }

                            // Overflow control
                            t = std::max(std::fabs(H(i, n - 1)), std::fabs(H(i, n)));
                            if ((eps * t) * t > 1)
                            {
                                for (int j = i; j <= n; j++)
                                {
                                    H(j, n - 1) /= t;
                                    H(j, n) /= t;
                                }
                            }
                        }
                    }
                }
            }

            viennacl::fast_copy(H.begin(), H.end(),  vcl_H);
            // viennacl::fast_copy(V.begin(), V.end(),  vcl_V);

            viennacl::matrix<SCALARTYPE, row_major, ALIGNMENT> tmp = V;

            V = viennacl::linalg::prod(trans(tmp), vcl_H);
        }

        template <typename SCALARTYPE, unsigned int ALIGNMENT>
        bool householder_twoside(
                            viennacl::matrix<SCALARTYPE, row_major, ALIGNMENT>& A,
                            viennacl::matrix<SCALARTYPE, row_major, ALIGNMENT>& Q,
                            viennacl::vector<SCALARTYPE, ALIGNMENT>& D,
                            vcl_size_t start)
        {
            viennacl::ocl::context & ctx = const_cast<viennacl::ocl::context &>(viennacl::traits::opencl_handle(A).context());

            if(start + 2 >= A.size1())
                return false;

            prepare_householder_vector(A, D, A.size1(), start + 1, start, start + 1, true);

            {
                viennacl::ocl::kernel& kernel = ctx.get_kernel(viennacl::linalg::opencl::kernels::svd<SCALARTYPE>::program_name(), SVD_HOUSEHOLDER_UPDATE_A_LEFT_KERNEL);

                viennacl::ocl::enqueue(kernel(
                                              A,
                                              D,
                                              static_cast<cl_uint>(start + 1),
                                              static_cast<cl_uint>(start),
                                              static_cast<cl_uint>(A.size1()),
                                              static_cast<cl_uint>(A.size2()),
                                              static_cast<cl_uint>(A.internal_size2()),
                                              viennacl::ocl::local_mem(static_cast<cl_uint>(128 * 4))
                                      ));
            }

            {
                viennacl::ocl::kernel& kernel = ctx.get_kernel(viennacl::linalg::opencl::kernels::svd<SCALARTYPE>::program_name(), SVD_HOUSEHOLDER_UPDATE_A_RIGHT_KERNEL);

                viennacl::ocl::enqueue(kernel(
                                              A,
                                              D,
                                              static_cast<cl_uint>(0),
                                              static_cast<cl_uint>(0),
                                              static_cast<cl_uint>(A.size1()),
                                              static_cast<cl_uint>(A.size2()),
                                              static_cast<cl_uint>(A.internal_size2()),
                                              viennacl::ocl::local_mem(static_cast<cl_uint>(128 * sizeof(SCALARTYPE)))
                                      ));
            }

            {
                viennacl::ocl::kernel& kernel = ctx.get_kernel(viennacl::linalg::opencl::kernels::svd<SCALARTYPE>::program_name(), SVD_HOUSEHOLDER_UPDATE_QL_KERNEL);

                viennacl::ocl::enqueue(kernel(
                                                Q,
                                                D,
                                                static_cast<cl_uint>(A.size1()),
                                                static_cast<cl_uint>(A.size2()),
                                                static_cast<cl_uint>(Q.internal_size2()),
                                                viennacl::ocl::local_mem(static_cast<cl_uint>(128 * sizeof(SCALARTYPE)))
                                            ));
            }

            return true;
        }

        template <typename SCALARTYPE, typename F, unsigned int ALIGNMENT>
        void tridiagonal_reduction(viennacl::matrix<SCALARTYPE, F, ALIGNMENT>& A,
                                    viennacl::matrix<SCALARTYPE, F, ALIGNMENT>& Q)
        {
            vcl_size_t sz = A.size1();

            viennacl::vector<SCALARTYPE> hh_vector(sz);

            for(vcl_size_t i = 0; i < sz; i++)
            {
                householder_twoside(A, Q, hh_vector, i);
            }

        }

        template <typename SCALARTYPE, typename F, unsigned int ALIGNMENT>
        void qr_method(viennacl::matrix<SCALARTYPE, F, ALIGNMENT> & A,
                       viennacl::matrix<SCALARTYPE, F, ALIGNMENT> & Q,
                       boost::numeric::ublas::vector<SCALARTYPE> & D,
                       boost::numeric::ublas::vector<SCALARTYPE> & E,
                       bool is_symmetric = true)
        {
            viennacl::ocl::context & ctx = const_cast<viennacl::ocl::context &>(viennacl::traits::opencl_handle(A).context());

            assert(A.size1() == A.size2() && bool("Input matrix must be square for QR method!"));

            D.resize(A.size1());
            E.resize(A.size1());

            viennacl::linalg::opencl::kernels::svd<SCALARTYPE>::init(ctx);

            Q = viennacl::identity_matrix<SCALARTYPE>(Q.size1(), ctx);

            // reduce to tridiagonal form
            detail::tridiagonal_reduction(A, Q);

            // pack diagonal and super-diagonal
            // ublas::vector<SCALARTYPE> D(A.size1()), E(A.size1());

            bidiag_pack(A, D, E);

            // find eigenvalues
            if(is_symmetric)
            {

                detail::tql2(Q, D, E);
                transpose(Q);
            }
            else
            {
                detail::hqr2(A, Q, D, E);
            }

            // std::cout << A << "\n";

            boost::numeric::ublas::matrix<float> eigen_values(A.size1(), A.size1());
            eigen_values.clear();

            for (vcl_size_t i = 0; i < A.size1(); i++)
            {
                if(std::fabs(E(i)) < EPS)
                {
                    eigen_values(i, i) = D(i);
                }
                else
                {
                    eigen_values(i, i) = D(i);
                    eigen_values(i, i + 1) = E(i);
                    eigen_values(i + 1, i) = -E(i);
                    eigen_values(i + 1, i + 1) = D(i);
                    i++;
                }
            }

            copy(eigen_values, A);
        }
    }


    template <typename SCALARTYPE, typename F, unsigned int ALIGNMENT>
    void qr_method_nsm(viennacl::matrix<SCALARTYPE, F, ALIGNMENT>& A,
                       viennacl::matrix<SCALARTYPE, F, ALIGNMENT>& Q,
                       boost::numeric::ublas::vector<SCALARTYPE>& D,
                       boost::numeric::ublas::vector<SCALARTYPE>& E
                      )
    {
        detail::qr_method(A, Q, D, E, false);
    }

    template <typename SCALARTYPE, typename F, unsigned int ALIGNMENT>
    void qr_method_sym(viennacl::matrix<SCALARTYPE, F, ALIGNMENT>& A,
                       viennacl::matrix<SCALARTYPE, F, ALIGNMENT>& Q,
                       boost::numeric::ublas::vector<SCALARTYPE>& D
                      )
    {
        boost::numeric::ublas::vector<SCALARTYPE> E(A.size1());

        detail::qr_method(A, Q, D, E, true);
    }

  }
}

#endif
