#ifndef VIENNACL_LINALG_LANCZOS_HPP_
#define VIENNACL_LINALG_LANCZOS_HPP_

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

/** @file viennacl/linalg/lanczos.hpp
*   @brief Generic interface for the Lanczos algorithm.
*
*   Contributed by Guenther Mader and Astrid Rupp.
*/

#include <cmath>
#include <vector>
#include "viennacl/vector.hpp"
#include "viennacl/compressed_matrix.hpp"
#include "viennacl/linalg/prod.hpp"
#include "viennacl/linalg/inner_prod.hpp"
#include "viennacl/linalg/norm_2.hpp"
#include "viennacl/io/matrix_market.hpp"
#include "viennacl/linalg/bisect.hpp"
#include <boost/random.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/matrix_expression.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/operation.hpp>
#include <boost/numeric/ublas/vector_expression.hpp>
#include <boost/numeric/ublas/io.hpp>

namespace viennacl
{
  namespace linalg
  {

    /** @brief A tag for the lanczos algorithm.
    */
    class lanczos_tag
    {
      public:

        enum
        {
          partial_reorthogonalization = 0,
          full_reorthogonalization,
          no_reorthogonalization
        };

        /** @brief The constructor
        *
        * @param factor                 Exponent of epsilon - tolerance for batches of Reorthogonalization
        * @param numeig                 Number of eigenvalues to be returned
        * @param met                    Method for Lanczos-Algorithm: 0 for partial Reorthogonalization, 1 for full Reorthogonalization and 2 for Lanczos without Reorthogonalization
        * @param krylov                 Maximum krylov-space size
        */

        lanczos_tag(double factor = 0.75,
                    vcl_size_t numeig = 10,
                    int met = 0,
                    vcl_size_t krylov = 100) : factor_(factor), num_eigenvalues_(numeig), method_(met), krylov_size_(krylov) {}

        /** @brief Sets the number of eigenvalues */
        void num_eigenvalues(int numeig){ num_eigenvalues_ = numeig; }

          /** @brief Returns the number of eigenvalues */
        vcl_size_t num_eigenvalues() const { return num_eigenvalues_; }

          /** @brief Sets the exponent of epsilon */
        void factor(double fct) { factor_ = fct; }

        /** @brief Returns the exponent */
        double factor() const { return factor_; }

        /** @brief Sets the size of the kylov space */
        void krylov_size(int max) { krylov_size_ = max; }

        /** @brief Returns the size of the kylov space */
        vcl_size_t  krylov_size() const { return krylov_size_; }

        /** @brief Sets the reorthogonalization method */
        void method(int met){ method_ = met; }

        /** @brief Returns the reorthogonalization method */
        int method() const { return method_; }


      private:
        double factor_;
        vcl_size_t num_eigenvalues_;
        int method_; // see enum defined above for possible values
        vcl_size_t krylov_size_;

    };


    namespace detail
    {
      /**
      *   @brief Implementation of the Lanczos PRO algorithm
      *
      *   @param A            The system matrix
      *   @param r            Random start vector
      *   @param size         Size of krylov-space
      *   @param tag          Lanczos_tag with several options for the algorithm
      *   @return             Returns the eigenvalues (number of eigenvalues equals size of krylov-space)
      */

      template< typename MatrixT, typename VectorT >
      std::vector<
              typename viennacl::result_of::cpu_value_type<typename MatrixT::value_type>::type
              >
      lanczosPRO (MatrixT const& A, VectorT & r, vcl_size_t size, lanczos_tag const & tag)
      {

        typedef typename viennacl::result_of::value_type<MatrixT>::type        ScalarType;
        typedef typename viennacl::result_of::cpu_value_type<ScalarType>::type    CPU_ScalarType;


        // generation of some random numbers, used for lanczos PRO algorithm
        boost::mt11213b mt;
        boost::normal_distribution<CPU_ScalarType> N(0, 1);
        boost::bernoulli_distribution<CPU_ScalarType> B(0.5);
        boost::triangle_distribution<CPU_ScalarType> T(-1, 0, 1);

        boost::variate_generator<boost::mt11213b&, boost::normal_distribution<CPU_ScalarType> >     get_N(mt, N);
        boost::variate_generator<boost::mt11213b&, boost::bernoulli_distribution<CPU_ScalarType> >  get_B(mt, B);
        boost::variate_generator<boost::mt11213b&, boost::triangle_distribution<CPU_ScalarType> >   get_T(mt, T);


        long i, j, k, index, retry, reorths;
        std::vector<long> l_bound(size/2), u_bound(size/2);
        bool second_step;
        CPU_ScalarType squ_eps, eta, temp, eps, retry_th;
        vcl_size_t n = r.size();
        std::vector< std::vector<CPU_ScalarType> > w(2, std::vector<CPU_ScalarType>(size));
        CPU_ScalarType cpu_beta;

        boost::numeric::ublas::vector<CPU_ScalarType> s(n);

        VectorT t(n);
        CPU_ScalarType inner_rt;
        ScalarType vcl_beta;
        ScalarType vcl_alpha;
        std::vector<CPU_ScalarType> alphas, betas;
        boost::numeric::ublas::matrix<CPU_ScalarType> Q(n, size);

        second_step = false;
        eps = std::numeric_limits<CPU_ScalarType>::epsilon();
        squ_eps = std::sqrt(eps);
        retry_th = 1e-2;
        eta = std::exp(std::log(eps) * tag.factor());
        reorths = 0;
        retry = 0;

        vcl_beta = viennacl::linalg::norm_2(r);

        r /= vcl_beta;

        detail::copy_vec_to_vec(r,s);
        boost::numeric::ublas::column(Q, 0) = s;

        VectorT u = viennacl::linalg::prod(A, r);
        vcl_alpha = viennacl::linalg::inner_prod(u, r);
        alphas.push_back(vcl_alpha);
        w[0][0] = 1;
        betas.push_back(vcl_beta);

        long batches = 0;
        for(i = 1;i < static_cast<long>(size); i++)
        {
          r = u - vcl_alpha * r;
          vcl_beta = viennacl::linalg::norm_2(r);

          betas.push_back(vcl_beta);
          r = r / vcl_beta;

          index = i % 2;
          w[index][i] = 1;
          k = (i + 1) % 2;
          w[index][0] = (betas[1] * w[k][1] + (alphas[0] - vcl_alpha) * w[k][0] - betas[i - 1] * w[index][0]) / vcl_beta + eps * 0.3 * get_N() * (betas[1] + vcl_beta);

          for(j = 1;j < i - 1;j++)
          {
                  w[index][j] = (betas[j + 1] * w[k][j + 1] + (alphas[j] - vcl_alpha) * w[k][j] + betas[j] * w[k][j - 1] - betas[i - 1] * w[index][j]) / vcl_beta + eps * 0.3 * get_N() * (betas[j + 1] + vcl_beta);
          }
          w[index][i - 1] = 0.6 * eps * n * get_N() * betas[1] / vcl_beta;

          if(second_step)
          {
            for(j = 0;j < batches;j++)
            {
              l_bound[j]++;
              u_bound[j]--;

              for(k = l_bound[j];k < u_bound[j];k++)
              {
                detail::copy_vec_to_vec(boost::numeric::ublas::column(Q, k), t);
                inner_rt = viennacl::linalg::inner_prod(r,t);
                r = r - inner_rt * t;
                w[index][k] = 1.5 * eps * get_N();
                reorths++;
              }
            }
            temp = viennacl::linalg::norm_2(r);
            r = r / temp;
            vcl_beta = vcl_beta * temp;
            second_step = false;
          }
          batches = 0;

          for(j = 0;j < i;j++)
          {
            if(std::fabs(w[index][j]) >= squ_eps)
            {
              detail::copy_vec_to_vec(boost::numeric::ublas::column(Q, j), t);
              inner_rt = viennacl::linalg::inner_prod(r,t);
              r = r - inner_rt * t;
              w[index][j] = 1.5 * eps * get_N();
              k = j - 1;
              reorths++;
              while(k >= 0 && std::fabs(w[index][k]) > eta)
              {
                detail::copy_vec_to_vec(boost::numeric::ublas::column(Q, k), t);
                inner_rt = viennacl::linalg::inner_prod(r,t);
                r = r - inner_rt * t;
                w[index][k] = 1.5 * eps * get_N();
                k--;
                reorths++;
              }
              l_bound[batches] = k + 1;
              k = j + 1;

              while(k < i && std::fabs(w[index][k]) > eta)
              {
                detail::copy_vec_to_vec(boost::numeric::ublas::column(Q, k), t);
                inner_rt = viennacl::linalg::inner_prod(r,t);
                r = r - inner_rt * t;
                w[index][k] = 1.5 * eps * get_N();
                k++;
                reorths++;
              }
              u_bound[batches] = k - 1;
              batches++;
              j = k;
            }
          }

          if(batches > 0)
          {
            temp = viennacl::linalg::norm_2(r);
            r = r / temp;
            vcl_beta = vcl_beta * temp;
            second_step = true;

            while(temp < retry_th)
            {
              for(j = 0;j < i;j++)
              {
                detail::copy_vec_to_vec(boost::numeric::ublas::column(Q, k), t);
                inner_rt = viennacl::linalg::inner_prod(r,t);
                r = r - inner_rt * t;
                reorths++;
              }
              retry++;
              temp = viennacl::linalg::norm_2(r);
              r = r / temp;
              vcl_beta = vcl_beta * temp;
            }
          }

          detail::copy_vec_to_vec(r,s);
          boost::numeric::ublas::column(Q, i) = s;

          cpu_beta = vcl_beta;
          s = - cpu_beta * boost::numeric::ublas::column(Q, i - 1);
          detail::copy_vec_to_vec(s, u);
          u += viennacl::linalg::prod(A, r);
          vcl_alpha = viennacl::linalg::inner_prod(u, r);
          alphas.push_back(vcl_alpha);
        }

        return bisect(alphas, betas);

      }


      /**
      *   @brief Implementation of the lanczos algorithm without reorthogonalization
      *
      *   @param A            The system matrix
      *   @param r            Random start vector
      *   @param size         Size of krylov-space
      *   @return             Returns the eigenvalues (number of eigenvalues equals size of krylov-space)
      */
      template< typename MatrixT, typename VectorT >
      std::vector<
              typename viennacl::result_of::cpu_value_type<typename MatrixT::value_type>::type
              >
      lanczos (MatrixT const& A, VectorT & r, vcl_size_t size, lanczos_tag)
      {

        typedef typename viennacl::result_of::value_type<MatrixT>::type        ScalarType;
        typedef typename viennacl::result_of::cpu_value_type<ScalarType>::type    CPU_ScalarType;

        ScalarType vcl_beta;
        ScalarType vcl_alpha;
        std::vector<CPU_ScalarType> alphas, betas;
        CPU_ScalarType norm;
        vcl_size_t n = r.size();
        VectorT u(n), t(n);
        boost::numeric::ublas::vector<CPU_ScalarType> s(r.size()), u_zero(n), q(n);
        boost::numeric::ublas::matrix<CPU_ScalarType> Q(n, size);

        u_zero = boost::numeric::ublas::zero_vector<CPU_ScalarType>(n);
        detail::copy_vec_to_vec(u_zero, u);
        norm = norm_2(r);

        for(vcl_size_t i = 0;i < size; i++)
        {
          r /= norm;
          vcl_beta = norm;

          detail::copy_vec_to_vec(r,s);
          boost::numeric::ublas::column(Q, i) = s;

          u += prod(A, r);
          vcl_alpha = inner_prod(u, r);
          r = u - vcl_alpha * r;
          norm = norm_2(r);

          q = boost::numeric::ublas::column(Q, i);
          detail::copy_vec_to_vec(q, t);

          u = - norm * t;
          alphas.push_back(vcl_alpha);
          betas.push_back(vcl_beta);
          s.clear();
        }

        return bisect(alphas, betas);
      }

      /**
      *   @brief Implementation of the Lanczos FRO algorithm
      *
      *   @param A            The system matrix
      *   @param r            Random start vector
      *   @param size         Size of krylov-space
      *   @return             Returns the eigenvalues (number of eigenvalues equals size of krylov-space)
      */
      template< typename MatrixT, typename VectorT >
      std::vector<
              typename viennacl::result_of::cpu_value_type<typename MatrixT::value_type>::type
              >
      lanczosFRO (MatrixT const& A, VectorT & r, vcl_size_t size, lanczos_tag)
      {

        typedef typename viennacl::result_of::value_type<MatrixT>::type        ScalarType;
        typedef typename viennacl::result_of::cpu_value_type<ScalarType>::type    CPU_ScalarType;

          CPU_ScalarType temp;
          CPU_ScalarType norm;
          ScalarType vcl_beta;
          ScalarType vcl_alpha;
          std::vector<CPU_ScalarType> alphas, betas;
          vcl_size_t n = r.size();
          VectorT u(n), t(n);
          ScalarType inner_rt;
          boost::numeric::ublas::vector<CPU_ScalarType> u_zero(n), s(r.size()), q(n);
          boost::numeric::ublas::matrix<CPU_ScalarType> Q(n, size);

          long reorths = 0;
          norm = norm_2(r);


          for(vcl_size_t i = 0; i < size; i++)
          {
            r /= norm;

            for(vcl_size_t j = 0; j < i; j++)
            {
              q = boost::numeric::ublas::column(Q, j);
              detail::copy_vec_to_vec(q, t);
              inner_rt = viennacl::linalg::inner_prod(r,t);
              r = r - inner_rt * t;
              reorths++;
            }
            temp = viennacl::linalg::norm_2(r);
            r = r / temp;
            vcl_beta = temp * norm;
            detail::copy_vec_to_vec(r,s);
            boost::numeric::ublas::column(Q, i) = s;

            u += viennacl::linalg::prod(A, r);
            vcl_alpha = viennacl::linalg::inner_prod(u, r);
            r = u - vcl_alpha * r;
            norm = viennacl::linalg::norm_2(r);
            q = boost::numeric::ublas::column(Q, i);
            detail::copy_vec_to_vec(q, t);
            u = - norm * t;
            alphas.push_back(vcl_alpha);
            betas.push_back(vcl_beta);
          }

          return bisect(alphas, betas);
      }

    } // end namespace detail

    /**
    *   @brief Implementation of the calculation of eigenvalues using lanczos
    *
    *   @param matrix        The system matrix
    *   @param tag           Tag with several options for the lanczos algorithm
    *   @return              Returns the n largest eigenvalues (n defined in the lanczos_tag)
    */
    template< typename MatrixT >
    std::vector< typename viennacl::result_of::cpu_value_type<typename MatrixT::value_type>::type >
    eig(MatrixT const & matrix, lanczos_tag const & tag)
    {
      typedef typename viennacl::result_of::value_type<MatrixT>::type           ScalarType;
      typedef typename viennacl::result_of::cpu_value_type<ScalarType>::type    CPU_ScalarType;
      typedef typename viennacl::result_of::vector_for_matrix<MatrixT>::type    VectorT;

      boost::mt11213b mt;
      boost::normal_distribution<CPU_ScalarType> N(0, 1);
      boost::bernoulli_distribution<CPU_ScalarType> B(0.5);
      boost::triangle_distribution<CPU_ScalarType> T(-1, 0, 1);

      boost::variate_generator<boost::mt11213b&, boost::normal_distribution<CPU_ScalarType> >     get_N(mt, N);
      boost::variate_generator<boost::mt11213b&, boost::bernoulli_distribution<CPU_ScalarType> >  get_B(mt, B);
      boost::variate_generator<boost::mt11213b&, boost::triangle_distribution<CPU_ScalarType> >   get_T(mt, T);

      std::vector<CPU_ScalarType> eigenvalues;
      vcl_size_t matrix_size = matrix.size1();
      VectorT r(matrix_size);
      std::vector<CPU_ScalarType> s(matrix_size);

      for(vcl_size_t i=0; i<s.size(); ++i)
        s[i] = 3.0 * get_B() + get_T() - 1.5;

      detail::copy_vec_to_vec(s,r);

      vcl_size_t size_krylov = (matrix_size < tag.krylov_size()) ? matrix_size
                                                                  : tag.krylov_size();

      switch(tag.method())
      {
        case lanczos_tag::partial_reorthogonalization:
          eigenvalues = detail::lanczosPRO(matrix, r, size_krylov, tag);
          break;
        case lanczos_tag::full_reorthogonalization:
          eigenvalues = detail::lanczosFRO(matrix, r, size_krylov, tag);
          break;
        case lanczos_tag::no_reorthogonalization:
          eigenvalues = detail::lanczos(matrix, r, size_krylov, tag);
          break;
      }

      std::vector<CPU_ScalarType> largest_eigenvalues;

      for(vcl_size_t i = 1; i<=tag.num_eigenvalues(); i++)
        largest_eigenvalues.push_back(eigenvalues[size_krylov-i]);


      return largest_eigenvalues;
    }




  } // end namespace linalg
} // end namespace viennacl
#endif
