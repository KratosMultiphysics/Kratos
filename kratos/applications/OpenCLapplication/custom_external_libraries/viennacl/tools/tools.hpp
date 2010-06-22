/* =======================================================================
   Copyright (c) 2010, Institute for Microelectronics, TU Vienna.
   http://www.iue.tuwien.ac.at
                             -----------------
                     ViennaCL - The Vienna Computing Library
                             -----------------
                            
   authors:    Karl Rupp                          rupp@iue.tuwien.ac.at
               Florian Rudolf                     flo.rudy+viennacl@gmail.com
               Josef Weinbub                      weinbub@iue.tuwien.ac.at

   license:    MIT (X11), see file LICENSE in the ViennaCL base directory

   file changelog: - May 28, 2010   New from scratch for first release
======================================================================= */

#ifndef _VIENNACL_TOOLS_HPP_
#define _VIENNACL_TOOLS_HPP_

#include <string>
#include <fstream>
#include <sstream>
#include "viennacl/forwards.h"

#include <vector>
#include <map>

namespace viennacl
{
  namespace tools
  {
    
    //supply suitable increment functions for the iterators:
    template <class SCALARTYPE, typename F, unsigned int ALIGNMENT>
    struct MATRIX_ITERATOR_INCREMENTER<viennacl::row_iteration, viennacl::matrix<SCALARTYPE, F, ALIGNMENT> >
    {
      static void apply(const viennacl::matrix<SCALARTYPE, F, ALIGNMENT> & mat, unsigned int & row, unsigned int & col) { ++row; }
    };

    template <class SCALARTYPE, typename F, unsigned int ALIGNMENT>
    struct MATRIX_ITERATOR_INCREMENTER<viennacl::col_iteration, viennacl::matrix<SCALARTYPE, F, ALIGNMENT> >
    {
      static void apply(const viennacl::matrix<SCALARTYPE, F, ALIGNMENT> & mat, unsigned int & row, unsigned int & col) { ++col; }
    };

    
    
    template <bool b, class T = void> 
    struct enable_if
    {
      typedef T   type;
    };

    template <class T> 
    struct enable_if<false, T> {};
    
    
    /** @brief A guard that checks whether the floating point type of GPU types is either float or double */
    template <typename T>
    struct CHECK_SCALAR_TEMPLATE_ARGUMENT
    {
        typedef typename T::ERROR_SCALAR_MUST_HAVE_TEMPLATE_ARGUMENT_FLOAT_OR_DOUBLE  ResultType;
    };
    
    template <>
    struct CHECK_SCALAR_TEMPLATE_ARGUMENT<float>
    {
        typedef float  ResultType;
    };
    
    template <>
    struct CHECK_SCALAR_TEMPLATE_ARGUMENT<double>
    {
        typedef double  ResultType;
    };

    
    /** @brief A const iterator for sparse matrices of type std::vector<std::map<unsigned int, SCALARTYPE> >
    *  
    *  The iterator behaves like ublas iterators. Attention: Iteration along first columns and then rows via .begin() is untested!
    *
    *  @tparam SCALARTYPE     either float or double
    *  @tparam is_iterator1   if true, this iterator iterates along increasing row indices, otherwise along increasiong column indices
    *  @tparam increment      if +1, this is a forward iterator, if -1 we have a reverse iterator
    */
    template <typename SCALARTYPE, bool is_iterator1, int increment>
    class const_sparse_matrix_adapted_iterator
    {
      typedef const_sparse_matrix_adapted_iterator<SCALARTYPE, is_iterator1, increment>    self_type;
      
      public:
        typedef self_type     iterator1;
        typedef self_type     iterator2;
        
        const_sparse_matrix_adapted_iterator(std::vector<std::map<unsigned int, SCALARTYPE> > const & mat, int i, int j)
         : _mat(mat), _i(i), _j(j)
        {
          if (i < static_cast<int>(_mat.size()) && j < static_cast<int>(_mat.size()))
            iter2 = _mat[i].begin();
          else
            iter2 = _mat[i].end();
        }
         
        SCALARTYPE operator*(void) const
        {
          if (is_iterator1)
          {
            typedef typename std::map<unsigned int, SCALARTYPE>::const_iterator  col_iterator;
            
            col_iterator colit = _mat[_i].find(_j);

            if (colit != _mat[_i].end())
              return colit->second;
            return 0.0;
          }
          else
            return iter2->second;
        }
        
        self_type & operator++(void)
        {
          if (is_iterator1)
            _i += increment;
          else
            ++iter2;
          return *this;
        }
        self_type & operator++(int) { self_type tmp = *this; ++(*this); return tmp; }
        
        bool operator==(self_type const & other) const
        {
          if (is_iterator1)
            return (_i == other._i);
          return (iter2 == other.iter2);
        }
        bool operator!=(self_type const & other) const { return !(*this == other); }
        
        unsigned int index1() const { return _i; }
        unsigned int index2() const
        { 
          if (is_iterator1)
            return 0;
          else
            return iter2->first;
        }
        
        const_sparse_matrix_adapted_iterator<SCALARTYPE, !is_iterator1, 1> begin() const
        {
          return const_sparse_matrix_adapted_iterator<SCALARTYPE, !is_iterator1, 1>(_mat, _i, iter2->first);
        }
        const_sparse_matrix_adapted_iterator<SCALARTYPE, !is_iterator1, 1> end() const
        {
          return const_sparse_matrix_adapted_iterator<SCALARTYPE, !is_iterator1, 1>(_mat, _i, _mat.size());
        }
        
      private:
        std::vector<std::map<unsigned int, SCALARTYPE> > const & _mat;
        typename std::map<unsigned int, SCALARTYPE>::const_iterator iter2;
        int _i;
        int _j;
    };
    
    /** @brief Adapts a constant sparse matrix type made up from std::vector<std::map<unsigned int, SCALARTYPE> > to basic ublas-compatibility.
    *
    *  @tparam SCALARTYPE   either float or double
    */
    template <typename SCALARTYPE>
    class const_sparse_matrix_adapter
    {
      public:
        typedef const_sparse_matrix_adapted_iterator<SCALARTYPE, true, 1>      const_iterator1;
        typedef const_sparse_matrix_adapted_iterator<SCALARTYPE, false, 1>     const_iterator2;

        typedef const_sparse_matrix_adapted_iterator<SCALARTYPE, true, -1>   const_reverse_iterator1;
        typedef SCALARTYPE    value_type;
        
        const_sparse_matrix_adapter(std::vector<std::map<unsigned int, SCALARTYPE> > const & mat) 
         : _mat(mat) {};
        
        unsigned int size1() const { return _mat.size(); }
        unsigned int size2() const { return _mat.size(); }  //we allow only square matrices
        
        const_iterator1 begin1() const { return const_iterator1(_mat, 0, 0); }
        const_iterator1 end1() const   { return const_iterator1(_mat, size1(), size2()); }

        const_reverse_iterator1 rbegin1() const { return const_reverse_iterator1(_mat, size1() - 1, 0); }
        const_reverse_iterator1 rend1() const   { return const_reverse_iterator1(_mat, -1, size2()); }

        const_iterator2 begin2() const { return const_iterator2(_mat, 0, 0); }
        const_iterator2 end2() const   { return const_iterator2(_mat, 0, size2()); }
        
        SCALARTYPE operator()(unsigned int i, unsigned int j) const
        {
          typedef typename std::map<unsigned int, SCALARTYPE>::const_iterator  col_iterator;
          
          col_iterator colit = _mat[i].find(j);

          if (colit != _mat[i].end())
            return colit->second;
          return 0.0;
        }

      private:
        std::vector<std::map<unsigned int, SCALARTYPE> > const & _mat;
    };
    
    
    /** @brief A non-const iterator for sparse matrices of type std::vector<std::map<unsigned int, SCALARTYPE> >
    *  
    *  The iterator behaves like ublas iterators. Attention: Iteration along first columns and then rows via .begin() is untested! Reverse iterators are missing!
    *
    *  @tparam SCALARTYPE     either float or double
    *  @tparam is_iterator1   if true, this iterator iterates along increasing row indices, otherwise along increasiong column indices
    */
    template <typename SCALARTYPE, bool is_iterator1>
    class sparse_matrix_adapted_iterator
    {
      typedef sparse_matrix_adapted_iterator<SCALARTYPE, is_iterator1>    self_type;
      
      public:
        typedef self_type     iterator1;
        typedef self_type     iterator2;
        
        sparse_matrix_adapted_iterator(std::vector<std::map<unsigned int, SCALARTYPE> > & mat, unsigned int i, unsigned int j)
         : _mat(mat), _i(i), _j(j)
        {
          if (i < _mat.size() && j < _mat.size())
            iter2 = _mat[i].begin();
          else
            iter2 = _mat[i].end();
        }
         
        SCALARTYPE & operator*(void)
        {
          if (is_iterator1)
          {
            return _mat[_i][_j];
          }
          else
            return iter2->second;
        }
        
        self_type & operator++(void)
        {
          if (is_iterator1)
            ++_i;
          else
            ++iter2;
          return *this;
        }
        self_type & operator++(int) { self_type tmp = *this; ++(*this); return tmp; }
        
        bool operator==(self_type const & other) const
        {
          if (is_iterator1)
            return (_i == other._i);
          return (iter2 == other.iter2);
        }
        bool operator!=(self_type const & other) const { return !(*this == other); }
        
        unsigned int index1() const { return _i; }
        unsigned int index2() const
        { 
          if (is_iterator1)
            return 0;
          else
            return iter2->first;
        }
        
        sparse_matrix_adapted_iterator<SCALARTYPE, !is_iterator1> begin() const
        {
          return sparse_matrix_adapted_iterator<SCALARTYPE, !is_iterator1>(_mat, _i, iter2->first);
        }
        sparse_matrix_adapted_iterator<SCALARTYPE, !is_iterator1> end() const
        {
          return sparse_matrix_adapted_iterator<SCALARTYPE, !is_iterator1>(_mat, _i, _mat.size());
        }
        
      private:
        std::vector<std::map<unsigned int, SCALARTYPE> > & _mat;
        typename std::map<unsigned int, SCALARTYPE>::iterator iter2;
        unsigned int _i;
        unsigned int _j;
    };
    
    
    
    /** @brief Adapts a non-const sparse matrix type made up from std::vector<std::map<unsigned int, SCALARTYPE> > to basic ublas-compatibility.
    *
    *  @tparam SCALARTYPE   either float or double
    */
    template <typename SCALARTYPE>
    class sparse_matrix_adapter : public const_sparse_matrix_adapter<SCALARTYPE>
    {
        typedef const_sparse_matrix_adapter<SCALARTYPE>   BaseType;
      public:
        typedef sparse_matrix_adapted_iterator<SCALARTYPE, true>      iterator1;
        typedef sparse_matrix_adapted_iterator<SCALARTYPE, false>     iterator2;
        
        sparse_matrix_adapter(std::vector<std::map<unsigned int, SCALARTYPE> > & mat) 
         : BaseType(mat), _mat(mat) { };
        
        iterator1 begin1() const { return iterator1(_mat, 0, 0); }
        iterator1 end1() const   { return iterator1(_mat, _mat.size(), _mat.size()); }

        iterator2 begin2() const { return iterator2(_mat, 0, 0); }
        iterator2 end2() const   { return iterator2(_mat, _mat.size(), _mat.size()); }
        
        SCALARTYPE & operator()(unsigned int i, unsigned int j) { return _mat[i][j]; }
        
        void resize(unsigned int i, unsigned int j, bool preserve = true)
        {
          if (i>0)
            _mat.resize(i);
          if (!preserve)
            clear();
        }
        
        void clear()
        {
          for (unsigned int i=0; i<_mat.size(); ++i)
            _mat[i].clear();
        }
        
      private:
        std::vector<std::map<unsigned int, SCALARTYPE> > & _mat;
    };
    
    
    /** @brief Reads a text from a file into a std::string
    *
    * @param filename   The filename
    * @return The text read from the file
    */
    std::string readTextFromFile(const std::string & filename)
    {
      std::ifstream f(filename.c_str());
      if (!f) return std::string();

      std::stringstream result;
      std::string tmp;
      while (getline(f, tmp))
        result << tmp << std::endl;

      return result.str();
    }

    /** @brief Replaces all occurances of a substring by another stringstream
    *
    * @param text   The string to search in
    * @param to_search  The substring to search for
    * @param to_replace The replacement for found substrings
    * @return The resulting string
    */
    std::string strReplace(const std::string & text, std::string to_search, std::string to_replace)
    {
      std::string::size_type pos=0;
        std::string result;
        std::string::size_type found;
        while( (found = text.find(to_search, pos)) != std::string::npos )
        {
        result.append(text.substr(pos,found-pos));
            result.append(to_replace);
            pos = found + to_search.length();
        }
        if (pos < text.length())
            result.append(text.substr(pos));
        return result;
    }

    /** @brief Rounds an integer to the next multiple of another integer
    *
    * @tparam INT_TYPE  The integer type
    * @param to_reach   The integer to be rounded up (ceil operation)
    * @param base       The base
    * @return The smallest multiple of 'base' such that to_reach <= base
    */
    template <class INT_TYPE>
    INT_TYPE roundUpToNextMultiple(INT_TYPE to_reach, INT_TYPE base)
    {
      if (to_reach % base == 0) return to_reach;
      return ((to_reach / base) + 1) * base;
    }
    
    
    /** @brief Create a double precision kernel out of a single precision kernel
    *
    * @param source   The source string
    * @return   The double precision kernel
    */
    std::string make_double_kernel(std::string const & source)
    {
      #ifdef VIENNACL_EXPERIMENTAL_DOUBLE_PRECISION_WITH_STREAM_SDK
        std::string result = "#pragma OPENCL EXTENSION cl_amd_fp64 : enable\n\n";
      #else
        std::string result = "#pragma OPENCL EXTENSION cl_khr_fp64 : enable\n\n";
      #endif
      result.append(strReplace(source, "float", "double"));
      return result;
    }

    /** @brief Removes the const qualifier from a type */
    template <typename T>
    struct CONST_REMOVER
    {
      typedef T   ResultType;
    };

    template <typename T>
    struct CONST_REMOVER<const T>
    {
      typedef T   ResultType;
    };


    /** @brief Extracts the vector type from one of the two arguments. Used for the vector_expression type.
    *
    * @tparam LHS   The left hand side operand of the vector_expression
    * @tparam RHS   The right hand side operand of the vector_expression
    */
    template <typename LHS, typename RHS>
    struct VECTOR_EXTRACTOR_IMPL
    {
      typedef typename LHS::ERROR_COULD_NOT_EXTRACT_VECTOR_INFORMATION_FROM_VECTOR_EXPRESSION  ResultType;
    };
    
    template <typename LHS, typename ScalarType, unsigned int A>
    struct VECTOR_EXTRACTOR_IMPL<LHS, viennacl::vector<ScalarType, A> >
    {
      typedef viennacl::vector<ScalarType, A>   ResultType;
    };

    template <typename RHS, typename ScalarType, unsigned int A>
    struct VECTOR_EXTRACTOR_IMPL<viennacl::vector<ScalarType, A>, RHS>
    {
      typedef viennacl::vector<ScalarType, A>   ResultType;
    };

    //resolve ambiguities for previous cases:
    template <typename ScalarType, unsigned int A>
    struct VECTOR_EXTRACTOR_IMPL<viennacl::vector<ScalarType, A>, viennacl::vector<ScalarType, A> >
    {
      typedef viennacl::vector<ScalarType, A>   ResultType;
    };

    template <typename LHS, typename RHS>
    struct VECTOR_EXTRACTOR
    {
      typedef typename VECTOR_EXTRACTOR_IMPL<typename CONST_REMOVER<LHS>::ResultType,
                                              typename CONST_REMOVER<RHS>::ResultType>::ResultType      ResultType;
    };

    /** @brief Deduces the size of the resulting vector represented by a vector_expression from the operands
    *
    * @tparam LHS   The left hand side operand
    * @tparam RHS   The right hand side operand
    * @tparam OP    The operation tag
    */
    template <typename LHS, typename RHS, typename OP>
    struct VECTOR_SIZE_DEDUCER
    {
      //Standard case: LHS is the vector type and carries the correct size
      static unsigned int size(LHS & lhs, RHS & rhs) { return lhs.size(); }
    };
    
    //special case: matrix-vector product: Return the number of rows of the matrix
    template <typename MatrixType, typename ScalarType, unsigned int A>
    struct VECTOR_SIZE_DEDUCER<MatrixType, viennacl::vector<ScalarType, A>, viennacl::op_prod>
    {
      static unsigned int size(MatrixType & lhs, viennacl::vector<ScalarType, A> & rhs) { return lhs.size1(); }
    };

    //special case: matrix-vector product: Return the number of rows of the matrix
    template <typename MatrixType, typename ScalarType, unsigned int A>
    struct VECTOR_SIZE_DEDUCER<MatrixType, const viennacl::vector<ScalarType, A>, viennacl::op_prod>
    {
      static unsigned int size(MatrixType & lhs, const viennacl::vector<ScalarType, A> & rhs) { return lhs.size1(); }
    };

    /** @brief Obtain the cpu scalar type from a type, including a GPU type like viennacl::scalar<T>
    *
    * @tparam T   Either a CPU scalar type or a GPU scalar type
    */
    template <typename T>
    struct CPU_SCALAR_TYPE_DEDUCER
    {
      typedef T       ResultType;
    };
    
    template <typename T>
    struct CPU_SCALAR_TYPE_DEDUCER< viennacl::scalar<T> >
    {
      typedef T       ResultType;
    };
    
  } //namespace tools
} //namespace viennacl
    

#endif
