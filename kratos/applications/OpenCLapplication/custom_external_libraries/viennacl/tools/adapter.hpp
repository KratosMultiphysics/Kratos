#ifndef VIENNACL_TOOLS_ADAPTER_HPP_
#define VIENNACL_TOOLS_ADAPTER_HPP_

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

/** @file viennacl/tools/adapter.hpp
    @brief Adapter classes for sparse matrices made of the STL type std::vector<std::map<SizeType, SCALARTYPE> >
*/

#include <string>
#include <fstream>
#include <sstream>
#include <assert.h>
#include "viennacl/forwards.h"

#include <vector>
#include <map>

namespace viennacl
{
  namespace tools
  {

    /** @brief A const iterator for sparse matrices of type std::vector<std::map<SizeType, SCALARTYPE> >
    *
    *  The iterator behaves like ublas iterators. Attention: Iteration along first columns and then rows via .begin() is untested!
    *
    *  @tparam SCALARTYPE     either float or double
    *  @tparam is_iterator1   if true, this iterator iterates along increasing row indices, otherwise along increasing column indices
    *  @tparam increment      if +1, this is a forward iterator, if -1 we have a reverse iterator
    */
    template <typename SCALARTYPE, typename SizeType, bool is_iterator1, bool is_forward>
    class const_sparse_matrix_adapted_iterator
    {
      typedef const_sparse_matrix_adapted_iterator<SCALARTYPE, SizeType, is_iterator1, is_forward>    self_type;

      public:
        typedef self_type     iterator1;
        typedef self_type     iterator2;
        typedef vcl_size_t   size_type;

        const_sparse_matrix_adapted_iterator(std::vector<std::map<SizeType, SCALARTYPE> > const & mat, int i, int j)
         : mat_(mat), i_(i), j_(j)
        {
          if (i < 0) //reverse iterator end
          {
            //iter2 = mat_[0].rend();  //reverse iterator end
          }
          else  //i_ is valid
          {
            if (j < 0)
            {
              //iter2 = mat_[i].rend();
            }
            else //j_ is valid
            {
              if (i_ < mat_.size() && mat_[i].size() > 0 )
              {
                //TODO: Start at entry j, not at the beginning
                if (static_cast<int>(mat_[i].rbegin()->first) < j)
                  iter2 = mat_[i].end();
                else
                  iter2 = mat_[i].begin();
              }
              else if (i_ < mat_.size() && mat_[i].size() == 0)
                iter2 = mat_[i].end();
              else //i is out of range -> end iterator requested
                iter2 = mat_.back().end(); //forward iterator end
            }
          }
        }

        SCALARTYPE operator*(void) const
        {
          if (is_iterator1)
          {
            typedef typename std::map<SizeType, SCALARTYPE>::const_iterator  col_iterator;

            col_iterator colit = mat_[i_].find(static_cast<unsigned int>(j_));

            if (colit != mat_[i_].end())
              return colit->second;
            return 0.0;
          }
          else
            return iter2->second;
        }

        self_type & operator++(void)
        {
          if (is_iterator1)
          {
            if (is_forward)
              ++i_;
            else
              --i_;
          }
          else
            ++iter2;
          return *this;
        }
        self_type operator++(int) { self_type tmp = *this; ++(*this); return tmp; }

        self_type operator+=(SizeType offset)
        {
          if (is_iterator1)
          {
            if (is_forward)
              i_ += offset;
            else
              i_ -= offset;
          }
          else
          {
            for (SizeType k=0; k<offset; ++k)
              ++iter2;  //Note: User must ensure that this is always valid...
          }
          return *this;
        }

        bool operator==(self_type const & other) const
        {
          return is_iterator1 ? (i_ == other.i_) : (iter2 == other.iter2);
        }

        bool operator!=(self_type const & other) const { return !(*this == other); }

        size_type index1() const { return i_; }
        size_type index2() const
        {
          if (is_iterator1)
            return 0;
          else
            return iter2->first;
        }

        const_sparse_matrix_adapted_iterator<SCALARTYPE, SizeType, !is_iterator1, true> begin() const
        {
          return const_sparse_matrix_adapted_iterator<SCALARTYPE, SizeType, !is_iterator1, true>(mat_, static_cast<int>(i_), 0);
        }
        const_sparse_matrix_adapted_iterator<SCALARTYPE, SizeType, !is_iterator1, true> end() const
        {
          int end_ = static_cast<int>(mat_[i_].size());
          if (end_ > 0)
            end_ = mat_[i_].rbegin()->first;
          return const_sparse_matrix_adapted_iterator<SCALARTYPE, SizeType, !is_iterator1, true>(mat_, static_cast<int>(i_), end_ + 1);
        }

      private:
        std::vector<std::map<SizeType, SCALARTYPE> > const & mat_;
        typename std::map<SizeType, SCALARTYPE>::const_iterator iter2;
        size_type i_;
        size_type j_;
    };

    /** @brief Adapts a constant sparse matrix type made up from std::vector<std::map<SizeType, SCALARTYPE> > to basic ublas-compatibility.
    *
    *  @tparam SCALARTYPE   either float or double
    */
    template <typename SCALARTYPE, typename SizeType = unsigned int>
    class const_sparse_matrix_adapter
    {
      public:
        typedef const_sparse_matrix_adapted_iterator<SCALARTYPE, SizeType, true, true>      const_iterator1;
        typedef const_sparse_matrix_adapted_iterator<SCALARTYPE, SizeType, false, true>     const_iterator2;

        typedef const_sparse_matrix_adapted_iterator<SCALARTYPE, SizeType, true, false>   const_reverse_iterator1;
        typedef SCALARTYPE    value_type;
        typedef vcl_size_t   size_type;

        const_sparse_matrix_adapter(std::vector<std::map<SizeType, SCALARTYPE> > const & mat)
         : mat_(mat), size1_(mat_.size()), size2_(mat_.size()) {}

        const_sparse_matrix_adapter(std::vector<std::map<SizeType, SCALARTYPE> > const & mat, size_type num_rows, size_type num_cols)
         : mat_(mat), size1_(num_rows), size2_(num_cols) {}

        size_type size1() const { return size1_; }
        size_type size2() const { return size2_; }

        const_iterator1 begin1() const { return const_iterator1(mat_, 0, 0); }
        const_iterator1 end1() const   { return const_iterator1(mat_, static_cast<int>(size1()), static_cast<int>(size2())); }

        const_reverse_iterator1 rbegin1() const { return const_reverse_iterator1(mat_, static_cast<int>(size1() - 1), 0); }
        const_reverse_iterator1 rend1() const   { return const_reverse_iterator1(mat_, -1, static_cast<int>(size2())); }

        const_iterator2 begin2() const { return const_iterator2(mat_, 0, 0); }
        const_iterator2 end2() const   { return const_iterator2(mat_, size1(), size2()); }

        SCALARTYPE operator()(SizeType i, SizeType j) const
        {
          typedef typename std::map<SizeType, SCALARTYPE>::const_iterator  col_iterator;

          col_iterator colit = mat_[i].find(j);

          if (colit != mat_[i].end())
            return colit->second;
          return 0.0;
        }

      private:
        std::vector<std::map<SizeType, SCALARTYPE> > const & mat_;
        size_type size1_;
        size_type size2_;
    };


    /** @brief A non-const iterator for sparse matrices of type std::vector<std::map<SizeType, SCALARTYPE> >
    *
    *  The iterator behaves like ublas iterators. Attention: Iteration along first columns and then rows via .begin() is untested! Reverse iterators are missing!
    *
    *  @tparam SCALARTYPE     either float or double
    *  @tparam is_iterator1   if true, this iterator iterates along increasing row indices, otherwise along increasiong column indices
    */
    template <typename SCALARTYPE, typename SizeType, bool is_iterator1>
    class sparse_matrix_adapted_iterator
    {
      typedef sparse_matrix_adapted_iterator<SCALARTYPE, SizeType, is_iterator1>    self_type;

      public:
        typedef self_type     iterator1;
        typedef self_type     iterator2;
        typedef vcl_size_t   size_type;

        sparse_matrix_adapted_iterator(std::vector<std::map<SizeType, SCALARTYPE> > & mat, int i, int j)
         : mat_(mat), i_(i), j_(j)
        {
          if (i < 0) //reverse iterator end
          {
            //iter2 = mat_[0].rend();  //reverse iterator end
          }
          else  //_i is valid
          {
            if (j < 0)
            {
              //iter2 = mat[i]_.rend();
            }
            else //_j is valid
            {
              if (i_ < mat_.size() && mat_[i].size() > 0 )
              {
                //TODO: Start at entry j, not at the beginning
                if (static_cast<int>(mat_[i].rbegin()->first) < j)
                  iter2 = mat_[i].end();
                else
                  iter2 = mat_[i].begin();
              }
              else if (i_ < mat_.size() && mat_[i].size() == 0)
                iter2 = mat_[i].end();
              else //i is out of range -> end iterator requested
                iter2 = mat_.back().end(); //forward iterator end
            }
          }
        }

        SCALARTYPE & operator*(void)
        {
          if (is_iterator1)
          {
            return mat_[i_][static_cast<SizeType>(j_)];
          }
          else
            return iter2->second;
        }

        self_type & operator++(void)
        {
          if (is_iterator1)
            ++i_;
          else
            ++iter2;
          return *this;
        }
        self_type operator++(int) { self_type tmp = *this; ++(*this); return tmp; }

        self_type operator+=(size_type offset)
        {
          if (is_iterator1)
            i_ += offset;
          else
          {
            for (size_type k=0; k<offset; ++k)
              ++iter2;  //Note: User must ensure that this is always valid...
          }
          return *this;
        }

        bool operator==(self_type const & other) const
        {
          if (is_iterator1)
            return (i_ == other.i_);
          return (iter2 == other.iter2);
        }
        bool operator!=(self_type const & other) const { return !(*this == other); }

        size_type index1() const { return i_; }
        size_type index2() const
        {
          if (is_iterator1)
            return 0;
          else
            return iter2->first;
        }

        sparse_matrix_adapted_iterator<SCALARTYPE, SizeType, !is_iterator1> begin() const
        {
          return sparse_matrix_adapted_iterator<SCALARTYPE, SizeType, !is_iterator1>(mat_, static_cast<int>(i_), 0);
        }
        sparse_matrix_adapted_iterator<SCALARTYPE, SizeType, !is_iterator1> end() const
        {
          int end_ = static_cast<int>(mat_[i_].size());
          if (end_ > 0)
            end_ = mat_[i_].rbegin()->first;
          return sparse_matrix_adapted_iterator<SCALARTYPE, SizeType, !is_iterator1>(mat_, static_cast<int>(i_), end_ + 1);
        }

      private:
        std::vector<std::map<SizeType, SCALARTYPE> > & mat_;
        typename std::map<SizeType, SCALARTYPE>::iterator iter2;
        size_type i_;
        size_type j_;
    };



    /** @brief Adapts a non-const sparse matrix type made up from std::vector<std::map<SizeType, SCALARTYPE> > to basic ublas-compatibility.
    *
    *  @tparam SCALARTYPE   either float or double
    */
    template <typename SCALARTYPE, typename SizeType = unsigned int>
    class sparse_matrix_adapter : public const_sparse_matrix_adapter<SCALARTYPE, SizeType>
    {
        typedef const_sparse_matrix_adapter<SCALARTYPE, SizeType>   BaseType;
      public:
        typedef sparse_matrix_adapted_iterator<SCALARTYPE, SizeType, true>      iterator1;
        typedef sparse_matrix_adapted_iterator<SCALARTYPE, SizeType, false>     iterator2;
        typedef const_sparse_matrix_adapted_iterator<SCALARTYPE, SizeType, true, true>      const_iterator1;
        typedef const_sparse_matrix_adapted_iterator<SCALARTYPE, SizeType, false, true>     const_iterator2;
        typedef SizeType                                              size_type;

        sparse_matrix_adapter(std::vector<std::map<SizeType, SCALARTYPE> > & mat)
         : BaseType(mat), mat_(mat), size1_(mat_.size()), size2_(mat_.size()) {}

        sparse_matrix_adapter(std::vector<std::map<SizeType, SCALARTYPE> > & mat,
                              vcl_size_t num_rows,
                              vcl_size_t num_cols)
         : BaseType(mat, num_rows, num_cols), mat_(mat), size1_(static_cast<size_type>(num_rows)), size2_(static_cast<size_type>(num_cols)) {}

        iterator1 begin1() { return iterator1(mat_, 0, 0); }
        iterator1 end1() { return iterator1(mat_, static_cast<int>(mat_.size()), static_cast<int>(mat_.back().size())); }

        const_iterator1 begin1() const { return const_iterator1(mat_, 0, 0); }
        const_iterator1 end1() const   { return const_iterator1(mat_, size1(), size2()); }

        iterator2 begin2() { return iterator2(mat_, 0, 0); }
        iterator2 end2() { return iterator2(mat_, mat_.size(), mat_.back().size()); }

        const_iterator2 begin2() const { return const_iterator2(mat_, 0, 0); }
        const_iterator2 end2() const   { return const_iterator2(mat_, size1(), size2()); }

        SCALARTYPE & operator()(vcl_size_t i, vcl_size_t j) { return mat_[i][static_cast<size_type>(j)]; }

        void resize(vcl_size_t i, vcl_size_t j, bool preserve = true)
        {
          if (i>0)
            mat_.resize(i);
          if (!preserve)
            clear();

          size1_ = static_cast<size_type>(i);
          size2_ = static_cast<size_type>(j);
        }

        void clear()
        {
          for (size_type i=0; i<mat_.size(); ++i)
            mat_[i].clear();
        }

        size_type size1() { return size1_; }
        size_type size1() const { return size1_; } //Note: Due to name hiding it is not sufficient to have it in the base class

        //assume a square matrix
        size_type size2() { return size2_; }
        size_type size2() const { return size2_; } //Note: Due to name hiding it is not sufficient to have it in the base class

      private:
        std::vector<std::map<SizeType, SCALARTYPE> > & mat_;
        size_type size1_;
        size_type size2_;
    };


  }
}
#endif
