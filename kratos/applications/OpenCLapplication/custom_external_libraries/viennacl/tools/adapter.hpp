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
======================================================================= */

#ifndef _VIENNACL_TOOLS_ADAPTER_HPP_
#define _VIENNACL_TOOLS_ADAPTER_HPP_

/** @file adapter.hpp
    @brief Adapter classes for sparse matrices made of the STL type std::vector<std::map<unsigned int, SCALARTYPE> >
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
    
    /** @brief A const iterator for sparse matrices of type std::vector<std::map<unsigned int, SCALARTYPE> >
    *  
    *  The iterator behaves like ublas iterators. Attention: Iteration along first columns and then rows via .begin() is untested!
    *
    *  @tparam SCALARTYPE     either float or double
    *  @tparam is_iterator1   if true, this iterator iterates along increasing row indices, otherwise along increasing column indices
    *  @tparam increment      if +1, this is a forward iterator, if -1 we have a reverse iterator
    */
    template <typename SCALARTYPE, bool is_iterator1, bool is_forward>
    class const_sparse_matrix_adapted_iterator
    {
      typedef const_sparse_matrix_adapted_iterator<SCALARTYPE, is_iterator1, is_forward>    self_type;
      
      public:
        typedef self_type     iterator1;
        typedef self_type     iterator2;
        
        const_sparse_matrix_adapted_iterator(std::vector<std::map<unsigned int, SCALARTYPE> > const & mat, int i, int j)
         : _mat(mat), _i(i), _j(j)
        {
          if (i < 0) //reverse iterator end
          {
           // iter2 = _mat[0].rend();  //reverse iterator end
          }
          else  //_i is valid
          {
            if (j < 0)
            {
              //iter2 = _mat[i].rend();
            }
            else //_j is valid
            {
              int mat_size = _mat.size();
              if (_i < mat_size && _j < mat_size )
              {
                //TODO: Start at entry j, not at the begin
                iter2 = _mat[i].begin();
              }
              else if (_i < mat_size  && _j >= mat_size )
                iter2 = _mat[i].end();
              else //i is out of range -> end iterator requested
                iter2 = _mat[_mat.size() - 1].end(); //forward iterator end
            }
          }
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
          {
            if (is_forward)
              ++_i;
            else
              --_i;
          }
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
        
        int index1() const { return _i; }
        int index2() const
        { 
          if (is_iterator1)
            return 0;
          else
            return iter2->first;
        }
        
        const_sparse_matrix_adapted_iterator<SCALARTYPE, !is_iterator1, true> begin() const
        {
          return const_sparse_matrix_adapted_iterator<SCALARTYPE, !is_iterator1, true>(_mat, _i, iter2->first);
        }
        const_sparse_matrix_adapted_iterator<SCALARTYPE, !is_iterator1, true> end() const
        {
          return const_sparse_matrix_adapted_iterator<SCALARTYPE, !is_iterator1, true>(_mat, _i, static_cast<unsigned int>(_mat.size()));
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
        typedef const_sparse_matrix_adapted_iterator<SCALARTYPE, true, true>      const_iterator1;
        typedef const_sparse_matrix_adapted_iterator<SCALARTYPE, false, true>     const_iterator2;

        typedef const_sparse_matrix_adapted_iterator<SCALARTYPE, true, false>   const_reverse_iterator1;
        typedef SCALARTYPE    value_type;
        
        const_sparse_matrix_adapter(std::vector<std::map<unsigned int, SCALARTYPE> > const & mat) 
         : _mat(mat) {};
        
        unsigned int size1() const { return static_cast<unsigned int>(_mat.size()); }
        unsigned int size2() const { return static_cast<unsigned int>(_mat.size()); }  //we allow only square matrices
        
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
        
        sparse_matrix_adapted_iterator(std::vector<std::map<unsigned int, SCALARTYPE> > & mat, int i, int j)
         : _mat(mat), _i(i), _j(j)
        {
          if (i < 0) //reverse iterator end
          {
            //iter2 = _mat[0].rend();  //reverse iterator end
          }
          else  //_i is valid
          {
            if (j < 0)
            {
              //iter2 = _mat[i].rend();
            }
            else //_j is valid
            {
              if (_i < _mat.size() && _j < _mat.size() )
              {
                //TODO: Start at entry j, not at the begin
                iter2 = _mat[i].begin();
              }
              else if (_i < _mat.size() && _j >= _mat.size())
                iter2 = _mat[i].end();
              else //i is out of range -> end iterator requested
                iter2 = _mat[_mat.size() - 1].end(); //forward iterator end
            }
          }
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
          return sparse_matrix_adapted_iterator<SCALARTYPE, !is_iterator1>(_mat, _i, static_cast<int>(_mat.size()));
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
        
        size_t size1() { return _mat.size(); }
        
        //assume a square matrix
        size_t size2() { return _mat.size(); }
        
      private:
        std::vector<std::map<unsigned int, SCALARTYPE> > & _mat;
    };
    

  }
}
#endif
