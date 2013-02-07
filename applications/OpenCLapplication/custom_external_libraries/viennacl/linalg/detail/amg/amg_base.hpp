#ifndef VIENNACL_LINALG_DETAIL_AMG_AMG_BASE_HPP_
#define VIENNACL_LINALG_DETAIL_AMG_AMG_BASE_HPP_

/* =========================================================================
   Copyright (c) 2010-2012, Institute for Microelectronics,
                            Institute for Analysis and Scientific Computing,
                            TU Wien.

                            -----------------
                  ViennaCL - The Vienna Computing Library
                            -----------------

   Project Head:    Karl Rupp                   rupp@iue.tuwien.ac.at
               
   (A list of authors and contributors can be found in the PDF manual)

   License:         MIT (X11), see file LICENSE in the base directory
============================================================================= */

/** @file amg_base.hpp
    @brief Helper classes and functions for the AMG preconditioner. Experimental.
    
    AMG code contributed by Markus Wagner
*/

#include <boost/numeric/ublas/vector.hpp>
#include <cmath>
#include <set>
#include <list>
#include <algorithm>

#include <map>
#ifdef _OPENMP
#include <omp.h>
#endif

#include "amg_debug.hpp"

#define VIENNACL_AMG_COARSE_RS 1
#define VIENNACL_AMG_COARSE_ONEPASS 2
#define VIENNACL_AMG_COARSE_RS0 3
#define VIENNACL_AMG_COARSE_RS3 4
#define VIENNACL_AMG_COARSE_AG 5
#define VIENNACL_AMG_INTERPOL_DIRECT 1
#define VIENNACL_AMG_INTERPOL_CLASSIC 2
#define VIENNACL_AMG_INTERPOL_AG 3
#define VIENNACL_AMG_INTERPOL_SA 4

namespace viennacl
{
  namespace linalg
  {
    namespace detail
    {
      namespace amg
      {
        /** @brief A tag for algebraic multigrid (AMG). Used to transport information from the user to the implementation.
        */
        class amg_tag
        {
          public:
            /** @brief The constructor.
            * @param coarse    Coarsening Routine (Default: VIENNACL_AMG_COARSE_CLASSIC)
            * @param interpol  Interpolation routine (Default: VIENNACL_AMG_INTERPOL_DIRECT)
            * @param threshold    Strength of dependence threshold for the coarsening process (Default: 0.25)
            * @param interpolweight  Interpolation parameter for SA interpolation and truncation parameter for direct+classical interpolation
            * @param jacobiweight  Weight of the weighted Jacobi smoother iteration step (Default: 1 = Regular Jacobi smoother)
            * @param presmooth    Number of presmoothing operations on every level (Default: 1)
            * @param postsmooth   Number of postsmoothing operations on every level (Default: 1)
            * @param coarselevels  Number of coarse levels that are constructed
            *      (Default: 0 = Optimize coarse levels for direct solver such that coarsest level has a maximum of COARSE_LIMIT points)
            *      (Note: Coarsening stops when number of coarse points = 0 and overwrites the parameter with actual number of coarse levels)
            */
            amg_tag(unsigned int coarse = 1,
                    unsigned int interpol = 1,
                    double threshold = 0.25,
                    double interpolweight = 0.2,
                    double jacobiweight = 1,
                    unsigned int presmooth = 1,
                    unsigned int postsmooth = 1,
                    unsigned int coarselevels = 0)
            : _coarse(coarse), _interpol(interpol),
              _threshold(threshold), _interpolweight(interpolweight), _jacobiweight(jacobiweight), 
              _presmooth(presmooth), _postsmooth(postsmooth), _coarselevels(coarselevels) {}; 

            // Getter-/Setter-Functions
            void set_coarse(unsigned int coarse) { if (coarse > 0) _coarse = coarse; }
            unsigned int get_coarse() const { return _coarse; }
            
            void set_interpol(unsigned int interpol) { if (interpol > 0) _interpol = interpol; }
            unsigned int get_interpol() const { return _interpol; }
            
            void set_threshold(double threshold) { if (threshold > 0 && threshold <= 1) _threshold = threshold; }
            double get_threshold() const{ return _threshold; }
            
            void set_as(double jacobiweight) { if (jacobiweight > 0 && jacobiweight <= 2) _jacobiweight = jacobiweight; }
            double get_interpolweight() const { return _interpolweight; }
            
            void set_interpolweight(double interpolweight) { if (interpolweight > 0 && interpolweight <= 2) _interpolweight = interpolweight; }
            double get_jacobiweight() const { return _jacobiweight; }
            
            void set_presmooth(int presmooth) { if (presmooth >= 0) _presmooth = presmooth; }
            unsigned int get_presmooth() const { return _presmooth; }
            
            void set_postsmooth(int postsmooth) { if (postsmooth >= 0) _postsmooth = postsmooth; }
            unsigned int get_postsmooth() const { return _postsmooth; }
            
            void set_coarselevels(int coarselevels)  { if (coarselevels >= 0) _coarselevels = coarselevels; }
            unsigned int get_coarselevels() const { return _coarselevels; }

          private:
            unsigned int _coarse, _interpol;
            double _threshold, _interpolweight, _jacobiweight;
            unsigned int _presmooth, _postsmooth, _coarselevels;
        };
        
        /** @brief A class for a scalar that can be written to the sparse matrix or sparse vector datatypes.
        *  @brief Values are only written to those datatypes if non-zero to optimize memory usage and performance.
        *  @brief Needed for the []- and ()-operators.
        */
        template <typename InternalType, typename IteratorType, typename ScalarType>
        class amg_nonzero_scalar
        {
          private:
            InternalType *_m;
            IteratorType _iter;
            unsigned int _i,_j;
            ScalarType _s;

          public:
            amg_nonzero_scalar();
            
            /** @brief The constructor.
            *  @param m    Pointer to the sparse vector/matrix the scalar will be written to
            *  @param iter    Iterator pointing to the respective element in the vector/matrix if available
            *  @param i    Row index scalar will be written to
            *  @param j    Col index scalar will be written to
            *  @param s    Value of the scalar (usually used as dummy here as it will be set by the assignment operator)
            */
            amg_nonzero_scalar(InternalType *m,
                              IteratorType & iter,
                              unsigned int i,
                              unsigned int j,
                              ScalarType s = 0): _m(m), _iter(iter), _i(i), _j(j), _s(s) {}
            
            /** @brief Assignment operator. Writes value into matrix at the given position.
            *  @param value  Value that will be written
            */
            ScalarType operator = (const ScalarType value)
            {
              _s = value;
              // Only write if scalar is nonzero
              if (_s == 0) return _s;
              // Write to _m using iterator _iter or indices (_i,_j)
              _m->addscalar (_iter,_i,_j,_s);
              return _s;
            }
            
            /** @brief Addition operator. Adds a constant.
            *  @param value  Value that will be written
            */
            ScalarType operator += (const ScalarType value)
            {
              // If zero is added, then no change necessary
              if (value == 0)
                return _s;
              
              _s += value;
              // Remove entry if resulting scalar is zero
              if (_s == 0)
              {
                _m->removescalar(_iter,_i);
                return _s;
              }
              //Write to _m using iterator _iter or indices (_i,_j)
              _m->addscalar (_iter,_i,_j,_s);
              return _s;
            }
            ScalarType operator ++ (int)
            {
              _s++;
              if (_s == 0)
                _m->removescalar(_iter,_i);
              _m->addscalar (_iter,_i,_j,_s);
              return _s;
            }
            ScalarType operator ++ ()
            {
              _s++;
              if (_s == 0)
                _m->removescalar(_iter,_i);
              _m->addscalar (_iter,_i,_j,_s);
              return _s;
            }
            operator ScalarType (void) { return _s;  }
        };
    
        /** @brief Defines an iterator for the sparse vector type.
        */
        template <typename InternalType>
        class amg_sparsevector_iterator
        {
          private:
            typedef amg_sparsevector_iterator<InternalType> self_type;  
            typedef typename InternalType::mapped_type ScalarType;
            
            InternalType & internal_vec;
            typename InternalType::iterator iter;
              
          public:
            
            /** @brief The constructor.
            *  @param vec    Internal sparse vector
            *  @param begin  Whether the iterator starts at the beginning or end of vec
            */
            amg_sparsevector_iterator(InternalType & vec, bool begin=true): internal_vec(vec)
            {
              if (begin)
                iter = internal_vec.begin();
              else
                iter = internal_vec.end();
            }
            
            bool operator == (self_type other)
            {
              if (iter == other.iter)
                return true;
              else
                return false;
            }
            bool operator != (self_type other)
            {
              if (iter != other.iter)
                return true;
              else
                return false;
            }
            
            self_type & operator ++ () const { iter++; return *this; }
            self_type & operator ++ () { iter++; return *this; }
            self_type & operator -- () const { iter--; return *this; }
            self_type & operator -- () { iter--; return *this; }  
            ScalarType & operator * () const { return (*iter).second; }
            ScalarType & operator * () { return (*iter).second; }
            unsigned int index() const { return (*iter).first; }
            unsigned int index() { return (*iter).first; }
        };
    
        /** @brief A class for the sparse vector type.
        */
        template <typename ScalarType>
        class amg_sparsevector
        {
          public:
            typedef ScalarType value_type;
      
          private:
            // A map is used internally which saves all non-zero elements with pairs of (index,value)
            typedef std::map<unsigned int,ScalarType> InternalType;
            typedef amg_sparsevector<ScalarType> self_type;
            typedef amg_nonzero_scalar<self_type,typename InternalType::iterator,ScalarType> NonzeroScalarType;
              
            // Size is only a dummy variable. Not needed for internal map structure but for compatible vector interface.
            unsigned int _size;
            InternalType internal_vector;
      
          public:
            typedef amg_sparsevector_iterator<InternalType> iterator;
            typedef typename InternalType::const_iterator const_iterator;
      
          public:
            /** @brief The constructor.
            *  @param size    Size of the vector
            */
            amg_sparsevector(unsigned int size = 0): _size(size)
            {
              internal_vector = InternalType();
            }
            
            void resize(unsigned int size) { _size = size; }
            unsigned int size() const { return _size;}
            
            // Returns number of non-zero entries in vector equal to the size of the underlying map.
            unsigned int internal_size() const { return internal_vector.size(); }
            // Delete underlying map.
            void clear() { internal_vector.clear();  }
            // Remove entry at position i.
            void remove(unsigned int i) { internal_vector.erase(i); }
            
            // Add s to the entry at position i
            void add (unsigned int i, ScalarType s)
            {
              typename InternalType::iterator iter = internal_vector.find(i);
              // If there is no entry at position i, add new entry at that position
              if (iter == internal_vector.end())
                addscalar(iter,i,i,s);
              else
              {
                s += (*iter).second;
                // If new value is zero, then erase the entry, otherwise write new value
                if (s != 0)
                  (*iter).second = s;
                else
                  internal_vector.erase(iter);
              }
            }
            
            // Write to the map. Is called from non-zero scalar type.
            template <typename IteratorType>
            void addscalar(IteratorType & iter, unsigned int i, unsigned int j, ScalarType s)
            {
              // Don't write if value is zero
              if (s == 0)
                return;
              
              // If entry is already present, overwrite value, otherwise make new entry
              if (iter != internal_vector.end())  
                (*iter).second = s;
              else
                internal_vector[i] = s;
            }
            
            // Remove value from the map. Is called from non-zero scalar type.
            template <typename IteratorType>
            void removescalar(IteratorType & iter, unsigned int i) { internal_vector.erase(iter); }   
            
            // Bracket operator. Returns non-zero scalar type with actual values of the respective entry which calls addscalar/removescalar after value is altered.
            NonzeroScalarType operator [] (unsigned int i)
            {
              typename InternalType::iterator it = internal_vector.find(i);
              // If value is already present then build non-zero scalar with actual value, otherwise 0.
              if (it != internal_vector.end())
                return NonzeroScalarType (this,it,i,i,(*it).second);
              else
                return NonzeroScalarType (this,it,i,i,0);  
            }
            
            // Use internal data structure directly for read-only access. No need to use non-zero scalar as no write access possible.
            ScalarType operator [] (unsigned int i) const
            {
              const_iterator it = internal_vector.find(i);
              
              if (it != internal_vector.end())
                return (*it).second;
              else
                return 0;
            }
            
            // Iterator functions.
            iterator begin() { return iterator(internal_vector); }
            const_iterator begin() const { return internal_vector.begin(); }
            iterator end() { return iterator(internal_vector,false); }
            const_iterator end() const { return internal_vector.end(); }
            
            // checks whether value at index i is nonzero. More efficient than doing [] == 0.
            bool isnonzero(unsigned int i) const { return internal_vector.find(i) != internal_vector.end();  }
            
            // Copies data into a ublas vector type.
            operator boost::numeric::ublas::vector<ScalarType> (void)
            {
              boost::numeric::ublas::vector<ScalarType> vec (_size);    
              for (iterator iter = begin(); iter != end(); ++iter)
                vec [iter.index()] = *iter;        
              return vec;
            } 
         };
    
        /** @brief A class for the sparse matrix type.
        *  Uses vector of maps as data structure for higher performance and lower memory usage.
        *  Uses similar interface as ublas::compressed_matrix.
        *  Can deal with transposed of matrix internally: Creation, Storage, Iterators, etc.
        */
        template <typename ScalarType>
        class amg_sparsematrix
        {
          public:
            typedef ScalarType value_type;
          private:
            typedef std::map<unsigned int,ScalarType> RowType;
            typedef std::vector<RowType> InternalType;
            typedef amg_sparsematrix<ScalarType> self_type;
            
            // Adapter is used for certain functionality, especially iterators.
            typedef typename viennacl::tools::sparse_matrix_adapter<ScalarType> AdapterType;
            typedef typename viennacl::tools::const_sparse_matrix_adapter<ScalarType> ConstAdapterType;
            
            // Non-zero scalar is used to write to the matrix.
            typedef amg_nonzero_scalar<self_type,typename RowType::iterator,ScalarType> NonzeroScalarType;

            // Holds matrix coefficients.
            InternalType internal_mat;
            // Holds matrix coefficient of transposed matrix if built. 
            // Note: Only internal_mat is written using operators and methods while internal_mat_trans is built from internal_mat using do_trans().
            InternalType internal_mat_trans;
            // Saves sizes.
            size_t s1, s2;
            
            // True if the transposed of the matrix is used (for calculations, iteration, etc.).
            bool transposed_mode;
            // True if the transposed is already built (saved in internal_mat_trans) and also up to date (no changes to internal_mat).
            bool transposed;
            
          public:          
            typedef typename AdapterType::iterator1 iterator1;
            typedef typename AdapterType::iterator2 iterator2;
            typedef typename ConstAdapterType::const_iterator1 const_iterator1;
            typedef typename ConstAdapterType::const_iterator2 const_iterator2;
            
            /** @brief Standard constructor. */
            amg_sparsematrix ()
            {
              transposed_mode = false;
              transposed = false;
            }
            
            /** @brief Constructor. Builds matrix of size (i,j).
              * @param i  Size of first dimension
              * @param j  Size of second dimension
              */
            amg_sparsematrix (unsigned int i, unsigned int j)
            {
              AdapterType a (internal_mat, i, j);
              a.resize (i,j,false);
              AdapterType a_trans (internal_mat_trans, j, i);
              a_trans.resize (j,i,false);
              s1 = i;
              s2 = j;
              a.clear();
              a_trans.clear();
              transposed_mode = false;
              transposed = false;
            }
            
            /** @brief Constructor. Builds matrix via std::vector<std::map> by copying memory
            * (Only necessary feature of this other matrix type is to have const iterators)
            * @param mat  Vector of maps
            */
            amg_sparsematrix (std::vector<std::map<unsigned int, ScalarType> > const & mat)
            {  
              AdapterType a (internal_mat, mat.size(), mat.size());
              AdapterType a_trans (internal_mat_trans, mat.size(), mat.size());
              a.resize(mat.size(), mat.size());
              a_trans.resize(mat.size(), mat.size());
              
              internal_mat = mat;  
              s1 = s2 = mat.size();
              
              transposed_mode = false;
              transposed = false;
            }
            
            /** @brief Constructor. Builds matrix via another matrix type.
              * (Only necessary feature of this other matrix type is to have const iterators)
              * @param mat  Matrix
              */
            template <typename MatrixType>
            amg_sparsematrix (MatrixType const & mat)
            {  
              AdapterType a (internal_mat, mat.size1(), mat.size2());
              AdapterType a_trans (internal_mat_trans, mat.size2(), mat.size1());
              a.resize(mat.size1(), mat.size2());
              a_trans.resize (mat.size2(), mat.size1());
              s1 = mat.size1();
              s2 = mat.size2();
              a.clear();
              a_trans.clear();
              
              for (typename MatrixType::const_iterator1 row_iter = mat.begin1(); row_iter != mat.end1(); ++row_iter)
              {
                for (typename MatrixType::const_iterator2 col_iter = row_iter.begin(); col_iter != row_iter.end(); ++col_iter)
                {
                  if (*col_iter != 0)
                  {
                    unsigned int x = col_iter.index1();
                    unsigned int y = col_iter.index2();
                    a (x,y) = *col_iter;
                    a_trans (y,x) = *col_iter;
                  }
                }
              }
              transposed_mode = false;
              transposed = true;
            }
                  
            // Build transposed of the current matrix.
            void do_trans()
            {
              // Do it only once if called in a parallel section
            #ifdef _OPENMP
              #pragma omp critical
            #endif
              { 
                // Only build transposed if it is not built or not up to date
                if (!transposed)
                {
                  // Mode has to be set to standard mode temporarily
                  bool save_mode = transposed_mode;
                  transposed_mode = false;
                  
                  for (iterator1 row_iter = begin1(); row_iter != end1(); ++row_iter)
                for (iterator2 col_iter = row_iter.begin(); col_iter != row_iter.end(); ++col_iter)
                  internal_mat_trans[col_iter.index2()][col_iter.index1()] = *col_iter;
                
                  transposed_mode = save_mode;
                  transposed = true;
                }
              }
            } //do_trans()
            
            // Set transposed mode (true=transposed, false=regular)
            void set_trans(bool mode)
            {
              transposed_mode = mode;
              if (mode)
                do_trans();
            }
            
            bool get_trans() const { return transposed_mode; }     
                  
            // Checks whether coefficient (i,j) is non-zero. More efficient than using (i,j) == 0.
            bool isnonzero (unsigned int i, unsigned int j) const
            {
              if (!transposed_mode)
              {
                if (internal_mat[i].find(j) != internal_mat[i].end())
                  return true;
                else
                  return false;
              }
              else
              {
                if (internal_mat[j].find(i) != internal_mat[j].end())
                  return true;
                else
                  return false;
              }
            } //isnonzero()
                
            // Add s to value at (i,j)
            void add (unsigned int i, unsigned int j, ScalarType s)
            {
              // If zero is added then do nothing.
              if (s == 0)
                return;
              
              typename RowType::iterator col_iter = internal_mat[i].find(j);
              // If there is no entry at position (i,j), then make new entry.
              if (col_iter == internal_mat[i].end())
                addscalar(col_iter,i,j,s);
              else
              {
                s += (*col_iter).second;
                // Update value and erase entry if value is zero.
                if (s != 0)
                  (*col_iter).second = s;
                else
                  internal_mat[i].erase(col_iter);
              }
              transposed = false;
            } //add()
            
            // Write to the internal data structure. Is called from non-zero scalar type.
            template <typename IteratorType>
            void addscalar(IteratorType & iter, unsigned int i, unsigned int j, ScalarType s)
            {
              // Don't write if value is zero
              if (s == 0)
                return;
              
              if (iter != internal_mat[i].end())  
                (*iter).second = s;
              else
                internal_mat[i][j] = s;
              
              transposed = false;
            }
            
            // Remove entry from internal data structure. Is called from non-zero scalar type.
            template <typename IteratorType>
            void removescalar(IteratorType & iter, unsigned int i)
            {
              internal_mat[i].erase(iter);
              transposed = false;
            }   
            
            // Return non-zero scalar at position (i,j). Value is written to the non-zero scalar and updated via addscalar()/removescalar().
            NonzeroScalarType operator()(unsigned int i, unsigned int j)
            {
              typename RowType::iterator iter;
              
              if (!transposed_mode)
              {
                iter = internal_mat[i].find(j);
                if (iter != internal_mat[i].end())
                  return NonzeroScalarType (this,iter,i,j,(*iter).second);
                else
                  return NonzeroScalarType (this,iter,i,j,0);
              }
              else
              {
                iter = internal_mat[j].find(i);
                if (iter != internal_mat[j].end())
                  return NonzeroScalarType (this,iter,j,i,(*iter).second);
                else
                  return NonzeroScalarType (this,iter,j,i,0);
              }
            }
            
            // For read-only access return the actual value directly. Non-zero datatype not needed as no write access possible.
            ScalarType operator()(unsigned int i, unsigned int j) const
            {
              typename RowType::const_iterator iter;
              
              if (!transposed_mode)
              {
                iter = internal_mat[i].find(j);
                if (iter != internal_mat[i].end())
                  return (*iter).second;
                else
                  return 0;
              }
              else
              {
                iter = internal_mat[j].find(i);
                if (iter != internal_mat[j].end())
                  return (*iter).second;
                else
                  return 0;
              }
            }
              
            void resize(unsigned int i, unsigned int j, bool preserve = true)
            {
              AdapterType a (internal_mat);
              a.resize(i,j,preserve);
              AdapterType a_trans (internal_mat_trans);
              a_trans.resize(j,i,preserve);
              s1 = i;
              s2 = j;
            }
            
            void clear() 
            {
              AdapterType a (internal_mat, s1, s2);
              a.clear();
              AdapterType a_trans (internal_mat_trans, s2, s1);
              a_trans.clear();
              transposed = true;
            }

            size_t size1()
            {
              if (!transposed_mode)
                return s1;
              else
                return s2;
            }
            
            size_t size1() const
            {
              if (!transposed_mode)
                return s1;
              else
                return s2;
            }
            
            
            size_t size2()
            {
              if (!transposed_mode)
                return s2;
              else
                return s1;
            }
            size_t size2() const
            {
              if (!transposed_mode)
                return s2;
              else
                return s1;
            }
            
            iterator1 begin1(bool trans = false)
            {
              if (!trans && !transposed_mode)
              {
                AdapterType a (internal_mat, s1, s2);
                return a.begin1();
              }
              else
              {
                do_trans();
                AdapterType a_trans (internal_mat_trans, s2, s1);
                return a_trans.begin1();
              }
            }
            
            iterator1 end1(bool trans = false)
            {
              if (!trans && !transposed_mode)
              {
                AdapterType a (internal_mat, s1, s2);
                return a.end1();
              }
              else
              {
                //do_trans();
                AdapterType a_trans (internal_mat_trans, s2, s1);
                return a_trans.end1();
              }
            }
            
            iterator2 begin2(bool trans = false)
            {
              if (!trans && !transposed_mode)
              {
                AdapterType a (internal_mat, s1, s2);
                return a.begin2();
              }
              else
              {
                do_trans();
                AdapterType a_trans (internal_mat_trans, s2, s1);
                return a_trans.begin2();
              }
            }
            
            iterator2 end2(bool trans = false)
            {
              if (!trans && !transposed_mode)
              {
                AdapterType a (internal_mat, s1, s2);
                return a.end2();
              }
              else
              {
                //do_trans();
                AdapterType a_trans (internal_mat_trans, s2, s1);
                return a_trans.end2();
              }
            }
            
            const_iterator1 begin1() const
            {
              // Const_iterator of transposed can only be used if transposed matrix is already built and up to date.
              assert((!transposed_mode || (transposed_mode && transposed)) && "Error: Cannot build const_iterator when transposed has not been built yet!");
                    ConstAdapterType a_const (internal_mat, s1, s2);
              return a_const.begin1();
            }
            
            const_iterator1 end1(bool trans = false) const
            {
              assert((!transposed_mode || (transposed_mode && transposed)) && "Error: Cannot build const_iterator when transposed has not been built yet!");
              ConstAdapterType a_const (internal_mat, trans ? s2 : s1, trans ? s1 : s2);
              return a_const.end1();
            }
            
            const_iterator2 begin2(bool trans = false) const
            {
              assert((!transposed_mode || (transposed_mode && transposed)) && "Error: Cannot build const_iterator when transposed has not been built yet!");
              ConstAdapterType a_const (internal_mat, trans ? s2 : s1, trans ? s1 : s2);
              return a_const.begin2();
            }
            
            const_iterator2 end2(bool trans = false) const
            {
              assert((!transposed_mode || (transposed_mode && transposed)) && "Error: Cannot build const_iterator when transposed has not been built yet!");
              ConstAdapterType a_const (internal_mat, trans ? s2 : s1, trans ? s1 : s2);
              return a_const.end2();
            }
            
            // Returns pointer to the internal data structure. Improves performance of copy operation to GPU.
            std::vector<std::map<unsigned int, ScalarType> > * get_internal_pointer()
            {    
              if (!transposed_mode)
                return &internal_mat;
              
              if (!transposed)
                do_trans();
              return &internal_mat_trans;
            }
            operator boost::numeric::ublas::compressed_matrix<ScalarType> (void)
            {
              boost::numeric::ublas::compressed_matrix<ScalarType> mat;
              mat.resize(size1(),size2(),false);
              mat.clear();
              
              for (iterator1 row_iter = begin1(); row_iter != end1(); ++row_iter)
                  for (iterator2 col_iter = row_iter.begin(); col_iter != row_iter.end(); ++col_iter)
                    mat (col_iter.index1(), col_iter.index2()) = *col_iter;
                  
              return mat;
            } 
            operator boost::numeric::ublas::matrix<ScalarType> (void)
            {
              boost::numeric::ublas::matrix<ScalarType> mat;
              mat.resize(size1(),size2(),false);
              mat.clear();
              
              for (iterator1 row_iter = begin1(); row_iter != end1(); ++row_iter)
                  for (iterator2 col_iter = row_iter.begin(); col_iter != row_iter.end(); ++col_iter)
                    mat (col_iter.index1(), col_iter.index2()) = *col_iter;
                  
              return mat;
            }
        };
          
        /** @brief A class for the AMG points.
        *   Saves point index and influence measure
        *  Holds information whether point is undecided, C or F point.
        *  Holds lists of points that are influenced by or influencing this point
        */
        class amg_point
        {
          private:
            typedef amg_sparsevector<amg_point*> ListType;
            
            unsigned int _index;
            unsigned int _influence;
            // Determines whether point is undecided.
            bool _undecided;
            // Determines wheter point is C point (true) or F point (false).
            bool _cpoint;
            unsigned int _coarse_index;
            // Index offset of parallel coarsening. In that case a point acts as if it had an index of _index-_offset and treats other points as if they had an index of index+_offset
            unsigned int _offset;
            // Aggregate the point belongs to.
            unsigned int _aggregate;
            
            // Holds all points influencing this point.
            ListType influencing_points;
            // Holds all points that are influenced by this point.
            ListType influenced_points;
      
          public:
            typedef ListType::iterator iterator;
            typedef ListType::const_iterator const_iterator;
            
            /** @brief The constructor.
            */
            amg_point (unsigned int index, unsigned int size): _index(index), _influence(0), _undecided(true), _cpoint(false), _coarse_index(0), _offset(0), _aggregate(0)
            {
              influencing_points = ListType(size);
              influenced_points = ListType(size);
            }
            
            void set_offset(unsigned int offset) { _offset = offset; }
            unsigned int get_offset() { return _offset; }
            void set_index(unsigned int index) { _index = index+_offset; }
            unsigned int get_index() const { return _index-_offset;  }
            unsigned int get_influence() const { return _influence;  }
            void set_aggregate(unsigned int aggregate) { _aggregate = aggregate; }
            unsigned int get_aggregate () { return _aggregate; }
            
            bool is_cpoint() const { return _cpoint && !_undecided;  }
            bool is_fpoint() const { return !_cpoint && !_undecided; }
            bool is_undecided() const { return _undecided; }
            
            // Returns number of influencing points
            unsigned int number_influencing() const  { return influencing_points.internal_size(); }
            // Returns true if *point is influencing this point
            bool is_influencing(amg_point* point) const { return influencing_points.isnonzero(point->get_index()+_offset); }
            // Add *point to influencing points
            void add_influencing_point(amg_point* point) { influencing_points[point->get_index()+_offset] = point;  }
            // Add *point to influenced points
            void add_influenced_point(amg_point* point) { influenced_points[point->get_index()+_offset] = point; }
            
            // Clear influencing points
            void clear_influencing() { influencing_points.clear(); }
            // Clear influenced points
            void clear_influenced() {influenced_points.clear(); }
            
            
            unsigned int get_coarse_index() const { return _coarse_index; }
            void set_coarse_index(unsigned int index) { _coarse_index = index; }
            
            // Calculates the initial influence measure equal to the number of influenced points.
            void calc_influence() { _influence = influenced_points.internal_size();  }
            
            // Add to influence measure.
            unsigned int add_influence(int add)
            {
              _influence += add;
              return _influence;
            }
            // Make this point C point. Only call via amg_pointvector.
            void make_cpoint() 
            { 
              _undecided = false;
              _cpoint = true; 
              _influence = 0;
            }
            // Make this point F point. Only call via amg_pointvector.
            void make_fpoint()
            {
              _undecided = false;
              _cpoint = false;
              _influence = 0;
            }
            // Switch point from F to C point. Only call via amg_pointvector.
            void switch_ftoc() { _cpoint = true; }  
            
            // Iterator handling for influencing and influenced points.
            iterator begin_influencing() { return influencing_points.begin(); }
            iterator end_influencing() { return influencing_points.end(); }
            const_iterator begin_influencing() const { return influencing_points.begin(); }
            const_iterator end_influencing() const { return influencing_points.end(); }
            iterator begin_influenced() { return influenced_points.begin();  }
            iterator end_influenced() { return influenced_points.end(); }
            const_iterator begin_influenced() const { return influenced_points.begin(); }
            const_iterator end_influenced() const { return influenced_points.end(); }
        };
        
        /** @brief Comparison class for the sorted set of points in amg_pointvector. Set is sorted by influence measure from lower to higher with the point-index as tie-breaker.
        */
        struct classcomp
        {
          // Function returns true if l comes before r in the ordering.
          bool operator() (amg_point* l, amg_point* r) const
          {
            // Map is sorted by influence number starting with the highest
            // If influence number is the same then lowest point index comes first
            return (l->get_influence() < r->get_influence() || (l->get_influence() == r->get_influence() && l->get_index() > r->get_index()));
          }
        };
      
        /** @brief A class for the AMG points.
        *  Holds pointers of type amg_point in a vector that can be accessed using [point-index].
        *  Additional list of pointers sorted by influence number and index to improve coarsening performance (see amg_coarse_classic_onepass() in amg_coarse.hpp)
        *  Constructs indices for C points on the coarse level, needed for interpolation.
        */
        class amg_pointvector
        {
          private:
            // Type for the sorted list
            typedef std::set<amg_point*,classcomp> ListType;
            // Type for the vector of pointers
            typedef std::vector<amg_point*> VectorType;
            VectorType pointvector;
            ListType pointlist;
            unsigned int _size;
            unsigned int c_points, f_points;
      
          public:
            typedef VectorType::iterator iterator;
            typedef VectorType::const_iterator const_iterator;
            
            /** @brief The constructor.
            *  @param size    Number of points
            */
            amg_pointvector(unsigned int size = 0): _size(size)
            {
              pointvector = VectorType(size);
              c_points = f_points = 0;
            }
            
            // Construct all the points dynamically and save pointers into vector.
            void init_points()
            {  
              for (unsigned int i=0; i<size(); ++i)
                pointvector[i] = new amg_point(i,size());
            }
            // Delete all the points.
            void delete_points()
            {  
              for (unsigned int i=0; i<size(); ++i)
                delete pointvector[i];
            }
            // Add point to the vector. Note: User has to make sure that point at point->get_index() does not exist yet, otherwise it will be overwritten!
            void add_point(amg_point *point)
            {
              pointvector[point->get_index()] = point;
              if (point->is_cpoint()) c_points++;
              else if (point->is_fpoint()) f_points++;
            }

            // Update C and F count for point *point.
            // Necessary if C and F points were constructed outside this data structure (e.g. by parallel coarsening RS0 or RS3).
            void update_cf(amg_point *point) 
            {
              if (point->is_cpoint()) c_points++;
              else if (point->is_fpoint()) f_points++;
            }
            // Clear the C and F point count.
            void clear_cf() { c_points = f_points = 0; }
            
            // Clear both point lists.
            void clear_influencelists()
            {
              for (iterator iter = pointvector.begin(); iter != pointvector.end(); ++iter)
              {
                (*iter)->clear_influencing();
                (*iter)->clear_influenced();
              }
            }
            
            amg_point* operator [] (unsigned int i) const { return pointvector[i]; }
            iterator begin() { return pointvector.begin(); }
            iterator end() { return pointvector.end(); }
            const_iterator begin() const { return pointvector.begin(); }
            const_iterator end() const { return pointvector.end(); }
            
            void resize(unsigned int size)
            {
              _size = size;
              pointvector = VectorType(size);
            }
            unsigned int size() const { return _size; }
            
            // Returns number of C points
            unsigned int get_cpoints() const { return c_points; }
            // Returns number of F points
            unsigned int get_fpoints() const { return f_points; }
            
            // Does the initial sorting of points into the list. Sorting is automatically done by the std::set data type.
            void sort()
            {
              for (iterator iter = begin(); iter != end(); ++iter)
                pointlist.insert(*iter);
            }
            // Returns the point with the highest influence measure
            amg_point* get_nextpoint()
            {
              // No points remaining? Return NULL such that coarsening will stop.
              if (pointlist.size() == 0)
                return NULL;
              // If point with highest influence measure (end of the list) has measure of zero, then no further C points can be constructed. Return NULL.
              if ((*(--pointlist.end()))->get_influence() == 0)
                return NULL;
              // Otherwise, return the point with highest influence measure located at the end of the list.
              else
                return *(--pointlist.end());
            }
            // Add "add" to influence measure for point *point in the sorted list.
            void add_influence(amg_point* point, unsigned int add)
            {
              ListType::iterator iter = pointlist.find(point);
              // If point is not in the list then stop.
              if (iter == pointlist.end()) return;
              
              // Save iterator and decrement
              ListType::iterator iter2 = iter;
              iter2--;
              
              // Point has to be erased first as changing the value does not re-order the std::set
              pointlist.erase(iter);
              point->add_influence(add);
              
              // Insert point back into the list. Using the iterator improves performance. The new position has to be at the same position or to the right of the old.
              pointlist.insert(iter2,point);
            }
            // Make *point to C point and remove from sorted list
            void make_cpoint(amg_point* point)
            {
              pointlist.erase(point);
              point->make_cpoint();
              c_points++;
            }
            // Make *point to F point and remove from sorted list
            void make_fpoint(amg_point* point)
            {
              pointlist.erase(point);
              point->make_fpoint();
              f_points++;
            }
            // Swich *point from F to C point
            void switch_ftoc(amg_point* point)
            { 
              point->switch_ftoc();
              c_points++;
              f_points--;
            }

            // Build vector of indices for C point on the coarse level.
            void build_index()
            {
              unsigned int count = 0;
              // Use simple counter for index creation.
              for (iterator iter = pointvector.begin(); iter != pointvector.end(); ++iter)
              {
                // Set index on coarse level using counter variable
                if ((*iter)->is_cpoint())
                {
                  (*iter)->set_coarse_index(count);
                  count++;
                }
              }
            }
            
            // Return information for debugging purposes
            template <typename MatrixType>
            void get_influence_matrix(MatrixType & mat) const
            {
              mat = MatrixType(size(),size());
              mat.clear();
              
              for (const_iterator row_iter = begin(); row_iter != end(); ++row_iter)
                for (amg_point::iterator col_iter = (*row_iter)->begin_influencing(); col_iter != (*row_iter)->end_influencing(); ++col_iter)
                  mat((*row_iter)->get_index(),(*col_iter)->get_index()) = true;  
            }
            template <typename VectorType>
            void get_influence(VectorType & vec) const
            {
              vec = VectorType(_size);
              vec.clear();
              
              for (const_iterator iter = begin(); iter != end(); ++iter)
                vec[(*iter)->get_index()] = (*iter)->get_influence();
            }
            template <typename VectorType>
            void get_sorting(VectorType & vec) const
            {
              vec = VectorType(pointlist.size());
              vec.clear();
              unsigned int i=0;
              
              for (ListType::const_iterator iter = pointlist.begin(); iter != pointlist.end(); ++iter)
              {
                vec[i] = (*iter)->get_index();
                i++;
              }
            }
            template <typename VectorType>
            void get_C(VectorType & vec) const
            {
              vec = VectorType(_size);
              vec.clear();
              
              for (const_iterator iter = begin(); iter != end(); ++iter)
              {
                if ((*iter)->is_cpoint())
                  vec[(*iter)->get_index()] = true;
              }
            }
            template <typename VectorType>
            void get_F(VectorType & vec) const
            {
              vec = VectorType(_size);
              vec.clear();
              
              for (const_iterator iter = begin(); iter != end(); ++iter)
              {
                if ((*iter)->is_fpoint())
                  vec[(*iter)->get_index()] = true;
              }
            }
            template <typename MatrixType>
            void get_Aggregates(MatrixType & mat) const
            {
              mat = MatrixType(_size,_size);
              mat.clear();
              
              for (const_iterator iter = begin(); iter != end(); ++iter)
              {
                if (!(*iter)->is_undecided())
                  mat((*iter)->get_aggregate(),(*iter)->get_index()) = true;
              }
            }
        };
        
        /** @brief A class for the matrix slicing for parallel coarsening schemes (RS0/RS3).
          * @brief Holds information on a per-processor basis and offers functionality to slice and join the data structures.
          */
        template <typename InternalType1, typename InternalType2>
        class amg_slicing
        {
            typedef typename InternalType1::value_type SparseMatrixType;
            typedef typename InternalType2::value_type PointVectorType;    
            
          public:
            // Data structures on a per-processor basis.
            boost::numeric::ublas::vector<InternalType1> A_slice;
            boost::numeric::ublas::vector<InternalType2> Pointvector_slice;
            // Holds the offsets showing the indices for which a new slice begins.
            boost::numeric::ublas::vector<boost::numeric::ublas::vector<unsigned int> > Offset;
            
            unsigned int _threads;
            unsigned int _levels;
            
            void init(unsigned int levels, unsigned int threads = 0)
            {
              // Either use the number of threads chosen by the user or the maximum number of threads available on the processor.
              if (threads == 0)
            #ifdef _OPENMP
                _threads = omp_get_num_procs();
            #else
              _threads = 1;
            #endif   
              else 
                _threads = threads;
              
              _levels = levels;
              
              A_slice.resize(_threads);
              Pointvector_slice.resize(_threads);
              // Offset has _threads+1 entries to also hold the total size
              Offset.resize(_threads+1);
              
              for (unsigned int i=0; i<_threads; ++i)
              {
                A_slice[i].resize(_levels);
                Pointvector_slice[i].resize(_levels);
                // Offset needs one more level for the build-up of the next offset
                Offset[i].resize(_levels+1);
              }
              Offset[_threads].resize(_levels+1);
            } //init()
            
            // Slice matrix A into as many parts as threads are used.
            void slice (unsigned int level, InternalType1 const & A, InternalType2 const & Pointvector)
            {
              // On the finest level, build a new slicing first.
              if (level == 0)
                slice_new (level, A);
              
              // On coarser levels use the same slicing as on the finest level (Points stay together on the same thread on all levels).
              // This is necessary as due to interpolation and galerkin product there only exist connections between points on the same thread on coarser levels.
              // Note: Offset is determined in amg_coarse_rs0() after fine level was built.
              slice_build (level, A, Pointvector);
            }

            // Join point data structure into Pointvector
            void join (unsigned int level, InternalType2 & Pointvector) const
            {
              typedef typename InternalType2::value_type PointVectorType;
              
              // Reset index offset of all points and update overall C and F point count
              Pointvector[level].clear_cf();
              for (typename PointVectorType::iterator iter = Pointvector[level].begin(); iter != Pointvector[level].end(); ++iter)
              {
                (*iter)->set_offset(0);
                Pointvector[level].update_cf(*iter);
              }
            }
              
          private:     
            /** @brief Slices mat into this->threads parts of (almost) equal size
            * @param level    Level for which slicing is requested
            * @param A     System matrix on all levels
            */
            void slice_new (unsigned int level, InternalType1 const & A)
            {  
              typedef typename SparseMatrixType::const_iterator1 ConstRowIterator;
              typedef typename SparseMatrixType::const_iterator2 ConstColIterator;
              
              // Determine index offset of all the slices (index of A[level] when the respective slice starts).
            #ifdef _OPENMP
              #pragma omp parallel for 
            #endif
              for (unsigned int i=0; i<=_threads; ++i)
              {
                // Offset of first piece is zero. Pieces 1,...,threads-1 have equal size while the last one might be greater.
                if (i == 0) Offset[i][level] = 0;
                else if (i == _threads) Offset[i][level] = A[level].size1();
                else Offset[i][level] = i*(A[level].size1()/_threads);
              }
            }   
            
            /** @brief Slices mat into pieces determined by this->Offset
            * @param level    Level to which Slices are saved
            * @param A     System matrix on all levels
            * @param Pointvector  Vector of points on all levels
            */
            void slice_build (unsigned int level, InternalType1 const & A, InternalType2 const & Pointvector)
            {
              typedef typename SparseMatrixType::const_iterator1 ConstRowIterator;
              typedef typename SparseMatrixType::const_iterator2 ConstColIterator;
              
              unsigned int x, y;
              amg_point *point;
              
            #ifdef _OPENMP
              #pragma omp parallel for private (x,y,point)
            #endif
              for (unsigned int i=0; i<_threads; ++i)
              {
                // Allocate space for the matrix slice and the pointvector.
                A_slice[i][level] = SparseMatrixType(Offset[i+1][level]-Offset[i][level],Offset[i+1][level]-Offset[i][level]);
                Pointvector_slice[i][level] = PointVectorType(Offset[i+1][level]-Offset[i][level]);
                
                // Iterate over the part that belongs to thread i (from Offset[i][level] to Offset[i+1][level]).
                ConstRowIterator row_iter = A[level].begin1();
                row_iter += Offset[i][level];
                x = row_iter.index1();
                    
                while (x < Offset[i+1][level] && row_iter != A[level].end1())
                {
                  // Set offset for point index and save point for the respective thread
                  point = Pointvector[level][x];
                  point->set_offset(Offset[i][level]);
                  Pointvector_slice[i][level].add_point(point);
                  
                  ConstColIterator col_iter = row_iter.begin();
                  y = col_iter.index2();
                  
                  // Save all coefficients from the matrix slice
                  while (y < Offset[i+1][level] && col_iter != row_iter.end())
                  {
                    if (y >= Offset[i][level])
                A_slice[i][level](x-Offset[i][level],y-Offset[i][level]) = *col_iter;
                    
                    ++col_iter;
                    y = col_iter.index2();
                  }
                  
                  ++row_iter;
                  x = row_iter.index1();
                }
              }
            }
        };  
        
        /** @brief Sparse matrix product. Calculates RES = A*B.
          * @param A    Left Matrix
          * @param B    Right Matrix
          * @param RES    Result Matrix
          */
        template <typename SparseMatrixType>
        void amg_mat_prod (SparseMatrixType & A, SparseMatrixType & B, SparseMatrixType & RES)
        {
          typedef typename SparseMatrixType::value_type ScalarType;
          typedef typename SparseMatrixType::iterator1 InternalRowIterator;
          typedef typename SparseMatrixType::iterator2 InternalColIterator;
          
          unsigned int x,y,z;
          ScalarType prod;
          RES = SparseMatrixType(A.size1(), B.size2());
          RES.clear();
          
    #ifdef _OPENMP
          #pragma omp parallel for private (x,y,z,prod) shared (A,B,RES)
    #endif
          for (x=0; x<A.size1(); ++x)
          {
            InternalRowIterator row_iter = A.begin1();
            row_iter += x;
            for (InternalColIterator col_iter = row_iter.begin(); col_iter != row_iter.end(); ++col_iter)
            {
              y = col_iter.index2(); 
              InternalRowIterator row_iter2 = B.begin1();
              row_iter2 += y;

              for(InternalColIterator col_iter2 = row_iter2.begin(); col_iter2 != row_iter2.end(); ++col_iter2)
              {
                z = col_iter2.index2();
                prod = *col_iter * *col_iter2;
                RES.add(x,z,prod);
              }
            }
          }
        }
        
        /** @brief Sparse Galerkin product: Calculates RES = trans(P)*A*P
          * @param A    Operator matrix (quadratic)
          * @param P    Prolongation/Interpolation matrix
          * @param RES    Result Matrix (Galerkin operator)
          */
        template <typename SparseMatrixType>
        void amg_galerkin_prod (SparseMatrixType & A, SparseMatrixType & P, SparseMatrixType & RES)
        {
          typedef typename SparseMatrixType::value_type ScalarType;
          typedef typename SparseMatrixType::iterator1 InternalRowIterator;
          typedef typename SparseMatrixType::iterator2 InternalColIterator;
          
          unsigned int x,y1,y2,z;
          amg_sparsevector<ScalarType> row;
          RES = SparseMatrixType(P.size2(), P.size2());
          RES.clear();
          
    #ifdef _OPENMP
          #pragma omp parallel for private (x,y1,y2,z,row) shared (A,P,RES)
    #endif      
          for (x=0; x<P.size2(); ++x)
          {
            row = amg_sparsevector<ScalarType>(A.size2());
            InternalRowIterator row_iter = P.begin1(true);
            row_iter += x;

            for (InternalColIterator col_iter = row_iter.begin(); col_iter != row_iter.end(); ++col_iter)
            {
              y1 = col_iter.index2(); 
              InternalRowIterator row_iter2 = A.begin1();
              row_iter2 += y1;
              
              for(InternalColIterator col_iter2 = row_iter2.begin(); col_iter2 != row_iter2.end(); ++col_iter2)
              {
                y2 = col_iter2.index2();
                row.add (y2, *col_iter * *col_iter2);
              }
            }
            for (typename amg_sparsevector<ScalarType>::iterator iter = row.begin(); iter != row.end(); ++iter)
            {
              y2 = iter.index();
              InternalRowIterator row_iter3 = P.begin1();
              row_iter3 += y2;
              
              for (InternalColIterator col_iter3 = row_iter3.begin(); col_iter3 != row_iter3.end(); ++col_iter3)
              {
                z = col_iter3.index2();
                RES.add (x, z, *col_iter3 * *iter);
              }
            }
          }
          
          #ifdef DEBUG
          std::cout << "Galerkin Operator: " << std::endl;
          printmatrix (RES);
          #endif
        }
        
        /** @brief Test triple-matrix product by comparing it to ublas functions. Very slow for large matrices!
          * @param A    Operator matrix (quadratic)
          * @param P    Prolongation/Interpolation matrix
          * @param A_i1    Result Matrix
          */
        template <typename SparseMatrixType>
        void test_triplematprod(SparseMatrixType & A, SparseMatrixType & P, SparseMatrixType  & A_i1)
        {
          typedef typename SparseMatrixType::value_type ScalarType;
          
          boost::numeric::ublas::compressed_matrix<ScalarType> A_temp (A.size1(), A.size2());
          A_temp = A;
          boost::numeric::ublas::compressed_matrix<ScalarType> P_temp (P.size1(), P.size2());
          P_temp = P;
          P.set_trans(true);
          boost::numeric::ublas::compressed_matrix<ScalarType> R_temp (P.size1(), P.size2());
          R_temp = P;
          P.set_trans(false);
          
          boost::numeric::ublas::compressed_matrix<ScalarType> RA (R_temp.size1(),A_temp.size2());
          RA = boost::numeric::ublas::prod(R_temp,A_temp);
          boost::numeric::ublas::compressed_matrix<ScalarType> RAP (RA.size1(),P_temp.size2());
          RAP = boost::numeric::ublas::prod(RA,P_temp);
          
          for (unsigned int x=0; x<RAP.size1(); ++x)
          {
            for (unsigned int y=0; y<RAP.size2(); ++y)
            {
              if (abs((ScalarType)RAP(x,y) - (ScalarType)A_i1(x,y)) > 0.0001)
                std::cout << x << " " << y << " " << RAP(x,y) << " " << A_i1(x,y) << std::endl;
            } 
          }
        }
        
        /** @brief Test if interpolation matrix makes sense. Only vanilla test though! Only checks if basic requirements are met!
          * @param A    Operator matrix (quadratic)
          * @param P    Prolongation/Interpolation matrix
          * @param Pointvector  Vector of points
          */
        template <typename SparseMatrixType, typename PointVectorType>
        void test_interpolation(SparseMatrixType & A, SparseMatrixType & P, PointVectorType & Pointvector)
        {
          for (unsigned int i=0; i<P.size1(); ++i)
          {
            if (Pointvector.is_cpoint(i))
            {
              bool set = false;
              for (unsigned int j=0; j<P.size2(); ++j)
              {
                if (P.isnonzero(i,j))
                {
                  if (P(i,j) != 1)
                    std::cout << "Error 1 in row " << i << std::endl;
                  if (P(i,j) == 1 && set)
                    std::cout << "Error 2 in row " << i << std::endl;
                  if (P(i,j) == 1 && !set)
                    set = true;
                }
              }
            }
            
            if (Pointvector.is_fpoint(i))
              for (unsigned int j=0; j<P.size2(); ++j)
              {
                if (P.isnonzero(i,j) && j> Pointvector.get_cpoints()-1)
                  std::cout << "Error 3 in row " << i << std::endl;
                if (P.isnonzero(i,j))
                {
                  bool set = false;
                  for (unsigned int k=0; k<P.size1(); ++k)
                  {
                    if (P.isnonzero(k,j))
                    {
                      if (Pointvector.is_cpoint(k) && P(k,j) == 1 && A.isnonzero(i,k))
                        set = true;      
                    }
                  }
                  if (!set)
                    std::cout << "Error 4 in row " << i << std::endl;
                }
              }
            }
        }
        
        
      } //namespace amg
    }
  }
}

#endif
