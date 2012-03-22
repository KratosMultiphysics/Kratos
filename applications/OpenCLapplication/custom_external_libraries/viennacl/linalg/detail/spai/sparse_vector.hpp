#ifndef VIENNACL_LINALG_DETAIL_SPAI_SPARSE_VECTOR_HPP
#define VIENNACL_LINALG_DETAIL_SPAI_SPARSE_VECTOR_HPP

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

/** @file viennacl/linalg/detail/spai/sparse_vector.hpp
    @brief Implementation of a helper sparse vector class for SPAI. Experimental in 1.2.x.
    
    SPAI code contributed by Nikolay Lukash
*/

#include <utility>
#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include <vector>
#include <math.h>
#include <map>
//local includes
//#include <omp.h>


namespace viennacl
{
  namespace linalg
  {
    namespace detail
    {
      namespace spai
      {
        
        /**
         * @brief Represents sparse vector based on std::map<unsigned int, ScalarType>
         */
        template <typename ScalarType>
        class sparse_vector{
        public:
            typedef typename std::map<unsigned int, ScalarType>::iterator iterator;
            typedef typename std::map<unsigned int, ScalarType>::const_iterator const_iterator;
            sparse_vector(){
            }
            
            /** @brief Set the index of the vector in the original matrix
             * May only be called once.
             */
            //getter
            ScalarType& operator[] (const unsigned int ind){
                return _v[ind];
                
            }
            
            void clear(){
                _v.clear();
            }
            
            const_iterator find(const unsigned int var) const{
                return _v.find(var);
            }
            
            iterator find(const unsigned int var){
                return _v.find(var);
            }
            
            const_iterator begin() const{
                return _v.begin();
            }
            
            const_iterator end() const{
                return _v.end();
            }
            
            
            iterator begin(){
                return _v.begin();
            }
            
            iterator end(){
                return _v.end();
            }
            
            
        private:
            unsigned int _size;
            std::map<unsigned int, ScalarType> _v;
        };
      }
    }
  }
}

#endif