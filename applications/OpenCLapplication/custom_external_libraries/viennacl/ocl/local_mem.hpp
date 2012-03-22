#ifndef VIENNACL_OCL_LOCAL_MEM_HPP_
#define VIENNACL_OCL_LOCAL_MEM_HPP_

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


/** @file local_mem.hpp
    @brief A local (shared) memory object for OpenCL
*/

namespace viennacl
{
  namespace ocl
  {
    /** @brief A class representing local (shared) OpenCL memory. Typically used as kernel argument */
    class local_mem
    {
      public:
        local_mem(unsigned int s) : size_(s) {}
        
        /** @brief Returns size in bytes */
        unsigned int size() const { return size_; }

        /** @brief Sets the size of the local memory in bytes */
        void size(unsigned int s) { size_ = s; }

      private:
        unsigned int size_;
    };
    
  }
}
#endif

