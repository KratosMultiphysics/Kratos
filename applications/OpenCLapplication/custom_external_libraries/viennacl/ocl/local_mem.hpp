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


#ifndef _VIENNACL_LOCAL_MEM_HPP_
#define _VIENNACL_LOCAL_MEM_HPP_

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

