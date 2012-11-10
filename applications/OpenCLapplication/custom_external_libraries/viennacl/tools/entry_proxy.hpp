#ifndef VIENNACL_TOOLS_ENTRY_PROXY_HPP_
#define VIENNACL_TOOLS_ENTRY_PROXY_HPP_

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

/** @file entry_proxy.hpp
    @brief A proxy class for entries in a vector
*/


#include "viennacl/forwards.h"
#include "viennacl/ocl/backend.hpp"
#include "viennacl/scalar.hpp"

namespace viennacl
{
    //proxy class for single vector entries (this is a slow operation!!)
    /**
    * @brief A proxy class for a single element of a vector or matrix. This proxy should not be noticed by end-users of the library.
    *
    * This proxy provides access to a single entry of a vector. If the element is assigned to a GPU object, no unnecessary transfers to the CPU and back to GPU are initiated.
    *
    * @tparam SCALARTYPE Either float or double
    */
    template <typename SCALARTYPE>
    class entry_proxy
    {
      public:
        /** @brief The constructor for the proxy class. Declared explicit to avoid any surprises created by the compiler.
        *
        * @param mem_offset The memory offset in multiples of sizeof(SCALARTYPE) relative to the memory pointed to by the handle
        * @param mem_handle A viennacl::ocl::handle for the memory buffer on the GPU.
        */
        explicit entry_proxy(unsigned int mem_offset, 
                             viennacl::ocl::handle<cl_mem> const & mem_handle) 
         : _index(mem_offset), _mem_handle(mem_handle) {};
        
         
        //operators:
        /** @brief Inplace addition of a CPU floating point value
        */
        entry_proxy & operator+=(SCALARTYPE value)
        {
          SCALARTYPE temp = read();
          temp += value; 
          write(temp);         
          return *this;
        }

        /** @brief Inplace subtraction of a CPU floating point value
        */
        entry_proxy &  operator-=(SCALARTYPE value)
        {
          SCALARTYPE temp = read();
          temp -= value; 
          write(temp);         
          return *this;
        }

        /** @brief Inplace multiplication with a CPU floating point value
        */
        entry_proxy &  operator*=(SCALARTYPE value)
        {
          SCALARTYPE temp = read();
          temp *= value; 
          write(temp);         
          return *this;
        }

        /** @brief Inplace division by a CPU floating point value
        */
        entry_proxy &  operator/=(SCALARTYPE value)
        {
          SCALARTYPE temp = read();
          temp /= value; 
          write(temp);         
          return *this;
        }

        /** @brief Assignment of a CPU floating point value
        */
        entry_proxy &  operator=(SCALARTYPE value)
        {
          write(value);
          return *this;
        }

        /** @brief Assignment of a GPU floating point value. Avoids unnecessary GPU->CPU->GPU transfers
        */
        entry_proxy & operator=(scalar<SCALARTYPE> const & value)
        {
          cl_int err = clEnqueueCopyBuffer(viennacl::ocl::get_queue().handle().get(), value.handle().get(), _mem_handle.get(), 0, sizeof(SCALARTYPE)*_index, sizeof(SCALARTYPE), 0, NULL, NULL);
          //assert(err == CL_SUCCESS);
          VIENNACL_ERR_CHECK(err);
          return *this;
        }

        /** @brief Assignment of another GPU value.
        */
        entry_proxy &  operator=(entry_proxy const & other)
        {
          cl_int err = clEnqueueCopyBuffer(viennacl::ocl::get_queue().handle().get(),
                                           other._mem_handle.get(), //src
                                           _mem_handle.get(),       //dest
                                           sizeof(SCALARTYPE) * other._index, //offset src
                                           sizeof(SCALARTYPE) * _index,       //offset dest
                                           sizeof(SCALARTYPE), 0, NULL, NULL);
          VIENNACL_ERR_CHECK(err);
          return *this;
        }

        //type conversion:
        // allows to write something like:
        //  double test = vector(4);
        /** @brief Conversion to a CPU floating point value.
        *
        *  This conversion allows to write something like
        *    double test = vector(4);
        *  However, one has to keep in mind that CPU<->GPU transfers are very slow compared to CPU<->CPU operations.
        */
        operator SCALARTYPE () const
        {
          SCALARTYPE temp = read();
          return temp;
        }
        
        /** @brief Returns the index of the represented element
        */
        unsigned int index() const { return _index; }
        
        /** @brief Returns the memory viennacl::ocl::handle
        */
        viennacl::ocl::handle<cl_mem> const & handle() const { return _mem_handle; }

      private:
        /** @brief Reads an element from the GPU to the CPU
        */
        SCALARTYPE read() const
        {
          SCALARTYPE temp;
          cl_int err;
          err = clEnqueueReadBuffer(viennacl::ocl::get_queue().handle().get(), _mem_handle.get(), CL_TRUE, sizeof(SCALARTYPE)*_index, sizeof(SCALARTYPE), &temp, 0, NULL, NULL);
          //assert(err == CL_SUCCESS);
          VIENNACL_ERR_CHECK(err);
          viennacl::ocl::get_queue().finish();
          return temp;
        }
        
        /** @brief Writes a floating point value to the GPU
        */
        void write(SCALARTYPE value)
        {
          cl_int err;
          err = clEnqueueWriteBuffer(viennacl::ocl::get_queue().handle().get(), _mem_handle.get(), CL_TRUE, sizeof(SCALARTYPE)*_index, sizeof(SCALARTYPE), &value, 0, NULL, NULL);
          //assert(err == CL_SUCCESS);
          VIENNACL_ERR_CHECK(err);
        }
        
        unsigned int _index;
        viennacl::ocl::handle<cl_mem> const & _mem_handle;
    }; //entry_proxy
    
}

#endif
