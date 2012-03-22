#ifndef VIENNACL_OCL_KERNEL_HPP_
#define VIENNACL_OCL_KERNEL_HPP_

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

/** @file kernel.hpp
    @brief Representation of an OpenCL kernel in ViennaCL.
*/

#ifdef __APPLE__
#include <OpenCL/cl.h>
#else
#include <CL/cl.h>
#endif

#include "viennacl/ocl/forwards.h"
#include "viennacl/ocl/backend.hpp"
#include "viennacl/ocl/handle.hpp"
#include "viennacl/ocl/program.hpp"
#include "viennacl/ocl/device.hpp"
#include "viennacl/ocl/local_mem.hpp"

namespace viennacl
{
  namespace ocl
  {
    
    /** @brief Represents an OpenCL kernel within ViennaCL */
    class kernel
    {
      template <typename KernelType>
      friend void enqueue(KernelType & k, viennacl::ocl::command_queue const & queue);
      
      
    public:
      kernel() : handle_(0)
      {
        #if defined(VIENNACL_DEBUG_ALL) || defined(VIENNACL_DEBUG_KERNEL)
        std::cout << "ViennaCL: Creating kernel object (default CTOR)" << std::endl;
        #endif
        set_work_size_defaults();
      }
      
      kernel(viennacl::ocl::handle<cl_program> const & prog, std::string const & name) 
       : handle_(0), program_(prog), name_(name), init_done_(false)
      {
        #if defined(VIENNACL_DEBUG_ALL) || defined(VIENNACL_DEBUG_KERNEL)
        std::cout << "ViennaCL: Creating kernel object (full CTOR)" << std::endl;
        #endif
        set_work_size_defaults();
      }
      
      kernel(kernel const & other) 
       : handle_(other.handle_), program_(other.program_), name_(other.name_), init_done_(other.init_done_)
      {
        #if defined(VIENNACL_DEBUG_ALL) || defined(VIENNACL_DEBUG_KERNEL)
        std::cout << "ViennaCL: Creating kernel object (Copy CTOR)" << std::endl;
        #endif
        local_work_size_[0] = other.local_work_size_[0];
        local_work_size_[1] = other.local_work_size_[1];
        
        global_work_size_[0] = other.global_work_size_[0];
        global_work_size_[1] = other.global_work_size_[1];
      }
      
      viennacl::ocl::kernel & operator=(const kernel & other)
      {
        #if defined(VIENNACL_DEBUG_ALL) || defined(VIENNACL_DEBUG_KERNEL)
        std::cout << "ViennaCL: Assigning kernel object" << std::endl;
        #endif
        handle_ = other.handle_;
        program_ = other.program_;
        name_ = other.name_;
        init_done_ = other.init_done_;
        local_work_size_[0] = other.local_work_size_[0];
        local_work_size_[1] = other.local_work_size_[1];
        global_work_size_[0] = other.global_work_size_[0];
        global_work_size_[1] = other.global_work_size_[1];
        return *this;
      }
      
      
      /** @brief Sets an unsigned integer argument at the provided position */
      void arg(unsigned int pos, cl_uint val)
      {
        init();
        #if defined(VIENNACL_DEBUG_ALL) || defined(VIENNACL_DEBUG_KERNEL)
        std::cout << "ViennaCL: Setting unsigned long kernel argument at pos " << pos << " for kernel " << name_ << std::endl;
        #endif
        cl_int err = clSetKernelArg(handle_.get(), pos, sizeof(cl_uint), (void*)&val);
        VIENNACL_ERR_CHECK(err);
      }

      /** @brief Sets a single precision floating point argument at the provided position */
      void arg(unsigned int pos, float val)
      {
        init();
        #if defined(VIENNACL_DEBUG_ALL) || defined(VIENNACL_DEBUG_KERNEL)
        std::cout << "ViennaCL: Setting floating point kernel argument at pos " << pos << " for kernel " << name_ << std::endl;
        #endif
        cl_int err = clSetKernelArg(handle_.get(), pos, sizeof(float), (void*)&val);
        VIENNACL_ERR_CHECK(err);
      }

      /** @brief Sets a double precision floating point argument at the provided position */
      void arg(unsigned int pos, double val)
      {
        init();
        #if defined(VIENNACL_DEBUG_ALL) || defined(VIENNACL_DEBUG_KERNEL)
        std::cout << "ViennaCL: Setting double precision kernel argument at pos " << pos << " for kernel " << name_ << std::endl;
        #endif
        cl_int err = clSetKernelArg(handle_.get(), pos, sizeof(double), (void*)&val);
        VIENNACL_ERR_CHECK(err);
      }

      //generic handling: call .handle() member
      /** @brief Sets an OpenCL memory object at the provided position */
      template<class VCL_TYPE>
      void arg(unsigned int pos, VCL_TYPE const & val)
      {
        init();
        #if defined(VIENNACL_DEBUG_ALL) || defined(VIENNACL_DEBUG_KERNEL)
        std::cout << "ViennaCL: Setting generic kernel argument at pos " << pos << " for kernel " << name_ << std::endl;
        #endif
        cl_mem temp = val.handle().get();
        cl_int err = clSetKernelArg(handle_.get(), pos, sizeof(cl_mem), (void*)&temp);
        VIENNACL_ERR_CHECK(err);
      }
      
      //forward handles directly:
      /** @brief Sets an OpenCL object at the provided position */
      template<class CL_TYPE>
      void arg(unsigned int pos, viennacl::ocl::handle<CL_TYPE> const & h)
      {
        //arg(pos, h);
        init();
        #if defined(VIENNACL_DEBUG_ALL) || defined(VIENNACL_DEBUG_KERNEL)
        std::cout << "ViennaCL: Setting handle kernel argument at pos " << pos << " for kernel " << name_ << std::endl;
        #endif
        CL_TYPE temp = h.get();
        cl_int err = clSetKernelArg(handle_.get(), pos, sizeof(CL_TYPE), (void*)&temp);
        VIENNACL_ERR_CHECK(err);
      }
      
      
      //local buffer argument:
      /** @brief Sets an OpenCL local memory object at the provided position */
      void arg(unsigned int pos, const local_mem & mem)
      {
        unsigned int size =  mem.size();
        init();
        #if defined(VIENNACL_DEBUG_ALL) || defined(VIENNACL_DEBUG_KERNEL)
        std::cout << "ViennaCL: Setting local memory kernel argument at pos " << pos << " for kernel " << name_ << std::endl;
        #endif
        cl_int err = clSetKernelArg(handle_.get(), pos, size, 0);
        VIENNACL_ERR_CHECK(err);
      }
      
      
      
      /** @brief Convenience function for setting one kernel parameter */
      template <typename T0>
      kernel & operator()(T0 const & t0)
      {
         arg(0, t0);
         return *this;
      }     

      /** @brief Convenience function for setting two kernel parameters */
      template <typename T0, typename T1>
      kernel & operator()(T0 const & t0, T1 const & t1)
      {
         arg(0, t0); arg(1, t1);
         return *this;
      }     

      /** @brief Convenience function for setting three kernel parameters */
      template <typename T0, typename T1, typename T2>
      kernel & operator()(T0 const & t0, T1 const & t1, T2 const & t2)
      {
         arg(0, t0); arg(1, t1); arg(2, t2);
         return *this;
      }     

      /** @brief Convenience function for setting four kernel parameters */
      template <typename T0, typename T1, typename T2, typename T3>
      kernel & operator()(T0 const & t0, T1 const & t1, T2 const & t2, T3 const & t3)
      {
         arg(0, t0); arg(1, t1); arg(2, t2); arg(3, t3);
         return *this;
      }     

      /** @brief Convenience function for setting five kernel parameters */
      template <typename T0, typename T1, typename T2, typename T3, typename T4>
      kernel & operator()(T0 const & t0, T1 const & t1, T2 const & t2, T3 const & t3, T4 const & t4)
      {
         arg(0, t0); arg(1, t1); arg(2, t2); arg(3, t3); arg(4, t4);
         return *this;
      }     

      /** @brief Convenience function for setting six kernel parameters */
      template <typename T0, typename T1, typename T2, typename T3, typename T4, typename T5>
      kernel & operator()(T0 const & t0, T1 const & t1, T2 const & t2, T3 const & t3, T4 const & t4, T5 const & t5)
      {
         arg(0, t0); arg(1, t1); arg(2, t2); arg(3, t3); arg(4, t4); arg(5, t5);
         return *this;
      }     

      /** @brief Convenience function for setting seven kernel parameters */
      template <typename T0, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6>
      kernel & operator()(T0 const & t0, T1 const & t1, T2 const & t2, T3 const & t3, T4 const & t4, T5 const & t5, T6 const & t6)
      {
         arg(0, t0); arg(1, t1); arg(2, t2); arg(3, t3); arg(4, t4); arg(5, t5); arg(6, t6);
         return *this;
      }     

      /** @brief Convenience function for setting eight kernel parameters */
      template <typename T0, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7>
      kernel & operator()(T0 const & t0, T1 const & t1, T2 const & t2, T3 const & t3, T4 const & t4, T5 const & t5, T6 const & t6, T7 const & t7)
      {
         arg(0, t0); arg(1, t1); arg(2, t2); arg(3, t3); arg(4, t4); arg(5, t5); arg(6, t6); arg(7, t7);
         return *this;
      }     

      /** @brief Convenience function for setting nine kernel parameters */
      template <typename T0, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8>
      kernel & operator()(T0 const & t0, T1 const & t1, T2 const & t2, T3 const & t3, T4 const & t4, T5 const & t5, T6 const & t6, T7 const & t7, T8 const & t8)
      {
         arg(0, t0); arg(1, t1); arg(2, t2); arg(3, t3); arg(4, t4); arg(5, t5); arg(6, t6); arg(7, t7); arg(8, t8);
         return *this;
      }     

      /** @brief Convenience function for setting ten kernel parameters */
      template <typename T0, typename T1, typename T2, typename T3, typename T4,
                typename T5, typename T6, typename T7, typename T8, typename T9>
      kernel & operator()(T0 const & t0, T1 const & t1, T2 const & t2, T3 const & t3, T4 const & t4,
                          T5 const & t5, T6 const & t6, T7 const & t7, T8 const & t8, T9 const & t9)
      {
         arg(0, t0); arg(1, t1); arg(2, t2); arg(3, t3); arg(4, t4); arg(5, t5); arg(6, t6); arg(7, t7); arg(8, t8); arg(9, t9);
         return *this;
      }     

      /** @brief Convenience function for setting eleven kernel parameters */
      template <typename T0, typename T1, typename T2, typename T3, typename T4, typename T5,
                typename T6, typename T7, typename T8, typename T9, typename T10>
      kernel & operator()(T0 const & t0, T1 const & t1, T2 const & t2, T3 const & t3, T4 const & t4, T5 const & t5,
                          T6 const & t6, T7 const & t7, T8 const & t8, T9 const & t9, T10 const & t10)
      {
         arg(0, t0); arg(1, t1); arg(2, t2); arg(3, t3); arg(4, t4); arg(5, t5); arg(6, t6); arg(7, t7); arg(8, t8); arg(9, t9); arg(10, t10);
         return *this;
      }     

      /** @brief Convenience function for setting twelve kernel parameters */
      template <typename T0, typename T1, typename T2, typename T3, typename T4, typename T5,
                typename T6, typename T7, typename T8, typename T9, typename T10, typename T11>
      kernel & operator()(T0 const & t0, T1 const & t1, T2 const & t2, T3 const & t3, T4 const & t4, T5 const & t5,
                          T6 const & t6, T7 const & t7, T8 const & t8, T9 const & t9, T10 const & t10, T11 const & t11)
      {
         arg(0, t0); arg(1, t1); arg(2, t2); arg(3, t3); arg(4, t4); arg(5, t5);
         arg(6, t6); arg(7, t7); arg(8, t8); arg(9, t9); arg(10, t10); arg(11, t11);
         return *this;
      }     

      /** @brief Convenience function for setting thirteen kernel parameters */
      template <typename T0, typename T1, typename T2, typename T3, typename T4, typename T5,
                typename T6, typename T7, typename T8, typename T9, typename T10, typename T11, typename T12>
      kernel & operator()(T0 const & t0, T1 const & t1, T2 const & t2, T3 const & t3, T4 const & t4, T5 const & t5,
                          T6 const & t6, T7 const & t7, T8 const & t8, T9 const & t9, T10 const & t10, T11 const & t11, T12 const & t12)
      {
         arg(0, t0); arg(1, t1); arg(2, t2); arg(3, t3); arg(4, t4); arg(5, t5);
         arg(6, t6); arg(7, t7); arg(8, t8); arg(9, t9); arg(10, t10); arg(11, t11); arg(12, t12);
         return *this;
      }     

      /** @brief Convenience function for setting fourteen kernel parameters */
      template <typename T0, typename T1, typename T2, typename T3, typename T4, typename T5,
                typename T6, typename T7, typename T8, typename T9, typename T10, typename T11,
                typename T12, typename T13>
      kernel & operator()(T0 const & t0, T1 const & t1, T2 const & t2, T3 const & t3, T4 const & t4, T5 const & t5,
                          T6 const & t6, T7 const & t7, T8 const & t8, T9 const & t9, T10 const & t10, T11 const & t11,
                          T12 const & t12, T13 const & t13)
      {
         arg(0, t0); arg(1, t1); arg(2, t2); arg(3, t3); arg(4, t4); arg(5, t5);
         arg(6, t6); arg(7, t7); arg(8, t8); arg(9, t9); arg(10, t10); arg(11, t11);
         arg(12, t12); arg(13, t13);
         return *this;
      }     

      /** @brief Convenience function for setting fifteen kernel parameters */
      template <typename T0, typename T1, typename T2, typename T3, typename T4, typename T5,
                typename T6, typename T7, typename T8, typename T9, typename T10, typename T11,
                typename T12, typename T13, typename T14>
      kernel & operator()(T0 const & t0, T1 const & t1, T2 const & t2, T3 const & t3, T4 const & t4, T5 const & t5,
                          T6 const & t6, T7 const & t7, T8 const & t8, T9 const & t9, T10 const & t10, T11 const & t11,
                          T12 const & t12, T13 const & t13, T14 const & t14)
      {
         arg(0, t0); arg(1, t1); arg(2, t2); arg(3, t3); arg(4, t4); arg(5, t5);
         arg(6, t6); arg(7, t7); arg(8, t8); arg(9, t9); arg(10, t10); arg(11, t11);
         arg(12, t12); arg(13, t13); arg(14, t14);
         return *this;
      }     

      /** @brief Convenience function for setting sixteen kernel parameters */
      template <typename T0, typename T1, typename T2, typename T3, typename T4, typename T5,
                typename T6, typename T7, typename T8, typename T9, typename T10, typename T11,
                typename T12, typename T13, typename T14, typename T15>
      kernel & operator()(T0 const & t0, T1 const & t1, T2 const & t2, T3 const & t3, T4 const & t4, T5 const & t5,
                          T6 const & t6, T7 const & t7, T8 const & t8, T9 const & t9, T10 const & t10, T11 const & t11,
                          T12 const & t12, T13 const & t13, T14 const & t14, T15 const & t15)
      {
         arg(0, t0); arg(1, t1); arg(2, t2); arg(3, t3); arg(4, t4); arg(5, t5);
         arg(6, t6); arg(7, t7); arg(8, t8); arg(9, t9); arg(10, t10); arg(11, t11);
         arg(12, t12); arg(13, t13); arg(14, t14); arg(15, t15);
         return *this;
      }     

      /** @brief Convenience function for setting seventeen kernel parameters */
      template <typename T0, typename T1, typename T2, typename T3, typename T4, typename T5,
                typename T6, typename T7, typename T8, typename T9, typename T10, typename T11,
                typename T12, typename T13, typename T14, typename T15, typename T16>
      kernel & operator()(T0 const & t0, T1 const & t1, T2 const & t2, T3 const & t3, T4 const & t4, T5 const & t5,
                          T6 const & t6, T7 const & t7, T8 const & t8, T9 const & t9, T10 const & t10, T11 const & t11,
                          T12 const & t12, T13 const & t13, T14 const & t14, T15 const & t15, T16 const & t16)
      {
         arg(0, t0); arg(1, t1); arg(2, t2); arg(3, t3); arg(4, t4); arg(5, t5);
         arg(6, t6); arg(7, t7); arg(8, t8); arg(9, t9); arg(10, t10); arg(11, t11);
         arg(12, t12); arg(13, t13); arg(14, t14); arg(15, t15); arg(16, t16);
         return *this;
      }     

      /** @brief Convenience function for setting eighteen kernel parameters */
      template <typename T0, typename T1, typename T2, typename T3, typename T4, typename T5,
                typename T6, typename T7, typename T8, typename T9, typename T10, typename T11,
                typename T12, typename T13, typename T14, typename T15, typename T16, typename T17>
      kernel & operator()(T0 const & t0, T1 const & t1, T2 const & t2, T3 const & t3, T4 const & t4, T5 const & t5,
                          T6 const & t6, T7 const & t7, T8 const & t8, T9 const & t9, T10 const & t10, T11 const & t11,
                          T12 const & t12, T13 const & t13, T14 const & t14, T15 const & t15, T16 const & t16, T17 const & t17)
      {
         arg(0, t0); arg(1, t1); arg(2, t2); arg(3, t3); arg(4, t4); arg(5, t5);
         arg(6, t6); arg(7, t7); arg(8, t8); arg(9, t9); arg(10, t10); arg(11, t11);
         arg(12, t12); arg(13, t13); arg(14, t14); arg(15, t15); arg(16, t16); arg(17, t17);
         return *this;
      }     

      /** @brief Convenience function for setting nineteen kernel parameters */
      template <typename T0, typename T1, typename T2, typename T3, typename T4, typename T5,
                typename T6, typename T7, typename T8, typename T9, typename T10, typename T11,
                typename T12, typename T13, typename T14, typename T15, typename T16, typename T17,
                typename T18>
      kernel & operator()(T0 const & t0, T1 const & t1, T2 const & t2, T3 const & t3, T4 const & t4, T5 const & t5,
                          T6 const & t6, T7 const & t7, T8 const & t8, T9 const & t9, T10 const & t10, T11 const & t11,
                          T12 const & t12, T13 const & t13, T14 const & t14, T15 const & t15, T16 const & t16, T17 const & t17,
                          T18 const & t18
                         )
      {
         arg(0, t0); arg(1, t1); arg(2, t2); arg(3, t3); arg(4, t4); arg(5, t5);
         arg(6, t6); arg(7, t7); arg(8, t8); arg(9, t9); arg(10, t10); arg(11, t11);
         arg(12, t12); arg(13, t13); arg(14, t14); arg(15, t15); arg(16, t16); arg(17, t17);
         arg(18, t18);
         return *this;
      }     

      /** @brief Convenience function for setting twenty kernel parameters */
      template <typename T0, typename T1, typename T2, typename T3, typename T4, typename T5,
                typename T6, typename T7, typename T8, typename T9, typename T10, typename T11,
                typename T12, typename T13, typename T14, typename T15, typename T16, typename T17,
                typename T18, typename T19>
      kernel & operator()(T0 const & t0, T1 const & t1, T2 const & t2, T3 const & t3, T4 const & t4, T5 const & t5,
                          T6 const & t6, T7 const & t7, T8 const & t8, T9 const & t9, T10 const & t10, T11 const & t11,
                          T12 const & t12, T13 const & t13, T14 const & t14, T15 const & t15, T16 const & t16, T17 const & t17,
                          T18 const & t18, T19 const & t19
                         )
      {
         arg(0, t0); arg(1, t1); arg(2, t2); arg(3, t3); arg(4, t4); arg(5, t5);
         arg(6, t6); arg(7, t7); arg(8, t8); arg(9, t9); arg(10, t10); arg(11, t11);
         arg(12, t12); arg(13, t13); arg(14, t14); arg(15, t15); arg(16, t16); arg(17, t17);
         arg(18, t18); arg(19, t19);
         return *this;
      }     

      /** @brief Convenience function for setting twentyone kernel parameters */
      template <typename T0, typename T1, typename T2, typename T3, typename T4, typename T5,
                typename T6, typename T7, typename T8, typename T9, typename T10, typename T11,
                typename T12, typename T13, typename T14, typename T15, typename T16, typename T17,
                typename T18, typename T19, typename T20>
      kernel & operator()(T0 const & t0, T1 const & t1, T2 const & t2, T3 const & t3, T4 const & t4, T5 const & t5,
                          T6 const & t6, T7 const & t7, T8 const & t8, T9 const & t9, T10 const & t10, T11 const & t11,
                          T12 const & t12, T13 const & t13, T14 const & t14, T15 const & t15, T16 const & t16, T17 const & t17,
                          T18 const & t18, T19 const & t19, T20 const & t20
                         )
      {
         arg(0, t0); arg(1, t1); arg(2, t2); arg(3, t3); arg(4, t4); arg(5, t5);
         arg(6, t6); arg(7, t7); arg(8, t8); arg(9, t9); arg(10, t10); arg(11, t11);
         arg(12, t12); arg(13, t13); arg(14, t14); arg(15, t15); arg(16, t16); arg(17, t17);
         arg(18, t18); arg(19, t19); arg(20, t20);
         return *this;
      }     

      /** @brief Convenience function for setting twentytwo kernel parameters */
      template <typename T0, typename T1, typename T2, typename T3, typename T4, typename T5,
                typename T6, typename T7, typename T8, typename T9, typename T10, typename T11,
                typename T12, typename T13, typename T14, typename T15, typename T16, typename T17,
                typename T18, typename T19, typename T20, typename T21>
      kernel & operator()(T0 const & t0, T1 const & t1, T2 const & t2, T3 const & t3, T4 const & t4, T5 const & t5,
                          T6 const & t6, T7 const & t7, T8 const & t8, T9 const & t9, T10 const & t10, T11 const & t11,
                          T12 const & t12, T13 const & t13, T14 const & t14, T15 const & t15, T16 const & t16, T17 const & t17,
                          T18 const & t18, T19 const & t19, T20 const & t20, T21 const & t21
                         )
      {
         arg(0, t0); arg(1, t1); arg(2, t2); arg(3, t3); arg(4, t4); arg(5, t5);
         arg(6, t6); arg(7, t7); arg(8, t8); arg(9, t9); arg(10, t10); arg(11, t11);
         arg(12, t12); arg(13, t13); arg(14, t14); arg(15, t15); arg(16, t16); arg(17, t17);
         arg(18, t18); arg(19, t19); arg(20, t20); arg(21, t21);
         return *this;
      }     

      /** @brief Convenience function for setting twentythree kernel parameters */
      template <typename T0, typename T1, typename T2, typename T3, typename T4, typename T5,
                typename T6, typename T7, typename T8, typename T9, typename T10, typename T11,
                typename T12, typename T13, typename T14, typename T15, typename T16, typename T17,
                typename T18, typename T19, typename T20, typename T21, typename T22>
      kernel & operator()(T0 const & t0, T1 const & t1, T2 const & t2, T3 const & t3, T4 const & t4, T5 const & t5,
                          T6 const & t6, T7 const & t7, T8 const & t8, T9 const & t9, T10 const & t10, T11 const & t11,
                          T12 const & t12, T13 const & t13, T14 const & t14, T15 const & t15, T16 const & t16, T17 const & t17,
                          T18 const & t18, T19 const & t19, T20 const & t20, T21 const & t21, T22 const & t22
                         )
      {
         arg(0, t0); arg(1, t1); arg(2, t2); arg(3, t3); arg(4, t4); arg(5, t5);
         arg(6, t6); arg(7, t7); arg(8, t8); arg(9, t9); arg(10, t10); arg(11, t11);
         arg(12, t12); arg(13, t13); arg(14, t14); arg(15, t15); arg(16, t16); arg(17, t17);
         arg(18, t18); arg(19, t19); arg(20, t20); arg(21, t21);  arg(22, t22);
         return *this;
      }     

      /** @brief Returns the local work size at the respective dimension
      *
      * @param index   Dimension index (currently either 0 or 1)
      */
      size_t local_work_size(int index = 0) const
      {
        assert(index == 0 || index == 1);
        return local_work_size_[index];
      }
      /** @brief Returns the global work size at the respective dimension
      *
      * @param index   Dimension index (currently either 0 or 1)
      */
      size_t global_work_size(int index = 0) const
      { 
        assert(index == 0 || index == 1);
        return global_work_size_[index];
      }

      /** @brief Sets the local work size at the respective dimension
      *
      * @param index   Dimension index (currently either 0 or 1)
      * @param s       The new local work size
      */
      void local_work_size(int index, size_t s)
      {
        #if defined(VIENNACL_DEBUG_ALL) || defined(VIENNACL_DEBUG_KERNEL)
        std::cout << "ViennaCL: Setting local work size to " << s << " at index " << index << " for kernel " << name_ << std::endl;
        #endif
        assert(index == 0 || index == 1);
        local_work_size_[index] = s;
      }
      /** @brief Sets the global work size at the respective dimension
      *
      * @param index   Dimension index (currently either 0 or 1)
      * @param s       The new global work size
      */
      void global_work_size(int index, size_t s)
      { 
        #if defined(VIENNACL_DEBUG_ALL) || defined(VIENNACL_DEBUG_KERNEL)
        std::cout << "ViennaCL: Setting global work size to " << s << " at index " << index << " for kernel " << name_ << std::endl;
        #endif
        assert(index == 0 || index == 1);
        global_work_size_[index] = s;
      }

      std::string const & name() const { return name_; }

      viennacl::ocl::handle<cl_kernel> const & handle() const { return handle_; }


    private:
      void create_kernel()
      {
        cl_int err;
        #if defined(VIENNACL_DEBUG_ALL) || defined(VIENNACL_DEBUG_KERNEL)
        std::cout << "ViennaCL: Building kernel " << name_ << std::endl;
        #endif
        handle_ = clCreateKernel(program_.get(), name_.c_str(), &err);
        
        if (err != CL_SUCCESS)
        {
          #if defined(VIENNACL_DEBUG_ALL) || defined(VIENNACL_DEBUG_KERNEL)
          std::cout << "ViennaCL: Could not create kernel '" << name_ << "'." << std::endl;
          #endif
          //std::cerr << "Could not build kernel '" << name_ << "'." << std::endl;
        }
        VIENNACL_ERR_CHECK(err);
      }

      void set_work_size_defaults()
      {
        if (viennacl::ocl::current_device().type() == CL_DEVICE_TYPE_GPU)
        {
          local_work_size_[0] = 128; local_work_size_[1] = 0;
          global_work_size_[0] = 128*128; global_work_size_[1] = 0;
        }
        else //assume CPU type:
        {
          //conservative assumption: one thread per CPU core:
          local_work_size_[0] = 1; local_work_size_[1] = 0;
          global_work_size_[0] = viennacl::ocl::current_device().max_compute_units(); global_work_size_[1] = 0;
        }
      }

      void init()
      {
        if (!init_done_)
        {
          create_kernel();
          init_done_ = true;
        }
      }
      
      viennacl::ocl::handle<cl_kernel> handle_;
      viennacl::ocl::handle<cl_program> program_;
      std::string name_;
      bool init_done_;
      size_t local_work_size_[2];
      size_t global_work_size_[2];
    };
    
  } //namespace ocl
} //namespace viennacl

#endif
