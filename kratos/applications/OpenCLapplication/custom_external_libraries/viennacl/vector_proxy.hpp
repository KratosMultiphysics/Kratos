#ifndef VIENNACL_VECTOR_PROXY_HPP_
#define VIENNACL_VECTOR_PROXY_HPP_

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

/** @file vector_proxy.hpp
    @brief Proxy classes for vectors.
*/

#include "viennacl/forwards.h"
#include "viennacl/range.hpp"
#include "viennacl/slice.hpp"
#include "viennacl/vector.hpp"
#include "viennacl/tools/entry_proxy.hpp"

namespace viennacl
{
  /** @brief Class for representing non-strided subvectors of a bigger vector x.
    *
    * In MATLAB notation, this could for example refer to the subvector x(3:8) of a vector x.
    */
  template <typename VectorType>
  class vector_range : public vector_base<typename VectorType::cpu_value_type>
  {
      typedef vector_range<VectorType>             self_type;
      typedef vector_base<typename VectorType::cpu_value_type> base_type;

    public:
      typedef typename VectorType::value_type      value_type;
      typedef range::size_type                     size_type;
      typedef range::difference_type               difference_type;
      typedef value_type                           reference;
      typedef const value_type &                   const_reference;
      typedef typename VectorType::const_iterator  const_iterator;
      typedef typename VectorType::iterator        iterator;

      typedef typename VectorType::cpu_value_type    cpu_value_type;

      static const int alignment = VectorType::alignment;

      vector_range(VectorType & v, range const & entry_range)
       : base_type(v.handle(), entry_range.size(), v.start() + v.stride() * entry_range.start(), v.stride()) {}


      using base_type::operator=;

  };



  /////////////////////////////////////////////////////////////
  ///////////////////////// CPU to GPU ////////////////////////
  /////////////////////////////////////////////////////////////

  template <typename VectorType, typename SCALARTYPE>
  void copy(const VectorType & cpu_vector,
            vector_range<vector<SCALARTYPE> > & gpu_vector_range )
  {
    assert(cpu_vector.end() - cpu_vector.begin() >= 0 && bool("Range must have nonnegative length!"));

    if (cpu_vector.end() - cpu_vector.begin() > 0)
    {
      //we require that the size of the gpu_vector is larger or equal to the cpu-size
      std::vector<SCALARTYPE> temp_buffer(cpu_vector.end() - cpu_vector.begin());
      std::copy(cpu_vector.begin(), cpu_vector.end(), temp_buffer.begin());
      viennacl::backend::memory_write(gpu_vector_range.handle(), sizeof(SCALARTYPE)*gpu_vector_range.start(), sizeof(SCALARTYPE)*temp_buffer.size(), &(temp_buffer[0]));
    }
  }


  /** @brief Transfer from a cpu vector to a gpu vector. Convenience wrapper for viennacl::linalg::fast_copy(cpu_vec.begin(), cpu_vec.end(), gpu_vec.begin());
  *
  * @param cpu_vec    A cpu vector. Type requirements: Iterator can be obtained via member function .begin() and .end()
  * @param gpu_vec    The gpu vector.
  */
  template <typename CPUVECTOR, typename VectorType>
  void fast_copy(const CPUVECTOR & cpu_vec, vector_range<VectorType> & gpu_vec)
  {
    viennacl::fast_copy(cpu_vec.begin(), cpu_vec.end(), gpu_vec.begin());
  }

  /////////////////////////////////////////////////////////////
  ///////////////////////// GPU to CPU ////////////////////////
  /////////////////////////////////////////////////////////////


  template <typename SCALARTYPE, typename VectorType>
  void copy(vector_range<vector<SCALARTYPE> > const & gpu_vector_range,
            VectorType & cpu_vector)
  {
    assert(cpu_vector.end() - cpu_vector.begin() >= 0 && bool("Range must have nonnegative length!"));

    if (cpu_vector.end() > cpu_vector.begin())
    {
      std::vector<SCALARTYPE> temp_buffer(cpu_vector.end() - cpu_vector.begin());
      viennacl::backend::memory_read(gpu_vector_range.handle(), sizeof(SCALARTYPE)*gpu_vector_range.start(), sizeof(SCALARTYPE)*temp_buffer.size(), &(temp_buffer[0]));

      //now copy entries to cpu_vec:
      std::copy(temp_buffer.begin(), temp_buffer.end(), cpu_vector.begin());
    }
  }


  /** @brief Transfer from a GPU vector range to a CPU vector. Convenience wrapper for viennacl::linalg::fast_copy(gpu_vec.begin(), gpu_vec.end(), cpu_vec.begin());
  *
  * @param gpu_vec    A gpu vector range.
  * @param cpu_vec    The cpu vector. Type requirements: Output iterator can be obtained via member function .begin()
  */
  template <typename VectorType, typename CPUVECTOR>
  void fast_copy(vector_range< VectorType > const & gpu_vec,
                 CPUVECTOR & cpu_vec )
  {
    viennacl::fast_copy(gpu_vec.begin(), gpu_vec.end(), cpu_vec.begin());
  }



  //
  // Convenience function
  //
  template <typename VectorType>
  vector_range<VectorType> project(VectorType & vec, viennacl::range const & r1)
  {
    return vector_range<VectorType>(vec, r1);
  }

  template <typename VectorType>
  vector_range<VectorType> project(viennacl::vector_range<VectorType> & vec, viennacl::range const & r1)
  {
    assert(r1.size() <= vec.size() && bool("Size of range invalid!"));
    return vector_range<VectorType>(vec, viennacl::range(vec.start() + r1.start(), vec.start() + r1.start() + r1.size()));
  }

//
//
//
/////////////////////////////// Slice /////////////////////////////////////////////
//
//
//



  /** @brief Class for representing strided subvectors of a bigger vector x.
    *
    * In MATLAB notation, this could for example refer to the subvector x(3:2:8) of a vector x.
    */
  template <typename VectorType>
  class vector_slice : public vector_base<typename VectorType::cpu_value_type>
  {
      typedef vector_slice<VectorType>             self_type;
      typedef vector_base<typename VectorType::cpu_value_type> base_type;

    public:
      typedef typename VectorType::value_type      value_type;
      typedef slice::size_type                     size_type;
      typedef slice::difference_type               difference_type;
      typedef value_type                           reference;
      typedef const value_type &                   const_reference;
      typedef typename VectorType::const_iterator  const_iterator;
      typedef typename VectorType::iterator        iterator;

      typedef typename VectorType::cpu_value_type  cpu_value_type;

      static const int alignment = VectorType::alignment;

      vector_slice(VectorType & v, slice const & entry_slice)
          : base_type(v.handle(), entry_slice.size(), v.start() + v.stride() * entry_slice.start(), v.stride() * entry_slice.stride()) {}


      using base_type::operator=;

  };


  /////////////////////////////////////////////////////////////
  ///////////////////////// CPU to GPU ////////////////////////
  /////////////////////////////////////////////////////////////

  template <typename VectorType, typename SCALARTYPE>
  void copy(const VectorType & cpu_vector,
            vector_slice<vector<SCALARTYPE> > & gpu_vector_slice )
  {
    if (cpu_vector.size() > 0)
    {
      std::vector<SCALARTYPE> temp_buffer(gpu_vector_slice.stride() * gpu_vector_slice.size());

      viennacl::backend::memory_read(gpu_vector_slice.handle(), sizeof(SCALARTYPE)*gpu_vector_slice.start(), sizeof(SCALARTYPE)*temp_buffer.size(), &(temp_buffer[0]));

      for (vcl_size_t i=0; i<cpu_vector.size(); ++i)
        temp_buffer[i * gpu_vector_slice.stride()] = cpu_vector[i];

      viennacl::backend::memory_write(gpu_vector_slice.handle(), sizeof(SCALARTYPE)*gpu_vector_slice.start(), sizeof(SCALARTYPE)*temp_buffer.size(), &(temp_buffer[0]));
    }
  }



  /////////////////////////////////////////////////////////////
  ///////////////////////// GPU to CPU ////////////////////////
  /////////////////////////////////////////////////////////////


  template <typename VectorType, typename SCALARTYPE>
  void copy(vector_slice<vector<SCALARTYPE> > const & gpu_vector_slice,
            VectorType & cpu_vector)
  {
    assert(gpu_vector_slice.end() - gpu_vector_slice.begin() >= 0 && bool("Range must have nonnegative length!"));

    if (gpu_vector_slice.end() - gpu_vector_slice.begin() > 0)
    {
      std::vector<SCALARTYPE> temp_buffer(gpu_vector_slice.stride() * gpu_vector_slice.size());
      viennacl::backend::memory_read(gpu_vector_slice.handle(), sizeof(SCALARTYPE)*gpu_vector_slice.start(), sizeof(SCALARTYPE)*temp_buffer.size(), &(temp_buffer[0]));

      for (vcl_size_t i=0; i<cpu_vector.size(); ++i)
        cpu_vector[i] = temp_buffer[i * gpu_vector_slice.stride()];
    }
  }





  //
  // Convenience functions
  //
  template <typename VectorType>
  vector_slice<VectorType> project(VectorType & vec, viennacl::slice const & s1)
  {
    assert(s1.size() <= vec.size() && bool("Size of slice larger than vector size!"));
    return vector_slice<VectorType>(vec, s1);
  }

  template <typename VectorType>
  vector_slice<VectorType> project(viennacl::vector_slice<VectorType> & vec, viennacl::slice const & s1)
  {
    assert(s1.size() <= vec.size() && bool("Size of slice larger than vector proxy!"));
    return vector_slice<VectorType>(vec, viennacl::slice(vec.start() + s1.start(), vec.stride() * s1.stride(), s1.size()));
  }

  // interaction with range and vector_range:

  template <typename VectorType>
  vector_slice<VectorType> project(viennacl::vector_slice<VectorType> & vec, viennacl::range const & r1)
  {
    assert(r1.size() <= vec.size() && bool("Size of slice larger than vector proxy!"));
    return vector_slice<VectorType>(vec, viennacl::slice(vec.start() + r1.start(), vec.stride(), r1.size()));
  }

  template <typename VectorType>
  vector_slice<VectorType> project(viennacl::vector_range<VectorType> & vec, viennacl::slice const & s1)
  {
    assert(s1.size() <= vec.size() && bool("Size of slice larger than vector proxy!"));
    return vector_slice<VectorType>(vec, viennacl::range(vec.start() + s1.start(), s1.stride(), s1.size()));
  }


}

#endif
