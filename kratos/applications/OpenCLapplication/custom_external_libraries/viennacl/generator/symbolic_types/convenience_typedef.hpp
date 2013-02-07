#ifndef VIENNACL_GENERATOR_SYMBOLIC_TYPES_CONVENIENCE_TYPEDEF_HPP
#define VIENNACL_GENERATOR_SYMBOLIC_TYPES_CONVENIENCE_TYPEDEF_HPP

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

/** @file convenience_typedef.hpp
 *  @brief Convenience typedefs for quick creation of symbolic types
 *
 *  Generator code contributed by Philippe Tillet
 */

#include "viennacl/generator/symbolic_types/symbolic_vector.hpp"
#include "viennacl/generator/symbolic_types/symbolic_matrix.hpp"
#include "viennacl/generator/symbolic_types/symbolic_scalars.hpp"

namespace viennacl 
{
  namespace generator
  {

    //Symbolic vectors : float
    
    typedef symbolic_vector<0,float,1> symv_0_f;
    typedef symbolic_vector<1,float,1> symv_1_f;
    typedef symbolic_vector<2,float,1> symv_2_f;
    typedef symbolic_vector<3,float,1> symv_3_f;
    typedef symbolic_vector<4,float,1> symv_4_f;
    typedef symbolic_vector<5,float,1> symv_5_f;
    typedef symbolic_vector<6,float,1> symv_6_f;

    typedef symbolic_vector<0,float,4> symv_0_f_4;
    typedef symbolic_vector<1,float,4> symv_1_f_4;
    typedef symbolic_vector<2,float,4> symv_2_f_4;
    typedef symbolic_vector<3,float,4> symv_3_f_4;
    typedef symbolic_vector<4,float,4> symv_4_f_4;
    typedef symbolic_vector<5,float,4> symv_5_f_4;
    typedef symbolic_vector<6,float,4> symv_6_f_4;

    typedef symbolic_vector<0,float,16> symv_0_f_16;
    typedef symbolic_vector<1,float,16> symv_1_f_16;
    typedef symbolic_vector<2,float,16> symv_2_f_16;
    typedef symbolic_vector<3,float,16> symv_3_f_16;
    typedef symbolic_vector<4,float,16> symv_4_f_16;
    typedef symbolic_vector<5,float,16> symv_5_f_16;
    typedef symbolic_vector<6,float,16> symv_6_f_16;


    //Symbolic vectors : double

    typedef symbolic_vector<0,double,1> symv_0_d;
    typedef symbolic_vector<1,double,1> symv_1_d;
    typedef symbolic_vector<2,double,1> symv_2_d;
    typedef symbolic_vector<3,double,1> symv_3_d;
    typedef symbolic_vector<4,double,1> symv_4_d;
    typedef symbolic_vector<5,double,1> symv_5_d;
    typedef symbolic_vector<6,double,1> symv_6_d;

    typedef symbolic_vector<0,double,4> symv_0_d_4;
    typedef symbolic_vector<1,double,4> symv_1_d_4;
    typedef symbolic_vector<2,double,4> symv_2_d_4;
    typedef symbolic_vector<3,double,4> symv_3_d_4;
    typedef symbolic_vector<4,double,4> symv_4_d_4;
    typedef symbolic_vector<5,double,4> symv_5_d_4;
    typedef symbolic_vector<6,double,4> symv_6_d_4;

    typedef symbolic_vector<0,double,16> symv_0_d_16;
    typedef symbolic_vector<1,double,16> symv_1_d_16;
    typedef symbolic_vector<2,double,16> symv_2_d_16;
    typedef symbolic_vector<3,double,16> symv_3_d_16;
    typedef symbolic_vector<4,double,16> symv_4_d_16;
    typedef symbolic_vector<5,double,16> symv_5_d_16;
    typedef symbolic_vector<6,double,16> symv_6_d_16;


    //Symbolic matrices : float

    typedef symbolic_matrix<0,float,viennacl::row_major,1> symm_0_f;
    typedef symbolic_matrix<1,float,viennacl::row_major,1> symm_1_f;
    typedef symbolic_matrix<2,float,viennacl::row_major,1> symm_2_f;
    typedef symbolic_matrix<3,float,viennacl::row_major,1> symm_3_f;
    typedef symbolic_matrix<4,float,viennacl::row_major,1> symm_4_f;
    typedef symbolic_matrix<5,float,viennacl::row_major,1> symm_5_f;
    typedef symbolic_matrix<6,float,viennacl::row_major,1> symm_6_f;

    typedef symbolic_matrix<0,float,viennacl::row_major,16> symm_0_f_r_16;
    typedef symbolic_matrix<1,float,viennacl::row_major,16> symm_1_f_r_16;
    typedef symbolic_matrix<2,float,viennacl::row_major,16> symm_2_f_r_16;
    typedef symbolic_matrix<3,float,viennacl::row_major,16> symm_3_f_r_16;
    typedef symbolic_matrix<4,float,viennacl::row_major,16> symm_4_f_r_16;
    typedef symbolic_matrix<5,float,viennacl::row_major,16> symm_5_f_r_16;
    typedef symbolic_matrix<6,float,viennacl::row_major,16> symm_6_f_r_16;


    //Symbolic matrices : double

    typedef symbolic_matrix<0,double,viennacl::row_major,1> symm_0_d;
    typedef symbolic_matrix<1,double,viennacl::row_major,1> symm_1_d;
    typedef symbolic_matrix<2,double,viennacl::row_major,1> symm_2_d;
    typedef symbolic_matrix<3,double,viennacl::row_major,1> symm_3_d;
    typedef symbolic_matrix<4,double,viennacl::row_major,1> symm_4_d;
    typedef symbolic_matrix<5,double,viennacl::row_major,1> symm_5_d;
    typedef symbolic_matrix<6,double,viennacl::row_major,1> symm_6_d;

    typedef symbolic_matrix<0,double,viennacl::row_major,16> symm_0_d_r_16;
    typedef symbolic_matrix<1,double,viennacl::row_major,16> symm_1_d_r_16;
    typedef symbolic_matrix<2,double,viennacl::row_major,16> symm_2_d_r_16;
    typedef symbolic_matrix<3,double,viennacl::row_major,16> symm_3_d_r_16;
    typedef symbolic_matrix<4,double,viennacl::row_major,16> symm_4_d_r_16;
    typedef symbolic_matrix<5,double,viennacl::row_major,16> symm_5_d_r_16;
    typedef symbolic_matrix<6,double,viennacl::row_major,16> symm_6_d_r_16;


    //CPU Symbolic scalar: float

    typedef cpu_symbolic_scalar<0,float> c_syms_0_f;
    typedef cpu_symbolic_scalar<1,float> c_syms_1_f;
    typedef cpu_symbolic_scalar<2,float> c_syms_2_f;
    typedef cpu_symbolic_scalar<3,float> c_syms_3_f;
    typedef cpu_symbolic_scalar<4,float> c_syms_4_f;
    typedef cpu_symbolic_scalar<5,float> c_syms_5_f;
    typedef cpu_symbolic_scalar<6,float> c_syms_6_f;


    //CPU Symbolic scalar: double

    typedef cpu_symbolic_scalar<0,double> c_syms_0_d;
    typedef cpu_symbolic_scalar<1,double> c_syms_1_d;
    typedef cpu_symbolic_scalar<2,double> c_syms_2_d;
    typedef cpu_symbolic_scalar<3,double> c_syms_3_d;
    typedef cpu_symbolic_scalar<4,double> c_syms_4_d;
    typedef cpu_symbolic_scalar<5,double> c_syms_5_d;
    typedef cpu_symbolic_scalar<6,double> c_syms_6_d;


    //GPU Symbolic scalar: float

    typedef gpu_symbolic_scalar<0,float> syms_0_f;
    typedef gpu_symbolic_scalar<1,float> syms_1_f;
    typedef gpu_symbolic_scalar<2,float> syms_2_f;
    typedef gpu_symbolic_scalar<3,float> syms_3_f;
    typedef gpu_symbolic_scalar<4,float> syms_4_f;
    typedef gpu_symbolic_scalar<5,float> syms_5_f;
    typedef gpu_symbolic_scalar<6,float> syms_6_f;


    //GPU Symbolic scalar: double

    typedef gpu_symbolic_scalar<0,double> syms_0_d;
    typedef gpu_symbolic_scalar<1,double> syms_1_d;
    typedef gpu_symbolic_scalar<2,double> syms_2_d;
    typedef gpu_symbolic_scalar<3,double> syms_3_d;
    typedef gpu_symbolic_scalar<4,double> syms_4_d;
    typedef gpu_symbolic_scalar<5,double> syms_5_d;
    typedef gpu_symbolic_scalar<6,double> syms_6_d;


  }
}

#endif


