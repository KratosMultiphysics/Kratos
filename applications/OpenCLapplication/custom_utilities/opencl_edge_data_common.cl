/*
==============================================================================
KratosOpenCLApplication
A library based on:
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi, Farshid Mossaiby
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
mossaiby@yahoo.com
- CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain


Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNERS.

The  above  copyright  notice  and  this permission  notice  shall  be
included in all copies or substantial portions of the Software.

THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT  SHALL THE AUTHORS OR COPYRIGHT HOLDERS  BE LIABLE FOR ANY
CLAIM, DAMAGES OR  OTHER LIABILITY, WHETHER IN AN  ACTION OF CONTRACT,
TORT  OR OTHERWISE, ARISING  FROM, OUT  OF OR  IN CONNECTION  WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

==============================================================================
*/

//
//   Project Name:        Kratos
//   Last Modified by:    $Author: mossaiby $
//   Date:                $Date: 2010-09-30 18:26:51 $
//   Revision:            $Revision: 1.00 $
//
//

//
// opencl_EdgeType_data_common.cl
//
// OpenCL kernels and functions used in opencl_EdgeType_data.h


//
// Macros to control the behavior of the code

#define USE_CONSERVATIVE_FORM_FOR_SCALAR_CONVECTION

#define USE_CONSERVATIVE_FORM_FOR_VECTOR_CONVECTION


//
// Enable OpenCL extensions

//
// Try enabling the AMD version

#pragma OPENCL EXTENSION cl_amd_fp64: enable

//
// If failed, try Khronos version

#ifndef cl_amd_fp64

	#pragma OPENCL EXTENSION cl_khr_fp64: enable

#endif


//
// OpenCL 1.0 adjustment

#if KRATOS_OCL_VERSION < 110

typedef double4 double3;

#endif


//
// Helper macros

#define KRATOS_OCL_LAPLACIANIJ_0_0(a)	((a).s0)
#define KRATOS_OCL_LAPLACIANIJ_0_1(a)	((a).s1)
#define KRATOS_OCL_LAPLACIANIJ_0_2(a)	((a).s2)
#define KRATOS_OCL_LAPLACIANIJ_1_0(a)	((a).s3)
#define KRATOS_OCL_LAPLACIANIJ_1_1(a)	((a).s4)
#define KRATOS_OCL_LAPLACIANIJ_1_2(a)	((a).s5)
#define KRATOS_OCL_LAPLACIANIJ_2_0(a)	((a).s6)
#define KRATOS_OCL_LAPLACIANIJ_2_1(a)	((a).s7)
#define KRATOS_OCL_LAPLACIANIJ_2_2(a)	((a).s8)

#define KRATOS_OCL_MASS(a)				((a).s9)

#define KRATOS_OCL_NI_DNJ_0(a)			((a).sa)
#define KRATOS_OCL_NI_DNJ_1(a)			((a).sb)
#define KRATOS_OCL_NI_DNJ_2(a)			((a).sc)

#define KRATOS_OCL_DNI_NJ_0(a)			((a).sd)
#define KRATOS_OCL_DNI_NJ_1(a)			((a).se)
#define KRATOS_OCL_DNI_NJ_2(a)			((a).sf)

#define KRATOS_OCL_COMP_0(a)			((a).s0)
#define KRATOS_OCL_COMP_1(a)			((a).s1)
#define KRATOS_OCL_COMP_2(a)			((a).s2)


//
// Fast math macros

// Currently these are not supported on GPUs

#ifdef __CPU__

	#define KRATOS_OCL_DIVIDE(x, y)			native_divide(x, y)
	#define KRATOS_OCL_RECIP(x)				native_recip(x)
	#define KRATOS_OCL_SQRT(x)				native_sqrt(x)

#else

	#define KRATOS_OCL_DIVIDE(x, y)			((x) / (y))
	#define KRATOS_OCL_RECIP(x)				(1.00 / (x))
	#define KRATOS_OCL_SQRT(x)				sqrt(x)

#endif

//
// Used types

typedef unsigned int IndexType;

typedef double3 VectorType;

typedef double ValueType;

typedef double16 EdgeType;


//
// OpenCL defines length() as length of the vector, so if we use double4 instead of double3, we have to take care of this

inline ValueType length3(double4 x)
{
	double4 t = x;
	t.s3 = 0.00;

	return KRATOS_OCL_SQRT(dot(t, t));
}

#if KRATOS_OCL_VERSION < 110
	#define KRATOS_OCL_LENGTH3(x)			length3(x)
	#define KRATOS_OCL_VECTOR3(x, y, z)		((VectorType)(x, y, z, 0.00))
#else
	#define KRATOS_OCL_LENGTH3(x)			length(x)
	#define KRATOS_OCL_VECTOR3(x, y, z)		((VectorType)(x, y, z))
#endif


// A dummy kernel for test
__kernel void Test(__global double *input, __global double *output, const double offset)
{
	const size_t id = get_global_id(0);
	const double iv = input[id];
	output[id] = 2.00 * sin(iv) * cos(iv) + offset;
}

//
// Edge specific kernels

inline void Add_Gp(const VectorType Ni_DNj, const VectorType DNi_Nj, VectorType *destination, const ValueType p_i, const ValueType p_j)
{
	 // destination[comp] -= Ni_DNj[comp] * p_j - DNi_Nj[comp] * p_i
	 *destination -= Ni_DNj * p_j - DNi_Nj * p_i;
}

inline void Sub_Gp(const VectorType Ni_DNj, const VectorType DNi_Nj, VectorType *destination, const ValueType p_i, const ValueType p_j)
{
	 // destination[comp] += Ni_DNj[comp] * p_j - DNi_Nj[comp] * p_i
	 *destination += Ni_DNj * p_j - DNi_Nj * p_i;
}

inline void Add_D_v(const VectorType Ni_DNj, ValueType *destination, const VectorType v_i, const VectorType v_j)
{
	// destination += Ni_DNj[comp] * (v_j[comp] - v_i[comp])
	*destination += dot(Ni_DNj, v_j - v_i);
}

inline void Sub_D_v(const VectorType Ni_DNj, ValueType *destination, const VectorType v_i, const VectorType v_j)
{
	// destination -= Ni_DNj[comp] * (v_j[comp] - v_i[comp])
	*destination -= dot(Ni_DNj, v_j - v_i);
}

inline void Add_grad_p(const VectorType Ni_DNj, VectorType *destination, const ValueType p_i, const ValueType p_j)
{
	// destination[comp] += Ni_DNj[comp] * (p_j - p_i)
	*destination += Ni_DNj * (p_j - p_i);
}

inline void Sub_grad_p(const VectorType Ni_DNj, VectorType *destination, const ValueType p_i, const ValueType p_j)
{
	// destination[comp] -= Ni_DNj[comp] * (p_j - p_i)
	*destination -= Ni_DNj * (p_j - p_i);
}

inline void Add_div_v(const VectorType Ni_DNj, const VectorType DNi_Nj, ValueType *destination, const VectorType v_i, const VectorType v_j)
{
	// destination -= Ni_DNj[comp]*v_j[comp] - DNi_Nj[comp]*v_i[comp]
	*destination -= dot(Ni_DNj, v_j) - dot(DNi_Nj, v_i);
}

inline void Sub_div_v(const VectorType Ni_DNj, const VectorType DNi_Nj, ValueType *destination, const VectorType v_i, const VectorType v_j)
{
	// destination += Ni_DNj[comp]*v_j[comp] - DNi_Nj[comp]*v_i[comp]
	*destination += dot(Ni_DNj, v_j) - dot(DNi_Nj, v_i);
}

inline void CalculateScalarLaplacian(const ValueType LaplacianIJ_0_0, const ValueType LaplacianIJ_1_1, const ValueType LaplacianIJ_2_2, ValueType *l_ij)
{
	// l_ij += LaplacianIJ(comp, comp)
	*l_ij = LaplacianIJ_0_0 + LaplacianIJ_1_1 + LaplacianIJ_2_2;
}

inline void Add_ConvectiveContribution(const VectorType Ni_DNj, VectorType *destination, const VectorType a_i, const VectorType U_i, const VectorType a_j, const VectorType U_j)
{

#ifdef USE_CONSERVATIVE_FORM_FOR_VECTOR_CONVECTION

	// temp += a_i[k_comp] * Ni_DNj[k_comp]
	// destination[l_comp] += temp * (U_j[l_comp] - U_i[l_comp])
	*destination += dot(a_i, Ni_DNj) * (U_j - U_i);

#else

	// aux_i += a_i[k_comp] * Ni_DNj[k_comp]
	// aux_j += a_j[k_comp] * Ni_DNj[k_comp]
	// destination[l_comp] += aux_j * U_j[l_comp] - aux_i * U_i[l_comp]
	*destination += dot(a_j * U_j - a_i * U_i, Ni_DNj);

#endif

}

inline void Sub_ConvectiveContribution(const VectorType Ni_DNj, VectorType *destination, const VectorType a_i, const VectorType U_i, const VectorType a_j, const VectorType U_j)
{

#ifdef USE_CONSERVATIVE_FORM_FOR_VECTOR_CONVECTION

	// temp += a_i[k_comp] * Ni_DNj[k_comp]
	// destination[l_comp] -= temp * (U_j[l_comp] - U_i[l_comp])
	*destination -= dot(a_i, Ni_DNj) * (U_j - U_i);

#else

	// aux_i += a_i[k_comp] * Ni_DNj[k_comp]
	// aux_j += a_j[k_comp] * Ni_DNj[k_comp]
	// destination[l_comp] -= aux_j * U_j[l_comp] - aux_i * U_i[l_comp]
	*destination -= dot(a_j * U_j - a_i * U_i, Ni_DNj);

#endif

}

inline void Add_ConvectiveContribution2(const VectorType Ni_DNj, ValueType *destination, const VectorType a_i, const ValueType phi_i, const VectorType a_j, const ValueType phi_j)
{

#ifdef USE_CONSERVATIVE_FORM_FOR_SCALAR_CONVECTION

	// temp += a_i[k_comp] * Ni_DNj[k_comp]
	*destination += dot(a_i, Ni_DNj) * (phi_j - phi_i);

#else

	// aux_i += a_i[k_comp] * Ni_DNj[k_comp]
	// aux_j += a_j[k_comp] * Ni_DNj[k_comp]
	*destination += dot(a_j * phi_j - a_i * phi_i, Ni_DNj);

#endif

}

inline void Sub_ConvectiveContribution2(const VectorType Ni_DNj, ValueType *destination, const VectorType a_i, const ValueType phi_i, const VectorType a_j, const ValueType phi_j)
{

#ifdef USE_CONSERVATIVE_FORM_FOR_SCALAR_CONVECTION

	// temp += a_i[k_comp] * Ni_DNj[k_comp]
	*destination -= dot(a_i, Ni_DNj) * (phi_j - phi_i);

#else

	// aux_i += a_i[k_comp] * Ni_DNj[k_comp]
	// aux_j += a_j[k_comp] * Ni_DNj[k_comp]
	*destination -= dot(a_j * phi_j - a_i * phi_i, Ni_DNj);

#endif

}

inline void CalculateConvectionStabilization_LOW(const VectorType LaplacianIJ_0, const VectorType LaplacianIJ_1, const VectorType LaplacianIJ_2, VectorType *stab_low, const VectorType a_i, const VectorType U_i, const VectorType a_j, const VectorType U_j)
{
	// conv_stab += a_i[k_comp] * a_i[m_comp] * LaplacianIJ(k_comp,m_comp)
	// stab_low[l_comp] = conv_stab * (U_j[l_comp] - U_i[l_comp])
	*stab_low = dot(a_i, KRATOS_OCL_VECTOR3(dot(a_i, LaplacianIJ_0), dot(a_i, LaplacianIJ_1), dot(a_i, LaplacianIJ_2))) * (U_j - U_i);
}

inline void CalculateConvectionStabilization_LOW2(const VectorType LaplacianIJ_0, const VectorType LaplacianIJ_1, const VectorType LaplacianIJ_2, ValueType *stab_low, const VectorType a_i, const ValueType phi_i, const VectorType a_j, const ValueType phi_j)
{
	// conv_stab += a_i[k_comp] * a_i[m_comp] * LaplacianIJ(k_comp,m_comp)
	*stab_low = dot(a_i, KRATOS_OCL_VECTOR3(dot(a_i, LaplacianIJ_0), dot(a_i, LaplacianIJ_1), dot(a_i, LaplacianIJ_2))) * (phi_j - phi_i);
}

inline void CalculateConvectionStabilization_HIGH(const VectorType Ni_DNj, VectorType *stab_high, const VectorType a_i, const VectorType pi_i, const VectorType a_j, const VectorType pi_j)
{

#ifdef USE_CONSERVATIVE_FORM_FOR_VECTOR_CONVECTION

	// temp += a_i[k_comp] * Ni_DNj[k_comp]
	// stab_high[l_comp] = -temp * (pi_j[l_comp] - pi_i[l_comp])
	*stab_high = dot(a_i, Ni_DNj) * (pi_i - pi_j);

#else

	// aux_i += a_i[k_comp] * Ni_DNj[k_comp]
	// aux_j += a_j[k_comp] * Ni_DNj[k_comp]
	// stab_high[l_comp] = -(aux_j * pi_j[l_comp] - aux_i * pi_i[l_comp])
	*stab_high = dot(a_i * pi_i - a_j * pi_j, Ni_DNj);

#endif

}

inline void CalculateConvectionStabilization_HIGH2(const VectorType Ni_DNj, ValueType *stab_high, const VectorType a_i, const ValueType pi_i, const VectorType a_j, const ValueType pi_j)
{

#ifdef USE_CONSERVATIVE_FORM_FOR_VECTOR_CONVECTION

	// temp += a_i[k_comp] * Ni_DNj[k_comp]
	*stab_high = dot(a_i, Ni_DNj) * (pi_i - pi_j);

#else

	// aux_i += a_i[k_comp] * Ni_DNj[k_comp]
	// aux_j += a_j[k_comp] * Ni_DNj[k_comp]
	*stab_high = dot(a_i * pi_i - a_j * pi_j, Ni_DNj);

#endif

}

inline void Add_StabContribution(VectorType *destination, const ValueType tau, const ValueType beta, const VectorType stab_low, const VectorType stab_high)
{
	// destination[l_comp] += tau * (stab_low[l_comp] - beta * stab_high[l_comp])
	*destination += tau * (stab_low - beta * stab_high);
}

inline void Sub_StabContribution(VectorType *destination, const ValueType tau, const ValueType beta, const VectorType stab_low, const VectorType stab_high)
{
	// destination[l_comp] -= tau * (stab_low[l_comp] - beta * stab_high[l_comp])
	*destination -= tau * (stab_low - beta * stab_high);
}

inline void Add_StabContribution2(ValueType *destination, const ValueType tau, const ValueType beta, const ValueType stab_low, const ValueType stab_high)
{
	*destination += tau * (stab_low - beta * stab_high);
}

inline void Sub_StabContribution2(ValueType *destination, const ValueType tau, const ValueType beta, const ValueType stab_low, const ValueType stab_high)
{
	*destination -= tau * (stab_low - beta * stab_high);
}

inline void Add_ViscousContribution(const ValueType LaplacianIJ_0_0, const ValueType LaplacianIJ_1_1, const ValueType LaplacianIJ_2_2, __global VectorType *destination, const VectorType U_i, const ValueType nu_i, const VectorType U_j, const ValueType nu_j)
{
	// L += LaplacianIJ(l_comp, l_comp)
	// destination[l_comp] += nu_i * L * (U_j[l_comp] - U_i[l_comp])
	*destination += nu_i * (LaplacianIJ_0_0 + LaplacianIJ_1_1 + LaplacianIJ_2_2) * (U_j - U_i);
}
