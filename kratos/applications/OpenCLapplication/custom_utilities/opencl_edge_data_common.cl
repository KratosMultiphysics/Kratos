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

#if __OPENCL_VERSION__ < 110

typedef double4 double3;

#endif


//
// Helper macros

#define KRATOS_OCL_LAPLACIANIJ_0_0(a)	((*a).s0)
#define KRATOS_OCL_LAPLACIANIJ_0_1(a)	((*a).s1)
#define KRATOS_OCL_LAPLACIANIJ_0_2(a)	((*a).s2)
#define KRATOS_OCL_LAPLACIANIJ_1_0(a)	((*a).s3)
#define KRATOS_OCL_LAPLACIANIJ_1_1(a)	((*a).s4)
#define KRATOS_OCL_LAPLACIANIJ_1_2(a)	((*a).s5)
#define KRATOS_OCL_LAPLACIANIJ_2_0(a)	((*a).s6)
#define KRATOS_OCL_LAPLACIANIJ_2_1(a)	((*a).s7)
#define KRATOS_OCL_LAPLACIANIJ_2_2(a)	((*a).s8)

#define KRATOS_OCL_MASS(a)				((*a).s9)

#define KRATOS_OCL_NI_DNJ_0(a)			((*a).sa)
#define KRATOS_OCL_NI_DNJ_1(a)			((*a).sb)
#define KRATOS_OCL_NI_DNJ_2(a)			((*a).sc)

#define KRATOS_OCL_DNI_NJ_0(a)			((*a).sd)
#define KRATOS_OCL_DNI_NJ_1(a)			((*a).se)
#define KRATOS_OCL_DNI_NJ_2(a)			((*a).sf)

#define KRATOS_OCL_COMP_0(a)			((*a).s0)
#define KRATOS_OCL_COMP_1(a)			((*a).s1)
#define KRATOS_OCL_COMP_2(a)			((*a).s2)


//
// Fast math macros

#define KRATOS_OCL_NATIVE_DIVIDE(x, y)	((x) / (y))		/* native_divide(x, y) */
#define KRATOS_OCL_NATIVE_RECIP(x)		(1.00 / (x))	/* native_recip(x) */
#define KRATOS_OCL_NATIVE_SQRT(x)		sqrt(x)			/* native_sqrt(x) */

//
// Used types

// TODO: Any better candidate? Or even a change in used fields to enable vector operations in kernels
typedef double16 EdgeType;

typedef unsigned int IndexType;

typedef double3 VectorType;

typedef double ValueType;


//
// OpenCL defines length() as length of the vector, so if we use double4 instead of double3, we have to take care of this

inline ValueType length3(double4 x)
{
	double4 t = x;
	t.s3 = 0.00;

	return KRATOS_OCL_NATIVE_SQRT(dot(t, t));
}

#if __OPENCL_VERSION__ < 110
	#define KRATOS_OCL_LENGTH3(x)		length3(x)
#else
	#define KRATOS_OCL_LENGTH3(x)		length(x)
#endif


// A dummy kernel for test
__kernel void Test(__global const double *input, __global double *output, const double offset)
{
	const size_t id = get_global_id(0);
	const double iv = input[id];
	output[id] = 2.00 * sin(iv) * cos(iv) + offset;
}

//
// Edge specific kernels

// TODO: Optimize kernels (vectorize, use built-in functions, avoid excessive use of __global data)
// TODO: Some __global attributes must be removed, depending on the usage, like Add_grad_p

void Add_Gp(__global EdgeType *a, __global VectorType *destination, const ValueType p_i, const ValueType p_j)
{
	 // destination[comp] -= Ni_DNj[comp] * p_j - DNi_Nj[comp] * p_i
	 KRATOS_OCL_COMP_0(destination) -= KRATOS_OCL_NI_DNJ_0(a) * p_j - KRATOS_OCL_DNI_NJ_0(a) * p_i;
	 KRATOS_OCL_COMP_1(destination) -= KRATOS_OCL_NI_DNJ_1(a) * p_j - KRATOS_OCL_DNI_NJ_1(a) * p_i;
	 KRATOS_OCL_COMP_2(destination) -= KRATOS_OCL_NI_DNJ_2(a) * p_j - KRATOS_OCL_DNI_NJ_2(a) * p_i;
}

void Sub_Gp(__global EdgeType *a, __global VectorType *destination, const ValueType p_i, const ValueType p_j)
{
	 // destination[comp] += Ni_DNj[comp] * p_j - DNi_Nj[comp] * p_i
	 KRATOS_OCL_COMP_0(destination) += KRATOS_OCL_NI_DNJ_0(a) * p_j - KRATOS_OCL_DNI_NJ_0(a) * p_i;
	 KRATOS_OCL_COMP_1(destination) += KRATOS_OCL_NI_DNJ_1(a) * p_j - KRATOS_OCL_DNI_NJ_1(a) * p_i;
	 KRATOS_OCL_COMP_2(destination) += KRATOS_OCL_NI_DNJ_2(a) * p_j - KRATOS_OCL_DNI_NJ_2(a) * p_i;
}

void Add_D_v(__global EdgeType *a, __global ValueType *destination, __global const VectorType *v_i, __global const VectorType *v_j)
{
	// destination += Ni_DNj[comp] * (v_j[comp] - v_i[comp])
	*destination +=
		KRATOS_OCL_NI_DNJ_0(a) * (KRATOS_OCL_COMP_0(v_j) - KRATOS_OCL_COMP_0(v_i)) +
		KRATOS_OCL_NI_DNJ_1(a) * (KRATOS_OCL_COMP_1(v_j) - KRATOS_OCL_COMP_1(v_i)) +
		KRATOS_OCL_NI_DNJ_2(a) * (KRATOS_OCL_COMP_2(v_j) - KRATOS_OCL_COMP_2(v_i));
}

void Sub_D_v(__global EdgeType *a, __global ValueType *destination, __global const VectorType *v_i, __global const VectorType *v_j)
{
	// destination -= Ni_DNj[comp] * (v_j[comp] - v_i[comp])
	*destination -=
		KRATOS_OCL_NI_DNJ_0(a) * (KRATOS_OCL_COMP_0(v_j) - KRATOS_OCL_COMP_0(v_i)) +
		KRATOS_OCL_NI_DNJ_1(a) * (KRATOS_OCL_COMP_1(v_j) - KRATOS_OCL_COMP_1(v_i)) +
		KRATOS_OCL_NI_DNJ_2(a) * (KRATOS_OCL_COMP_2(v_j) - KRATOS_OCL_COMP_2(v_i));
}

void Add_grad_p(__global EdgeType *a, VectorType *destination, const ValueType p_i, const ValueType p_j)
{
	// destination[comp] += Ni_DNj[comp] * (p_j - p_i)
	ValueType dp = p_j - p_i;

	KRATOS_OCL_COMP_0(destination) += KRATOS_OCL_NI_DNJ_0(a) * dp;
	KRATOS_OCL_COMP_1(destination) += KRATOS_OCL_NI_DNJ_1(a) * dp;
	KRATOS_OCL_COMP_2(destination) += KRATOS_OCL_NI_DNJ_2(a) * dp;
}

void Sub_grad_p(__global EdgeType *a, __global VectorType *destination, const ValueType p_i, const ValueType p_j)
{
	// destination[comp] -= Ni_DNj[comp] * (p_j - p_i)
	ValueType dp = p_j - p_i;

	KRATOS_OCL_COMP_0(destination) -= KRATOS_OCL_NI_DNJ_0(a) * dp;
	KRATOS_OCL_COMP_1(destination) -= KRATOS_OCL_NI_DNJ_1(a) * dp;
	KRATOS_OCL_COMP_2(destination) -= KRATOS_OCL_NI_DNJ_2(a) * dp;
}

void Add_div_v(__global EdgeType *a, __global ValueType *destination, __global const VectorType *v_i, __global const VectorType *v_j)
{
	// destination -= Ni_DNj[comp]*v_j[comp] - DNi_Nj[comp]*v_i[comp]
	*destination -=
		KRATOS_OCL_NI_DNJ_0(a) * KRATOS_OCL_COMP_0(v_j) - KRATOS_OCL_DNI_NJ_0(a) * KRATOS_OCL_COMP_0(v_i) +
		KRATOS_OCL_NI_DNJ_1(a) * KRATOS_OCL_COMP_1(v_j) - KRATOS_OCL_DNI_NJ_1(a) * KRATOS_OCL_COMP_1(v_i) +
		KRATOS_OCL_NI_DNJ_2(a) * KRATOS_OCL_COMP_2(v_j) - KRATOS_OCL_DNI_NJ_2(a) * KRATOS_OCL_COMP_2(v_i);
}

void Sub_div_v(__global EdgeType *a, __global ValueType *destination, __global const VectorType *v_i, __global const VectorType *v_j)
{
	// destination += Ni_DNj[comp]*v_j[comp] - DNi_Nj[comp]*v_i[comp]
	*destination +=
		KRATOS_OCL_NI_DNJ_0(a) * KRATOS_OCL_COMP_0(v_j) - KRATOS_OCL_DNI_NJ_0(a) * KRATOS_OCL_COMP_0(v_i) +
		KRATOS_OCL_NI_DNJ_1(a) * KRATOS_OCL_COMP_1(v_j) - KRATOS_OCL_DNI_NJ_1(a) * KRATOS_OCL_COMP_1(v_i) +
		KRATOS_OCL_NI_DNJ_2(a) * KRATOS_OCL_COMP_2(v_j) - KRATOS_OCL_DNI_NJ_2(a) * KRATOS_OCL_COMP_2(v_i);
}

void CalculateScalarLaplacian(__global EdgeType *a, ValueType *l_ij)
{
	// l_ij += LaplacianIJ(comp, comp)
	*l_ij = KRATOS_OCL_LAPLACIANIJ_0_0(a) + KRATOS_OCL_LAPLACIANIJ_1_1(a) + KRATOS_OCL_LAPLACIANIJ_2_2(a);
}

void Add_ConvectiveContribution(__global EdgeType *a, __global VectorType *destination,
	__global const VectorType *a_i, __global const VectorType *U_i,
	__global const VectorType *a_j, __global const VectorType *U_j)
{

#ifdef USE_CONSERVATIVE_FORM_FOR_VECTOR_CONVECTION

	// temp += a_i[k_comp] * Ni_DNj[k_comp]
	ValueType temp =
		KRATOS_OCL_COMP_0(a_i) * KRATOS_OCL_NI_DNJ_0(a) +
		KRATOS_OCL_COMP_1(a_i) * KRATOS_OCL_NI_DNJ_1(a) +
		KRATOS_OCL_COMP_2(a_i) * KRATOS_OCL_NI_DNJ_2(a);

	// destination[l_comp] += temp * (U_j[l_comp] - U_i[l_comp])
	KRATOS_OCL_COMP_0(destination) += temp * (KRATOS_OCL_COMP_0(U_j) - KRATOS_OCL_COMP_0(U_i));
	KRATOS_OCL_COMP_1(destination) += temp * (KRATOS_OCL_COMP_1(U_j) - KRATOS_OCL_COMP_1(U_i));
	KRATOS_OCL_COMP_2(destination) += temp * (KRATOS_OCL_COMP_2(U_j) - KRATOS_OCL_COMP_2(U_i));

#else

	// aux_i += a_i[k_comp] * Ni_DNj[k_comp]
	ValueType aux_i =
		KRATOS_OCL_COMP_0(a_i) * KRATOS_OCL_NI_DNJ_0(a) +
		KRATOS_OCL_COMP_1(a_i) * KRATOS_OCL_NI_DNJ_1(a) +
		KRATOS_OCL_COMP_2(a_i) * KRATOS_OCL_NI_DNJ_2(a);

	// aux_j += a_j[k_comp] * Ni_DNj[k_comp]
	ValueType aux_j =
		KRATOS_OCL_COMP_0(a_j) * KRATOS_OCL_NI_DNJ_0(a) +
		KRATOS_OCL_COMP_1(a_j) * KRATOS_OCL_NI_DNJ_1(a) +
		KRATOS_OCL_COMP_2(a_j) * KRATOS_OCL_NI_DNJ_2(a);

	// destination[l_comp] += aux_j * U_j[l_comp] - aux_i * U_i[l_comp]
	KRATOS_OCL_COMP_0(destination) += aux_j * KRATOS_OCL_COMP_0(U_j) - aux_i * KRATOS_OCL_COMP_0(U_i);
	KRATOS_OCL_COMP_1(destination) += aux_j * KRATOS_OCL_COMP_1(U_j) - aux_i * KRATOS_OCL_COMP_1(U_i);
	KRATOS_OCL_COMP_2(destination) += aux_j * KRATOS_OCL_COMP_2(U_j) - aux_i * KRATOS_OCL_COMP_2(U_i);

#endif

}

void Sub_ConvectiveContribution(__global EdgeType *a, __global VectorType *destination,
	__global const VectorType *a_i, __global const VectorType *U_i,
	__global const VectorType *a_j, __global const VectorType *U_j)
{

#ifdef USE_CONSERVATIVE_FORM_FOR_VECTOR_CONVECTION

	// temp += a_i[k_comp] * Ni_DNj[k_comp]
	ValueType temp =
		KRATOS_OCL_COMP_0(a_i) * KRATOS_OCL_NI_DNJ_0(a) +
		KRATOS_OCL_COMP_1(a_i) * KRATOS_OCL_NI_DNJ_1(a) +
		KRATOS_OCL_COMP_2(a_i) * KRATOS_OCL_NI_DNJ_2(a);

	// destination[l_comp] -= temp * (U_j[l_comp] - U_i[l_comp])
	KRATOS_OCL_COMP_0(destination) -= temp * (KRATOS_OCL_COMP_0(U_j) - KRATOS_OCL_COMP_0(U_i));
	KRATOS_OCL_COMP_1(destination) -= temp * (KRATOS_OCL_COMP_1(U_j) - KRATOS_OCL_COMP_1(U_i));
	KRATOS_OCL_COMP_2(destination) -= temp * (KRATOS_OCL_COMP_2(U_j) - KRATOS_OCL_COMP_2(U_i));

#else

	// aux_i += a_i[k_comp] * Ni_DNj[k_comp]
	ValueType aux_i =
		KRATOS_OCL_COMP_0(a_i) * KRATOS_OCL_NI_DNJ_0(a) +
		KRATOS_OCL_COMP_1(a_i) * KRATOS_OCL_NI_DNJ_1(a) +
		KRATOS_OCL_COMP_2(a_i) * KRATOS_OCL_NI_DNJ_2(a);

	// aux_j += a_j[k_comp] * Ni_DNj[k_comp]
	ValueType aux_j =
		KRATOS_OCL_COMP_0(a_j) * KRATOS_OCL_NI_DNJ_0(a) +
		KRATOS_OCL_COMP_1(a_j) * KRATOS_OCL_NI_DNJ_1(a) +
		KRATOS_OCL_COMP_2(a_j) * KRATOS_OCL_NI_DNJ_2(a);

	// destination[l_comp] -= aux_j * U_j[l_comp] - aux_i * U_i[l_comp]
	KRATOS_OCL_COMP_0(destination) -= aux_j * KRATOS_OCL_COMP_0(U_j) - aux_i * KRATOS_OCL_COMP_0(U_i);
	KRATOS_OCL_COMP_1(destination) -= aux_j * KRATOS_OCL_COMP_1(U_j) - aux_i * KRATOS_OCL_COMP_1(U_i);
	KRATOS_OCL_COMP_2(destination) -= aux_j * KRATOS_OCL_COMP_2(U_j) - aux_i * KRATOS_OCL_COMP_2(U_i);

#endif

}

void Add_ConvectiveContribution2(__global EdgeType *a, __global ValueType *destination,
	__global const VectorType *a_i, const ValueType phi_i,
	__global const VectorType *a_j, const ValueType phi_j)
{

#ifdef USE_CONSERVATIVE_FORM_FOR_SCALAR_CONVECTION

	// temp += a_i[k_comp] * Ni_DNj[k_comp]
	ValueType temp =
		KRATOS_OCL_COMP_0(a_i) * KRATOS_OCL_NI_DNJ_0(a) +
		KRATOS_OCL_COMP_1(a_i) * KRATOS_OCL_NI_DNJ_1(a) +
		KRATOS_OCL_COMP_2(a_i) * KRATOS_OCL_NI_DNJ_2(a);

	*destination += temp * (phi_j - phi_i);

#else

	// aux_i += a_i[k_comp] * Ni_DNj[k_comp]
	ValueType aux_i =
		KRATOS_OCL_COMP_0(a_i) * KRATOS_OCL_NI_DNJ_0(a) +
		KRATOS_OCL_COMP_1(a_i) * KRATOS_OCL_NI_DNJ_1(a) +
		KRATOS_OCL_COMP_2(a_i) * KRATOS_OCL_NI_DNJ_2(a);

	// aux_j += a_j[k_comp] * Ni_DNj[k_comp]
	ValueType aux_j =
		KRATOS_OCL_COMP_0(a_j) * KRATOS_OCL_NI_DNJ_0(a) +
		KRATOS_OCL_COMP_1(a_j) * KRATOS_OCL_NI_DNJ_1(a) +
		KRATOS_OCL_COMP_2(a_j) * KRATOS_OCL_NI_DNJ_2(a);

	*destination += aux_j * phi_j - aux_i * phi_i;

#endif

}

void Sub_ConvectiveContribution2(__global EdgeType *a, ValueType *destination,
	const VectorType *a_i, const ValueType phi_i,
	const VectorType *a_j, const ValueType phi_j)
{

#ifdef USE_CONSERVATIVE_FORM_FOR_SCALAR_CONVECTION

	// temp += a_i[k_comp] * Ni_DNj[k_comp]
	ValueType temp =
		KRATOS_OCL_COMP_0(a_i) * KRATOS_OCL_NI_DNJ_0(a) +
		KRATOS_OCL_COMP_1(a_i) * KRATOS_OCL_NI_DNJ_1(a) +
		KRATOS_OCL_COMP_2(a_i) * KRATOS_OCL_NI_DNJ_2(a);

	*destination -= temp * (phi_j - phi_i);

#else

	// aux_i += a_i[k_comp] * Ni_DNj[k_comp]
	ValueType aux_i =
		KRATOS_OCL_COMP_0(a_i) * KRATOS_OCL_NI_DNJ_0(a) +
		KRATOS_OCL_COMP_1(a_i) * KRATOS_OCL_NI_DNJ_1(a) +
		KRATOS_OCL_COMP_2(a_i) * KRATOS_OCL_NI_DNJ_2(a);

	// aux_j += a_j[k_comp] * Ni_DNj[k_comp]
	ValueType aux_j =
		KRATOS_OCL_COMP_0(a_j) * KRATOS_OCL_NI_DNJ_0(a) +
		KRATOS_OCL_COMP_1(a_j) * KRATOS_OCL_NI_DNJ_1(a) +
		KRATOS_OCL_COMP_2(a_j) * KRATOS_OCL_NI_DNJ_2(a);

	*destination -= aux_j * phi_j - aux_i * phi_i;

#endif

}

void CalculateConvectionStabilization_LOW(__global EdgeType *a, __global VectorType *stab_low,
	__global const VectorType *a_i, __global const VectorType *U_i,
	__global const VectorType *a_j, __global const VectorType *U_j)
{
	// conv_stab += a_i[k_comp] * a_i[m_comp] * LaplacianIJ(k_comp,m_comp)
	ValueType conv_stab =
		KRATOS_OCL_COMP_0(a_i) * KRATOS_OCL_COMP_0(a_i) * KRATOS_OCL_LAPLACIANIJ_0_0(a) +
		KRATOS_OCL_COMP_0(a_i) * KRATOS_OCL_COMP_1(a_i) * KRATOS_OCL_LAPLACIANIJ_0_1(a) +
		KRATOS_OCL_COMP_0(a_i) * KRATOS_OCL_COMP_2(a_i) * KRATOS_OCL_LAPLACIANIJ_0_2(a) +
		KRATOS_OCL_COMP_1(a_i) * KRATOS_OCL_COMP_0(a_i) * KRATOS_OCL_LAPLACIANIJ_1_0(a) +
		KRATOS_OCL_COMP_1(a_i) * KRATOS_OCL_COMP_1(a_i) * KRATOS_OCL_LAPLACIANIJ_1_1(a) +
		KRATOS_OCL_COMP_1(a_i) * KRATOS_OCL_COMP_2(a_i) * KRATOS_OCL_LAPLACIANIJ_1_2(a) +
		KRATOS_OCL_COMP_2(a_i) * KRATOS_OCL_COMP_0(a_i) * KRATOS_OCL_LAPLACIANIJ_2_0(a) +
		KRATOS_OCL_COMP_2(a_i) * KRATOS_OCL_COMP_1(a_i) * KRATOS_OCL_LAPLACIANIJ_2_1(a) +
		KRATOS_OCL_COMP_2(a_i) * KRATOS_OCL_COMP_2(a_i) * KRATOS_OCL_LAPLACIANIJ_2_2(a);

	// stab_low[l_comp] = conv_stab * (U_j[l_comp] - U_i[l_comp])
	KRATOS_OCL_COMP_0(stab_low) = conv_stab * (KRATOS_OCL_COMP_0(U_j) - KRATOS_OCL_COMP_0(U_j));
	KRATOS_OCL_COMP_1(stab_low) = conv_stab * (KRATOS_OCL_COMP_1(U_j) - KRATOS_OCL_COMP_1(U_j));
	KRATOS_OCL_COMP_2(stab_low) = conv_stab * (KRATOS_OCL_COMP_2(U_j) - KRATOS_OCL_COMP_2(U_j));
}

void CalculateConvectionStabilization_LOW2(__global EdgeType *a, ValueType *stab_low,
	const VectorType *a_i, const ValueType phi_i,
	const VectorType *a_j, const ValueType phi_j)
{
	// conv_stab += a_i[k_comp] * a_i[m_comp] * LaplacianIJ(k_comp,m_comp)
	ValueType conv_stab =
		KRATOS_OCL_COMP_0(a_i) * KRATOS_OCL_COMP_0(a_i) * KRATOS_OCL_LAPLACIANIJ_0_0(a) +
		KRATOS_OCL_COMP_0(a_i) * KRATOS_OCL_COMP_1(a_i) * KRATOS_OCL_LAPLACIANIJ_0_1(a) +
		KRATOS_OCL_COMP_0(a_i) * KRATOS_OCL_COMP_2(a_i) * KRATOS_OCL_LAPLACIANIJ_0_2(a) +
		KRATOS_OCL_COMP_1(a_i) * KRATOS_OCL_COMP_0(a_i) * KRATOS_OCL_LAPLACIANIJ_1_0(a) +
		KRATOS_OCL_COMP_1(a_i) * KRATOS_OCL_COMP_1(a_i) * KRATOS_OCL_LAPLACIANIJ_1_1(a) +
		KRATOS_OCL_COMP_1(a_i) * KRATOS_OCL_COMP_2(a_i) * KRATOS_OCL_LAPLACIANIJ_1_2(a) +
		KRATOS_OCL_COMP_2(a_i) * KRATOS_OCL_COMP_0(a_i) * KRATOS_OCL_LAPLACIANIJ_2_0(a) +
		KRATOS_OCL_COMP_2(a_i) * KRATOS_OCL_COMP_1(a_i) * KRATOS_OCL_LAPLACIANIJ_2_1(a) +
		KRATOS_OCL_COMP_2(a_i) * KRATOS_OCL_COMP_2(a_i) * KRATOS_OCL_LAPLACIANIJ_2_2(a);

	*stab_low = conv_stab * (phi_j - phi_i);
}

void CalculateConvectionStabilization_HIGH(__global EdgeType *a, __global VectorType *stab_high,
	__global const VectorType *a_i, __global const VectorType *pi_i,
	__global const VectorType *a_j, __global const VectorType *pi_j)
{

#ifdef USE_CONSERVATIVE_FORM_FOR_VECTOR_CONVECTION

	// temp += a_i[k_comp] * Ni_DNj[k_comp]
	ValueType temp =
		KRATOS_OCL_COMP_0(a_i) * KRATOS_OCL_NI_DNJ_0(a) +
		KRATOS_OCL_COMP_1(a_i) * KRATOS_OCL_NI_DNJ_1(a) +
		KRATOS_OCL_COMP_2(a_i) * KRATOS_OCL_NI_DNJ_2(a);

	// stab_high[l_comp] = -temp * (pi_j[l_comp] - pi_i[l_comp])
	KRATOS_OCL_COMP_0(stab_high) = -temp * (KRATOS_OCL_COMP_0(pi_j) - KRATOS_OCL_COMP_0(pi_i)); // check if the minus sign is correct
	KRATOS_OCL_COMP_1(stab_high) = -temp * (KRATOS_OCL_COMP_1(pi_j) - KRATOS_OCL_COMP_1(pi_i)); // check if the minus sign is correct
	KRATOS_OCL_COMP_2(stab_high) = -temp * (KRATOS_OCL_COMP_2(pi_j) - KRATOS_OCL_COMP_2(pi_i)); // check if the minus sign is correct

#else

	// aux_i += a_i[k_comp] * Ni_DNj[k_comp]
	ValueType aux_i =
		KRATOS_OCL_COMP_0(a_i) * KRATOS_OCL_NI_DNJ_0(a) +
		KRATOS_OCL_COMP_1(a_i) * KRATOS_OCL_NI_DNJ_1(a) +
		KRATOS_OCL_COMP_2(a_i) * KRATOS_OCL_NI_DNJ_2(a);

	// aux_j += a_j[k_comp] * Ni_DNj[k_comp]
	ValueType aux_j =
		KRATOS_OCL_COMP_0(a_j) * KRATOS_OCL_NI_DNJ_0(a) +
		KRATOS_OCL_COMP_1(a_j) * KRATOS_OCL_NI_DNJ_1(a) +
		KRATOS_OCL_COMP_2(a_j) * KRATOS_OCL_NI_DNJ_2(a);

	// stab_high[l_comp] =  -(aux_j * pi_j[l_comp] - aux_i * pi_i[l_comp])
	KRATOS_OCL_COMP_0(stab_high) = -(aux_j * KRATOS_OCL_COMP_0(pi_j) - aux_i * KRATOS_OCL_COMP_0(pi_i));
	KRATOS_OCL_COMP_1(stab_high) = -(aux_j * KRATOS_OCL_COMP_1(pi_j) - aux_i * KRATOS_OCL_COMP_1(pi_i));
	KRATOS_OCL_COMP_2(stab_high) = -(aux_j * KRATOS_OCL_COMP_2(pi_j) - aux_i * KRATOS_OCL_COMP_2(pi_i));

#endif

}

void CalculateConvectionStabilization_HIGH2(__global EdgeType *a, ValueType *stab_high,
	const VectorType *a_i, const ValueType pi_i,
	const VectorType *a_j, const ValueType pi_j)
{

#ifdef USE_CONSERVATIVE_FORM_FOR_VECTOR_CONVECTION

	// temp += a_i[k_comp] * Ni_DNj[k_comp]
	ValueType temp =
		KRATOS_OCL_COMP_0(a_i) * KRATOS_OCL_NI_DNJ_0(a) +
		KRATOS_OCL_COMP_1(a_i) * KRATOS_OCL_NI_DNJ_1(a) +
		KRATOS_OCL_COMP_2(a_i) * KRATOS_OCL_NI_DNJ_2(a);

	*stab_high  = -temp * (pi_j - pi_i); // check if the minus sign is correct

#else

	// aux_i += a_i[k_comp] * Ni_DNj[k_comp]
	ValueType aux_i =
		KRATOS_OCL_COMP_0(a_i) * KRATOS_OCL_NI_DNJ_0(a) +
		KRATOS_OCL_COMP_1(a_i) * KRATOS_OCL_NI_DNJ_1(a) +
		KRATOS_OCL_COMP_2(a_i) * KRATOS_OCL_NI_DNJ_2(a);

	// aux_j += a_j[k_comp] * Ni_DNj[k_comp]
	ValueType aux_j =
		KRATOS_OCL_COMP_0(a_j) * KRATOS_OCL_NI_DNJ_0(a) +
		KRATOS_OCL_COMP_1(a_j) * KRATOS_OCL_NI_DNJ_1(a) +
		KRATOS_OCL_COMP_2(a_j) * KRATOS_OCL_NI_DNJ_2(a);

	*stab_high = -(aux_j * pi_j - aux_i * pi_i);

#endif

}

void Add_StabContribution(__global EdgeType *a, __global VectorType *destination,
	const ValueType tau, const ValueType beta,
	__global const VectorType *stab_low, __global const VectorType *stab_high)
{
	// destination[l_comp] += tau * (stab_low[l_comp] - beta * stab_high[l_comp])
	KRATOS_OCL_COMP_0(destination) += tau * (KRATOS_OCL_COMP_0(stab_low) - beta * KRATOS_OCL_COMP_0(stab_high));
	KRATOS_OCL_COMP_1(destination) += tau * (KRATOS_OCL_COMP_1(stab_low) - beta * KRATOS_OCL_COMP_1(stab_high));
	KRATOS_OCL_COMP_2(destination) += tau * (KRATOS_OCL_COMP_2(stab_low) - beta * KRATOS_OCL_COMP_2(stab_high));
}

void Sub_StabContribution(__global EdgeType *a, __global VectorType *destination,
	const ValueType tau, const ValueType beta,
	__global const VectorType *stab_low, __global const VectorType *stab_high)
{
	// destination[l_comp] -= tau * (stab_low[l_comp] - beta * stab_high[l_comp])
	KRATOS_OCL_COMP_0(destination) -= tau * (KRATOS_OCL_COMP_0(stab_low) - beta * KRATOS_OCL_COMP_0(stab_high));
	KRATOS_OCL_COMP_1(destination) -= tau * (KRATOS_OCL_COMP_1(stab_low) - beta * KRATOS_OCL_COMP_1(stab_high));
	KRATOS_OCL_COMP_2(destination) -= tau * (KRATOS_OCL_COMP_2(stab_low) - beta * KRATOS_OCL_COMP_2(stab_high));
}

void Add_StabContribution2(__global EdgeType *a, __global ValueType *destination,
	const ValueType tau, const ValueType beta,
	const ValueType stab_low, const ValueType stab_high)
{
	*destination += tau * (stab_low - beta * stab_high);
}

void Sub_StabContribution2(__global EdgeType *a, ValueType *destination,
	const ValueType tau, const ValueType beta,
	const ValueType stab_low, const ValueType stab_high)
{
	*destination -= tau * (stab_low - beta * stab_high);
}

void Add_ViscousContribution(__global EdgeType *a, __global VectorType *destination,
	__global const VectorType *U_i, const ValueType nu_i,
	__global const VectorType *U_j, const ValueType nu_j)
{
	// L += LaplacianIJ(l_comp, l_comp)
	ValueType L = KRATOS_OCL_LAPLACIANIJ_0_0(a) + KRATOS_OCL_LAPLACIANIJ_1_1(a) + KRATOS_OCL_LAPLACIANIJ_2_2(a);

	ValueType nu_avg = 0.5 * (nu_i + nu_j);

	// destination[l_comp] += nu_i * L * (U_j[l_comp] - U_i[l_comp])
	KRATOS_OCL_COMP_0(destination) += nu_i * L * (KRATOS_OCL_COMP_0(U_j) - KRATOS_OCL_COMP_0(U_i));
	KRATOS_OCL_COMP_1(destination) += nu_i * L * (KRATOS_OCL_COMP_1(U_j) - KRATOS_OCL_COMP_1(U_i));
	KRATOS_OCL_COMP_2(destination) += nu_i * L * (KRATOS_OCL_COMP_2(U_j) - KRATOS_OCL_COMP_2(U_i));
}

