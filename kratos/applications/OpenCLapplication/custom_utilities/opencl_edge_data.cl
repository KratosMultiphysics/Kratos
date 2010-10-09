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
// opencl_edge_data.cl
//
// OpenCL kernels and functions used in opencl_edge_data.h

// Enable OpenCL extension
#pragma OPENCL EXTENSION cl_amd_fp64: enable

// We suggest defining the following macro
#define USE_CONSERVATIVE_FORM_FOR_SCALAR_CONVECTION

// We suggest defining the following macro
#define USE_CONSERVATIVE_FORM_FOR_VECTOR_CONVECTION

// Some useful macros, will be renamed if not consistent
#define KRATOS_OCL_LAPLACIANIJ_0_0(a)	(*a).s0
#define KRATOS_OCL_LAPLACIANIJ_0_1(a)	(*a).s1
#define KRATOS_OCL_LAPLACIANIJ_0_2(a)	(*a).s2
#define KRATOS_OCL_LAPLACIANIJ_1_0(a)	(*a).s3
#define KRATOS_OCL_LAPLACIANIJ_1_1(a)	(*a).s4
#define KRATOS_OCL_LAPLACIANIJ_1_2(a)	(*a).s5
#define KRATOS_OCL_LAPLACIANIJ_2_0(a)	(*a).s6
#define KRATOS_OCL_LAPLACIANIJ_2_1(a)	(*a).s7
#define KRATOS_OCL_LAPLACIANIJ_2_2(a)	(*a).s8

#define KRATOS_OCL_MASS(a)				(*a).s9

#define KRATOS_OCL_NI_DNJ_0(a)			(*a).sa
#define KRATOS_OCL_NI_DNJ_1(a)			(*a).sb
#define KRATOS_OCL_NI_DNJ_2(a)			(*a).sc

#define KRATOS_OCL_DNI_NJ_0(a)			(*a).sd
#define KRATOS_OCL_DNI_NJ_1(a)			(*a).se
#define KRATOS_OCL_DNI_NJ_2(a)			(*a).sf

#define KRATOS_OCL_COMP(a, n)			(*a).s[n]

#define KRATOS_OCL_COMP_0(a)			(*a).x
#define KRATOS_OCL_COMP_1(a)			(*a).y
#define KRATOS_OCL_COMP_2(a)			(*a).z


// A dummy kernel for test
__kernel void Test(__global const double *input, __global double *output, const double offset)
{
	const size_t id = get_global_id(0);
	const double iv = input[id];
	output[id] = 2.00 * sin(iv) * cos(iv) + offset;
}

// Common kernels

__kernel void SetToZero(__global double *Vector, const unsigned int n)
{
	// Get work item index
	const size_t id = get_global_id(0);

	// Check if we are in the range
	if (id < n)
	{
		Vector[id] = 0.00;
	}
}

__kernel void Add_Minv_value(__global double *DestinationVector, __global const double *Origin1Vector, const double Value, __global const double *MinvVector, __global const double *OriginVector, const unsigned int n)
{
	// Get work item index
	const size_t id = get_global_id(0);

	// Check if we are in the range
	if (id < n)
	{
		DestinationVector[id] = Origin1Vector[id] + Value * MinvVector[id] * OriginVector[id];
	}
}

// Edge specific kernels

void Add_Gp(double16 *a, double4 *destination, const double p_i, const double p_j)
{
	 // destination[comp] -= Ni_DNj[comp] * p_j - DNi_Nj[comp] * p_i
	 KRATOS_OCL_COMP_0(destination) -= KRATOS_OCL_NI_DNJ_0(a) * p_j - KRATOS_OCL_DNI_NJ_0(a) * p_i;
	 KRATOS_OCL_COMP_1(destination) -= KRATOS_OCL_NI_DNJ_1(a) * p_j - KRATOS_OCL_DNI_NJ_1(a) * p_i;
	 KRATOS_OCL_COMP_2(destination) -= KRATOS_OCL_NI_DNJ_2(a) * p_j - KRATOS_OCL_DNI_NJ_2(a) * p_i;
}

void Sub_Gp(double16 *a, double4 *destination, const double p_i, const double p_j)
{
	 // destination[comp] += Ni_DNj[comp] * p_j - DNi_Nj[comp] * p_i
	 KRATOS_OCL_COMP_0(destination) += KRATOS_OCL_NI_DNJ_0(a) * p_j - KRATOS_OCL_DNI_NJ_0(a) * p_i;
	 KRATOS_OCL_COMP_1(destination) += KRATOS_OCL_NI_DNJ_1(a) * p_j - KRATOS_OCL_DNI_NJ_1(a) * p_i;
	 KRATOS_OCL_COMP_2(destination) += KRATOS_OCL_NI_DNJ_2(a) * p_j - KRATOS_OCL_DNI_NJ_2(a) * p_i;
}

void Add_D_v(double16 *a, double *destination, const double4 *v_i, const double4 *v_j)
{
	// destination += Ni_DNj[comp] * (v_j[comp] - v_i[comp])
	*destination +=
		KRATOS_OCL_NI_DNJ_0(a) * (KRATOS_OCL_COMP_0(v_j) - KRATOS_OCL_COMP_0(v_i)) +
		KRATOS_OCL_NI_DNJ_1(a) * (KRATOS_OCL_COMP_1(v_j) - KRATOS_OCL_COMP_1(v_i)) +
		KRATOS_OCL_NI_DNJ_2(a) * (KRATOS_OCL_COMP_2(v_j) - KRATOS_OCL_COMP_2(v_i));
}

void Sub_D_v(double16 *a, double *destination, const double4 *v_i, const double4 *v_j)
{
	// destination -= Ni_DNj[comp] * (v_j[comp] - v_i[comp])
	*destination -=
		KRATOS_OCL_NI_DNJ_0(a) * (KRATOS_OCL_COMP_0(v_j) - KRATOS_OCL_COMP_0(v_i)) +
		KRATOS_OCL_NI_DNJ_1(a) * (KRATOS_OCL_COMP_1(v_j) - KRATOS_OCL_COMP_1(v_i)) +
		KRATOS_OCL_NI_DNJ_2(a) * (KRATOS_OCL_COMP_2(v_j) - KRATOS_OCL_COMP_2(v_i));
}

void Add_grad_p(double16 *a, double4 *destination, const double p_i, const double p_j)
{
	// destination[comp] += Ni_DNj[comp] * (p_j - p_i)
	double dp = p_j-p_i;

	KRATOS_OCL_COMP_0(destination) += KRATOS_OCL_NI_DNJ_0(a) * dp;
	KRATOS_OCL_COMP_1(destination) += KRATOS_OCL_NI_DNJ_1(a) * dp;
	KRATOS_OCL_COMP_2(destination) += KRATOS_OCL_NI_DNJ_2(a) * dp;
}

void Sub_grad_p(double16 *a, double4 *destination, const double p_i, const double p_j)
{
	// destination[comp] -= Ni_DNj[comp] * (p_j - p_i)
	double dp = p_j-p_i;

	KRATOS_OCL_COMP_0(destination) -= KRATOS_OCL_NI_DNJ_0(a) * dp;
	KRATOS_OCL_COMP_1(destination) -= KRATOS_OCL_NI_DNJ_1(a) * dp;
	KRATOS_OCL_COMP_2(destination) -= KRATOS_OCL_NI_DNJ_2(a) * dp;
}

void Add_div_v(double16 *a, double *destination, const double4 *v_i, const double4 *v_j)
{
	// destination -= Ni_DNj[comp]*v_j[comp] - DNi_Nj[comp]*v_i[comp]
	*destination -=
		KRATOS_OCL_NI_DNJ_0(a) * KRATOS_OCL_COMP_0(v_j) - KRATOS_OCL_DNI_NJ_0(a) * KRATOS_OCL_COMP_0(v_i) +
		KRATOS_OCL_NI_DNJ_1(a) * KRATOS_OCL_COMP_1(v_j) - KRATOS_OCL_DNI_NJ_1(a) * KRATOS_OCL_COMP_1(v_i) +
		KRATOS_OCL_NI_DNJ_2(a) * KRATOS_OCL_COMP_2(v_j) - KRATOS_OCL_DNI_NJ_2(a) * KRATOS_OCL_COMP_2(v_i);
}

void Sub_div_v(double16 *a, double *destination, const double4 *v_i, const double4 *v_j)
{
	// destination += Ni_DNj[comp]*v_j[comp] - DNi_Nj[comp]*v_i[comp]
	*destination +=
		KRATOS_OCL_NI_DNJ_0(a) * KRATOS_OCL_COMP_0(v_j) - KRATOS_OCL_DNI_NJ_0(a) * KRATOS_OCL_COMP_0(v_i) +
		KRATOS_OCL_NI_DNJ_1(a) * KRATOS_OCL_COMP_1(v_j) - KRATOS_OCL_DNI_NJ_1(a) * KRATOS_OCL_COMP_1(v_i) +
		KRATOS_OCL_NI_DNJ_2(a) * KRATOS_OCL_COMP_2(v_j) - KRATOS_OCL_DNI_NJ_2(a) * KRATOS_OCL_COMP_2(v_i);
}

void CalculateScalarLaplacian(double16 *a, double *l_ij)
{
	// l_ij += LaplacianIJ(comp, comp)
	*l_ij = KRATOS_OCL_LAPLACIANIJ_0_0(a) + KRATOS_OCL_LAPLACIANIJ_1_1(a) + KRATOS_OCL_LAPLACIANIJ_2_2(a);
}

void Add_ConvectiveContribution(double16 *a, double4 *destination,
	const double4 *a_i, const double4 *U_i,
	const double4 *a_j, const double4 *U_j)
{

#ifdef USE_CONSERVATIVE_FORM_FOR_VECTOR_CONVECTION

	// temp += a_i[k_comp] * Ni_DNj[k_comp]
	double temp =
		KRATOS_OCL_COMP_0(a_i) * KRATOS_OCL_NI_DNJ_0(a) +
		KRATOS_OCL_COMP_1(a_i) * KRATOS_OCL_NI_DNJ_1(a) +
		KRATOS_OCL_COMP_2(a_i) * KRATOS_OCL_NI_DNJ_2(a);

	// destination[l_comp] += temp * (U_j[l_comp] - U_i[l_comp])
	KRATOS_OCL_COMP_0(destination) += temp * (KRATOS_OCL_COMP_0(U_j) - KRATOS_OCL_COMP_0(U_i));
	KRATOS_OCL_COMP_1(destination) += temp * (KRATOS_OCL_COMP_1(U_j) - KRATOS_OCL_COMP_1(U_i));
	KRATOS_OCL_COMP_2(destination) += temp * (KRATOS_OCL_COMP_2(U_j) - KRATOS_OCL_COMP_2(U_i));

#else

	// aux_i += a_i[k_comp] * Ni_DNj[k_comp]
	double aux_i =
		KRATOS_OCL_COMP_0(a_i) * KRATOS_OCL_NI_DNJ_0(a) +
		KRATOS_OCL_COMP_1(a_i) * KRATOS_OCL_NI_DNJ_1(a) +
		KRATOS_OCL_COMP_2(a_i) * KRATOS_OCL_NI_DNJ_2(a);

	// aux_j += a_j[k_comp] * Ni_DNj[k_comp]
	double aux_j =
		KRATOS_OCL_COMP_0(a_j) * KRATOS_OCL_NI_DNJ_0(a) +
		KRATOS_OCL_COMP_1(a_j) * KRATOS_OCL_NI_DNJ_1(a) +
		KRATOS_OCL_COMP_2(a_j) * KRATOS_OCL_NI_DNJ_2(a);

	// destination[l_comp] += aux_j * U_j[l_comp] - aux_i * U_i[l_comp]
	KRATOS_OCL_COMP_0(destination) += aux_j * KRATOS_OCL_COMP_0(U_j) - aux_i * KRATOS_OCL_COMP_0(U_i);
	KRATOS_OCL_COMP_1(destination) += aux_j * KRATOS_OCL_COMP_1(U_j) - aux_i * KRATOS_OCL_COMP_1(U_i);
	KRATOS_OCL_COMP_2(destination) += aux_j * KRATOS_OCL_COMP_2(U_j) - aux_i * KRATOS_OCL_COMP_2(U_i);

#endif

}

void Sub_ConvectiveContribution(double16 *a, double4 *destination,
	const double4 *a_i, const double4 *U_i,
	const double4 *a_j, const double4 *U_j)
{

#ifdef USE_CONSERVATIVE_FORM_FOR_VECTOR_CONVECTION

	// temp += a_i[k_comp] * Ni_DNj[k_comp]
	double temp =
		KRATOS_OCL_COMP_0(a_i) * KRATOS_OCL_NI_DNJ_0(a) +
		KRATOS_OCL_COMP_1(a_i) * KRATOS_OCL_NI_DNJ_1(a) +
		KRATOS_OCL_COMP_2(a_i) * KRATOS_OCL_NI_DNJ_2(a);

	// destination[l_comp] -= temp * (U_j[l_comp] - U_i[l_comp])
	KRATOS_OCL_COMP_0(destination) -= temp * (KRATOS_OCL_COMP_0(U_j) - KRATOS_OCL_COMP_0(U_i));
	KRATOS_OCL_COMP_1(destination) -= temp * (KRATOS_OCL_COMP_1(U_j) - KRATOS_OCL_COMP_1(U_i));
	KRATOS_OCL_COMP_2(destination) -= temp * (KRATOS_OCL_COMP_2(U_j) - KRATOS_OCL_COMP_2(U_i));

#else

	// aux_i += a_i[k_comp] * Ni_DNj[k_comp]
	double aux_i =
		KRATOS_OCL_COMP_0(a_i) * KRATOS_OCL_NI_DNJ_0(a) +
		KRATOS_OCL_COMP_1(a_i) * KRATOS_OCL_NI_DNJ_1(a) +
		KRATOS_OCL_COMP_2(a_i) * KRATOS_OCL_NI_DNJ_2(a);

	// aux_j += a_j[k_comp] * Ni_DNj[k_comp]
	double aux_j =
		KRATOS_OCL_COMP_0(a_j) * KRATOS_OCL_NI_DNJ_0(a) +
		KRATOS_OCL_COMP_1(a_j) * KRATOS_OCL_NI_DNJ_1(a) +
		KRATOS_OCL_COMP_2(a_j) * KRATOS_OCL_NI_DNJ_2(a);

	// destination[l_comp] -= aux_j * U_j[l_comp] - aux_i * U_i[l_comp]
	KRATOS_OCL_COMP_0(destination) -= aux_j * KRATOS_OCL_COMP_0(U_j) - aux_i * KRATOS_OCL_COMP_0(U_i);
	KRATOS_OCL_COMP_1(destination) -= aux_j * KRATOS_OCL_COMP_1(U_j) - aux_i * KRATOS_OCL_COMP_1(U_i);
	KRATOS_OCL_COMP_2(destination) -= aux_j * KRATOS_OCL_COMP_2(U_j) - aux_i * KRATOS_OCL_COMP_2(U_i);

#endif

}

void Add_ConvectiveContribution2(double16 *a, double *destination,
	const double4 *a_i, const double phi_i,
	const double4 *a_j, const double phi_j)
{

#ifdef USE_CONSERVATIVE_FORM_FOR_SCALAR_CONVECTION

	// temp += a_i[k_comp] * Ni_DNj[k_comp]
	double temp =
		KRATOS_OCL_COMP_0(a_i) * KRATOS_OCL_NI_DNJ_0(a) +
		KRATOS_OCL_COMP_1(a_i) * KRATOS_OCL_NI_DNJ_1(a) +
		KRATOS_OCL_COMP_2(a_i) * KRATOS_OCL_NI_DNJ_2(a);

	*destination += temp * (phi_j - phi_i);

#else

	// aux_i += a_i[k_comp] * Ni_DNj[k_comp]
	double aux_i =
		KRATOS_OCL_COMP_0(a_i) * KRATOS_OCL_NI_DNJ_0(a) +
		KRATOS_OCL_COMP_1(a_i) * KRATOS_OCL_NI_DNJ_1(a) +
		KRATOS_OCL_COMP_2(a_i) * KRATOS_OCL_NI_DNJ_2(a);

	// aux_j += a_j[k_comp] * Ni_DNj[k_comp]
	double aux_j =
		KRATOS_OCL_COMP_0(a_j) * KRATOS_OCL_NI_DNJ_0(a) +
		KRATOS_OCL_COMP_1(a_j) * KRATOS_OCL_NI_DNJ_1(a) +
		KRATOS_OCL_COMP_2(a_j) * KRATOS_OCL_NI_DNJ_2(a);

	*destination += aux_j * phi_j - aux_i * phi_i;

#endif

}

void Sub_ConvectiveContribution2(double16 *a, double *destination,
	const double4 *a_i, const double phi_i,
	const double4 *a_j, const double phi_j)
{

#ifdef USE_CONSERVATIVE_FORM_FOR_SCALAR_CONVECTION

	// temp += a_i[k_comp] * Ni_DNj[k_comp]
	double temp =
		KRATOS_OCL_COMP_0(a_i) * KRATOS_OCL_NI_DNJ_0(a) +
		KRATOS_OCL_COMP_1(a_i) * KRATOS_OCL_NI_DNJ_1(a) +
		KRATOS_OCL_COMP_2(a_i) * KRATOS_OCL_NI_DNJ_2(a);

	*destination -= temp * (phi_j - phi_i);

#else

	// aux_i += a_i[k_comp] * Ni_DNj[k_comp]
	double aux_i =
		KRATOS_OCL_COMP_0(a_i) * KRATOS_OCL_NI_DNJ_0(a) +
		KRATOS_OCL_COMP_1(a_i) * KRATOS_OCL_NI_DNJ_1(a) +
		KRATOS_OCL_COMP_2(a_i) * KRATOS_OCL_NI_DNJ_2(a);

	// aux_j += a_j[k_comp] * Ni_DNj[k_comp]
	double aux_j =
		KRATOS_OCL_COMP_0(a_j) * KRATOS_OCL_NI_DNJ_0(a) +
		KRATOS_OCL_COMP_1(a_j) * KRATOS_OCL_NI_DNJ_1(a) +
		KRATOS_OCL_COMP_2(a_j) * KRATOS_OCL_NI_DNJ_2(a);

	*destination -= aux_j * phi_j - aux_i * phi_i;

#endif

}

void CalculateConvectionStabilization_LOW(double16 *a, double4 *stab_low,
	const double4 *a_i, const double4 *U_i,
	const double4 *a_j, const double4 *U_j)
{
	// conv_stab += a_i[k_comp] * a_i[m_comp] * LaplacianIJ(k_comp,m_comp)
	double conv_stab =
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

void CalculateConvectionStabilization_LOW2(double16 *a, double *stab_low,
	const double4 *a_i, const double phi_i,
	const double4 *a_j, const double phi_j)
{
	// conv_stab += a_i[k_comp] * a_i[m_comp] * LaplacianIJ(k_comp,m_comp)
	double conv_stab =
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

void CalculateConvectionStabilization_HIGH(double16 *a, double4 *stab_high,
	const double4 *a_i, const double4 *pi_i,
	const double4 *a_j, const double4 *pi_j)
{

#ifdef USE_CONSERVATIVE_FORM_FOR_VECTOR_CONVECTION

	// temp += a_i[k_comp] * Ni_DNj[k_comp]
	double temp =
		KRATOS_OCL_COMP_0(a_i) * KRATOS_OCL_NI_DNJ_0(a) +
		KRATOS_OCL_COMP_1(a_i) * KRATOS_OCL_NI_DNJ_1(a) +
		KRATOS_OCL_COMP_2(a_i) * KRATOS_OCL_NI_DNJ_2(a);

	// stab_high[l_comp] = -temp * (pi_j[l_comp] - pi_i[l_comp])
	KRATOS_OCL_COMP_0(stab_high) = -temp * (KRATOS_OCL_COMP_0(pi_j) - KRATOS_OCL_COMP_0(pi_i)); // check if the minus sign is correct
	KRATOS_OCL_COMP_1(stab_high) = -temp * (KRATOS_OCL_COMP_1(pi_j) - KRATOS_OCL_COMP_1(pi_i)); // check if the minus sign is correct
	KRATOS_OCL_COMP_2(stab_high) = -temp * (KRATOS_OCL_COMP_2(pi_j) - KRATOS_OCL_COMP_2(pi_i)); // check if the minus sign is correct

#else

	// aux_i += a_i[k_comp] * Ni_DNj[k_comp]
	double aux_i =
		KRATOS_OCL_COMP_0(a_i) * KRATOS_OCL_NI_DNJ_0(a) +
		KRATOS_OCL_COMP_1(a_i) * KRATOS_OCL_NI_DNJ_1(a) +
		KRATOS_OCL_COMP_2(a_i) * KRATOS_OCL_NI_DNJ_2(a);

	// aux_j += a_j[k_comp] * Ni_DNj[k_comp]
	double aux_j =
		KRATOS_OCL_COMP_0(a_j) * KRATOS_OCL_NI_DNJ_0(a) +
		KRATOS_OCL_COMP_1(a_j) * KRATOS_OCL_NI_DNJ_1(a) +
		KRATOS_OCL_COMP_2(a_j) * KRATOS_OCL_NI_DNJ_2(a);

	// stab_high[l_comp] =  -(aux_j * pi_j[l_comp] - aux_i * pi_i[l_comp])
	KRATOS_OCL_COMP_0(stab_high) = -(aux_j * KRATOS_OCL_COMP_0(pi_j) - aux_i * KRATOS_OCL_COMP_0(pi_i));
	KRATOS_OCL_COMP_1(stab_high) = -(aux_j * KRATOS_OCL_COMP_1(pi_j) - aux_i * KRATOS_OCL_COMP_1(pi_i));
	KRATOS_OCL_COMP_2(stab_high) = -(aux_j * KRATOS_OCL_COMP_2(pi_j) - aux_i * KRATOS_OCL_COMP_2(pi_i));

#endif

}

void CalculateConvectionStabilization_HIGH2(double16 *a, double *stab_high,
	const double4 *a_i, const double pi_i,
	const double4 *a_j, const double pi_j)
{

#ifdef USE_CONSERVATIVE_FORM_FOR_VECTOR_CONVECTION

	// temp += a_i[k_comp] * Ni_DNj[k_comp]
	double temp =
		KRATOS_OCL_COMP_0(a_i) * KRATOS_OCL_NI_DNJ_0(a) +
		KRATOS_OCL_COMP_1(a_i) * KRATOS_OCL_NI_DNJ_1(a) +
		KRATOS_OCL_COMP_2(a_i) * KRATOS_OCL_NI_DNJ_2(a);

	*stab_high  = -temp * (pi_j - pi_i); // check if the minus sign is correct

#else

	// aux_i += a_i[k_comp] * Ni_DNj[k_comp]
	double aux_i =
		KRATOS_OCL_COMP_0(a_i) * KRATOS_OCL_NI_DNJ_0(a) +
		KRATOS_OCL_COMP_1(a_i) * KRATOS_OCL_NI_DNJ_1(a) +
		KRATOS_OCL_COMP_2(a_i) * KRATOS_OCL_NI_DNJ_2(a);

	// aux_j += a_j[k_comp] * Ni_DNj[k_comp]
	double aux_j =
		KRATOS_OCL_COMP_0(a_j) * KRATOS_OCL_NI_DNJ_0(a) +
		KRATOS_OCL_COMP_1(a_j) * KRATOS_OCL_NI_DNJ_1(a) +
		KRATOS_OCL_COMP_2(a_j) * KRATOS_OCL_NI_DNJ_2(a);

	*stab_high = -(aux_j * pi_j- aux_i * pi_i);

#endif

}

void Add_StabContribution(double16 *a, double4 *destination,
	const double tau, const double beta,
	const double4 *stab_low, const double4 *stab_high)
{
	// destination[l_comp] += tau * (stab_low[l_comp] - beta * stab_high[l_comp])
	KRATOS_OCL_COMP_0(destination) += tau * (KRATOS_OCL_COMP_0(stab_low) - beta * KRATOS_OCL_COMP_0(stab_high));
	KRATOS_OCL_COMP_1(destination) += tau * (KRATOS_OCL_COMP_1(stab_low) - beta * KRATOS_OCL_COMP_1(stab_high));
	KRATOS_OCL_COMP_2(destination) += tau * (KRATOS_OCL_COMP_2(stab_low) - beta * KRATOS_OCL_COMP_2(stab_high));
}

void Sub_StabContribution(double16 *a, double4 *destination,
	const double tau, const double beta,
	const double4 *stab_low, const double4 *stab_high)
{
	// destination[l_comp] -= tau * (stab_low[l_comp] - beta * stab_high[l_comp])
	KRATOS_OCL_COMP_0(destination) -= tau * (KRATOS_OCL_COMP_0(stab_low) - beta * KRATOS_OCL_COMP_0(stab_high));
	KRATOS_OCL_COMP_1(destination) -= tau * (KRATOS_OCL_COMP_1(stab_low) - beta * KRATOS_OCL_COMP_1(stab_high));
	KRATOS_OCL_COMP_2(destination) -= tau * (KRATOS_OCL_COMP_2(stab_low) - beta * KRATOS_OCL_COMP_2(stab_high));
}

void Add_StabContribution2(double16 *a, double *destination,
	const double tau, const double beta,
	const double stab_low, const double stab_high)
{
	*destination += tau * (stab_low - beta * stab_high);
}

void Sub_StabContribution2(double16 *a, double *destination,
	const double tau, const double beta,
	const double stab_low, const double stab_high)
{
	*destination -= tau * (stab_low - beta * stab_high);
}

void Add_ViscousContribution(double16 *a, double4 *destination,
	const double4 *U_i, const double nu_i,
	const double4 *U_j, const double nu_j)
{
	// L += LaplacianIJ(l_comp, l_comp)
	double L = KRATOS_OCL_LAPLACIANIJ_0_0(a) + KRATOS_OCL_LAPLACIANIJ_1_1(a) + KRATOS_OCL_LAPLACIANIJ_2_2(a);

	double nu_avg = 0.5 * (nu_i + nu_j);

	// destination[l_comp] += nu_i * L * (U_j[l_comp] - U_i[l_comp])
	KRATOS_OCL_COMP_0(destination) += nu_i * L * (KRATOS_OCL_COMP_0(U_j) - KRATOS_OCL_COMP_0(U_i));
	KRATOS_OCL_COMP_1(destination) += nu_i * L * (KRATOS_OCL_COMP_1(U_j) - KRATOS_OCL_COMP_1(U_i));
	KRATOS_OCL_COMP_2(destination) += nu_i * L * (KRATOS_OCL_COMP_2(U_j) - KRATOS_OCL_COMP_2(U_i));
}
