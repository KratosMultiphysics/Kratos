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
__kernel void Test(__global double *input, __global double *output, double offset)
{
	size_t id = get_global_id(0);
	output[id] = 2.00 * sin(input[id]) * cos(input[id]) + offset;
}

//
// Array3
//
// Helper function to access a 2D array using its 1D memory with _Dim2 = 3

inline __global double *Array3(__global double *_Array, unsigned int _i, unsigned int _j)
{
	return _Array + _i * 3 + _j;
}

//
// CalculateAdvectiveVelocity
//
// OpenCL version of PureConvectionEdgeBased::CalculateAdvectiveVelocity

__kernel void CalculateAdvectiveVelocity(__global const double *mUn, __global const double *mUn1, __global double *mA, const double coefficient, const unsigned int n_node)
{
	// Get work item index
	unsigned int i_node = get_global_id(0);

	// Check if we are in the range
	if (i_node > n_node)
	{
		return;
	}

	//
	// Original code:
	//
	// array_1d <double, TDim> &a_i = mA[i_node];
	// const array_1d<double, TDim>&  Un_i = mUn[i_node];
	// const array_1d<double, TDim>&  Un1_i = mUn1[i_node];
	//
	// for (unsigned int k_comp = 0; k_comp < TDim; k_comp++)
	//   a_i[k_comp] = coefficient * Un1_i[k_comp] + (1.0 - coefficient)* Un_i[k_comp];

	*Array3(mA, i_node, 0) = coefficient * *Array3(mUn1, i_node, 0) + (1.0 - coefficient) * *Array3(mUn, i_node, 0);
	*Array3(mA, i_node, 1) = coefficient * *Array3(mUn1, i_node, 1) + (1.0 - coefficient) * *Array3(mUn, i_node, 1);
	*Array3(mA, i_node, 2) = coefficient * *Array3(mUn1, i_node, 2) + (1.0 - coefficient) * *Array3(mUn, i_node, 2);
}
