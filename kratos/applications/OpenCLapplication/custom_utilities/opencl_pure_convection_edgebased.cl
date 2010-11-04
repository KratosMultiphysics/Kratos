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


#include "opencl_edge_data_common.cl"

/*

Currently used functions:

Add_grad_p
Sub_ConvectiveContribution2
CalculateConvectionStabilization_LOW2
CalculateConvectionStabilization_HIGH2
Sub_StabContribution2
CalculateScalarLaplacian

*/

// TODO: Optimize kernels (vectorize, use built-in functions, avoid excessive use of __global data)

//
// CalculateAdvectiveVelocity
//
// OpenCL version of PureConvectionEdgeBased::CalculateAdvectiveVelocity

__kernel void CalculateAdvectiveVelocity(__global VectorType *mUn, __global VectorType *mUn1, __global VectorType *mA, const ValueType coefficient, const IndexType n_nodes)
{
	// Get work item index
	__private const size_t i_node = get_global_id(0);

	// Check if we are in the range
	if (i_node < n_nodes)
	{
		// mA[i_node] = coefficient * mUn1[i_node] + (1.00 - coefficient) * mUn[i_node];
		mA[i_node] = mix(mUn[i_node], mUn1[i_node], coefficient);
	}
}

//
// Solve1
//
// Part of Solve()

__kernel void Solve1(__global ValueType *Hmin, __global VectorType *A, __global ValueType *Tau, const ValueType time_inv, const IndexType n_nodes)
{
	// Get work item index
	__private const size_t i_node = get_global_id(0);

	// Check if we are in the range
	if (i_node < n_nodes)
	{
		// Tau[i_node] = 1.00 / (2.00 * length(A[i_node]) / Hmin[i_node] + 0.01 * time_inv);
		Tau[i_node] = KRATOS_OCL_RECIP(KRATOS_OCL_DIVIDE(2.00 * KRATOS_OCL_LENGTH3(A[i_node]), Hmin[i_node]) + 0.01 * time_inv);
	}
}

//
// CalculateRHS1
//
// Part of CalculateRHS

__kernel void CalculateRHS1(__global VectorType *Pi, __global ValueType *phi, __global IndexType *RowStartIndex, __global IndexType *ColumnIndex, __global EdgeType *EdgeValues, __global ValueType *InvertedMass, const IndexType n_nodes)
{
	// Get work item index
	__private const size_t i_node = get_global_id(0);

	// Check if we are in the range
	if (i_node < n_nodes)
	{
		__private VectorType Temp_Pi_i_node = 0.00;
		__private ValueType Temp_Phi_i_node = phi[i_node];

		for (IndexType csr_index = RowStartIndex[i_node]; csr_index != RowStartIndex[i_node + 1]; csr_index++)
		{
			EdgeType CurrentEdge = EdgeValues[csr_index];
			Add_grad_p(KRATOS_OCL_VECTOR3(KRATOS_OCL_NI_DNJ_0(CurrentEdge), KRATOS_OCL_NI_DNJ_1(CurrentEdge), KRATOS_OCL_NI_DNJ_2(CurrentEdge)), &Temp_Pi_i_node, Temp_Phi_i_node, phi[ColumnIndex[csr_index]]);
		}

		// Apply inverted mass matrix
		Temp_Pi_i_node *= InvertedMass[i_node];
		Pi[i_node] = Temp_Pi_i_node;
	}
}

//
// CalculateRHS2
//
// Part of CalculateRHS

__kernel void CalculateRHS2(__global VectorType *Pi, __global ValueType *phi, __global IndexType *RowStartIndex, __global IndexType *ColumnIndex, __global VectorType *x, __global ValueType *Beta, const IndexType n_nodes)
{
	// Get work item index
	__private const size_t i_node = get_global_id(0);

	// Check if we are in the range
	if (i_node < n_nodes)
	{
		__private ValueType Temp_Beta_i_node = 0.00;

		__private ValueType h = 0.00;
		__private ValueType n = 0.00;

		__private VectorType Temp_x_i_node = x[i_node];
		__private VectorType Temp_Pi_i_node = Pi[i_node];
		__private ValueType Temp_Phi_i_node = phi[i_node];

		for (IndexType csr_index = RowStartIndex[i_node]; csr_index != RowStartIndex[i_node + 1]; csr_index++)
		{
			__private IndexType j_neighbour = ColumnIndex[csr_index];
			__private VectorType dir = x[j_neighbour] - Temp_x_i_node;

			h += KRATOS_OCL_LENGTH3(dir);
			n += 1.00;

			__private ValueType Temp_1 = fabs(phi[j_neighbour] - Temp_Phi_i_node);

			Temp_Beta_i_node += fabs(KRATOS_OCL_DIVIDE(Temp_1 - fabs(0.5 * dot(dir, Temp_Pi_i_node + Pi[j_neighbour])), Temp_1 + 1e-6));
		}

		Temp_Beta_i_node /= n;
		h /= n;

		Beta[i_node] = Temp_Beta_i_node > 1.00 ? h : h * Temp_Beta_i_node;
	}
}

//
// CalculateRHS3
//
// Part of CalculateRHS

__kernel void CalculateRHS3(__global VectorType *Pi, __global ValueType *phi, __global IndexType *RowStartIndex, __global IndexType *ColumnIndex, __global EdgeType *EdgeValues, __global VectorType *convective_velocity, __global ValueType *Beta, __global ValueType *rhs, __global ValueType *Tau, const IndexType n_nodes)
{
	// Get work item index
	__private const size_t i_node = get_global_id(0);

	// Check if we are in the range
	if (i_node < n_nodes)
	{
		__private VectorType Temp_convective_velocity_i_node = convective_velocity[i_node];
		__private ValueType Temp_Phi_i_node = phi[i_node];
		__private ValueType Temp_Tau_i_node = Tau[i_node];
		__private ValueType Temp_Beta_i_node = Beta[i_node];
		__private ValueType Temp_rhs_i_node = 0.00;

		__private ValueType pi_i = dot(Pi[i_node], Temp_convective_velocity_i_node);

		for (IndexType csr_index = RowStartIndex[i_node]; csr_index != RowStartIndex[i_node + 1]; csr_index++)
		{
			__private IndexType j_neighbour = ColumnIndex[csr_index];

			__private VectorType Temp_convective_velocity_j_neighbour = convective_velocity[j_neighbour];
			__private ValueType Temp_Phi_j_neighbour = phi[j_neighbour];

			EdgeType CurrentEdge = EdgeValues[csr_index];

			// Convection operator
			Sub_ConvectiveContribution2(KRATOS_OCL_VECTOR3(KRATOS_OCL_NI_DNJ_0(CurrentEdge), KRATOS_OCL_NI_DNJ_1(CurrentEdge), KRATOS_OCL_NI_DNJ_2(CurrentEdge)), &Temp_rhs_i_node, Temp_convective_velocity_i_node, Temp_Phi_i_node, Temp_convective_velocity_j_neighbour, Temp_Phi_j_neighbour);

			__private ValueType stab_low;
			__private ValueType stab_high;

			// Calculate stabilization part
			CalculateConvectionStabilization_LOW2(
				KRATOS_OCL_VECTOR3(KRATOS_OCL_LAPLACIANIJ_0_0(CurrentEdge), KRATOS_OCL_LAPLACIANIJ_0_1(CurrentEdge), KRATOS_OCL_LAPLACIANIJ_0_2(CurrentEdge)),
				KRATOS_OCL_VECTOR3(KRATOS_OCL_LAPLACIANIJ_1_0(CurrentEdge), KRATOS_OCL_LAPLACIANIJ_1_1(CurrentEdge), KRATOS_OCL_LAPLACIANIJ_1_2(CurrentEdge)),
				KRATOS_OCL_VECTOR3(KRATOS_OCL_LAPLACIANIJ_2_0(CurrentEdge), KRATOS_OCL_LAPLACIANIJ_2_1(CurrentEdge), KRATOS_OCL_LAPLACIANIJ_2_2(CurrentEdge)),
				&stab_low, Temp_convective_velocity_i_node, Temp_Phi_i_node, Temp_convective_velocity_j_neighbour, Temp_Phi_j_neighbour);

			CalculateConvectionStabilization_HIGH2(KRATOS_OCL_VECTOR3(KRATOS_OCL_NI_DNJ_0(CurrentEdge), KRATOS_OCL_NI_DNJ_1(CurrentEdge), KRATOS_OCL_NI_DNJ_2(CurrentEdge)), &stab_high, Temp_convective_velocity_i_node, Temp_Phi_i_node, Temp_convective_velocity_j_neighbour, Temp_Phi_j_neighbour);

			Sub_StabContribution2(0, &Temp_rhs_i_node, Temp_Tau_i_node, 1.00, stab_low, stab_high);


			__private ValueType laplacian_ij;

			CalculateScalarLaplacian(KRATOS_OCL_LAPLACIANIJ_0_0(CurrentEdge), KRATOS_OCL_LAPLACIANIJ_1_1(CurrentEdge), KRATOS_OCL_LAPLACIANIJ_2_2(CurrentEdge), &laplacian_ij);

			Temp_rhs_i_node -= 0.35 * laplacian_ij * (Temp_Phi_j_neighbour - Temp_Phi_i_node) * Temp_Beta_i_node * KRATOS_OCL_LENGTH3(Temp_convective_velocity_i_node);
		}

		rhs[i_node] = Temp_rhs_i_node;
	}
}
