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


// TODO: Optimize kernels (vectorize, use built-in functions, avoid excessive use of __global data)

//
// CalculateAdvectiveVelocity
//
// OpenCL version of PureConvectionEdgeBased::CalculateAdvectiveVelocity

__kernel void CalculateAdvectiveVelocity(__global const VectorType *mUn, __global const VectorType *mUn1, __global VectorType *mA, const ValueType coefficient, const IndexType n_nodes)
{
	// Get work item index
	IndexType i_node = get_global_id(0);

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

__kernel void Solve1(__global const ValueType *Hmin, __global const VectorType *A, __global ValueType *Tau, const ValueType time_inv, const IndexType n_nodes)
{
	// Get work item index
	IndexType i_node = get_global_id(0);

	// Check if we are in the range
	if (i_node < n_nodes)
	{
		// Tau[i_node] = 1.00 / (2.00 * length(A[i_node]) / Hmin[i_node] + 0.01 * time_inv);
		Tau[i_node] = KRATOS_OCL_NATIVE_RECIP(KRATOS_OCL_NATIVE_DIVIDE(2.00 * KRATOS_OCL_LENGTH3(A[i_node]), Hmin[i_node]) + 0.01 * time_inv);
	}
}

//
// CalculateRHS1
//
// Part of CalculateRHS

__kernel void CalculateRHS1(__global VectorType *Pi, __global const ValueType *phi, __global const IndexType *RowStartIndex, __global const IndexType *ColumnIndex, __global const EdgeType *EdgeValues, __global const ValueType *InvertedMass, const IndexType n_nodes)
{
	// Get work item index
	IndexType i_node = get_global_id(0);

	// Check if we are in the range
	if (i_node < n_nodes)
	{
		VectorType Temp_Pi_i_node = 0.00;
		ValueType Temp_Phi_i_node = phi[i_node];

		for (IndexType csr_index = RowStartIndex[i_node]; csr_index != RowStartIndex[i_node + 1]; csr_index++)
		{
			IndexType j_neighbour = ColumnIndex[csr_index];

			Add_grad_p(&EdgeValues[csr_index], &Temp_Pi_i_node, Temp_Phi_i_node, phi[j_neighbour]);
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

__kernel void CalculateRHS2(__global VectorType *Pi, __global const ValueType *phi, __global const IndexType *RowStartIndex, __global const IndexType *ColumnIndex, __global const VectorType *x, __global ValueType *Beta, const IndexType n_nodes)
{
	// Get work item index
	IndexType i_node = get_global_id(0);

	// Check if we are in the range
	if (i_node < n_nodes)
	{
		VectorType dir;

		ValueType Temp_Beta_i_node = 0.00;

		ValueType n = 0.00;
		ValueType h = 0.00;

		VectorType Temp_x_i_node = x[i_node];
		VectorType Temp_Pi_i_node = Pi[i_node];
		ValueType Temp_Phi_i_node = phi[i_node];

		for (IndexType csr_index = RowStartIndex[i_node]; csr_index != RowStartIndex[i_node + 1]; csr_index++)
		{
			IndexType j_neighbour = ColumnIndex[csr_index];

			ValueType Temp_Phi_j_neighbour = phi[j_neighbour];

			dir = x[j_neighbour] - Temp_x_i_node;

			ValueType proj = 0.5 * dot(dir, Temp_Pi_i_node + Pi[j_neighbour]);

			ValueType numerator = fabs(fabs(Temp_Phi_j_neighbour - Temp_Phi_i_node) - fabs(proj));
			ValueType denominator = fabs(fabs(Temp_Phi_j_neighbour - Temp_Phi_i_node) + 1e-6);

			ValueType beta = KRATOS_OCL_NATIVE_DIVIDE(numerator, denominator);

			Temp_Beta_i_node += beta;

			n += 1.0;
			h += KRATOS_OCL_LENGTH3(dir);
		}

		Temp_Beta_i_node /= n;
		h /= n;

		if (Temp_Beta_i_node > 1.00)
		{
			Temp_Beta_i_node = 1.00;
		}

		Temp_Beta_i_node *= h;
		Beta[i_node] = Temp_Beta_i_node;
	}
}

//
// CalculateRHS3
//
// Part of CalculateRHS

__kernel void CalculateRHS3(__global VectorType *Pi, __global const ValueType *phi, __global const IndexType *RowStartIndex, __global const IndexType *ColumnIndex, __global const EdgeType *EdgeValues, __global const VectorType *convective_velocity, __global const ValueType *Beta, __global ValueType *rhs, __global const ValueType *Tau, const IndexType n_nodes)
{
	// Get work item index
	IndexType i_node = get_global_id(0);

	// Check if we are in the range
	if (i_node < n_nodes)
	{
		ValueType stab_low;
		ValueType stab_high;

		VectorType Temp_convective_velocity_i_node = convective_velocity[i_node];
		ValueType Temp_convective_velocity_i_node_length = KRATOS_OCL_LENGTH3(Temp_convective_velocity_i_node);
		ValueType Temp_Phi_i_node = phi[i_node];
		ValueType Temp_Tau_i_node = Tau[i_node];
		ValueType Temp_Beta_i_node = Beta[i_node];
		ValueType Temp_rhs_i_node = 0.00;

		ValueType pi_i = dot(Pi[i_node], Temp_convective_velocity_i_node);

		for (IndexType csr_index = RowStartIndex[i_node]; csr_index != RowStartIndex[i_node + 1]; csr_index++)
		{
			IndexType j_neighbour = ColumnIndex[csr_index];

			VectorType Temp_convective_velocity_j_neighbour = convective_velocity[j_neighbour];
			ValueType Temp_Phi_j_neighbour = phi[j_neighbour];

			ValueType pi_j = dot(Pi[j_neighbour], Temp_convective_velocity_i_node);

			// Convection operator
			Sub_ConvectiveContribution2(&EdgeValues[csr_index], &Temp_rhs_i_node, &Temp_convective_velocity_i_node, Temp_Phi_i_node, &Temp_convective_velocity_j_neighbour, Temp_Phi_j_neighbour);

			// Calculate stabilization part
			CalculateConvectionStabilization_LOW2(&EdgeValues[csr_index], &stab_low, &Temp_convective_velocity_i_node, Temp_Phi_i_node, &Temp_convective_velocity_j_neighbour, Temp_Phi_j_neighbour);

			CalculateConvectionStabilization_HIGH2(&EdgeValues[csr_index], &stab_high, &Temp_convective_velocity_i_node, Temp_Phi_i_node, &Temp_convective_velocity_j_neighbour, Temp_Phi_j_neighbour);

			Sub_StabContribution2(&EdgeValues[csr_index], &Temp_rhs_i_node, Temp_Tau_i_node, 1.00, stab_low, stab_high);


			ValueType laplacian_ij = 0.00;

			CalculateScalarLaplacian(&EdgeValues[csr_index], &laplacian_ij);

			ValueType capturing = laplacian_ij * (Temp_Phi_j_neighbour - Temp_Phi_i_node);

			Temp_rhs_i_node -= 0.35 * capturing * Temp_Beta_i_node * Temp_convective_velocity_i_node_length;
			rhs[i_node] = Temp_rhs_i_node;
		}
	}
}
