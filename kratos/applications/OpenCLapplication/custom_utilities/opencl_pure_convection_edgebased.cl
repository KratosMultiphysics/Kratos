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
// OpenCL kernels and functions used in opencl_pure_convection_edgebased.h


#include "opencl_edge_data_common.cl"


//
// CalculateAdvectiveVelocity
//
// OpenCL version of PureConvectionEdgeBased::CalculateAdvectiveVelocity

__kernel void CalculateAdvectiveVelocity(__global VectorType *mUn, __global VectorType *mUn1, __global VectorType *mA, const ValueType coefficient, const IndexType n_nodes)
{
	// Get work item index
	const size_t i_node = get_global_id(0);

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
	const size_t i_node = get_global_id(0);

	// Check if we are in the range
	if (i_node < n_nodes)
	{
		// TODO: Use mad()
		// Tau[i_node] = 1.00 / (2.00 * length(A[i_node]) / Hmin[i_node] + 0.01 * time_inv);
		Tau[i_node] = KRATOS_OCL_RECIP(KRATOS_OCL_DIVIDE(2.00 * KRATOS_OCL_LENGTH3(A[i_node]), Hmin[i_node]) + 0.01 * time_inv);
	}
}

//
// CalculateRHS1
//
// Part of CalculateRHS

__kernel void CalculateRHS1(__global VectorType *Pi, __global ValueType *phi, __global IndexType *RowStartIndex, __global IndexType *ColumnIndex, __read_only image2d_t EdgeValues, __global ValueType *InvertedMass, const IndexType n_nodes, __local IndexType *Bounds)
{
	// Get work item index
	const size_t i_node = get_global_id(0);
	const size_t i_thread = get_local_id(0);

	// Reading for loop bounds

	if (i_thread == 0)
	{
		Bounds[0] = RowStartIndex[i_node];
	}

	barrier(CLK_LOCAL_MEM_FENCE);

	Bounds[i_thread + 1] = RowStartIndex[i_node + 1];

	barrier(CLK_LOCAL_MEM_FENCE);

	// Check if we are in the range
	if (i_node < n_nodes)
	{
		VectorType Temp_Pi_i_node = 0.00;
		ValueType Temp_Phi_i_node = phi[i_node];

		for (IndexType csr_index = Bounds[i_thread]; csr_index != Bounds[i_thread + 1]; csr_index++)
		{
			//EdgeType CurrentEdge = EdgeValues[csr_index];
			EdgeType CurrentEdge = ReadDouble16FromDouble16Image(EdgeValues, csr_index);

			//VectorType Ni_DNj = KRATOS_OCL_VECTOR3(KRATOS_OCL_NI_DNJ_0(EdgeValues[csr_index]), KRATOS_OCL_NI_DNJ_1(EdgeValues[csr_index]), KRATOS_OCL_NI_DNJ_2(EdgeValues[csr_index]));
			VectorType Ni_DNj = KRATOS_OCL_VECTOR3(KRATOS_OCL_NI_DNJ_0(CurrentEdge), KRATOS_OCL_NI_DNJ_1(CurrentEdge), KRATOS_OCL_NI_DNJ_2(CurrentEdge));

			Add_grad_p(Ni_DNj, &Temp_Pi_i_node, Temp_Phi_i_node, phi[ColumnIndex[csr_index]]);
		}

		// Apply inverted mass matrix
		Pi[i_node] = Temp_Pi_i_node * InvertedMass[i_node];
	}
}

//
// CalculateRHS2
//
// Part of CalculateRHS

__kernel void CalculateRHS2(__global VectorType *Pi, __global ValueType *phi, __global IndexType *RowStartIndex, __global IndexType *ColumnIndex, __global VectorType *x, __global ValueType *Beta, const IndexType n_nodes, __local IndexType *Bounds)
{
	// Get work item index
	const size_t i_node = get_global_id(0);
	const size_t i_thread = get_local_id(0);

	// Reading for loop bounds

	if (i_thread == 0)
	{
		Bounds[0] = RowStartIndex[i_node];
	}

	barrier(CLK_LOCAL_MEM_FENCE);

	Bounds[i_thread + 1] = RowStartIndex[i_node + 1];

	barrier(CLK_LOCAL_MEM_FENCE);

	// Check if we are in the range
	if (i_node < n_nodes)
	{
		ValueType Temp_Beta_i_node = 0.00;

		ValueType h = 0.00;
		ValueType n = 0.00;

		VectorType Temp_x_i_node = x[i_node];
		VectorType Temp_Pi_i_node = Pi[i_node];
		ValueType Temp_Phi_i_node = phi[i_node];

		for (IndexType csr_index = Bounds[i_thread]; csr_index != Bounds[i_thread + 1]; csr_index++)
		{
			IndexType j_neighbour = ColumnIndex[csr_index];
			VectorType dir = x[j_neighbour] - Temp_x_i_node;

			h += KRATOS_OCL_LENGTH3(dir);
			n += 1.00;

			ValueType Temp_1 = fabs(phi[j_neighbour] - Temp_Phi_i_node);

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

__kernel void CalculateRHS3(__global VectorType *Pi, __global ValueType *phi, __global IndexType *RowStartIndex, __global IndexType *ColumnIndex, __read_only image2d_t EdgeValues, __global VectorType *convective_velocity, __global ValueType *Beta, __global ValueType *rhs, __global ValueType *Tau, const IndexType n_nodes, __local IndexType *Bounds)
{
	// Get work item index
	const size_t i_node = get_global_id(0);
	const size_t i_thread = get_local_id(0);

	// Reading for loop bounds

	if (i_thread == 0)
	{
		Bounds[0] = RowStartIndex[i_node];
	}

	barrier(CLK_LOCAL_MEM_FENCE);

	Bounds[i_thread + 1] = RowStartIndex[i_node + 1];

	barrier(CLK_LOCAL_MEM_FENCE);

	// Check if we are in the range
	if (i_node < n_nodes)
	{
		VectorType Temp_convective_velocity_i_node = convective_velocity[i_node];
		ValueType Temp_Phi_i_node = phi[i_node];
		ValueType Temp_Tau_i_node = Tau[i_node];
		ValueType Temp_Beta_i_node = Beta[i_node];
		ValueType Temp_rhs_i_node = 0.00;

		ValueType pi_i = dot(Pi[i_node], Temp_convective_velocity_i_node);

		for (IndexType csr_index = Bounds[i_thread]; csr_index != Bounds[i_thread + 1]; csr_index++)
		{
			IndexType j_neighbour = ColumnIndex[csr_index];

			ValueType Temp_Phi_j_neighbour = phi[j_neighbour];

			ValueType stab_low;
			ValueType stab_high;

			//EdgeType CurrentEdge = EdgeValues[csr_index];
			EdgeType CurrentEdge = ReadDouble16FromDouble16Image(EdgeValues, csr_index);

			//VectorType Ni_DNj = KRATOS_OCL_VECTOR3(KRATOS_OCL_NI_DNJ_0(EdgeValues[csr_index]), KRATOS_OCL_NI_DNJ_1(EdgeValues[csr_index]), KRATOS_OCL_NI_DNJ_2(EdgeValues[csr_index]));
			VectorType Ni_DNj = KRATOS_OCL_VECTOR3(KRATOS_OCL_NI_DNJ_0(CurrentEdge), KRATOS_OCL_NI_DNJ_1(CurrentEdge), KRATOS_OCL_NI_DNJ_2(CurrentEdge));
			//ReadVectorFromDouble16Image(EdgeValues, csr_index, 5, &Ni_DNj);

			Sub_ConvectiveContribution2(Ni_DNj, &Temp_rhs_i_node, Temp_convective_velocity_i_node, Temp_Phi_i_node, Temp_Phi_j_neighbour);
			CalculateConvectionStabilization_HIGH2(Ni_DNj, &stab_high, Temp_convective_velocity_i_node, Temp_Phi_i_node, Temp_Phi_j_neighbour);

			//VectorType Lij0 = KRATOS_OCL_VECTOR3(KRATOS_OCL_LAPLACIANIJ_0_0(EdgeValues[csr_index]), KRATOS_OCL_LAPLACIANIJ_0_1(EdgeValues[csr_index]), KRATOS_OCL_LAPLACIANIJ_0_2(EdgeValues[csr_index]));
			//VectorType Lij1 = KRATOS_OCL_VECTOR3(KRATOS_OCL_LAPLACIANIJ_1_0(EdgeValues[csr_index]), KRATOS_OCL_LAPLACIANIJ_1_1(EdgeValues[csr_index]), KRATOS_OCL_LAPLACIANIJ_1_2(EdgeValues[csr_index]));
			//VectorType Lij2 = KRATOS_OCL_VECTOR3(KRATOS_OCL_LAPLACIANIJ_2_0(EdgeValues[csr_index]), KRATOS_OCL_LAPLACIANIJ_2_1(EdgeValues[csr_index]), KRATOS_OCL_LAPLACIANIJ_2_2(EdgeValues[csr_index]));
			VectorType Lij0 = KRATOS_OCL_VECTOR3(KRATOS_OCL_LAPLACIANIJ_0_0(CurrentEdge), KRATOS_OCL_LAPLACIANIJ_0_1(CurrentEdge), KRATOS_OCL_LAPLACIANIJ_0_2(CurrentEdge));
			VectorType Lij1 = KRATOS_OCL_VECTOR3(KRATOS_OCL_LAPLACIANIJ_1_0(CurrentEdge), KRATOS_OCL_LAPLACIANIJ_1_1(CurrentEdge), KRATOS_OCL_LAPLACIANIJ_1_2(CurrentEdge));
			VectorType Lij2 = KRATOS_OCL_VECTOR3(KRATOS_OCL_LAPLACIANIJ_2_0(CurrentEdge), KRATOS_OCL_LAPLACIANIJ_2_1(CurrentEdge), KRATOS_OCL_LAPLACIANIJ_2_2(CurrentEdge));
			//Read3VectorFromDouble16Image(EdgeValues, csr_index, 0, &Lij0, &Lij1, &Lij2);

			ValueType laplacian_ij;
			CalculateScalarLaplacian(Lij0.x, Lij1.y, Lij2.z, &laplacian_ij);

			CalculateConvectionStabilization_LOW2(Lij0, Lij1, Lij2, &stab_low, Temp_convective_velocity_i_node, Temp_Phi_i_node, Temp_Phi_j_neighbour);
			Sub_StabContribution2(&Temp_rhs_i_node, Temp_Tau_i_node, 1.00, stab_low, stab_high);

			Temp_rhs_i_node -= 0.35 * laplacian_ij * (Temp_Phi_j_neighbour - Temp_Phi_i_node) * Temp_Beta_i_node * KRATOS_OCL_LENGTH3(Temp_convective_velocity_i_node);
		}

		rhs[i_node] = Temp_rhs_i_node;
	}
}
