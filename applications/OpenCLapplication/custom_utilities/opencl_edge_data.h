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


#if !defined(KRATOS_OPENCL_EDGE_DATA_H_INCLUDED)
#define KRATOS_OPENCL_EDGE_DATA_H_INCLUDED


// System includes
#include <string>
#include <iostream>
#include <algorithm>

// External includes


// Project includes
//#include "includes/define.h"
//#include "spaces/ublas_space.h"
#include "includes/model_part.h"
//#include "includes/node.h"
//#include "geometries/geometry.h"
#include "utilities/geometry_utilities.h"
//#include "incompressible_fluid_application/custom_utilities/edge_data.h"
#include "opencl_interface.h"

// Some useful macros, will be renamed if not consistent
#define KRATOS_OCL_LAPLACIANIJ_0_0(a)	a.s[0]
#define KRATOS_OCL_LAPLACIANIJ_0_1(a)	a.s[1]
#define KRATOS_OCL_LAPLACIANIJ_0_2(a)	a.s[2]
#define KRATOS_OCL_LAPLACIANIJ_1_0(a)	a.s[3]
#define KRATOS_OCL_LAPLACIANIJ_1_1(a)	a.s[4]
#define KRATOS_OCL_LAPLACIANIJ_1_2(a)	a.s[5]
#define KRATOS_OCL_LAPLACIANIJ_2_0(a)	a.s[6]
#define KRATOS_OCL_LAPLACIANIJ_2_1(a)	a.s[7]
#define KRATOS_OCL_LAPLACIANIJ_2_2(a)	a.s[8]

#define KRATOS_OCL_MASS(a)				a.s[9]

#define KRATOS_OCL_NI_DNJ_0(a)			a.s[10]
#define KRATOS_OCL_NI_DNJ_1(a)			a.s[11]
#define KRATOS_OCL_NI_DNJ_2(a)			a.s[12]

#define KRATOS_OCL_DNI_NJ_0(a)			a.s[13]
#define KRATOS_OCL_DNI_NJ_1(a)			a.s[14]
#define KRATOS_OCL_DNI_NJ_2(a)			a.s[15]

#define KRATOS_OCL_COMP_0(a)			a.s[0]
#define KRATOS_OCL_COMP_1(a)			a.s[1]
#define KRATOS_OCL_COMP_2(a)			a.s[2]


namespace Kratos
{

	//
	// AllocateArray
	//
	// Helper function to allocate an array
	// Pass the address of the array variable

	template <typename _Type> void AllocateArray(_Type **_Array, unsigned int _Size)
	{
		// Allocate memory
		*_Array = new _Type[_Size];
	}

	//
	// FreeArray
	//
	// Helper function to free an array
	// Pass the address of the array variable

	template <typename _Type> void FreeArray(_Type **_Array)
	{
		// Free memory
		delete [] *_Array;
	}

	//
	// Sqr
	//
	// Helper function to return square of a number

	inline double Sqr(double _X)
	{
		return _X * _X;
	}

	//
	// Norm2_3
	//
	// Helper function to calculate norm 2 of a vector with _Dim2 = 3

	inline double Norm2_3(cl_double3 X)
	{
		return sqrt(Sqr(KRATOS_OCL_COMP_0(X)) + Sqr(KRATOS_OCL_COMP_1(X)) + Sqr(KRATOS_OCL_COMP_2(X)));
	}

	//
	// OpenCLMatrixContainer
	//
	// A helper class for edge based mathods

	class OpenCLMatrixContainer
	{
		public:

			//
			// Used types

			typedef cl_double16 *EdgesVectorType;
			typedef cl_uint *IndicesVectorType;
			typedef cl_double3 *CalcVectorType;
			typedef cl_double *ValuesVectorType;

			//typedef UblasSpace <double, CompressedMatrix, Vector> SparseSpaceType;
			//typedef MatrixContainer <3, SparseSpaceType> MatrixContainerType;

			//
			// OpenCLMatrixContainer
			//
			// Constructor

			OpenCLMatrixContainer(cl_device_type device_type): mDeviceGroup(device_type, true)
			{
				// Loading OpenCL program
				// TODO: Add optimization flags here
				mpOpenCLEdgeData = mDeviceGroup.BuildProgramFromFile("opencl_edge_data.cl");

				// Register kernels
				mkAdd_Minv_value1 = mDeviceGroup.RegisterKernel(mpOpenCLEdgeData, "Add_Minv_value1");
				mkAdd_Minv_value3 = mDeviceGroup.RegisterKernel(mpOpenCLEdgeData, "Add_Minv_value3");
				mkSetToZero = mDeviceGroup.RegisterKernel(mpOpenCLEdgeData, "SetToZero");
			}

			//
			// ~OpenCLMatrixContainer
			//
			// Destructor

			~OpenCLMatrixContainer()
			{
				// Nothing to do!
			}

			//
			// Functions to return private values

			inline unsigned int GetNumberEdges()
			{
				return mNumberEdges;
			}

			inline EdgesVectorType GetEdgeValues()
			{
				return mNonzeroEdgeValues;
			}

			inline IndicesVectorType GetColumnIndex()
			{
				return mColumnIndex;
			}

			inline IndicesVectorType GetRowStartIndex()
			{
				return mRowStartIndex;
			}

			inline ValuesVectorType GetLumpedMass()
			{
				return mLumpedMassMatrix;
			}

			inline ValuesVectorType GetInvertedMass()
			{
				return mInvertedMassMatrix;
			}

			inline CalcVectorType GetDiagGradient()
			{
				return mDiagGradientMatrix;
			}

			inline ValuesVectorType GetHmin()
			{
				return mHmin;
			}

			// OpenCL stuff

			inline OpenCL::DeviceGroup &GetDeviceGroup()
			{
				return mDeviceGroup;
			}

			inline cl_uint GetEdgeValuesBuffer()
			{
				return mbNonzeroEdgeValues;
			}

			inline cl_uint GetColumnIndexBuffer()
			{
				return mbColumnIndex;
			}

			inline cl_uint GetRowStartIndexBuffer()
			{
				return mbRowStartIndex;
			}

			inline cl_uint GetLumpedMassBuffer()
			{
				return mbLumpedMassMatrix;
			}

			inline cl_uint GetInvertedMassBuffer()
			{
				return mbInvertedMassMatrix;
			}

			inline cl_uint GetDiagGradientBuffer()
			{
				return mbDiagGradientMatrix;
			}

			inline cl_uint GetHminBuffer()
			{
				return mbHmin;
			}

			//
			// ConstructCSRVector
			//
			// Function to size and initialize the vector of CSR tuples

			void ConstructCSRVector(ModelPart &model_part)
			{
				KRATOS_TRY

				// Remark: No colouring algorithm is used here (symmetry is neglected)
				//         respectively edge ij is considered different from edge ji

				// Defining the number of nodes and edges
				unsigned int n_nodes = model_part.Nodes().size();
				mNumberEdges = 0;

				// Counter to assign and get global nodal index
				unsigned int i_node = 0;

				// Counting the edges connecting the nodes
				for (ModelPart::NodesContainerType::iterator node_it = model_part.NodesBegin(); node_it != model_part.NodesEnd(); node_it++)
				{
					// Counting neighbours of each node
					mNumberEdges += (node_it -> GetValue(NEIGHBOUR_NODES)).size();

					// Assigning global index to each node
					node_it -> FastGetSolutionStepValue(AUX_INDEX) = static_cast <double> (i_node++);
				}

				// Error message in case number of nodes does not coincide with number of indices
				if (i_node != n_nodes)
				{
					KRATOS_WATCH("ERROR - Highest nodal index doesn't coincide with number of nodes!");
				}

				// Allocating memory for block of CSR data
				// TODO: Can this be optimized?
				AllocateArray(&mNonzeroEdgeValues, mNumberEdges);

				AllocateArray(&mColumnIndex, mNumberEdges);
				AllocateArray(&mRowStartIndex, n_nodes + 1);

				AllocateArray(&mLumpedMassMatrix, n_nodes);
				AllocateArray(&mInvertedMassMatrix, n_nodes);
				AllocateArray(&mHmin, n_nodes);

				AllocateArray(&mDiagGradientMatrix, n_nodes);

				// Allocating buffers on OpenCL device
				// TODO: Should these be read-only?
				// TODO: Try using Page-Locked memory, mapping, one chunk of memory, etc.
				mbNonzeroEdgeValues = mDeviceGroup.CreateBuffer(mNumberEdges * sizeof(cl_double16), CL_MEM_READ_ONLY);

				mbColumnIndex = mDeviceGroup.CreateBuffer(mNumberEdges * sizeof(cl_uint), CL_MEM_READ_ONLY);
				mbRowStartIndex = mDeviceGroup.CreateBuffer((n_nodes + 1) * sizeof(cl_uint), CL_MEM_READ_ONLY);

				mbLumpedMassMatrix = mDeviceGroup.CreateBuffer(n_nodes * sizeof(cl_double), CL_MEM_READ_ONLY);
				mbInvertedMassMatrix = mDeviceGroup.CreateBuffer(n_nodes * sizeof(cl_double), CL_MEM_READ_ONLY);
				mbHmin = mDeviceGroup.CreateBuffer(n_nodes * sizeof(cl_double), CL_MEM_READ_ONLY);

				mbDiagGradientMatrix = mDeviceGroup.CreateBuffer(n_nodes * sizeof(cl_double3), CL_MEM_READ_ONLY);

				// Initializing the CSR vector

				// Temporary variable as the row start index of a node depends on the number of neighbours of the previous one
				unsigned int row_start_temp = 0;

				// Main loop over all nodes
				for (ModelPart::NodesContainerType::iterator node_it = model_part.NodesBegin(); node_it != model_part.NodesEnd(); node_it++)
				{
					// Getting the global index of the node
					i_node = static_cast <unsigned int> (node_it -> FastGetSolutionStepValue(AUX_INDEX));

					// Determining its neighbours
					WeakPointerVector < Node<3> > &neighb_nodes = node_it -> GetValue(NEIGHBOUR_NODES);

					// Number of neighbours of node i determines row start index for the following node
					unsigned int n_neighbours = neighb_nodes.size();

					// Reserving memory for work array
					std::vector<unsigned int> work_array;
					work_array.reserve(n_neighbours);

					// Nested loop over the neighbouring nodes
					for (WeakPointerVector < Node<3> >::iterator neighb_it = neighb_nodes.begin(); neighb_it != neighb_nodes.end(); neighb_it++)
					{
						// Getting global index of the neighbouring node
						work_array.push_back(static_cast <unsigned int> (neighb_it -> FastGetSolutionStepValue(AUX_INDEX)));
					}

					// Reordering neighbours following their global indices
					std::sort(work_array.begin(), work_array.end());

					// Setting current row start index
					mRowStartIndex[i_node] = row_start_temp;

					// Nested loop over the by now ordered neighbours
					for (unsigned int counter = 0; counter < n_neighbours; counter++)
					{
						// Getting global index of the neighbouring node
						unsigned int j_neighbour = work_array[counter];

						// Calculating CSR index
						unsigned int csr_index = mRowStartIndex[i_node] + counter;

						// Saving column index j of the original matrix
						mColumnIndex[csr_index] = j_neighbour;

						// Initializing the CSR vector entries with zero
						KRATOS_OCL_LAPLACIANIJ_0_0(mNonzeroEdgeValues[csr_index]) = 0.00;
						KRATOS_OCL_LAPLACIANIJ_0_1(mNonzeroEdgeValues[csr_index]) = 0.00;
						KRATOS_OCL_LAPLACIANIJ_0_2(mNonzeroEdgeValues[csr_index]) = 0.00;
						KRATOS_OCL_LAPLACIANIJ_1_0(mNonzeroEdgeValues[csr_index]) = 0.00;
						KRATOS_OCL_LAPLACIANIJ_1_1(mNonzeroEdgeValues[csr_index]) = 0.00;
						KRATOS_OCL_LAPLACIANIJ_1_2(mNonzeroEdgeValues[csr_index]) = 0.00;
						KRATOS_OCL_LAPLACIANIJ_2_0(mNonzeroEdgeValues[csr_index]) = 0.00;
						KRATOS_OCL_LAPLACIANIJ_2_1(mNonzeroEdgeValues[csr_index]) = 0.00;
						KRATOS_OCL_LAPLACIANIJ_2_2(mNonzeroEdgeValues[csr_index]) = 0.00;

						KRATOS_OCL_MASS(mNonzeroEdgeValues[csr_index]) = 0.00;

						KRATOS_OCL_NI_DNJ_0(mNonzeroEdgeValues[csr_index]) = 0.00;
						KRATOS_OCL_NI_DNJ_1(mNonzeroEdgeValues[csr_index]) = 0.00;
						KRATOS_OCL_NI_DNJ_2(mNonzeroEdgeValues[csr_index]) = 0.00;

						KRATOS_OCL_DNI_NJ_0(mNonzeroEdgeValues[csr_index]) = 0.00;
						KRATOS_OCL_DNI_NJ_1(mNonzeroEdgeValues[csr_index]) = 0.00;
						KRATOS_OCL_DNI_NJ_2(mNonzeroEdgeValues[csr_index]) = 0.00;
					}

					// Preparing row start index for next node
					row_start_temp += n_neighbours;
				}

				// Adding last entry (necessary for abort criterion of loops)
				mRowStartIndex[n_nodes] = mNumberEdges;

				// Initializing node based values

				// Lumped mass matrix (elements Mi)
				for (i_node = 0; i_node < n_nodes; i_node++)
				{
					mLumpedMassMatrix[i_node] = 0.00;
				}

				// Set the heights to a huge number
				for (i_node = 0; i_node<n_nodes; i_node++)
				{
					mHmin[i_node] = 1e10;
				}

				// Diagonal of gradient matrix (elements Gii)
				for (i_node = 0; i_node < n_nodes; i_node++)
				{
					KRATOS_OCL_COMP_0(mDiagGradientMatrix[i_node]) = 0.00;
					KRATOS_OCL_COMP_1(mDiagGradientMatrix[i_node]) = 0.00;
					KRATOS_OCL_COMP_2(mDiagGradientMatrix[i_node]) = 0.00;
				}

				KRATOS_CATCH("")
			}

			//
			// BuildCSRData
			//
			// Function to precalculate CSR data

			void BuildCSRData(ModelPart &model_part)
			{
				KRATOS_TRY

				// Pre-calculating CSR data

				//Defining temporary local variables for elementwise addition

				// Shape functions
				array_1d <double, 3 + 1> N;

				// Shape function derivatives
				boost::numeric::ublas::bounded_matrix <double, 3 + 1, 3> dN_dx;

				// Volume
				double volume;

				// Weighting factor
				double weighting_factor = 1.0 / static_cast <double> (3 + 1);

				// Elemental matrices
				boost::numeric::ublas::bounded_matrix <double, 3 + 1, 3 + 1> mass_consistent;

				array_1d <double, 3 + 1> mass_lumped;

				// Global indices of elemental nodes
				array_1d <unsigned int, 3 + 1> nodal_indices;

				array_1d <double, 3 + 1> heights;

				// Loop over all elements
				for (ModelPart::ElementsContainerType::iterator elem_it = model_part.ElementsBegin(); elem_it != model_part.ElementsEnd(); elem_it++)
				{
					// Local element-wise calculations

					// Getting geometry data of the element
					GeometryUtils::CalculateGeometryData(elem_it -> GetGeometry(), dN_dx, N, volume);

					// Calculate lenght of the heights of the element
					for (unsigned int ie_node = 0; ie_node <= 3; ie_node++)
					{
						heights[ie_node] = dN_dx(ie_node, 0) * dN_dx(ie_node, 0);

						for (unsigned int comp = 1; comp < 3; comp++)
						{
							heights[ie_node] += dN_dx(ie_node, comp) * dN_dx(ie_node, comp);
						}

						heights[ie_node] = 1.0 / sqrt(heights[ie_node]);
					}


					// Setting up elemental mass matrices
					CalculateMassMatrix(mass_consistent, volume);

					noalias(mass_lumped) = ZeroVector(3 + 1);

					for (unsigned int ie_node = 0; ie_node <= 3; ie_node++)
					{
						for (unsigned int je_node = 0; je_node <= 3; je_node++)
						{
							mass_lumped[ie_node] += mass_consistent(ie_node, je_node);
						}
					}

					// Corresponding to Ni * dOmega respectively Nj * dOmega
					double weighted_volume = volume * weighting_factor;

					// Assembling global data structure

					// Loop over the nodes of the element to determine their global indices
					for (unsigned int ie_node = 0; ie_node <= 3; ie_node++)
					{
						nodal_indices[ie_node] = static_cast <unsigned int> (elem_it -> GetGeometry()[ie_node].FastGetSolutionStepValue(AUX_INDEX));
					}

					// Assembling global "edge matrices" by adding local contributions
					for (unsigned int ie_node = 0; ie_node <= 3; ie_node++)
					{
						// Check the heights and change the value if minimal is found
						if (mHmin[nodal_indices[ie_node]] > heights[ie_node])
							mHmin[nodal_indices[ie_node]] = heights[ie_node];

						for (unsigned int je_node = 0; je_node <= 3; je_node++)
						{
							// Remark: There is no edge linking node i with itself!
							// Diagonal terms
							if (ie_node != je_node)
							{
								// Calculating CSR index from global index
								unsigned int csr_index = GetCSRIndex(nodal_indices[ie_node], nodal_indices[je_node]);

								// Assigning precalculated element data to the referring edges

								// Contribution to edge mass
								KRATOS_OCL_MASS(mNonzeroEdgeValues[csr_index]) += mass_consistent(ie_node, je_node);


								// Contribution to edge laplacian
								KRATOS_OCL_LAPLACIANIJ_0_0(mNonzeroEdgeValues[csr_index]) += dN_dx(ie_node, 0) * dN_dx(je_node, 0) * volume;
								KRATOS_OCL_LAPLACIANIJ_0_1(mNonzeroEdgeValues[csr_index]) += dN_dx(ie_node, 0) * dN_dx(je_node, 1) * volume;
								KRATOS_OCL_LAPLACIANIJ_0_2(mNonzeroEdgeValues[csr_index]) += dN_dx(ie_node, 0) * dN_dx(je_node, 2) * volume;
								KRATOS_OCL_LAPLACIANIJ_1_0(mNonzeroEdgeValues[csr_index]) += dN_dx(ie_node, 1) * dN_dx(je_node, 0) * volume;
								KRATOS_OCL_LAPLACIANIJ_1_1(mNonzeroEdgeValues[csr_index]) += dN_dx(ie_node, 1) * dN_dx(je_node, 1) * volume;
								KRATOS_OCL_LAPLACIANIJ_1_2(mNonzeroEdgeValues[csr_index]) += dN_dx(ie_node, 1) * dN_dx(je_node, 2) * volume;
								KRATOS_OCL_LAPLACIANIJ_2_0(mNonzeroEdgeValues[csr_index]) += dN_dx(ie_node, 2) * dN_dx(je_node, 0) * volume;
								KRATOS_OCL_LAPLACIANIJ_2_1(mNonzeroEdgeValues[csr_index]) += dN_dx(ie_node, 2) * dN_dx(je_node, 1) * volume;
								KRATOS_OCL_LAPLACIANIJ_2_2(mNonzeroEdgeValues[csr_index]) += dN_dx(ie_node, 2) * dN_dx(je_node, 2) * volume;

								// Contribution to edge gradient
								KRATOS_OCL_NI_DNJ_0(mNonzeroEdgeValues[csr_index]) += dN_dx(je_node, 0) * weighted_volume;
								KRATOS_OCL_NI_DNJ_1(mNonzeroEdgeValues[csr_index]) += dN_dx(je_node, 1) * weighted_volume;
								KRATOS_OCL_NI_DNJ_2(mNonzeroEdgeValues[csr_index]) += dN_dx(je_node, 2) * weighted_volume;

								// Contribution to transposed edge gradient
								KRATOS_OCL_DNI_NJ_0(mNonzeroEdgeValues[csr_index]) += dN_dx(ie_node, 0) * weighted_volume;
								KRATOS_OCL_DNI_NJ_1(mNonzeroEdgeValues[csr_index]) += dN_dx(ie_node, 1) * weighted_volume;
								KRATOS_OCL_DNI_NJ_2(mNonzeroEdgeValues[csr_index]) += dN_dx(ie_node, 2) * weighted_volume;
							}
						}
					}

					// Assembling node based vectors
					for (unsigned int ie_node = 0; ie_node <= 3; ie_node++)
					{
						// Diagonal of the global lumped mass matrix
						mLumpedMassMatrix[nodal_indices[ie_node]] += mass_lumped[ie_node];
					}

					for (unsigned int ie_node = 0; ie_node <= 3; ie_node++)
					{
						// Diagonal of the global gradient matrix
						KRATOS_OCL_COMP_0(mDiagGradientMatrix[nodal_indices[ie_node]]) += dN_dx(ie_node, 0) * weighted_volume;
						KRATOS_OCL_COMP_1(mDiagGradientMatrix[nodal_indices[ie_node]]) += dN_dx(ie_node, 1) * weighted_volume;
						KRATOS_OCL_COMP_2(mDiagGradientMatrix[nodal_indices[ie_node]]) += dN_dx(ie_node, 2) * weighted_volume;
					}
				}

				unsigned int n_nodes = model_part.Nodes().size();

				// Copy mass matrix to inverted mass matrix
				for (unsigned int inode = 0; inode < n_nodes; inode++)
				{
					mInvertedMassMatrix[inode] = mLumpedMassMatrix[inode];
				}

				// Calculating inverted mass matrix (this requires syncronization for MPI paraellelism
				for (unsigned int inode = 0; inode < n_nodes; inode++)
				{
					mInvertedMassMatrix[inode] = 1.00 / mInvertedMassMatrix[inode];
				}

				// Copy data to OpenCL device
				// TODO: Single device code
				mDeviceGroup.CopyBuffer(mbNonzeroEdgeValues, OpenCL::HostToDevice, OpenCL::VoidPList(1, mNonzeroEdgeValues));

				mDeviceGroup.CopyBuffer(mbColumnIndex, OpenCL::HostToDevice, OpenCL::VoidPList(1, mColumnIndex));
				mDeviceGroup.CopyBuffer(mbRowStartIndex, OpenCL::HostToDevice, OpenCL::VoidPList(1, mRowStartIndex));

				mDeviceGroup.CopyBuffer(mbLumpedMassMatrix, OpenCL::HostToDevice, OpenCL::VoidPList(1, mLumpedMassMatrix));
				mDeviceGroup.CopyBuffer(mbInvertedMassMatrix, OpenCL::HostToDevice, OpenCL::VoidPList(1, mInvertedMassMatrix));
				mDeviceGroup.CopyBuffer(mbHmin, OpenCL::HostToDevice, OpenCL::VoidPList(1, mHmin));

				mDeviceGroup.CopyBuffer(mbDiagGradientMatrix, OpenCL::HostToDevice, OpenCL::VoidPList(1, mDiagGradientMatrix));

				KRATOS_CATCH("")
			}

			//
			// GetCSRIndex
			//
			// Function to calculate CSR index of edge ij

			unsigned int GetCSRIndex(unsigned int NodeI, unsigned int NeighbourJ)
			{
				KRATOS_TRY

				// Index indicating data position of edge ij
				unsigned int csr_index;

				// Searching for coincidence of stored column index and neighbour index j
				for (csr_index = mRowStartIndex[NodeI]; csr_index != mRowStartIndex[NodeI + 1]; csr_index++)
				{
					if (mColumnIndex[csr_index] == NeighbourJ)
					{
						break;
					}
				}

				// Returning CSR index of edge ij
				return csr_index;

				KRATOS_CATCH("")
			}

			//
			// FillCoordinatesFromDatabase
			//
			// Function to access database

			void FillCoordinatesFromDatabase(CalcVectorType rDestination, ModelPart::NodesContainerType &rNodes, cl_uint _BufferIndex)
			{
				KRATOS_TRY

				// Loop over all nodes
				for (ModelPart::NodesContainerType::iterator node_it = rNodes.begin(); node_it != rNodes.end(); node_it++)
				{
					// Get the global index of node i
					unsigned int i_node = static_cast <unsigned int> (node_it -> FastGetSolutionStepValue(AUX_INDEX));

					// Save value in the destination vector
					KRATOS_OCL_COMP_0(rDestination[i_node]) = (*node_it)[0];
					KRATOS_OCL_COMP_1(rDestination[i_node]) = (*node_it)[1];
					KRATOS_OCL_COMP_2(rDestination[i_node]) = (*node_it)[2];
				}

				// TODO: Single device code
				mDeviceGroup.CopyBuffer(_BufferIndex, OpenCL::HostToDevice, OpenCL::VoidPList(1, rDestination));

				KRATOS_CATCH("");
			}

			//
			// FillVectorFromDatabase
			//
			// Function to access database

			void FillVectorFromDatabase(Variable <array_1d <double, 3> > &rVariable, CalcVectorType rDestination, ModelPart::NodesContainerType &rNodes, cl_uint _BufferIndex)
			{
				KRATOS_TRY

				// Loop over all nodes
				for (ModelPart::NodesContainerType::iterator node_it = rNodes.begin(); node_it != rNodes.end(); node_it++)
				{
					// Get the global index of node i
					unsigned int i_node = static_cast <unsigned int> (node_it -> FastGetSolutionStepValue(AUX_INDEX));

					// Get the requested value in vector form
					array_1d <double, 3> &vector = node_it -> FastGetSolutionStepValue(rVariable);

					// Save value in the destination vector
					KRATOS_OCL_COMP_0(rDestination[i_node]) = vector[0];
					KRATOS_OCL_COMP_1(rDestination[i_node]) = vector[1];
					KRATOS_OCL_COMP_2(rDestination[i_node]) = vector[2];
				}

				// TODO: Single device code
				mDeviceGroup.CopyBuffer(_BufferIndex, OpenCL::HostToDevice, OpenCL::VoidPList(1, rDestination));

				KRATOS_CATCH("");
			}

			//
			// FillOldVectorFromDatabase
			//
			// Function to access database

			void FillOldVectorFromDatabase(Variable <array_1d <double, 3> > &rVariable, CalcVectorType rDestination, ModelPart::NodesContainerType &rNodes, cl_uint _BufferIndex)
			{
				KRATOS_TRY

				// Loop over alle nodes
				for (ModelPart::NodesContainerType::iterator node_it = rNodes.begin(); node_it != rNodes.end(); node_it++)
				{
					// Get the global index of node i
					unsigned int i_node = static_cast <unsigned int> (node_it -> FastGetSolutionStepValue(AUX_INDEX));

					// Get the requested value in vector form
					array_1d <double, 3> &vector = node_it -> FastGetSolutionStepValue(rVariable, 1);

					// Save value in the destination vector
					KRATOS_OCL_COMP_0(rDestination[i_node]) = vector[0];
					KRATOS_OCL_COMP_1(rDestination[i_node]) = vector[1];
					KRATOS_OCL_COMP_2(rDestination[i_node]) = vector[2];
				}

				// TODO: Single device code
				mDeviceGroup.CopyBuffer(_BufferIndex, OpenCL::HostToDevice, OpenCL::VoidPList(1, rDestination));

				KRATOS_CATCH("");
			}

			//
			// FillScalarFromDatabase
			//
			// Function to access database

			void FillScalarFromDatabase(Variable <double> &rVariable, ValuesVectorType rDestination, ModelPart::NodesContainerType &rNodes, cl_uint _BufferIndex)
			{
				KRATOS_TRY

				// Loop over all nodes
				for (ModelPart::NodesContainerType::iterator node_it = rNodes.begin(); node_it != rNodes.end(); node_it++)
				{
					// Get the global index of node i
					unsigned int i_node = static_cast <unsigned int> (node_it -> FastGetSolutionStepValue(AUX_INDEX));

					// Get the requested scalar value
					double &scalar = node_it -> FastGetSolutionStepValue(rVariable);

					// Save value in the destination vector
					rDestination[i_node] = scalar;
				}

				// TODO: Single device code
				mDeviceGroup.CopyBuffer(_BufferIndex, OpenCL::HostToDevice, OpenCL::VoidPList(1, rDestination));

				KRATOS_CATCH("");
			}

			//
			// FillOldScalarFromDatabase
			//
			// Function to access database

			void FillOldScalarFromDatabase(Variable <double> &rVariable, ValuesVectorType rDestination, ModelPart::NodesContainerType &rNodes, cl_uint _BufferIndex)
			{
				KRATOS_TRY

				// Loop over all nodes
				for (ModelPart::NodesContainerType::iterator node_it = rNodes.begin(); node_it != rNodes.end(); node_it++)
				{
					// Get the global index of node i
					unsigned int i_node = static_cast <unsigned int> (node_it -> FastGetSolutionStepValue(AUX_INDEX));

					// Get the requested scalar value
					double &scalar = node_it -> FastGetSolutionStepValue(rVariable, 1);

					// Save value in the destination vector
					rDestination[i_node] = scalar;
				}

				// TODO: Single device code
				mDeviceGroup.CopyBuffer(_BufferIndex, OpenCL::HostToDevice, OpenCL::VoidPList(1, rDestination));

				KRATOS_CATCH("");
			}

			//
			// WriteVectorToDatabase
			//
			// Function to access database

			void WriteVectorToDatabase(Variable <array_1d <double, 3> > &rVariable, CalcVectorType rOrigin, ModelPart::NodesContainerType &rNodes, cl_uint _BufferIndex)
			{
				KRATOS_TRY

				// TODO: Single device code
				mDeviceGroup.CopyBuffer(_BufferIndex, OpenCL::DeviceToHost, OpenCL::VoidPList(1, rOrigin));

				// Loop over all nodes
				for (ModelPart::NodesContainerType::iterator node_it = rNodes.begin(); node_it != rNodes.end(); node_it++)
				{
					// Get the global index of node i
					unsigned int i_node = static_cast <unsigned int> (node_it -> FastGetSolutionStepValue(AUX_INDEX));

					// Get reference of destination
					array_1d <double, 3> &vector = node_it -> FastGetSolutionStepValue(rVariable);

					// Save vector in database
					vector[0] = KRATOS_OCL_COMP_0(rOrigin[i_node]);
					vector[1] = KRATOS_OCL_COMP_1(rOrigin[i_node]);
					vector[2] = KRATOS_OCL_COMP_2(rOrigin[i_node]);
				}

				KRATOS_CATCH("");
			}

			//
			// WriteScalarToDatabase
			//
			// Function to access database

			void WriteScalarToDatabase(Variable <double> &rVariable, ValuesVectorType rOrigin, ModelPart::NodesContainerType &rNodes, cl_uint _BufferIndex)
			{
				KRATOS_TRY

				// TODO: Single device code
				mDeviceGroup.CopyBuffer(_BufferIndex, OpenCL::DeviceToHost, OpenCL::VoidPList(1, rOrigin));

				// Loop over all nodes
				for (ModelPart::NodesContainerType::iterator node_it = rNodes.begin(); node_it != rNodes.end(); node_it++)
				{
					// Get the global index of node i
					unsigned int i_node = static_cast <unsigned int> (node_it -> FastGetSolutionStepValue(AUX_INDEX));

					// Get reference of destination
					double &scalar = node_it -> FastGetSolutionStepValue(rVariable);

					// Save scalar in database
					scalar = rOrigin[i_node];
				}

				KRATOS_CATCH("");
			}

			//
			// Add_Minv_value1
			//
			// destination = origin1 + value * Minv * origin
			// cl_double version

			void Add_Minv_value1(cl_uint DestinationBufferIndex, cl_uint Origin1BufferIndex, const cl_double Value, cl_uint MinvBufferIndex, cl_uint OriginBufferIndex)
			{
				// TODO: Check if buffers are of the same size
				// TODO: Single device code

				cl_uint n = mDeviceGroup.BufferLengths[OriginBufferIndex][0] / sizeof(cl_double);

				// Setting arguments
				mDeviceGroup.SetBufferAsKernelArg(mkAdd_Minv_value1, 0, DestinationBufferIndex);
				mDeviceGroup.SetBufferAsKernelArg(mkAdd_Minv_value1, 1, Origin1BufferIndex);
				mDeviceGroup.SetKernelArg(mkAdd_Minv_value1, 2, Value);
				mDeviceGroup.SetBufferAsKernelArg(mkAdd_Minv_value1, 3, MinvBufferIndex);
				mDeviceGroup.SetBufferAsKernelArg(mkAdd_Minv_value1, 4, OriginBufferIndex);
				mDeviceGroup.SetKernelArg(mkAdd_Minv_value1, 5, n);

				// Execute OpenCL kernel
				mDeviceGroup.ExecuteKernel(0, mkAdd_Minv_value1, n);
			}

			//
			// Add_Minv_value3
			//
			// destination = origin1 + value * Minv * origin
			// cl_double3 version

			void Add_Minv_value3(cl_uint DestinationBufferIndex, cl_uint Origin1BufferIndex, const cl_double Value, cl_uint MinvBufferIndex, cl_uint OriginBufferIndex)
			{
				// TODO: Check if buffers are of the same size
				// TODO: Single device code

				cl_uint n = mDeviceGroup.BufferLengths[OriginBufferIndex][0] / sizeof(cl_double3);

				// Setting arguments
				mDeviceGroup.SetBufferAsKernelArg(mkAdd_Minv_value3, 0, DestinationBufferIndex);
				mDeviceGroup.SetBufferAsKernelArg(mkAdd_Minv_value3, 1, Origin1BufferIndex);
				mDeviceGroup.SetKernelArg(mkAdd_Minv_value3, 2, Value);
				mDeviceGroup.SetBufferAsKernelArg(mkAdd_Minv_value3, 3, MinvBufferIndex);
				mDeviceGroup.SetBufferAsKernelArg(mkAdd_Minv_value3, 4, OriginBufferIndex);
				mDeviceGroup.SetKernelArg(mkAdd_Minv_value3, 5, n);

				// Execute OpenCL kernel
				mDeviceGroup.ExecuteKernel(0, mkAdd_Minv_value3, n);
			}

			//
			// SetToZero
			//
			// Zeros a buffer on the OpenCL device group

			void SetToZero(cl_uint BufferIndex)
			{
				// TODO: Single device code
				cl_uint n = mDeviceGroup.BufferLengths[BufferIndex][0] / sizeof(double);

				// Setting arguments
				mDeviceGroup.SetBufferAsKernelArg(mkSetToZero, 0, BufferIndex);
				mDeviceGroup.SetKernelArg(mkSetToZero, 1, n);

				// Execute OpenCL kernel
				mDeviceGroup.ExecuteKernel(0, mkSetToZero, n);
			}

			//
			// AssignVectorToVector
			//
			// Copies a vector to another

			void AssignVectorToVector(cl_uint SourceBufferIndex, cl_uint DestinationBufferIndex)
			{
				// TODO: Single device code
				mDeviceGroup.CopyBufferToBuffer(SourceBufferIndex, DestinationBufferIndex);
			}

			//
			// Clear
			//
			// Frees allocated memory

			void Clear()
			{
				KRATOS_TRY

				// TODO: OpenCL buffers are not freed, so user cannot use Clear() and start over
				std::cout <<
					"Warning!" << std::endl <<
					"OpenCL buffers are not freed, so user cannot use Clear() and start over." << std::endl;

				// Free memory
				FreeArray(&mNonzeroEdgeValues);
				FreeArray(&mColumnIndex);
				FreeArray(&mRowStartIndex);
				FreeArray(&mLumpedMassMatrix);
				FreeArray(&mInvertedMassMatrix);
				FreeArray(&mDiagGradientMatrix);
				FreeArray(&mHmin);

				KRATOS_CATCH("")
			}

			//
			// CalculateMassMatrix
			//
			// Function to set up elemental mass matrices

			void CalculateMassMatrix(boost::numeric::ublas::bounded_matrix <double, 3, 3> &mass_consistent, double volume)
			{
				for (unsigned int i_node = 0; i_node <= 3; i_node++)
				{
					// Diagonal terms
					mass_consistent(i_node, i_node) = 0.16666666666666666667 * volume; // 1/6

					// Non-diagonal terms
					double temp = 0.08333333333333333333 * volume; // 1/12

					for (unsigned int j_neighbour = i_node + 1; j_neighbour <= 3; j_neighbour++)
					{
						// Taking advantage of symmetry
						mass_consistent(i_node, j_neighbour) = temp;
						mass_consistent(j_neighbour, i_node) = temp;
					}
				}
			}

			//
			// CalculateMassMatrix
			//
			// Function to set up elemental mass matrices

			void CalculateMassMatrix(boost::numeric::ublas::bounded_matrix<double, 4, 4> &mass_consistent, double volume)
			{
				for (unsigned int i_node = 0; i_node <= 3; i_node++)
				{
					// Diagonal terms
					mass_consistent(i_node, i_node) = 0.1 * volume;

					// Non-diagonal terms
					double temp = 0.05 * volume;

					for (unsigned int j_neighbour = i_node + 1; j_neighbour <= 3; j_neighbour++)
					{
						// Taking advantage of symmetry
						mass_consistent(i_node, j_neighbour) = temp;
						mass_consistent(j_neighbour, i_node) = temp;
					}
				}
			}

		private:

			// Number of edges
			unsigned int mNumberEdges;

			// OpenCL stuff
			OpenCL::DeviceGroup mDeviceGroup;

			// OpenCL program and kernels
			cl_uint mpOpenCLEdgeData, mkAdd_Minv_value1, mkAdd_Minv_value3, mkSetToZero;

			// OpenCL buffers
			cl_uint mbNonzeroEdgeValues, mbColumnIndex, mbRowStartIndex, mbLumpedMassMatrix, mbInvertedMassMatrix, mbDiagGradientMatrix, mbHmin;

			// CSR data vector for storage of the G, L and consistent M components of edge ij
			EdgesVectorType mNonzeroEdgeValues;

			// Vector to store column indices of nonzero matrix elements for each row
			IndicesVectorType mColumnIndex;

			// Index vector to access the start of matrix row i in the column vector
			IndicesVectorType mRowStartIndex;

			// Inverse of the mass matrix ... for parallel calculation each subdomain should contain this correctly calculated (including contributions of the neighbours)
			ValuesVectorType mInvertedMassMatrix;

			// Minimum height around one node
			ValuesVectorType mHmin;

			// Lumped mass matrix (separately stored due to lack of diagonal elements of the consistent mass matrix)
			ValuesVectorType mLumpedMassMatrix;

			// Diagonal of the gradient matrix (separately stored due to special calculations)
			CalcVectorType mDiagGradientMatrix;

	};

} //namespace Kratos

#endif //KRATOS_OPENCL_EDGE_DATA_H_INCLUDED defined
