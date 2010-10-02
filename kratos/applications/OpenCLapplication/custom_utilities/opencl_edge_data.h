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
#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/node.h"
//#include "geometries/geometry.h"
#include "utilities/geometry_utilities.h"
#include "incompressible_fluid_application.h"

// Some useful macros, will be renamed if not consistent
#define KRATOS_OCL_LAPLACIANIJ_0_0(a)	a.s0
#define KRATOS_OCL_LAPLACIANIJ_0_1(a)	a.s1
#define KRATOS_OCL_LAPLACIANIJ_0_2(a)	a.s2
#define KRATOS_OCL_LAPLACIANIJ_1_0(a)	a.s3
#define KRATOS_OCL_LAPLACIANIJ_1_1(a)	a.s4
#define KRATOS_OCL_LAPLACIANIJ_1_2(a)	a.s5
#define KRATOS_OCL_LAPLACIANIJ_2_0(a)	a.s6
#define KRATOS_OCL_LAPLACIANIJ_2_1(a)	a.s7
#define KRATOS_OCL_LAPLACIANIJ_2_2(a)	a.s8

#define KRATOS_OCL_MASS(a)				a.s9

#define KRATOS_OCL_NI_DNJ_0(a)			a.sa
#define KRATOS_OCL_NI_DNJ_1(a)			a.sb
#define KRATOS_OCL_NI_DNJ_2(a)			a.sc

#define KRATOS_OCL_DNI_NJ_0(a)			a.sd
#define KRATOS_OCL_DNI_NJ_1(a)			a.se
#define KRATOS_OCL_DNI_NJ_2(a)			a.sf

#define KRATOS_OCL_COMP(a, n)			a.s[n]

#define KRATOS_OCL_COMP_0(a)			a.x
#define KRATOS_OCL_COMP_1(a)			a.y
#define KRATOS_OCL_COMP_2(a)			a.z

namespace Kratos
{

	void FillOpenCLEdge(cl_double16& a, const EdgesStructureType& b)
	{
		KRATOS_OCL_LAPLACIANIJ_0_0(a) = b.LaplacianIJ(0, 0);
		KRATOS_OCL_LAPLACIANIJ_0_1(a) = b.LaplacianIJ(0, 1);
		KRATOS_OCL_LAPLACIANIJ_0_2(a) = b.LaplacianIJ(0, 2);
		KRATOS_OCL_LAPLACIANIJ_1_0(a) = b.LaplacianIJ(1, 0);
		KRATOS_OCL_LAPLACIANIJ_1_1(a) = b.LaplacianIJ(1, 1);
		KRATOS_OCL_LAPLACIANIJ_1_2(a) = b.LaplacianIJ(1, 2);
		KRATOS_OCL_LAPLACIANIJ_2_0(a) = b.LaplacianIJ(2, 0);
		KRATOS_OCL_LAPLACIANIJ_2_1(a) = b.LaplacianIJ(2, 1);
		KRATOS_OCL_LAPLACIANIJ_2_2(a) = b.LaplacianIJ(2, 2);

		KRATOS_OCL_MASS(a) = b.Mass;

		KRATOS_OCL_NI_DNJ_0(a) = b.Ni_DNj[0];
		KRATOS_OCL_NI_DNJ_1(a) = b.Ni_DNj[1];
		KRATOS_OCL_NI_DNJ_2(a) = b.Ni_DNj[2];

		KRATOS_OCL_DNI_NJ_0(a) = b.DNi_Nj[0];
		KRATOS_OCL_DNI_NJ_1(a) = b.DNi_Nj[1];
		KRATOS_OCL_DNI_NJ_2(a) = b.DNi_Nj[2];
	}

	//
	// OpenCLMatrixContainer
	//
	// A helper class for edge based mathods

	class OpenCLMatrixContainer
	{
		public:

			//
			// OpenCLMatrixContainer
			//
			// Constructor

			OpenCLMatrixContainer(MatrixContainer _matrix_container): matrix_container(_matrix_container)
			{
				// Nothing to do!
			}

			//
			// ~OpenCLMatrixContainer
			//
			// Destructor

			~OpenCLMatrixContainer()
			{
				// Nothing to do!
			}

		private:

			MatrixContainer matrix_container;
	}


////////////////// Up to here!

	//class definition of matrices using CSR format
	class OpenCLMatrixContainer
	{
		public:
			//name for the self defined structure
			typedef EdgesStructureType<TDim> CSR_Tuple;
			typedef std::vector<CSR_Tuple> EdgesVectorType;
			//name for row start and column index vectors
			typedef std::vector<unsigned int> IndicesVectorType;
			//names for separately stored node based values
			typedef std::vector<double> ValuesVectorType;
			typedef std::vector< array_1d<double,TDim> > CalcVectorType;

			//constructor and destructor
			MatrixContainer(){};
			~MatrixContainer(){};

			//functions to return private values
			inline unsigned int GetNumberEdges(){return mNumberEdges;}
			inline EdgesVectorType& GetEdgeValues(){return mNonzeroEdgeValues;}
			inline IndicesVectorType& GetColumnIndex(){return mColumnIndex;}
			inline IndicesVectorType& GetRowStartIndex(){return mRowStartIndex;}
			inline ValuesVectorType& GetLumpedMass(){return mLumpedMassMatrix;}
			inline ValuesVectorType& GetInvertedMass(){return mInvertedMassMatrix;}
			inline CalcVectorType& GetDiagGradient(){return mDiagGradientMatrix;}
			inline ValuesVectorType& GetHmin(){return mHmin;}

			//********************************************************
			//function to size and initialize the vector of CSR tuples
			void ConstructCSRVector(ModelPart& model_part)
			{
			KRATOS_TRY

				//SIZE OF CSR VECTOR

				//defining the number of nodes and edges
				unsigned int n_nodes = model_part.Nodes().size();
				//remark: no colouring algorithm is used here (symmetry is neglected)
				//        respectively edge ij is considered different from edge ji
				mNumberEdges = 0;
				//counter to assign and get global nodal index
				unsigned int i_node = 0;

				//counting the edges connecting the nodes
				for (typename ModelPart::NodesContainerType::iterator node_it=model_part.NodesBegin(); node_it!=model_part.NodesEnd(); node_it++)
				{
					//counting neighbours of each node
					mNumberEdges += (node_it->GetValue(NEIGHBOUR_NODES)).size();
					//DIAGONAL TERMS
					//mNumberEdges++;

					//assigning global index to each node
					node_it->FastGetSolutionStepValue(AUX_INDEX) = static_cast<double>(i_node++);
				}
				//error message in case number of nodes does not coincide with number of indices
				if (i_node != n_nodes)
					KRATOS_WATCH("ERROR - Highest nodal index doesn't coincide with number of nodes!");

				//allocating memory for block of CSR data
				mNonzeroEdgeValues.resize(mNumberEdges);
				mColumnIndex.resize(mNumberEdges);
				mRowStartIndex.resize(n_nodes+1);
				mLumpedMassMatrix.resize(n_nodes);
				mInvertedMassMatrix.resize(n_nodes);
				mDiagGradientMatrix.resize(n_nodes);
				mHmin.resize(n_nodes);

				//INITIALIZING OF THE CSR VECTOR

				//temporary variable as the row start index of a node depends on the number of neighbours of the previous one
				unsigned int row_start_temp = 0;
				//main loop over all nodes
				for (typename ModelPart::NodesContainerType::iterator node_it=model_part.NodesBegin(); node_it!=model_part.NodesEnd(); node_it++)
				{
					//getting the global index of the node
					i_node = static_cast<unsigned int>(node_it->FastGetSolutionStepValue(AUX_INDEX));
					//determining its neighbours
					WeakPointerVector< Node<3> >& neighb_nodes = node_it->GetValue(NEIGHBOUR_NODES);
					//number of neighbours of node i determines row start index for the following node
					unsigned int n_neighbours = neighb_nodes.size();
					//DIAGONAL TERMS
					//n_neighbours++;

					//reserving memory for work array
					std::vector<unsigned int> work_array;
					work_array.reserve(n_neighbours);
					//DIAGONAL TERMS
					//work_array.push_back(i_node);

					//nested loop over the neighbouring nodes
					for (WeakPointerVector< Node<3> >::iterator neighb_it=neighb_nodes.begin(); neighb_it!=neighb_nodes.end(); neighb_it++)
					{
						//getting global index of the neighbouring node
						work_array.push_back(static_cast<unsigned int>(neighb_it->FastGetSolutionStepValue(AUX_INDEX)));
					}
					//reordering neighbours following their global indices
					std::sort(work_array.begin(),work_array.end());

					//setting current row start index
					mRowStartIndex[i_node] = row_start_temp;
					//nested loop over the by now ordered neighbours
					for (unsigned int counter = 0; counter < n_neighbours; counter++)
					{
						//getting global index of the neighbouring node
						unsigned int j_neighbour = work_array[counter];
						//calculating CSR index
						unsigned int csr_index = mRowStartIndex[i_node]+counter;

						//saving column index j of the original matrix
						mColumnIndex[csr_index] = j_neighbour;
						//initializing the CSR vector entries with zero
						mNonzeroEdgeValues[csr_index].Mass = 0.0;

						//mNonzeroEdgeValues[csr_index].Laplacian = 0.0;
						noalias(mNonzeroEdgeValues[csr_index].LaplacianIJ) = ZeroMatrix(TDim,TDim);
						noalias(mNonzeroEdgeValues[csr_index].Ni_DNj) = ZeroVector(TDim);
						//TRANSPOSED GRADIENT
						noalias(mNonzeroEdgeValues[csr_index].DNi_Nj) = ZeroVector(TDim);
					}
					//preparing row start index for next node
					row_start_temp += n_neighbours;
				}
				//adding last entry (necessary for abort criterion of loops)
				mRowStartIndex[n_nodes] = mNumberEdges;

				//INITIALIZING NODE BASED VALUES

				//lumped mass matrix (elements Mi)
				for (i_node=0; i_node<n_nodes; i_node++)
					mLumpedMassMatrix[i_node] = 0.0;

				//set the heights to a huge number
				for (i_node=0; i_node<n_nodes; i_node++)
					mHmin[i_node] = 1e10;

				//diagonal of gradient matrix (elements Gii)
				for (i_node=0; i_node<n_nodes; i_node++)
					noalias(mDiagGradientMatrix[i_node]) = ZeroVector(TDim);

			KRATOS_CATCH("")
			}

			//*********************************
			//function to precalculate CSR data
			void BuildCSRData(ModelPart& model_part)
			{
			KRATOS_TRY

				//PRECALCULATING CSR DATA

				//defining temporary local variables for elementwise addition
				//shape functions
				array_1d<double, TDim+1> N;
				//shape function derivatives
				boost::numeric::ublas::bounded_matrix <double, TDim+1,TDim> dN_dx;
				//volume
				double volume;
				//weighting factor
				double weighting_factor = 1.0 / static_cast<double>(TDim+1);
				//elemental matrices
				boost::numeric::ublas::bounded_matrix <double, TDim+1,TDim+1> mass_consistent;
				//boost::numeric::ublas::bounded_matrix <double, TDim+1,TDim+1> laplacian;
				array_1d<double, TDim+1> mass_lumped;
				//global indices of elemental nodes
				array_1d<unsigned int, TDim+1> nodal_indices;

				array_1d<double, TDim+1> heights;


				//loop over all elements
				for (typename ModelPart::ElementsContainerType::iterator elem_it=model_part.ElementsBegin(); elem_it!=model_part.ElementsEnd(); elem_it++)
				{
					//LOCAL ELEMENTWISE CALCULATIONS

					//getting geometry data of the element
					GeometryUtils::CalculateGeometryData(elem_it->GetGeometry(), dN_dx, N, volume);

					//calculate lenght of the heights of the element
					for (unsigned int ie_node=0; ie_node<=TDim; ie_node++)
					{
						heights[ie_node] = dN_dx(ie_node,0)*dN_dx(ie_node,0);
						for (unsigned int comp=1; comp<TDim; comp++)
						{
							heights[ie_node] += dN_dx(ie_node,comp)*dN_dx(ie_node,comp);
						}
						heights[ie_node] = 1.0/sqrt(heights[ie_node]);
// KRATOS_WATCH(heights);
					}


					//setting up elemental mass matrices
					CalculateMassMatrix(mass_consistent, volume);
					noalias(mass_lumped) = ZeroVector(TDim+1);
					for (unsigned int ie_node=0; ie_node<=TDim; ie_node++)
					{
						for (unsigned int je_node=0; je_node<=TDim; je_node++)
						{
							//mass_consistent(ie_node,je_node) = N(ie_node) * N(je_node) * volume;
							mass_lumped[ie_node] += mass_consistent(ie_node,je_node);
						}
						//mass_lumped[ie_node] = volume * N[ie_node];
					}
					/*OLD DATA STRUCTURE
					//calculating elemental laplacian matrix
					noalias(laplacian) = ZeroMatrix(TDim+1,TDim+1);
					for (unsigned int ie_node=0; ie_node<=TDim; ie_node++)
						for (unsigned int je_node=ie_node+1; je_node<=TDim; je_node++)
							//componentwise multiplication
							for (unsigned int component=0; component<TDim; component++)
							{
								//taking advantage of symmetry
								double temp = dN_dx(ie_node,component) * dN_dx(je_node,component) * volume;
								laplacian(ie_node,je_node) += temp;
								laplacian(je_node,ie_node) += temp;
							}

					//multiply gradient with volume referring to each gauss point
					dN_dx *= (volume / double(TDim+1));*/
					//(corresponding to Ni * dOmega respectively Nj * dOmega)
					double weighted_volume = volume * weighting_factor;

					//ASSEMBLING GLOBAL DATA STRUCTURE

					//loop over the nodes of the element to determine their global indices
					for (unsigned int ie_node=0; ie_node<=TDim; ie_node++)
						nodal_indices[ie_node] = static_cast<unsigned int>(elem_it->GetGeometry()[ie_node].FastGetSolutionStepValue(AUX_INDEX));

					//assembling global "edge matrices" by adding local contributions
					for (unsigned int ie_node=0; ie_node<=TDim; ie_node++)
					{
						//check the heights and change the value if minimal is found
						if( mHmin[ nodal_indices[ie_node] ] > heights[ie_node])
							mHmin[ nodal_indices[ie_node] ] = heights[ie_node];

						for (unsigned int je_node=0; je_node<=TDim; je_node++)
						{
							//remark: there is no edge linking node i with itself!
							//DIAGONAL TERMS
							if (ie_node != je_node)
							{
								//calculating CSR index from global index
								unsigned int csr_index = GetCSRIndex(nodal_indices[ie_node], nodal_indices[je_node]);

								//assigning precalculated element data to the referring edges
								//contribution to edge mass
								mNonzeroEdgeValues[csr_index].Mass += mass_consistent(ie_node,je_node);


								//contribution to edge laplacian
								/*OLD DATA STRUCTURE
								mNonzeroEdgeValues[csr_index].Laplacian = laplacian(ie_node,je_node);*/
								boost::numeric::ublas::bounded_matrix <double,TDim,TDim>& laplacian = mNonzeroEdgeValues[csr_index].LaplacianIJ;
								for (unsigned int l_comp=0; l_comp<TDim; l_comp++)
									for (unsigned int k_comp=0; k_comp<TDim; k_comp++)
										laplacian(l_comp,k_comp) += dN_dx(ie_node,l_comp) * dN_dx(je_node,k_comp) * volume;

								//contribution to edge gradient
								array_1d<double, TDim>&  gradient = mNonzeroEdgeValues[csr_index].Ni_DNj;
								for (unsigned int l_comp=0; l_comp<TDim; l_comp++)
									//gradient[l_comp] += dN_dx(je_node,l_comp);
									gradient[l_comp] += dN_dx(je_node,l_comp) * weighted_volume;
								//TRANSPOSED GRADIENT
								//contribution to transposed edge gradient
								array_1d<double, TDim>& transp_gradient = mNonzeroEdgeValues[csr_index].DNi_Nj;
								for (unsigned int l_comp=0; l_comp<TDim; l_comp++)
									//transp_gradient[l_comp] += dN_dx(ie_node,l_comp);
									transp_gradient[l_comp] += dN_dx(ie_node,l_comp) * weighted_volume;
							}
						}
					}

					//assembling node based vectors
					for (unsigned int ie_node=0; ie_node<=TDim; ie_node++)
						//diagonal of the global lumped mass matrix
						mLumpedMassMatrix[nodal_indices[ie_node]] += mass_lumped[ie_node];
					for (unsigned int ie_node=0; ie_node<=TDim; ie_node++)
					{
						//diagonal of the global gradient matrix
						array_1d<double, TDim>&  gradient = mDiagGradientMatrix[nodal_indices[ie_node]];
						for (unsigned int component=0; component<TDim; component++)
							//gradient[component] += dN_dx(ie_node,component);
							gradient[component] += dN_dx(ie_node,component) * weighted_volume;
					}
				}

				//copy mass matrix to inverted mass matrix
				for(unsigned int inode=0; inode<mLumpedMassMatrix.size(); inode++)
				{
					mInvertedMassMatrix[inode] = mLumpedMassMatrix[inode];
				}

				//perform MPI syncronization between the domains

				//calculating inverted mass matrix (this requires syncronization for MPI paraellelism
				for(unsigned int inode=0; inode<mInvertedMassMatrix.size(); inode++)
				{
					mInvertedMassMatrix[inode] = 1.0/mInvertedMassMatrix[inode];
				}


			KRATOS_CATCH("")
			}

			//******************************************
			//function to calculate CSR index of edge ij
			unsigned int GetCSRIndex(unsigned int NodeI, unsigned int NeighbourJ)
			{
			KRATOS_TRY

				//index indicating data position of edge ij
				unsigned int csr_index;
				//searching for coincidence of stored column index and neighbour index j
				for (csr_index=mRowStartIndex[NodeI]; csr_index!=mRowStartIndex[NodeI+1]; csr_index++)
					if (mColumnIndex[csr_index] == NeighbourJ)
						break;

				//returning CSR index of edge ij
				return csr_index;

			KRATOS_CATCH("")
			}

			//***********************************************
			//function to get pointer to CSR tuple of edge ij
			CSR_Tuple* GetTuplePointer(unsigned int NodeI, unsigned int NeighbourJ)
			{
			KRATOS_TRY

				//index indicating data position of edge ij
				unsigned int csr_index;
				//searching for coincidence of stored column index and neighbour index j
				for (csr_index=mRowStartIndex[NodeI]; csr_index!=mRowStartIndex[NodeI+1]; csr_index++)
					if (mColumnIndex[csr_index] == NeighbourJ)
						break;

				//returning pointer to CSR tuple of edge ij
				return &mNonzeroEdgeValues[csr_index];

			KRATOS_CATCH("")
			}

			//*******************************
			//function to free dynamic memory
			void Clear()
			{
			KRATOS_TRY

				mNonzeroEdgeValues.clear();
				mColumnIndex.clear();
				mRowStartIndex.clear();
				mInvertedMassMatrix.clear();
				mLumpedMassMatrix.clear();
				mDiagGradientMatrix.clear();
				mHmin.clear();

			KRATOS_CATCH("")
			}
			//****************************
			//functions to access database
			//(note that this is already thought for parallel;
			// for a single processor this could be done in a faster way)
			void FillCoordinatesFromDatabase( CalcVectorType& rDestination, ModelPart::NodesContainerType& rNodes)
			{

				KRATOS_TRY

				//loop over alle nodes
						for (typename ModelPart::NodesContainerType::iterator node_it=rNodes.begin(); node_it!=rNodes.end(); node_it++)
				{
					//get the global index of node i
					unsigned int i_node = static_cast<unsigned int>(node_it->FastGetSolutionStepValue(AUX_INDEX));


					//save value in the destination vector
					for(unsigned int component = 0; component < TDim; component++)
						(rDestination[i_node])[component] = (*node_it)[component];
				}

				KRATOS_CATCH("");
			}

			//****************************
			//functions to access database
			//(note that this is already thought for parallel;
			// for a single processor this could be done in a faster way)
			void FillVectorFromDatabase(Variable<array_1d<double,3> >& rVariable, CalcVectorType& rDestination, ModelPart::NodesContainerType& rNodes)
			{

				KRATOS_TRY

				//loop over alle nodes
						for (typename ModelPart::NodesContainerType::iterator node_it=rNodes.begin(); node_it!=rNodes.end(); node_it++)
				{
					//get the global index of node i
					unsigned int i_node = static_cast<unsigned int>(node_it->FastGetSolutionStepValue(AUX_INDEX));

					//get the requested value in vector form
					array_1d<double,3>& vector = node_it->FastGetSolutionStepValue(rVariable);
					//save value in the destination vector
					for(unsigned int component = 0; component < TDim; component++)
						(rDestination[i_node])[component] = vector[component];
				}

				KRATOS_CATCH("");
			}
			void FillOldVectorFromDatabase(Variable<array_1d<double,3> >& rVariable, CalcVectorType& rDestination, ModelPart::NodesContainerType& rNodes)
			{

				KRATOS_TRY

				//loop over alle nodes
						for (typename ModelPart::NodesContainerType::iterator node_it=rNodes.begin(); node_it!=rNodes.end(); node_it++)
				{
					//get the global index of node i
					unsigned int i_node = static_cast<unsigned int>(node_it->FastGetSolutionStepValue(AUX_INDEX));

					//get the requested value in vector form
					array_1d<double,3>& vector = node_it->FastGetSolutionStepValue(rVariable,1);
					//save value in the destination vector
					for(unsigned int component = 0; component < TDim; component++)
						(rDestination[i_node])[component] = vector[component];
				}

				KRATOS_CATCH("");
			}

			void FillScalarFromDatabase(Variable<double>& rVariable, ValuesVectorType& rDestination, ModelPart::NodesContainerType& rNodes)
			{
				KRATOS_TRY

				//loop over all nodes
						for (typename ModelPart::NodesContainerType::iterator node_it=rNodes.begin(); node_it!=rNodes.end(); node_it++)
				{
					//get the global index of node i
					unsigned int i_node = static_cast<unsigned int>(node_it->FastGetSolutionStepValue(AUX_INDEX));

					//get the requested scalar value
					double& scalar = node_it->FastGetSolutionStepValue(rVariable);
					//save value in the destination vector
					rDestination[i_node] = scalar;
				}

				KRATOS_CATCH("");
			}
			void FillOldScalarFromDatabase(Variable<double>& rVariable, ValuesVectorType& rDestination, ModelPart::NodesContainerType& rNodes)
			{
				KRATOS_TRY

				//loop over all nodes
						for (typename ModelPart::NodesContainerType::iterator node_it=rNodes.begin(); node_it!=rNodes.end(); node_it++)
				{
					//get the global index of node i
					unsigned int i_node = static_cast<unsigned int>(node_it->FastGetSolutionStepValue(AUX_INDEX));

					//get the requested scalar value
					double& scalar = node_it->FastGetSolutionStepValue(rVariable,1);
					//save value in the destination vector
					rDestination[i_node] = scalar;
				}

				KRATOS_CATCH("");
			}

			void WriteVectorToDatabase(Variable<array_1d<double,3> >& rVariable, CalcVectorType& rOrigin, ModelPart::NodesContainerType& rNodes)
			{
				KRATOS_TRY

				//loop over alle nodes
						for (typename ModelPart::NodesContainerType::iterator node_it=rNodes.begin(); node_it!=rNodes.end(); node_it++)
				{
					//get the global index of node i
					unsigned int i_node = static_cast<unsigned int>(node_it->FastGetSolutionStepValue(AUX_INDEX));

					//get reference of destination
					array_1d<double,3>& vector = node_it->FastGetSolutionStepValue(rVariable);
					//save vector in database
					for(unsigned int component = 0; component < TDim; component++)
						vector[component] = (rOrigin[i_node])[component];
				}

				KRATOS_CATCH("");
			}

			void WriteScalarToDatabase(Variable<double>& rVariable, ValuesVectorType& rOrigin, ModelPart::NodesContainerType& rNodes)
			{
				KRATOS_TRY

				//loop over all nodes
						for (typename ModelPart::NodesContainerType::iterator node_it=rNodes.begin(); node_it!=rNodes.end(); node_it++)
				{
					//get the global index of node i
					unsigned int i_node = static_cast<unsigned int>(node_it->FastGetSolutionStepValue(AUX_INDEX));

					//get reference of destination
					double& scalar = node_it->FastGetSolutionStepValue(rVariable);
					//save scalar in database
					scalar = rOrigin[i_node];
				}

				KRATOS_CATCH("");
			}

			//*********************************************************************
			//destination = origin1 + value * Minv*origin
			void Add_Minv_value(
					    CalcVectorType& destination,
				            const CalcVectorType& origin1,
						const double value,
						const ValuesVectorType& Minv_vec,
						const CalcVectorType& origin
					   )
			{
				KRATOS_TRY
						int loop_size = destination.size();
				#pragma omp parallel for
				for (int i_node = 0; i_node < loop_size; i_node++)
				{
					array_1d<double, TDim>& dest = destination[i_node];
					const double m_inv = Minv_vec[i_node];
					const array_1d<double, TDim>& origin_vec1 = origin1[i_node];
					const array_1d<double, TDim>& origin_value = origin[i_node];

					double temp = value * m_inv;
					for(unsigned int comp = 0; comp < TDim; comp++)
						dest[comp] = origin_vec1[comp] + temp * origin_value[comp];
				}


				KRATOS_CATCH("")
			}

			void Add_Minv_value(
					    ValuesVectorType& destination,
				            const ValuesVectorType& origin1,
						const double value,
						const ValuesVectorType& Minv_vec,
						const ValuesVectorType& origin
					   )
			{
				KRATOS_TRY
						int loop_size = destination.size();
				#pragma omp parallel for
				for (int i_node = 0; i_node < loop_size; i_node++)
				{
					double& dest = destination[i_node];
					const double m_inv = Minv_vec[i_node];
					const double& origin_vec1 = origin1[i_node];
					const double& origin_value = origin[i_node];

					double temp = value * m_inv;
					dest = origin_vec1 + temp * origin_value;
				}


				KRATOS_CATCH("")
			}

			//**********************************************************************
			void SetToZero(	CalcVectorType& data_vector)
			{
				int loop_size = data_vector.size();
				#pragma omp parallel for
				for (int i_node = 0; i_node < loop_size; i_node++)
				{
					array_1d<double,TDim>& aaa = data_vector[i_node];
					for (unsigned int comp = 0; comp < TDim; comp++)
						aaa[comp] = 0.0;
				}
			}

			void SetToZero(	ValuesVectorType& data_vector)
			{
				int loop_size = data_vector.size();
				#pragma omp parallel for
				for (int i_node = 0; i_node < loop_size; i_node++)
				{
					 data_vector[i_node]= 0.0;;
				}
			}

			//**********************************************************************
			void AssignVectorToVector(	const CalcVectorType& origin,
					CalcVectorType& destination
						 )
			{
				int loop_size = origin.size();
				#pragma omp parallel for
				for (int i_node = 0; i_node < loop_size; i_node++)
				{
					const array_1d<double,TDim>& orig = origin[i_node];
					array_1d<double,TDim>& dest = destination[i_node];
					for (unsigned int comp = 0; comp < TDim; comp++)
						dest[comp] = orig[comp];
				}
			}


			void AssignVectorToVector(	const ValuesVectorType& origin,
					ValuesVectorType& destination
						 )
			{
				int loop_size = origin.size();
				#pragma omp parallel for
				for (int i_node = 0; i_node < loop_size; i_node++)
				{
					destination[i_node] = origin[i_node] ;
				}
			}



		private:
			//number of edges
			unsigned int mNumberEdges;

			//CSR data vector for storage of the G, L and consistent M components of edge ij
			EdgesVectorType mNonzeroEdgeValues;

			//vector to store column indices of nonzero matrix elements for each row
			IndicesVectorType mColumnIndex;

			//index vector to access the start of matrix row i in the column vector
			IndicesVectorType mRowStartIndex;

			//inverse of the mass matrix ... for parallel calculation each subdomain should contain this correctly calculated (including contributions of the neighbours)
			ValuesVectorType mInvertedMassMatrix;

			//minimum height around one node
			ValuesVectorType mHmin;

			//lumped mass matrix (separately stored due to lack of diagonal elements of the consistent mass matrix)
			ValuesVectorType mLumpedMassMatrix;
			//diagonal of the gradient matrix (separately stored due to special calculations)
			CalcVectorType mDiagGradientMatrix;

			//*******************************************
			//functions to set up elemental mass matrices
			void CalculateMassMatrix(boost::numeric::ublas::bounded_matrix<double,3,3>& mass_consistent, double volume)
			{
				for (unsigned int i_node=0; i_node<=TDim; i_node++)
				{
					//diagonal terms
					mass_consistent(i_node,i_node) = 0.16666666666666666667 * volume; //1/6
					//non-diagonal terms
					double temp = 0.08333333333333333333 * volume; // 1/12
					for(unsigned int j_neighbour=i_node+1; j_neighbour<=TDim; j_neighbour++)
					{
						//taking advantage of symmetry
						mass_consistent(i_node,j_neighbour) = temp;
						mass_consistent(j_neighbour,i_node) = temp;
					}
				}
			}
			void CalculateMassMatrix(boost::numeric::ublas::bounded_matrix<double,4,4>& mass_consistent, double volume)
			{
				for (unsigned int i_node=0; i_node<=TDim; i_node++)
				{
					//diagonal terms
					mass_consistent(i_node,i_node) = 0.1 * volume;
					//non-diagonal terms
					double temp = 0.05 * volume;
					for(unsigned int j_neighbour=i_node+1; j_neighbour<=TDim; j_neighbour++)
					{
						//taking advantage of symmetry
						mass_consistent(i_node,j_neighbour) = temp;
						mass_consistent(j_neighbour,i_node) = temp;
					}
				}
			}

	};

} //namespace Kratos

#endif //KRATOS_OPENCL_EDGE_DATA_H_INCLUDED defined


