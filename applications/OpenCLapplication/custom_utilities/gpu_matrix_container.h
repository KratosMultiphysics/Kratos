/*
==============================================================================
KratosOpenCLApplication 
A library based on:
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi
pooyan@cimne.upc.edu 
rrossi@cimne.upc.edu
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
//   Last Modified by:    $Author: antonia $
//   Date:                $Date: 2009-01-14 08:26:51 $
//   Revision:            $Revision: 1.11 $
//
//


#if !defined(KRATOS_EDGE_DATA_H_INCLUDED )
#define  KRATOS_EDGE_DATA_H_INCLUDED


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
#include "OpenCLapplication.h"


namespace Kratos
{			
	class GpuGpuMatrixContainer3D
	{
		public:
// 			//name for the self defined structure
// 			typedef EdgesStructureType<TDim> CSR_Tuple;
// 			typedef std::vector<CSR_Tuple> EdgesVectorType;
// 			//name for row start and column index vectors
// 			typedef std::vector<unsigned int> IndicesVectorType;
// 			//names for separately stored node based values
// 			typedef std::vector<double> ValuesVectorType;
// 			typedef std::vector< array_1d<double,TDim> > CalcVectorType;

			//constructor and destructor
			GpuMatrixContainer(/*const MatrixContainer& cpu_matrix_container*/)
			{
				//here you should make the copy
			};
			~GpuMatrixContainer(){};

/*			//functions to return private values
			inline unsigned int GetNumberEdges(){return mNumberEdges;}
			inline EdgesVectorType& GetEdgeValues(){return mNonzeroEdgeValues;}
			inline IndicesVectorType& GetColumnIndex(){return mColumnIndex;}
			inline IndicesVectorType& GetRowStartIndex(){return mRowStartIndex;}
			inline ValuesVectorType& GetLumpedMass(){return mLumpedMassMatrix;}
			inline ValuesVectorType& GetInvertedMass(){return mInvertedMassMatrix;}
			inline CalcVectorType& GetDiagGradient(){return mDiagGradientMatrix;}
			inline ValuesVectorType& GetHmin(){return mHmin;}

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
*/
			

	};

} //namespace Kratos

#endif //KRATOS_EDGE_DATA_H_INCLUDED defined


