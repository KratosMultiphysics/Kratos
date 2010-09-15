/*
==============================================================================
KratosPFEMApplication 
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
//   Last Modified by:    $Author: rrossi $
//   Date:                $Date: 2009-01-13 15:39:56 $
//   Revision:            $Revision: 1.3 $
//
//


#if !defined(KRATOS_PURE_CONVECTION_EDGEBASED_SOLVER_H_INCLUDED)
#define  KRATOS_PURE_CONVECTION_EDGEBASED_SOLVER_H_INCLUDED

  #define SPLIT_OSS
//   #define SYMM_PRESS


// System includes
#include <string>
#include <iostream>
#include <algorithm>

// #include <omp.h>


// External includes 


// Project includes
#include "includes/define.h"
#include "includes/model_part.h"
#include "includes/node.h"
//#include "geometries/geometry.h"
#include "utilities/geometry_utilities.h"
#include "incompressible_fluid_application.h"


namespace Kratos
{
					
				
	template<unsigned int TDim, class MatrixContainer, class TSparseSpace, class TLinearSolver>
	class PureConvectionEdgeBased
	{
		public:
			//name for the self defined structure
			typedef EdgesStructureType<TDim> CSR_Tuple;
			typedef std::vector<CSR_Tuple> EdgesVectorType;

			//name for row start and column index vectors
			typedef std::vector<unsigned int> IndicesVectorType;
			//defining matrix type for test calculations
			typedef std::vector< array_1d<double, TDim> > CalcVectorType;
			//defining type for local storage of nodal values
			typedef std::vector<double> ValuesVectorType;

			//defining types for matrix operations
			typedef typename TSparseSpace::MatrixType TSystemMatrixType;
			typedef typename TSparseSpace::VectorType TSystemVectorType;

			//constructor and destructor
			PureConvectionEdgeBased(MatrixContainer& mr_matrix_container,
				    ModelPart& mr_model_part
				    )
			: mr_matrix_container(mr_matrix_container),mr_model_part(mr_model_part)
			{};
			~PureConvectionEdgeBased(){};

			//***********************************
			//function to initialize fluid solver
			void Initialize(
				       )
			{
			KRATOS_TRY
				
			
				//get number of nodes
				unsigned int n_nodes = mr_model_part.Nodes().size();
			//unsigned int n_edges = mr_matrix_container.GetNumberEdges();

				//size data vectors
				mWork.resize(n_nodes);
				mPi.resize(n_nodes);
				mUn.resize(n_nodes);
				mUn1.resize(n_nodes);
				mphi_n.resize(n_nodes);
				mphi_n1.resize(n_nodes);
				mA.resize(n_nodes);
				mHmin.resize(n_nodes);
				mTau.resize(n_nodes);
				mBeta.resize(n_nodes);
				mx.resize(n_nodes);

				//read variables from Kratos
				mr_matrix_container.FillVectorFromDatabase(VELOCITY, mUn1, mr_model_part.Nodes());
				mr_matrix_container.FillOldVectorFromDatabase(VELOCITY, mUn, mr_model_part.Nodes());

				mr_matrix_container.FillScalarFromDatabase(DISTANCE, mphi_n1, mr_model_part.Nodes());
				mr_matrix_container.FillOldScalarFromDatabase(DISTANCE, mphi_n, mr_model_part.Nodes());

				mr_matrix_container.FillCoordinatesFromDatabase(mx, mr_model_part.Nodes());

				//set flag for first time step
				mFirstStep = true;

				ValuesVectorType& aaa = mr_matrix_container.GetHmin();
				for (unsigned int i_node = 0; i_node < n_nodes; i_node++)
				{
					mHmin[i_node] = aaa[i_node];
				}



			KRATOS_CATCH("")
			}

			//***************************************
			//function to set adequate time step size
			void ComputeTimeStep(double CFLNumber)
			{
			KRATOS_TRY

				//local variable for time step size
				double delta_t = 1e10;
			
				//getting value of current velocity 
				mr_matrix_container.FillVectorFromDatabase(VELOCITY, mUn1, mr_model_part.Nodes());
				
				//loop over all nodes
				double n_nodes = mUn1.size();
				for (unsigned int i_node = 0; i_node < n_nodes; i_node++)
				{
					//use CFL condition to compute time step size
					double delta_t_i = CFLNumber * mHmin[i_node] / norm_2(mUn1[i_node] ) ;
					//choose the overall minimum of delta_t_i
					if (delta_t_i < delta_t)
						delta_t = delta_t_i;
				}
				
				//perform MPI syncronization of the dt (minimum should be kept) 

				//write time step size to Kratos
				ProcessInfo& CurrentProcessInfo = mr_model_part.GetProcessInfo();
				CurrentProcessInfo[DELTA_TIME] = delta_t;

			KRATOS_CATCH("")
			}

			//**********************************************************************************
			//function to solve fluid equations - fractional step 1: compute fractional momentum
			void Solve()
			{
			KRATOS_TRY

				//PREREQUISITES

				//variables for node based data handling
				ModelPart::NodesContainerType& rNodes = mr_model_part.Nodes();
				int n_nodes = rNodes.size();

				//storage of nodal values in local variables
				ValuesVectorType rhs;
				rhs.resize(n_nodes);

				//read variables from Kratos
				mr_matrix_container.FillVectorFromDatabase(VELOCITY, mUn1, mr_model_part.Nodes());
				mr_matrix_container.FillOldVectorFromDatabase(VELOCITY, mUn, mr_model_part.Nodes());

				mr_matrix_container.FillScalarFromDatabase(DISTANCE, mphi_n1, mr_model_part.Nodes());
				mr_matrix_container.FillOldScalarFromDatabase(DISTANCE, mphi_n, mr_model_part.Nodes());
				

				//read time step size from Kratos
				ProcessInfo& CurrentProcessInfo = mr_model_part.GetProcessInfo();
				double delta_t = CurrentProcessInfo[DELTA_TIME];
				
				//compute advective velocity - area average of the current velocity
				double coefficient = 1;
				CalculateAdvectiveVelocity(mUn,mUn1,mA, coefficient);

				//compute intrinsic time
				double time_inv = 1.0/delta_t;

				#pragma omp parallel for firstprivate(time_inv)	
				for (int i_node = 0; i_node < n_nodes; i_node++)
				{
					double& h_i = mHmin[i_node];
					array_1d<double, TDim>& a_i = mA[i_node];

					double vel_norm = norm_2(a_i);
					mTau[i_node] = 1.0 / (2.0 * vel_norm/h_i + 0.01*time_inv );
				}
				mr_matrix_container.AssignVectorToVector(mphi_n, mWork); //mWork = mphi_n
				
				//first step of Runge Kutta
// 				mr_matrix_container.AssignVectorToVector(mphi_n,mphi_n1); //mphi_n1 = mphi_n
				mr_matrix_container.SetToZero(rhs);
				CalculateRHS( mphi_n1,mA,rhs);
	
				mr_matrix_container.Add_Minv_value(mWork,mWork,     delta_t/6.0    , mr_matrix_container.GetInvertedMass(), rhs);
				mr_matrix_container.Add_Minv_value(mphi_n1, mphi_n,   0.5*delta_t  , mr_matrix_container.GetInvertedMass(), rhs);

				//second step
				mr_matrix_container.SetToZero(rhs);
				CalculateRHS(mphi_n1,mA,rhs);
				mr_matrix_container.Add_Minv_value(mWork,mWork,  delta_t/3.0    , mr_matrix_container.GetInvertedMass(), rhs);
				mr_matrix_container.Add_Minv_value(mphi_n1, mphi_n,   0.5*delta_t  , mr_matrix_container.GetInvertedMass(),rhs);

				//third step
				CalculateAdvectiveVelocity(mUn, mUn1,mA, coefficient);
				mr_matrix_container.SetToZero(rhs);
				CalculateRHS( mphi_n1,mA,rhs);
				mr_matrix_container.Add_Minv_value(mWork,mWork,     delta_t/3.0    , mr_matrix_container.GetInvertedMass(), rhs);
				mr_matrix_container.Add_Minv_value(mphi_n1, mphi_n,     delta_t    , mr_matrix_container.GetInvertedMass(), rhs);

				//fourth step
				CalculateAdvectiveVelocity(mUn, mUn1,mA, coefficient);
				mr_matrix_container.SetToZero(rhs);
				CalculateRHS( mphi_n1,mA,rhs );
				mr_matrix_container.Add_Minv_value(mWork,mWork,   delta_t/6.0 , mr_matrix_container.GetInvertedMass(), rhs);
				//compute right-hand side
				mr_matrix_container.AssignVectorToVector(mWork,mphi_n1);

				mr_matrix_container.WriteScalarToDatabase(DISTANCE, mphi_n1, mr_model_part.Nodes());



			KRATOS_CATCH("")
			}

					

			//*********************************************************************
			//function to calculate right-hand side of fractional momentum equation
			void CalculateRHS(
					const ValuesVectorType& mphi, 
					const CalcVectorType& convective_velocity, 
					ValuesVectorType& rhs)
			{
				KRATOS_TRY

				int n_nodes = mphi.size();
				
				//calculating the convective projection
				#pragma omp parallel for 
				for (int i_node = 0; i_node < n_nodes; i_node++)
				{
					 array_1d<double, TDim>& pi_i = mPi[i_node];  
					 
					//setting to zero the projection
					for (unsigned int l_comp = 0; l_comp < TDim; l_comp++)
						pi_i[l_comp] = 0.0;
					
// 					double& pi_i = mPi[i_node];  
					const double& phi_i = mphi[i_node];

					const array_1d<double, TDim> a_i = convective_velocity[i_node];
					
					//loop to all the edges surrounding node I
					for (unsigned int csr_index=mr_matrix_container.GetRowStartIndex()[i_node]; csr_index!=mr_matrix_container.GetRowStartIndex()[i_node+1]; csr_index++)
					{
						unsigned int j_neighbour = mr_matrix_container.GetColumnIndex()[csr_index];

						const array_1d<double, TDim>& a_j = convective_velocity[j_neighbour];
						const double& phi_j = mphi[j_neighbour];
						
						
						CSR_Tuple& edge_ij = mr_matrix_container.GetEdgeValues()[csr_index];
						
//  						edge_ij.Add_ConvectiveContribution(pi_i,a_i,phi_i,a_j,phi_j);						
						edge_ij.Add_grad_p(pi_i,phi_i,phi_j);
						
					}

					//apply inverted mass matrix
					const double m_inv = mr_matrix_container.GetInvertedMass()[i_node];
// 					pi_i *= m_inv;

					for (unsigned int l_comp = 0; l_comp < TDim; l_comp++)
						pi_i[l_comp] *= m_inv;
// KRATOS_WATCH(pi_i);
				}
				
				//calculating limitor
				array_1d<double, TDim> dir;
				#pragma omp parallel for private(dir)
				for (int i_node = 0; i_node < n_nodes; i_node++)
				{
					const array_1d<double, TDim>& x_i = mx[i_node];
					const array_1d<double, TDim>& proj_i = mPi[i_node];
					const double& p_i = mphi[i_node];
					double& beta_i = mBeta[i_node];
					beta_i = 0.0;
					double n = 0.0;
					double h=0.0;

					    for (unsigned int csr_index = mr_matrix_container.GetRowStartIndex()[i_node]; csr_index != mr_matrix_container.GetRowStartIndex()[i_node + 1]; csr_index++)
					    {
						  unsigned int j_neighbour = mr_matrix_container.GetColumnIndex()[csr_index];
						    CSR_Tuple& edge_ij = mr_matrix_container.GetEdgeValues()[csr_index];
						    const double& p_j = mphi[j_neighbour];
						    
						    const array_1d<double, TDim>& x_j = mx[j_neighbour];
						    
						    for (unsigned int l_comp = 0; l_comp < TDim; l_comp++)
							    dir[l_comp] = x_j[l_comp] - x_i[l_comp];
						     						    
						    double lenght = dir[0]*dir[0];
 						    for (unsigned int l_comp = 1; l_comp < TDim; l_comp++)
 							    lenght += dir[l_comp]*dir[l_comp];
 						    lenght = sqrt(lenght);
						
						    const array_1d<double, TDim>& proj_j = mPi[j_neighbour];

						    double proj = 0.0;
						    for (unsigned int comp = 0; comp < TDim; comp++)
							proj += 0.5*dir[comp]*(proj_i[comp]+proj_j[comp]);
// 							proj += dir[comp]*pi_i[comp];
						    
						    double numerator = fabs( fabs(p_j - p_i) - fabs(proj) );
						    double denom = fabs( fabs(p_j - p_i)  + 1e-6);
						    
 						    double beta = numerator/denom;

						    beta_i += beta;
						    n += 1.0;
						    h += lenght;
						    
/*						    if(beta_i < beta)
							beta_i = beta;*/
					    }
					   
					beta_i /= n;
					h /=n;
					
					if(beta_i > 1.0)
							beta_i = 1.0;
					
					beta_i*=h;
//  KRATOS_WATCH(beta_i);
				}
				
				//perform MPI syncronization

				//calculating the RHS
				double stab_low;
				double stab_high;
				array_1d<double, TDim> aux;
				#pragma omp parallel for private(stab_low,stab_high,aux) 
				for ( int i_node = 0; i_node < n_nodes; i_node++)
				{
					double& rhs_i = rhs[i_node];
					const double& phi_i = mphi[i_node];
					const double& beta_i = mBeta[i_node];
					const array_1d<double, TDim>& a_i = convective_velocity[i_node]; 
					const array_1d<double, TDim>& proj_i = mPi[i_node];
// 					const array_1d<double, TDim>& x_i = mx[i_node];

					double pi_i = proj_i[0]*a_i[0];
					for (unsigned int l_comp = 1; l_comp < TDim; l_comp++)
						pi_i += proj_i[l_comp]*a_i[l_comp];
					//double& h_i = mHmin[i_node];
					
					double norm_a = a_i[0]*a_i[0];
					for (unsigned int l_comp = 1; l_comp < TDim; l_comp++)
						norm_a += a_i[l_comp]*a_i[l_comp];
					norm_a = sqrt(norm_a);

					
					//initializing with the external forces (e.g. gravity)
					rhs_i = 0.0;

					//loop to all the edges surrounding node I
					for (unsigned int csr_index=mr_matrix_container.GetRowStartIndex()[i_node]; csr_index!=mr_matrix_container.GetRowStartIndex()[i_node+1]; csr_index++)
					{
						unsigned int j_neighbour = mr_matrix_container.GetColumnIndex()[csr_index];
											
// 						const array_1d<double, TDim>& x_j = mx[j_neighbour];

						//double& rhs_j = rhs[j_neighbour];
						const double& phi_j = mphi[j_neighbour];
						const array_1d<double, TDim>& a_j = convective_velocity[j_neighbour]; 

						const array_1d<double, TDim>& proj_j = mPi[j_neighbour];
						double pi_j = proj_j[0]*a_i[0];
						for (unsigned int l_comp = 1; l_comp < TDim; l_comp++)
						  pi_j += proj_j[l_comp]*a_i[l_comp];
						
						//double& h_j = mHmin[j_neighbour];
						
						CSR_Tuple& edge_ij = mr_matrix_container.GetEdgeValues()[csr_index];

						//convection operator						
						edge_ij.Sub_ConvectiveContribution(rhs_i,a_i,phi_i,a_j,phi_j);
						
						//calculate stabilization part
						edge_ij.CalculateConvectionStabilization_LOW( stab_low,a_i,phi_i,a_j,phi_j);
												
						double edge_tau = mTau[i_node];
						
						
						
						edge_ij.CalculateConvectionStabilization_HIGH( stab_high,a_i,pi_i,a_j,pi_j);
						
 						edge_ij.Sub_StabContribution( rhs_i, edge_tau, 1.0, stab_low, stab_high);
/*						
						double laplacian_ij=0.0;
						edge_ij.CalculateScalarLaplacian( laplacian_ij );
						double capturing= laplacian_ij * (phi_j - phi_i);
  						rhs_i-= 0.1*0.5*capturing*beta_i*norm_a*mHmin[i_node];*/
						
						
						
// 						 for (unsigned int l_comp = 0; l_comp < TDim; l_comp++)
// 							    dir[l_comp] = x_j[l_comp] - x_i[l_comp];
// 						    
// 						    double proj = 0.0;
// 						    for (unsigned int comp = 0; comp < TDim; comp++)
// 							proj += 0.5*dir[comp]*(proj_i[comp]+proj_j[comp]);
// // 							proj += dir[comp]*proj_i[comp];
// 						    
// 						    double numerator = fabs( fabs(pi_j - pi_i) - fabs(proj) );
// 						    double denom = fabs( fabs(pi_j - pi_i) + fabs(proj) + 1e-6);
// 						    
//  						    double beta = numerator/denom;
// 						 
// 						 if(beta > 1.0)
// 							beta = 1.0;
						 
						 double beta = beta_i;
						 
						 double coeff = 0.35; //=0.5*0.5;
						 double laplacian_ij=0.0;
						edge_ij.CalculateScalarLaplacian( laplacian_ij );
						double capturing= laplacian_ij * (phi_j - phi_i);
   						rhs_i-= coeff*capturing*beta*norm_a;
												
// 						double aaa = 0.0;
// 						for (unsigned int k_comp = 0; k_comp < TDim; k_comp++)
// 							for (unsigned int m_comp = 0; m_comp < TDim; m_comp++)
// 								aaa += a_i[k_comp] * a_i[m_comp] * edge_ij.LaplacianIJ(k_comp,m_comp);
// 						aaa/=(norm_a*norm_a);
// 						double capturing2 = aaa * (phi_j - phi_i);
// 						
//     						rhs_i-= coeff*(capturing - capturing2)*beta*norm_a;
						//

					}

// KRATOS_WATCH(rhs_i);

				}



				KRATOS_CATCH("")
			}


			void CalculateAdvectiveVelocity(
					const CalcVectorType&  mUn,
					const CalcVectorType&  mUn1,
					CalcVectorType&   mA, 
					double coefficient)
			{
				int n_nodes = mUn1.size();
				#pragma omp parallel for 
				for (int i_node = 0; i_node < n_nodes; i_node++)
				{	
					
					//reference for advective velocity of node i
					array_1d<double, TDim>&  a_i 	= mA[i_node];
					const array_1d<double, TDim>&  Un_i 	= mUn[i_node];
					const array_1d<double, TDim>&  Un1_i 	= mUn1[i_node];
	
					for (unsigned int k_comp = 0; k_comp < TDim; k_comp++)
						a_i[k_comp] = coefficient * Un1_i[k_comp] + (1.0 - coefficient)* Un_i[k_comp];

				}


			}
						






			//*******************************
			//function to free dynamic memory
			void Clear()
			{
			KRATOS_TRY

				mWork.clear();
				mPi.clear();
				mUn.clear();
				mUn1.clear();
 				mA.clear();
				mphi_n.clear();
				mphi_n1.clear();
 				mHmin.clear();
 				mTau.clear();
				mBeta.clear();
				mx.clear();

			KRATOS_CATCH("")
			}
			


		private:
			MatrixContainer& mr_matrix_container;
			ModelPart& mr_model_part;
			
   			bool msmooth_convective_velocity;
			bool minclude_shock_capturing;
			
			//nodal values
			//velocity vector U at time steps n and n+1
			CalcVectorType  mUn1,mUn;
			CalcVectorType  mPi;
			//pressure vector p at time steps n and n+1
			CalcVectorType mx;
			ValuesVectorType mBeta;
			ValuesVectorType mWork;
			ValuesVectorType mphi_n, mphi_n1; //variable to be convected

			//advective velocity vector
 			CalcVectorType mA;

			//minimum length of the edges surrounding edges surrounding each nodal point
			ValuesVectorType mHmin;

			//flag for first time step
			bool mFirstStep;

			//intrinsic time step size
			ValuesVectorType mTau;
	
	};
} //namespace Kratos

#endif //KRATOS_PURE_CONVECTION_EDGEBASED_SOLVER_H_INCLUDED defined


