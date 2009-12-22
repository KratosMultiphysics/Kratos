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
//   Last Modified by:    $Author: antonia $
//   Date:                $Date: 2009-01-14 16:24:38 $
//   Revision:            $Revision: 1.11 $
//
//


#if !defined(KRATOS_LEVELSET_FLUID_SOLVER_H_INCLUDED)
#define  KRATOS_LEVELSET_FLUID_SOLVER_H_INCLUDED

  #define SPLIT_OSS
//        #define SYMM_PRESS


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
	class LevelSetFluidSolver
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
			LevelSetFluidSolver(MatrixContainer& mr_matrix_container,
				    ModelPart& mr_model_part,
				    bool include_shock_capturing,
				    bool smooth_convective_velocity
				    )
			: mr_matrix_container(mr_matrix_container),mr_model_part(mr_model_part)
			{
				//options
				minclude_shock_capturing = include_shock_capturing;
				msmooth_convective_velocity = smooth_convective_velocity;	
			
			};
			~LevelSetFluidSolver(){};

			//***********************************
			//function to initialize fluid solver
			void Initialize(
				       )
			{
			KRATOS_TRY
				
			
				//get number of nodes
				unsigned int n_nodes = mr_model_part.Nodes().size();
				unsigned int n_edges = mr_matrix_container.GetNumberEdges();
				//size data vectors
				mWork.resize(n_nodes);
				mvel_n.resize(n_nodes);
				mvel_n1.resize(n_nodes);
				mInitMom.resize(n_nodes);
				mCurrMom.resize(n_nodes);
				mPn.resize(n_nodes);
				mPn1.resize(n_nodes);
				mViscosity.resize(n_nodes);
				mRho.resize(n_nodes);
				mRhoOld.resize(n_nodes);
				mC2inv.resize(n_nodes);
				mA.resize(n_nodes);
				mHmin.resize(n_nodes);
				mHavg.resize(n_nodes);
				mNodalFlag.resize(n_nodes);
				mdistances.resize(n_nodes);
				mEps.resize(n_nodes);
				mEpsOld.resize(n_nodes);
				mD.resize(n_nodes);

				mTauPressure.resize(n_nodes);
				mTauConvection.resize(n_nodes);
				
				mPi.resize(n_nodes);
				mXi.resize(n_nodes);
				mBodyForce.resize(n_nodes);
				mDrag.resize(n_nodes);
				mx.resize(n_nodes);

				mCp.resize(n_nodes);

				mMach.resize(n_nodes);

				mEdgeDimensions.resize(n_edges);
				mBeta.resize(n_edges);
				for (unsigned int csr_index = 0; csr_index < n_edges; csr_index++)
					mBeta[csr_index] = 1.0;

				ValuesVectorType external_pressure;
				external_pressure.resize(n_nodes);

				//read velocity and pressure data from Kratos
				mr_matrix_container.FillVectorFromDatabase(VELOCITY, mvel_n, mr_model_part.Nodes());
				mr_matrix_container.FillScalarFromDatabase(EXTERNAL_PRESSURE, external_pressure, mr_model_part.Nodes());
				mr_matrix_container.FillScalarFromDatabase(IS_BOUNDARY, mNodalFlag, mr_model_part.Nodes());
				mr_matrix_container.FillScalarFromDatabase(DENSITY, mRho, mr_model_part.Nodes());
				mr_matrix_container.FillScalarFromDatabase(PRESSURE, mPn1, mr_model_part.Nodes());
				mr_matrix_container.FillOldScalarFromDatabase(PRESSURE, mPn, mr_model_part.Nodes());
				mr_matrix_container.FillVectorFromDatabase(VELOCITY, mvel_n1, mr_model_part.Nodes());
				//mr_matrix_container.FillCoordinatesFromDatabase(mx, mr_model_part.Nodes());
				//set flag for first time step
				mFirstStep = true;

				//loop to categorize boundary nodes
				for (unsigned int i_node = 0; i_node < n_nodes; i_node++)
				{
					//differentiate between types of boundary condition
					switch (static_cast<unsigned int>(mNodalFlag[i_node]))
					{
						case 1:
							//velocity inlet
							mVelocityInletList.push_back(i_node);
							mVelocityInlet.push_back(mvel_n[i_node]);
							mDensityInlet.push_back(mRho[i_node]);
							mDissipationList.push_back(i_node);
							break;
						case 2:
							//no-slip condition
							mNoSlipBoundaryList.push_back(i_node);
							break;
						case 3:
							//slip condition
							mSlipBoundaryList.push_back(i_node);
							break;
						case 4:
							//mixed condition (slip and pressure node)
							mPressureOutletList.push_back(i_node);
							mPressureOutlet.push_back(external_pressure[i_node]);
							mSlipBoundaryList.push_back(i_node);
							mDissipationList.push_back(i_node);
							break;
						case 5:
							//pressure outlet
							mPressureOutletList.push_back(i_node);
							mPressureOutlet.push_back(external_pressure[i_node]);
							mDissipationList.push_back(i_node);
							break;
					}
				}

				//print number of nodes corresponding to the different types of boundary conditions
				KRATOS_WATCH(mVelocityInletList.size())
				KRATOS_WATCH(mDensityInlet.size())
				KRATOS_WATCH(mPressureOutletList.size())
				KRATOS_WATCH(mSlipBoundaryList.size())
				KRATOS_WATCH(mNoSlipBoundaryList.size())
				KRATOS_WATCH(mDissipationList.size())


				//determine number of edges and entries
				unsigned int n_nonzero_entries = 2 * n_edges + n_nodes;
				//allocate memory for variables
				mL.resize(n_nodes,n_nodes,n_nonzero_entries);

				//loop over all nodes
				for (unsigned int i_node = 0; i_node < n_nodes; i_node++)
				{
					//flag for considering diagonal matrix elements
					bool flag = 0;

					//loop over all neighbours
					for (unsigned int csr_index=mr_matrix_container.GetRowStartIndex()[i_node]; csr_index!=mr_matrix_container.GetRowStartIndex()[i_node+1]; csr_index++)
					{
						//get global index of neighbouring node j
						unsigned int j_neighbour = mr_matrix_container.GetColumnIndex()[csr_index];
						//define matrix structure row by row (the order does matter!)
						if ((j_neighbour > i_node) && (flag == 0))
						{
							//add diagonal/nodal contribution
							mL.push_back(i_node, i_node, 0.0);
							flag = 1;
						}
						//add non-diagonal/edge contribution
						mL.push_back(i_node, j_neighbour, 0.0);
					}
					//if diagonal element is the last non-zero element of the row
					if (flag == 0)
						mL.push_back(i_node, i_node, 0.0);
				}

				//compute area normals
				CalculateNormals(mr_model_part.Conditions());
//  				WriteVectorToDatabase(NORMAL, mPressureNormal, mr_model_part.Nodes());
				mr_matrix_container.WriteVectorToDatabase(NORMAL, mSlipNormal, mr_model_part.Nodes());

				//compute minimum length of the surrounding edges
				CalculateEdgeLengths(mr_model_part.Nodes());

				//prepare initial momentum for first time step
				for (unsigned int i_node = 0; i_node < n_nodes; i_node++)
				{
					double& rho_i = mRho[i_node];
					array_1d<double, TDim>& u_i = mvel_n1[i_node];
					array_1d<double, TDim>& U_i = mInitMom[i_node];

					//compute initial momentum for iteration of step 1
					for (unsigned int component = 0; component < TDim; component++)
						U_i[component] = rho_i * u_i[component];
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
			
				//getting value of current velocity and of viscosity
				mr_matrix_container.FillVectorFromDatabase(VELOCITY, mvel_n1, mr_model_part.Nodes());
				mr_matrix_container.FillScalarFromDatabase(VISCOSITY, mViscosity, mr_model_part.Nodes());
				mr_matrix_container.FillCoordinatesFromDatabase(mx, mr_model_part.Nodes());
//*******************				
				//loop over all nodes
				double n_nodes = mvel_n1.size();
				for (unsigned int i_node = 0; i_node < n_nodes; i_node++)
				{
					array_1d<double, TDim>& v_i = mvel_n1[i_node];
// 					KRATOS_WATCH(v_i);
					array_1d<double, TDim>& x_i = mx[i_node];
// 					KRATOS_WATCH(x_i);

					//use CFL condition to compute time step size
					double delta_t_i = CFLNumber * 1.0 / (norm_2(v_i)/mHmin[i_node] + 2.0 * mViscosity[i_node]/(mHmin[i_node]*mHmin[i_node]) );
					
					//considering the most restrictive case of neighbor's velocities with similar direction but opposite sense.
					//loop over all neighbours
					for (unsigned int csr_index=mr_matrix_container.GetRowStartIndex()[i_node]; csr_index!=mr_matrix_container.GetRowStartIndex()[i_node+1]; csr_index++)
					{
						//get global index of neighbouring node j
						unsigned int j_neighbour = mr_matrix_container.GetColumnIndex()[csr_index];
						
						array_1d<double, TDim>& v_j = mvel_n1[j_neighbour];
						array_1d<double, TDim>& x_j = mx[j_neighbour];
						array_1d<double, TDim> edge_dir = ZeroVector(TDim);
// 						KRATOS_WATCH(x_j);
// 						KRATOS_WATCH(v_j);
						//Calculate edge direction
						edge_dir[0] = x_j[0]-x_i[0];
						edge_dir[1] = x_j[1]-x_i[1];
						edge_dir[2] = x_j[2]-x_i[2];
						
// 						KRATOS_WATCH(edge_dir);

						double aux = norm_2(edge_dir);
// 						KRATOS_WATCH(aux);

						if (aux == 0.0)
						{	edge_dir = ZeroVector(TDim);}
						else
						{	//normalized edge direction
							edge_dir /= norm_2(edge_dir);
// 							KRATOS_WATCH(edge_dir);
						}
						//int aux = inner_prod(v_i,v_j);
						double v_i_par = inner_prod(v_i, edge_dir);
						double v_j_par = inner_prod(v_j, edge_dir);
// 						KRATOS_WATCH(v_i_par);
// 						KRATOS_WATCH(v_j_par);

						if ((v_i_par >= 0.0 && v_j_par <= 0.0) || (v_i_par <= 0.0 && v_j_par >= 0.0))
						{	
							double delta_t_j = CFLNumber / ((fabs(v_i_par) + fabs(v_j_par))/mHmin[i_node] + 2.0 * mViscosity[i_node]/(mHmin[i_node]*mHmin[i_node]));
// 							KRATOS_WATCH(delta_t_j);
// 							KRATOS_WATCH(delta_t_i);
							if (delta_t_j < delta_t_i)
								delta_t_i = delta_t_j;
						}
					}
					//choose the overall minimum of delta_t_i
					if (delta_t_i < delta_t)
						delta_t = delta_t_i;
				}
//*******************
				//perform MPI syncronization of the dt (minimum should be kept) 

				//write time step size to Kratos
				ProcessInfo& CurrentProcessInfo = mr_model_part.GetProcessInfo();
				CurrentProcessInfo[DELTA_TIME] = delta_t;

			KRATOS_CATCH("")
			}

			//**********************************************************************************
			//function to solve fluid equations - fractional step 1: compute fractional momentum
			Vector SolveStep1()
			{
			KRATOS_TRY

				//PREREQUISITES

				//variables for node based data handling
				ModelPart::NodesContainerType& rNodes = mr_model_part.Nodes();
				int n_nodes = rNodes.size();
				//storage of nodal values in local variables
				CalcVectorType rhs;
				rhs.resize(n_nodes);

				//read velocity and pressure data from Kratos
				mr_matrix_container.FillVectorFromDatabase(VELOCITY, mvel_n1, rNodes);
				mr_matrix_container.FillOldVectorFromDatabase(VELOCITY, mvel_n, rNodes);
				
				mr_matrix_container.FillScalarFromDatabase(PRESSURE, mPn1, rNodes);
				mr_matrix_container.FillOldScalarFromDatabase(PRESSURE, mPn, rNodes);
				
				mr_matrix_container.FillScalarFromDatabase(DENSITY, mRho, rNodes);
				mr_matrix_container.FillOldScalarFromDatabase(DENSITY, mRhoOld, rNodes);
				
				mr_matrix_container.FillVectorFromDatabase(BODY_FORCE, mBodyForce, rNodes);
				mr_matrix_container.FillScalarFromDatabase(VISCOSITY, mViscosity, rNodes);

				mr_matrix_container.FillScalarFromDatabase(DISTANCE, mdistances, mr_model_part.Nodes());

				mr_matrix_container.FillScalarFromDatabase(POROSITY, mEps, rNodes);
				mr_matrix_container.FillOldScalarFromDatabase(POROSITY, mEpsOld, rNodes);

				//read time step size from Kratos
				ProcessInfo& CurrentProcessInfo = mr_model_part.GetProcessInfo();
				double delta_t = CurrentProcessInfo[DELTA_TIME];

				#pragma omp parallel for 	
				for ( int i_node = 0; i_node < n_nodes; i_node++)
				{
					// -> mCurrMom
					//compute the momentum at the current step -> mCurrMom
					double& rho_i = mRho[i_node];

					const array_1d<double, TDim>& u_i = mvel_n1[i_node];
					array_1d<double, TDim>& U_i = mCurrMom[i_node];
					for (unsigned int comp = 0; comp < TDim; comp++)
						U_i[comp] = rho_i * u_i[comp];
					
					// -> mInitMom
					double& rho_i_old = mRhoOld[i_node];

					//compute the momentum at the beginning of the step
					const array_1d<double, TDim>& u_i_old = mvel_n[i_node];
					array_1d<double, TDim>& U_i_old = mInitMom[i_node];
					for (unsigned int comp = 0; comp < TDim; comp++)
						U_i_old[comp] = rho_i_old * u_i_old[comp];
					
					//compute volumetric body force
					array_1d<double, TDim>& f_i = mBodyForce[i_node];
					for (unsigned int comp = 0; comp < TDim; comp++)
						f_i[comp] *= rho_i;
				}


				DivideByPorosity(mCurrMom, mCurrMom, mEps);
				DivideByPorosity(mInitMom, mInitMom, mEpsOld);

				DivideByPorosity(mvel_n, mvel_n, mEpsOld);
				DivideByPorosity(mvel_n1, mvel_n1, mEps);

				//compute advective velocity - area average of the current velocity
				CalculateAdvectiveVelocity(mvel_n1, mA, msmooth_convective_velocity);

				//compute intrinsic time
				double time_inv = 1.0/delta_t;
//     time_inv = 0.0;
				#pragma omp parallel for firstprivate(time_inv)	
				for (int i_node = 0; i_node < n_nodes; i_node++)
				{
// 					double& h_i = mHavg[i_node];
					double& h_i = mHmin[i_node];
					array_1d<double, TDim>& a_i = mA[i_node];
					const double nu_i = mViscosity[i_node];

//					mTau[i_node] = 1.0 / (0.5 * norm_2(a_i)/h_i + time_inv);
 					double vel_norm = norm_2(a_i);
//   					mTauPressure[i_node] = 1.0 / (2.0 * vel_norm/h_i + 0.01*time_inv + nu_i /(h_i*h_i) );
 					mTauPressure[i_node] = delta_t;
					mTauConvection[i_node] = 1.0 / (2.0 * vel_norm/h_i +  0.01*time_inv + nu_i /(h_i*h_i) );
					
					if (mTauPressure[i_node] < delta_t)
						mTauPressure[i_node] =  delta_t;
					else if(mTauPressure[i_node] > 100.0*delta_t)
						mTauPressure[i_node] = 100.0*delta_t;
				}

				//compute pressure switch
				if (mFirstStep == false)
					if(minclude_shock_capturing == true)
						ComputeMonotonicityPreserving();
				
				mr_matrix_container.AssignVectorToVector(mInitMom,mWork); //mWork = mvel_n NO!!!->  mWork = mU_iold
				
				//first step of Runge Kutta
				mr_matrix_container.AssignVectorToVector(mvel_n,mvel_n1); //mvel_n1 = mvel_n
				mr_matrix_container.AssignVectorToVector(mInitMom,mCurrMom); 

// double start_prod = omp_get_wtime();				
				CalculateAdvectiveVelocity(mvel_n1, mA, msmooth_convective_velocity);
				mr_matrix_container.SetToZero(rhs);
				CalculateRHS( mCurrMom, mPn1, mA, mBodyForce, mViscosity, rhs);
/*double norma=0.0;
for (int i_node = 0; i_node < n_nodes; i_node++)
	for (int kkk = 0; kkk < TDim; kkk++)
		norma += rhs[i_node][kkk]*rhs[i_node][kkk];
KRATOS_WATCH(norma);*/
				
				mr_matrix_container.Add_Minv_value(mWork,mWork,     delta_t/6.0    , mr_matrix_container.GetInvertedMass(), rhs);
				mr_matrix_container.Add_Minv_value(mCurrMom,mInitMom,   0.5*delta_t  , mr_matrix_container.GetInvertedMass(), rhs);
				ApplyVelocityBC(mCurrMom);
/*mr_matrix_container.WriteVectorToDatabase(CONV_PROJ, mA, rNodes);
mr_matrix_container.WriteScalarToDatabase(TEMPERATURE, mTauConvection, rNodes);*/
				
				//second step
				CalculateVelocity(mvel_n1,mCurrMom,mRho);
				CalculateAdvectiveVelocity( mvel_n1, mA, msmooth_convective_velocity);
				mr_matrix_container.SetToZero(rhs);
				CalculateRHS( mCurrMom, mPn1, mA, mBodyForce,mViscosity, rhs );
				mr_matrix_container.Add_Minv_value(mWork,mWork,  delta_t/3.0    , mr_matrix_container.GetInvertedMass(), rhs);
				mr_matrix_container.Add_Minv_value(mCurrMom,mInitMom,   0.5*delta_t  , mr_matrix_container.GetInvertedMass(),rhs);
				ApplyVelocityBC(mCurrMom);
				
				//third step
				CalculateVelocity(mvel_n1,mCurrMom,mRho);
				CalculateAdvectiveVelocity( mvel_n1, mA, msmooth_convective_velocity);
				mr_matrix_container.SetToZero(rhs);
				CalculateRHS( mCurrMom, mPn1, mA, mBodyForce,mViscosity, rhs);
				mr_matrix_container.Add_Minv_value(mWork,mWork,     delta_t/3.0    , mr_matrix_container.GetInvertedMass(), rhs);
				mr_matrix_container.Add_Minv_value(mCurrMom,mInitMom,     delta_t    , mr_matrix_container.GetInvertedMass(), rhs);
				ApplyVelocityBC(mCurrMom);
				
				//fourth step
				CalculateVelocity(mvel_n1,mCurrMom,mRho);
				CalculateAdvectiveVelocity( mvel_n1, mA, msmooth_convective_velocity);
				mr_matrix_container.SetToZero(rhs);
				CalculateRHS( mCurrMom, mPn1, mA, mBodyForce,mViscosity, rhs );
				mr_matrix_container.Add_Minv_value(mWork,mWork,   delta_t/6.0 , mr_matrix_container.GetInvertedMass(), rhs);
				ApplyVelocityBC(mCurrMom);

				//compute right-hand side
				mr_matrix_container.AssignVectorToVector(mWork,mCurrMom);
 				ApplyVelocityBC(mCurrMom);


// 				//compute ratio for iteration
 				Vector stop_criteria(TDim);
				noalias(stop_criteria) = ZeroVector(TDim);
//  				stop_criteria[0] = 0.0;
//  				stop_criteria[1] = 0.0;

				return stop_criteria;

			KRATOS_CATCH("")
			}

					

			//*********************************************************************
			//function to calculate right-hand side of fractional momentum equation
			void CalculateRHS(
					const CalcVectorType& momentum, 
					const ValuesVectorType& pressure, 
					const CalcVectorType& convective_velocity, 
					const CalcVectorType& body_force, 
     					const ValuesVectorType& viscosity,
					CalcVectorType& rhs)
			{
				KRATOS_TRY

				int n_nodes = momentum.size();
				
				//calculating the convective projection
				#pragma omp parallel for 
				for (int i_node = 0; i_node < n_nodes; i_node++)
				{
					array_1d<double, TDim>& pi_i = mPi[i_node];  //******************

					//setting to zero
					for (unsigned int l_comp = 0; l_comp < TDim; l_comp++)
						pi_i[l_comp] = 0.0;
					
					const array_1d<double, TDim>& a_i = convective_velocity[i_node];
					const array_1d<double, TDim>& U_i = momentum[i_node];
					
					//const double& p_i = pressure[i_node];

					for (unsigned int csr_index=mr_matrix_container.GetRowStartIndex()[i_node]; csr_index!=mr_matrix_container.GetRowStartIndex()[i_node+1]; csr_index++)
					{
						unsigned int j_neighbour = mr_matrix_container.GetColumnIndex()[csr_index];
						const array_1d<double, TDim>& a_j = convective_velocity[j_neighbour];
						const array_1d<double, TDim>& U_j = momentum[j_neighbour];
						//const double& p_j = pressure[j_neighbour];
						
						
						CSR_Tuple& edge_ij = mr_matrix_container.GetEdgeValues()[csr_index];
						
 						edge_ij.Add_ConvectiveContribution(pi_i,a_i,U_i,a_j,U_j);
						
//         						edge_ij.Add_grad_p(pi_i,p_i,p_j);
//     						edge_ij.Sub_grad_p(pi_i,p_i,p_j);
					}
					const double m_inv = mr_matrix_container.GetInvertedMass()[i_node];
					
					for (unsigned int l_comp = 0; l_comp < TDim; l_comp++)
						pi_i[l_comp] *= m_inv;
				}
				

				
				
				//perform MPI syncronization

				//calculating the RHS
				array_1d<double,TDim> stab_low;
				array_1d<double,TDim> stab_high;
				#pragma omp parallel for private(stab_low,stab_high) 
				for ( int i_node = 0; i_node < n_nodes; i_node++)
				{
					double dist = mdistances[i_node];
					if (dist < 0.0) //node is inside domain ---- if outside do nothing
					{
						array_1d<double, TDim>& rhs_i = rhs[i_node];
						const array_1d<double, TDim>& f_i = body_force[i_node];
						const array_1d<double, TDim>& a_i = convective_velocity[i_node]; 
						const array_1d<double, TDim>& U_i = momentum[i_node];
						const array_1d<double, TDim>& pi_i = mPi[i_node];
						const double& p_i = pressure[i_node];
						const double& nu_i = viscosity[i_node];
						//double& h_i = mHmin[i_node];
	
						//initializing with the external forces (e.g. gravity)
						double& m_i = mr_matrix_container.GetLumpedMass()[i_node];
						for (unsigned int comp = 0; comp < TDim; comp++)
							rhs_i[comp] = m_i * f_i[comp];

						//porous contribution
						double eps = mEps[i_node];
						double d = mD[i_node]; //diameter of the particle	
						double kinv = 150.0*(1.0-eps)*(1.0-eps)/(eps*eps*eps*d*d);
						double norm_u_2 = 0.0;
						for (unsigned int comp = 0; comp < TDim; comp++)
							norm_u_2 = a_i[comp]*a_i[comp];
// 							norm_u_2 = U_i[comp]*U_i[comp];

						//CORRECTED Term
						double nonlin_term = kinv * nu_i * eps + 1.75 * sqrt(norm_u_2 *  kinv / (eps * 150.0));
						//ERROR IN WRITING THE NON LINEAR TERM//
// 						double nonlin_term = kinv * nu_i * eps + 1.75 * norm_u_2  * sqrt(kinv / ( eps * 150.0));
						for (unsigned int comp = 0; comp < TDim; comp++)
							rhs_i[comp] -= m_i * nonlin_term * U_i[comp];

						
	
						//convective term
						for (unsigned int csr_index=mr_matrix_container.GetRowStartIndex()[i_node]; csr_index!=mr_matrix_container.GetRowStartIndex()[i_node+1]; csr_index++)
						{
							unsigned int j_neighbour = mr_matrix_container.GetColumnIndex()[csr_index];
	
							const array_1d<double, TDim>& a_j = convective_velocity[j_neighbour];
							const array_1d<double, TDim>& U_j = momentum[j_neighbour];
							const array_1d<double, TDim>& pi_j = mPi[j_neighbour];
							const double& p_j = pressure[j_neighbour];
							const double& nu_j = viscosity[j_neighbour];
							
							CSR_Tuple& edge_ij = mr_matrix_container.GetEdgeValues()[csr_index];
							
							edge_ij.Sub_ConvectiveContribution(rhs_i,a_i,U_i,a_j,U_j);
							
							//take care! we miss including a B.C.  for the external pressure
							edge_ij.Add_Gp(rhs_i,p_i,p_j);
							
							edge_ij.Sub_ViscousContribution(rhs_i,U_i,nu_i,U_j,nu_j);
							
							//add stabilization
// 							edge_ij.CalculateConvectionStabilization_LOW( stab_low,a_i,U_i,p_i,a_j,U_j,p_j);
							edge_ij.CalculateConvectionStabilization_LOW( stab_low,a_i,U_i,a_j,U_j);
														
							double edge_tau = mTauConvection[i_node];
								
							edge_ij.CalculateConvectionStabilization_HIGH( stab_high,a_i,pi_i,a_j,pi_j);
							
							double beta = mBeta[csr_index];
							edge_ij.Sub_StabContribution( rhs_i, edge_tau, beta, stab_low, stab_high);
	
						}
					}

				}


				//boundary integrals --> finishing the calculation of the pressure gradient
				int loop_size1 = mPressureOutletList.size();
				#pragma omp parallel for 
				for (int i_pressure = 0; i_pressure < loop_size1; i_pressure++)
				{
					unsigned int i_node = mPressureOutletList[i_pressure];
					array_1d<double, TDim>& rhs_i = rhs[i_node];
					const double& p_ext_i = mPressureOutlet[i_pressure];
					const array_1d<double, TDim>& an_i = mPressureNormal[i_node];

 					for (unsigned int comp = 0; comp < TDim; comp++)
 						rhs_i[comp] -= an_i[comp] *  p_ext_i;
					
					
// const array_1d<double, TDim>& U_i = momentum[i_node];
// const array_1d<double, TDim>& a_i = convective_velocity[i_node]; 
// double temp = 0.0;
// double scalar_prod = 0.0;
// for (unsigned int comp = 0; comp < TDim; comp++)
// {
// 	scalar_prod += an_i[comp] *  U_i[comp];
// 	temp += an_i[comp] *  an_i[comp];
// }
// temp = sqrt(temp);
// for (unsigned int comp = 0; comp < TDim; comp++)
// // 	rhs_i[comp] -= U_i[comp] * temp;
// //  	rhs_i[comp] -= an_i[comp] *  scalar_prod / temp;
//   	rhs_i[comp] -= a_i[comp] *  scalar_prod / temp;

					
				}

				KRATOS_CATCH("")
			}

			//*************************************************************************
			//function to solve fluid equations - fractional step 2: calculate pressure
			void SolveStep2(typename TLinearSolver::Pointer pLinearSolver)
			{
			KRATOS_TRY

				//PREREQUISITES

				//allocate memory for variables
				ModelPart::NodesContainerType& rNodes = mr_model_part.Nodes();
				int n_nodes = rNodes.size();
				//unknown and right-hand side vector
				TSystemVectorType dp, rhs;
				dp.resize(n_nodes);
				rhs.resize(n_nodes);
				array_1d<double, TDim> dU_i, dU_j, work_array;
				//read time step size from Kratos
				ProcessInfo& CurrentProcessInfo = mr_model_part.GetProcessInfo();
				double delta_t = CurrentProcessInfo[DELTA_TIME];

#ifdef _OPENMP
				double time_inv = 0.0; //1.0/delta_t;
#endif

#ifdef SPLIT_OSS
				
// 				#pragma omp parallel for firstprivate(time_inv), private(work_array)
				for (int i_node = 0; i_node < n_nodes; i_node++)
				{
					array_1d<double, TDim>& xi_i = mXi[i_node];
					for (unsigned int comp = 0; comp < TDim; comp++)
						xi_i[comp] = 0.0;
					
					const double& p_i = mPn1[i_node];
	
					for (unsigned int csr_index=mr_matrix_container.GetRowStartIndex()[i_node]; 		csr_index!=mr_matrix_container.GetRowStartIndex()[i_node+1]; csr_index++)
					{
						//get global index of neighbouring node j
						unsigned int j_neighbour = mr_matrix_container.GetColumnIndex()[csr_index];
						const double& p_j = mPn1[j_neighbour];

						//projection of pressure gradients
						CSR_Tuple& edge_ij = mr_matrix_container.GetEdgeValues()[csr_index];
						
   						edge_ij.Add_grad_p(xi_i,p_i,p_j);
// // //   						edge_ij.Sub_grad_p(xi_i,p_i,p_j);
					}
	
					const double& m_inv = mr_matrix_container.GetInvertedMass()[i_node];
					for (unsigned int l_comp = 0; l_comp < TDim; l_comp++)
						xi_i[l_comp] *= m_inv;
				}
#endif

				//loop over all nodes
				#pragma omp parallel for firstprivate(time_inv)
				for (int i_node = 0; i_node < n_nodes; i_node++)
				{
					double& rhs_i = rhs[i_node];
					rhs_i = 0.0;
					double& p_i = mPn1[i_node];
					double& eps_i = mEps[i_node];
					array_1d<double, TDim>& U_i_curr = mCurrMom[i_node];
					//array_1d<double, TDim>& a_i = mA[i_node];
					double& rho_i = mRho[i_node];

#ifdef SPLIT_OSS 
					array_1d<double, TDim>& xi_i = mXi[i_node];
#else
					
					array_1d<double, TDim>& pi_i = mPi[i_node];
#endif
					//const double& h_i = mHavg[i_node];
					double l_ii = 0.0;

					//loop over all neighbours
					for (unsigned int csr_index=mr_matrix_container.GetRowStartIndex()[i_node]; csr_index!=mr_matrix_container.GetRowStartIndex()[i_node+1]; csr_index++)
					{
						unsigned int j_neighbour = mr_matrix_container.GetColumnIndex()[csr_index];
						double& p_j = mPn1[j_neighbour];
						double& eps_j = mEps[j_neighbour];
						array_1d<double, TDim>& U_j_curr = mCurrMom[j_neighbour];
						//array_1d<double, TDim>& a_j = mA[j_neighbour];
#ifdef SPLIT_OSS 
						array_1d<double, TDim>& xi_j = mXi[j_neighbour];
#else
 						array_1d<double, TDim>& pi_j = mPi[j_neighbour];
#endif

						//const double& h_j = mHavg[j_neighbour];

						CSR_Tuple& edge_ij = mr_matrix_container.GetEdgeValues()[csr_index];
						
#ifdef SYMM_PRESS
						double edge_tau = 0.5*( mTauPressure[i_node] + mTauPressure[j_neighbour]);
#else
						double edge_tau = mTauPressure[i_node];
#endif
//     						double edge_tau = CalculateEdgeTau(time_inv,h_i,a_i,h_j,a_j); 
//   						
						
						

						//compute laplacian operator
						double sum_l_ikjk;
						edge_ij.CalculateScalarLaplacian(sum_l_ikjk);
						double sum_l_ikjk_onlystab = sum_l_ikjk * (edge_tau);
						sum_l_ikjk *= (delta_t + edge_tau);

						//assemble right-hand side
						//pressure contribution
						rhs_i -= sum_l_ikjk_onlystab * (p_j - p_i);
						
						//other part of the residual
#if !defined(SPLIT_OSS)
						array_1d<double, TDim>& a_j = mA[j_neighbour];
						boost::numeric::ublas::bounded_matrix<double,TDim,TDim>& L = edge_ij.LaplacianIJ;
						for(unsigned int i = 0; i<TDim; i++)
							for(unsigned int j = 0; j<TDim; j++)
								rhs_i -= edge_tau * a_i[j] * L(i,j) * (U_j_curr[j] - U_i_curr[j]);
#endif
							
						
						//calculating the divergence of the fract vel
       						edge_ij.Sub_D_v(rhs_i,U_i_curr * eps_i,U_j_curr * eps_j);
//       						edge_ij.Sub_D_v(rhs_i,a_i*rho_i,a_j*rho_i);

						//high order stabilizing term
						double temp = 0.0;
#ifdef SPLIT_OSS 
// 						edge_ij.Add_div_v(temp,mTauPressure[i_node]*xi_i,mTauPressure[j_neighbour]*xi_j);
 						edge_ij.Add_div_v(temp,xi_i,xi_j);
#else
       						edge_ij.Add_div_v(temp,pi_i,pi_j);
#endif
						temp *= mBeta[csr_index];
						
   						rhs_i += edge_tau * temp;
// 						rhs_i += temp;

						//assemble laplacian matrix
						mL(i_node, j_neighbour) = sum_l_ikjk;
						l_ii -= sum_l_ikjk;
					}
					mL(i_node, i_node) = l_ii;

					//add density variation contribution
					const double& rho_i_old = mRhoOld[i_node];
					const double& m_i = mr_matrix_container.GetLumpedMass()[i_node];
					rhs_i -= m_i * (rho_i - rho_i_old)/delta_t;
					
					//add mass contribution for compressible flows
/*					double& m_i = mr_matrix_container.GetLumpedMass()[i_node];
					mL(i_node, i_node) += mC2inv[i_node] * m_i / delta_t;*/
				}
				
				//find the max diagonal term
				double max_diag = 0.0;
				for (int i_node = 0; i_node < n_nodes; i_node++)
				{
					double L_diag = mL(i_node, i_node);
					if(fabs(L_diag) > fabs(max_diag)) max_diag = L_diag;
				}
				
				
				

				//respect pressure boundary conditions by penalization
				double huge = max_diag * 1e30;
				for (unsigned int i_pressure = 0; i_pressure < mPressureOutletList.size(); i_pressure++)
				{
					unsigned int i_node = mPressureOutletList[i_pressure];
					mL(i_node, i_node) = huge;
					rhs[i_node] = 0.0;
				}
				
				//modification for level_set
				mr_matrix_container.FillScalarFromDatabase(DISTANCE, mdistances, mr_model_part.Nodes());
				
				//selecting nodes for fixing pressure
// 				std::vector< unsigned int > aux(mdistances.size());
// 				for (unsigned int i_dist = 0; i_dist < mdistances.size(); i_dist++)
// 					aux[i_dist] = 0;
// 				for (unsigned int i_dist = 0; i_dist < mdistances.size(); i_dist++)
// 				{
// 					if(mdistances[i_dist] > 0)
// 					{
// 						aux[i_dist] = 1;
// /*						for (unsigned int csr_index=mr_matrix_container.GetRowStartIndex()[i_dist]; csr_index!=mr_matrix_container.GetRowStartIndex()[i_dist+1]; csr_index++)
// 						{
// 							unsigned int j_neighbour = mr_matrix_container.GetColumnIndex()[csr_index];
// 							aux[j_neighbour] = 1;
// 							
// 						}*/
// 					}
// 					
// 				}
		
				for (unsigned int i_dist = 0; i_dist < mdistances.size(); i_dist++)
				{
// 					if(aux[i_dist] != 0)
					if(mdistances[i_dist] > 0)
					{
// 						mPn1[i_dist] = 0.0;
						mL(i_dist, i_dist) = huge;
						rhs[i_dist] = 0.0;
					}
				}

				//set starting vector for iterative solvers
				for (int i_node = 0; i_node < n_nodes; i_node++)
					dp[i_node] = 0.0;

				//solve linear equation system L dp = rhs
				pLinearSolver->Solve(mL,dp,rhs);
				KRATOS_WATCH(*pLinearSolver)

				//update pressure
				for (int i_node = 0; i_node < n_nodes; i_node++)
					mPn1[i_node] += dp[i_node];	
				for (unsigned int i_pressure = 0; i_pressure < mPressureOutletList.size(); i_pressure++)
				{
					unsigned int i_node = mPressureOutletList[i_pressure];
					mPn1[i_node] = mPressureOutlet[i_pressure];
				}

				//calculate density variation from pressure variation
// 				for (unsigned int i_node = 0; i_node < n_nodes; i_node++)
// 					mRho[i_node] = mRhoOld[i_node] + dp[i_node] * mC2inv[i_node];
// 				for (unsigned int i_density = 0; i_density < mDensityInlet.size(); i_density++)
// 				{
// 					unsigned int i_node = mVelocityInletList[i_density];
// 					mRho[i_node] = mDensityInlet[i_density];
// 				}

				//write pressure and density to Kratos
				mr_matrix_container.WriteScalarToDatabase(PRESSURE, mPn1, rNodes);
// 				mr_matrix_container.WriteScalarToDatabase(DENSITY, mRho, rNodes);

			KRATOS_CATCH("")
			}

			//**********************************************************************************
			//function to solve fluid equations - fractional step 3: correct fractional momentum
			void SolveStep3()
			{
			KRATOS_TRY
				//get number of nodes
				ModelPart::NodesContainerType& rNodes = mr_model_part.Nodes();

				int n_nodes = rNodes.size();
//CANCELLAAAAAAA è necessario??! Non lo sto riempendo con nulla....e ad ogni passo di tempo è nuovo....
				mr_matrix_container.FillVectorFromDatabase(SEEPAGE_DRAG, mDrag, rNodes);
				//CORRECT FRACTIONAL MOMENTUM

				//define work array
				array_1d<double, TDim> correction;
				//read time step size from Kratos
				ProcessInfo& CurrentProcessInfo = mr_model_part.GetProcessInfo();
				double delta_t = CurrentProcessInfo[DELTA_TIME];

				//compute end of step momentum
				#pragma omp parallel for private(correction) firstprivate(delta_t)
				for (int i_node = 0; i_node < n_nodes; i_node++)
				{
					double dist = mdistances[i_node];
					if (dist < 0.0) //node is inside domain ---- if outside do nothing
					{
						array_1d<double, TDim>& U_i_curr = mCurrMom[i_node];
						double delta_p_i = mPn1[i_node] - mPn[i_node];
						const double m_inv = mr_matrix_container.GetInvertedMass()[i_node];
						
						//setting to zero
						for (unsigned int l_comp = 0; l_comp < TDim; l_comp++)
							correction[l_comp] = 0.0;
						
						//compute edge contributions dt*M^(-1)Gp
						for (unsigned int csr_index=mr_matrix_container.GetRowStartIndex()[i_node]; csr_index!=mr_matrix_container.GetRowStartIndex()[i_node+1]; csr_index++)
						{
							unsigned int j_neighbour = mr_matrix_container.GetColumnIndex()[csr_index];
							double delta_p_j = mPn1[j_neighbour] - mPn[j_neighbour];
							
							CSR_Tuple& edge_ij = mr_matrix_container.GetEdgeValues()[csr_index];
											
 							edge_ij.Add_Gp(correction,delta_p_i,delta_p_j);
// 							edge_ij.Sub_Gp(correction,delta_p_i,delta_p_j);
						}
						//compute prefactor
						double coefficient = delta_t * m_inv;
						
						//correct fractional momentum
						for (unsigned int comp = 0; comp < TDim; comp++)
							U_i_curr[comp] += coefficient * correction[comp];
					}
				}

				ApplyVelocityBC(mCurrMom);
				
				CalculateVelocity(mvel_n1,mCurrMom,mRho);

				MultiplyByPorosity(mvel_n1, mvel_n1, mEps);
				


				//write velocity of time step n+1 to Kratos
				mr_matrix_container.WriteVectorToDatabase(VELOCITY, mvel_n1, rNodes);
				
				CalculateDrag(mA, mCurrMom, mDrag, mViscosity);
				//CALCULATE THE DRAG MATRIX TO PASS TO THE SOLID PART
				mr_matrix_container.WriteVectorToDatabase(SEEPAGE_DRAG, mDrag, rNodes);

			KRATOS_CATCH("")
			}

			//************************************
			//function to calculate speed of sound
			void SolveStep4(ModelPart::NodesContainerType& rNodes)
			{
			KRATOS_TRY

				//get number of nodes
				int n_nodes = mC2inv.size();

				//compute speed of sound using equation of state
				#pragma omp parallel for 
				for (int i_node = 0; i_node < n_nodes; i_node++)
				{
					double& rho_i = mRho[i_node];
					double p_i_abs = mPn1[i_node];

					mC2inv[i_node] = rho_i / (mGamma * p_i_abs);
				}

			KRATOS_CATCH("")
			}
			
			//************************************
			void ApplyVelocityBC(CalcVectorType& MomentumArray)
			{
				KRATOS_TRY
			
				//velocity inlet
				int inlet_size = mVelocityInletList.size();
				#pragma omp parallel for schedule(static)
				for (int i_velocity = 0; i_velocity < inlet_size; i_velocity++)
				{
					unsigned int i_node = mVelocityInletList[i_velocity];
					array_1d<double, TDim>& u_i = mVelocityInlet[i_velocity];
					double& rho_i = mDensityInlet[i_velocity];
					array_1d<double, TDim>& U_i = MomentumArray[i_node];
	
					for (unsigned int comp = 0; comp < TDim; comp++)
						U_i[comp] = rho_i * u_i[comp];
				}
				
				//slip condition
				int slip_size = mSlipBoundaryList.size();
				#pragma omp parallel for 
				for (int i_slip = 0; i_slip < slip_size; i_slip++)
				{
					unsigned int i_node = mSlipBoundaryList[i_slip];
					array_1d<double, TDim>& U_i = MomentumArray[i_node];
					array_1d<double, TDim>& an_i = mSlipNormal[i_node];
					double projection_length = 0.0;
					double normalization = 0.0;
					for (unsigned int comp = 0; comp < TDim; comp++)
					{
						projection_length += U_i[comp] * an_i[comp];
						normalization += an_i[comp] * an_i[comp];
					}
					projection_length /= normalization;
						//tangential momentum as difference between original and normal momentum
					for (unsigned int comp = 0; comp < TDim; comp++)
						U_i[comp] -= projection_length * an_i[comp];
				}
				
				//no-slip condition
				int no_slip_size = mNoSlipBoundaryList.size();
				#pragma omp parallel for 
				for (int i_noslip = 0; i_noslip < no_slip_size; i_noslip++)
				{
					unsigned int i_node = mNoSlipBoundaryList[i_noslip];
					array_1d<double, TDim>& U_i = MomentumArray[i_node];
					noalias(U_i) = ZeroVector(TDim);
				}
				
				KRATOS_CATCH("")
			}

			//********************************
			//function to compute coefficients
			void ExtrapolateVelocities(unsigned int extrapolation_layers)
			{
			KRATOS_TRY	

				typedef Node<3> PointType;
				typedef PointerVector<PointType >           PointVector;
				typedef PointVector::iterator 	PointIterator;

				//reset is visited flag
				for( ModelPart::NodesContainerType::iterator inode = mr_model_part.NodesBegin();
								inode != mr_model_part.NodesEnd();
								inode++)	
				{
					inode->GetValue(IS_VISITED) = 0;
				}

			
				//generate a container with the layers to be extrapolated
				std::vector< PointVector > layers(extrapolation_layers);

				//detect the nodes inside the fluid surface
				for( ModelPart::NodesContainerType::iterator inode = mr_model_part.NodesBegin();
								inode != mr_model_part.NodesEnd();
								inode++)	
				{
					if( inode->FastGetSolutionStepValue(DISTANCE) <= 0.0) //candidates are only the ones inside the fluid domain
					{
						WeakPointerVector< Node<3> >& neighb_nodes = inode->GetValue(NEIGHBOUR_NODES); 
						for( WeakPointerVector< Node<3> >::iterator i =	neighb_nodes.begin(); i != neighb_nodes.end(); i++) 
						{ 
							if(i->FastGetSolutionStepValue(DISTANCE) > 0) //add the node as free surface if one of its neighb is outside
							{
								if( inode->GetValue(IS_VISITED) == 0)
								{
									layers[0].push_back( *(inode.base() ) );	
									inode->GetValue(IS_VISITED) = 1;
								}
							}
						} 
					}
				}

// 				//reset is visited flag
// 				for( ModelPart::NodesContainerType::iterator inode = mr_model_part.NodesBegin();
// 								inode != mr_model_part.NodesEnd();
// 								inode++)	
// 				{
// 					inode->GetValue(IS_VISITED) = 0;
// 				}

				//fill the following layers by neighbour relationships
				//each layer fills the following
				for(unsigned int il = 0; il<extrapolation_layers-1; il++)
				{
					for( PointIterator iii=(layers[il]).begin(); iii!=(layers[il]).end(); iii++)
					{
						WeakPointerVector< Node<3> >& neighb_nodes = iii->GetValue(NEIGHBOUR_NODES); 
						for(WeakPointerVector< Node<3> >::iterator jjj=neighb_nodes.begin(); jjj !=neighb_nodes.end(); jjj++) //destination = origin1 + value * Minv*origin
						{ 
		
							if( jjj->FastGetSolutionStepValue(DISTANCE) > 0 &&
							jjj->GetValue(IS_VISITED) == 0.0 )
							{
								layers[il+1].push_back( Node<3>::Pointer( *(jjj.base() ) ) );
								jjj->GetValue(IS_VISITED) = double(il+2.0);
							}
						}
					}
				}

				//perform extrapolation layer by layer by making an average 
				//of the neighbours of lower order
				array_1d<double,3> aux;	
				for(unsigned int il = 1; il<extrapolation_layers; il++)
				{
					for( PointIterator iii=layers[il].begin(); iii!=layers[il].end(); iii++)
					{
// 						noalias(aux) = ZeroVector(3);
// 						double dist_min = 10000000000.0;
// 
// 						array_1d<double,3>& coords_I = iii->Coordinates();
// 						
// 						WeakPointerVector< Node<3> >& neighb_nodes = iii->GetValue(NEIGHBOUR_NODES); 
// 						for(WeakPointerVector< Node<3> >::iterator j=neighb_nodes.begin(); 	j !=neighb_nodes.end(); j++) 
// 						{ 
// 							if(j->GetValue(IS_VISITED) < il+1) //if it is on the next layer
// 							{
// 								array_1d<double,3>& coords_J = j->Coordinates();
// 
// 								double dist = 0.0;
// 								for (unsigned int comp = 0; comp < TDim; comp++)
// 									dist += pow(coords_I[comp]-coords_J[comp],2);
// 								
// 								if(dist < dist_min)
// 								{
// 									dist_min = dist;
// 									noalias( iii->FastGetSolutionStepValue(VELOCITY) ) = j->FastGetSolutionStepValue(VELOCITY);
// 								}
// 								
// 							}
// 						} 


						//extrapolate the average velocity
						noalias(aux) = ZeroVector(3);
						double avg_number = 0.0;
						
						WeakPointerVector< Node<3> >& neighb_nodes = iii->GetValue(NEIGHBOUR_NODES); 
						for(WeakPointerVector< Node<3> >::iterator i=neighb_nodes.begin(); 	i !=neighb_nodes.end(); i++) 
						{ 
							if(i->GetValue(IS_VISITED) < il+1 && i->GetValue(IS_VISITED))
							{
								noalias(aux) += i->FastGetSolutionStepValue(VELOCITY);
								avg_number += 1.0;
							}
						} 
						if(avg_number != 0.0)
						aux /= avg_number;
						noalias( iii->FastGetSolutionStepValue(VELOCITY) ) = aux;
// 						noalias( iii->FastGetSolutionStepValue(VELOCITY,1) ) = aux;
					}
				}

// 				mr_matrix_container.FillVectorFromDatabase(VELOCITY, mvel_n1, mr_model_part.Nodes());
// 				mr_matrix_container.FillScalarFromDatabase(DISTANCE, mdistances, mr_model_part.Nodes());
// 
// 				unsigned int n_nodes = mPn1.size();
// 
// 				//pressure coefficient
// // 				#pragma omp parallel for 
// 				for (int i_node = 0; i_node < n_nodes; i_node++)
// 				{
// 					const double dist_i = mdistances[i_node];
// 					
// 
// 					if( dist_i > 0.0)
// 					{
// 						double nn = 0.0;
// 						
// 						array_1d<double, TDim>& vel_i = mvel_n1[i_node];
// 
// 						for (unsigned int comp = 0; comp < TDim; comp++)
// 							vel_i[comp] = 0.0;						
// 												
// 						//compute edge contributions dt*M^(-1)Gp
// 						for (unsigned int csr_index=mr_matrix_container.GetRowStartIndex()[i_node]; csr_index!=mr_matrix_container.GetRowStartIndex()[i_node+1]; csr_index++)
// 						{
// 							unsigned int j_neighbour = mr_matrix_container.GetColumnIndex()[csr_index];
// 
// 							const double dist_j =  mdistances[j_neighbour];
// 
// 							if(dist_j <= 0.0)
// 							{
// 								const array_1d<double, TDim>& vel_j = mvel_n1[j_neighbour];
// 								
// 								for (unsigned int comp = 0; comp < TDim; comp++)
// 									vel_i[comp] += vel_j[comp];
// 
// 								nn += 1.0;
// 								
// 							}
// 						}
// 
// 						if(nn> 1e-6) //it should be either 0 1 .. N
// 						{
// // std::cout << "inode= " << i_node << "nn = " << nn << std::endl;
// 
// 							double inv_nn = 1.0/nn;
// 							for (unsigned int comp = 0; comp < TDim; comp++)
// 									vel_i[comp] *= inv_nn;
// KRATOS_WATCH(vel_i);
// 						}
// 							
// 					}
// 				}
// 
// 				
// 
// // 				ApplyVelocityBC(mCurrMom);
// 
// 				
// 				//write velocity of time step n+1 to Kratos
// 				mr_matrix_container.WriteVectorToDatabase(VELOCITY, mvel_n1, mr_model_part.Nodes());

				

			KRATOS_CATCH("")
			}

			void ChangeSignToDistance()
			{
				KRATOS_TRY	

				for( ModelPart::NodesContainerType::iterator inode = mr_model_part.NodesBegin();
								inode != mr_model_part.NodesEnd();
								inode++)	
				{
					double dist = inode->FastGetSolutionStepValue(DISTANCE);
					inode->FastGetSolutionStepValue(DISTANCE) = -dist;
				}

				KRATOS_CATCH("")
			}

			void MarkNodesByDistance(double min, double max )
			{
				KRATOS_TRY	

				for( ModelPart::NodesContainerType::iterator inode = mr_model_part.NodesBegin();
								inode != mr_model_part.NodesEnd();
								inode++)	
				{
					double dist = inode->FastGetSolutionStepValue(DISTANCE);
					if(dist > min && dist < max)
						inode->GetValue(IS_VISITED) = 1;
					else
						inode->GetValue(IS_VISITED) = 0;
				}

				KRATOS_CATCH("")
			}

			void SaveScalarVariableToOldStep(Variable<double>& rVar)
			{
				KRATOS_TRY	

				for( ModelPart::NodesContainerType::iterator inode = mr_model_part.NodesBegin();
								inode != mr_model_part.NodesEnd();
								inode++)	
				{
					inode->FastGetSolutionStepValue(rVar,1) = inode->FastGetSolutionStepValue(rVar);
				}

				KRATOS_CATCH("")
			}

			void MarkExternalAndMixedNodes( )
			{
				KRATOS_TRY	

				for( ModelPart::NodesContainerType::iterator inode = mr_model_part.NodesBegin();
								inode != mr_model_part.NodesEnd();
								inode++)	
				{
					inode->GetValue(IS_VISITED) = 0;
				}

				//detect the nodes inside the fluid surface
				for( ModelPart::NodesContainerType::iterator inode = mr_model_part.NodesBegin();
								inode != mr_model_part.NodesEnd();
								inode++)	
				{
					if( inode->FastGetSolutionStepValue(DISTANCE) > 0.0) //candidates are only the ones inside the fluid domain
					{
						inode->GetValue(IS_VISITED) = 1;
						WeakPointerVector< Node<3> >& neighb_nodes = inode->GetValue(NEIGHBOUR_NODES); 
						for( WeakPointerVector< Node<3> >::iterator i =	neighb_nodes.begin(); i != neighb_nodes.end(); i++) 
						{ 
							i->GetValue(IS_VISITED) = 1;
						} 
					}
				}
				KRATOS_CATCH("")
			}

			void MarkInternalAndMixedNodes( )
			{
				KRATOS_TRY	

				for( ModelPart::NodesContainerType::iterator inode = mr_model_part.NodesBegin();
								inode != mr_model_part.NodesEnd();
								inode++)	
				{
					inode->GetValue(IS_VISITED) = 0;
				}

				//detect the nodes inside the fluid surface
				for( ModelPart::NodesContainerType::iterator inode = mr_model_part.NodesBegin();
								inode != mr_model_part.NodesEnd();
								inode++)	
				{
					if( inode->FastGetSolutionStepValue(DISTANCE) <= 0.0) //candidates are only the ones inside the fluid domain
					{
						inode->GetValue(IS_VISITED) = 1;
						WeakPointerVector< Node<3> >& neighb_nodes = inode->GetValue(NEIGHBOUR_NODES); 
						for( WeakPointerVector< Node<3> >::iterator i =	neighb_nodes.begin(); i != neighb_nodes.end(); i++) 
						{ 
							i->GetValue(IS_VISITED) = 1;
						} 
					}
				}
				KRATOS_CATCH("")
			}

			void CalculateVariablesDistribution(double rho_dense, double rho_light, double nu_dense, double nu_light, double eps, const array_1d<double,3>& body_force)
			{
				KRATOS_TRY	

				for( ModelPart::NodesContainerType::iterator inode = mr_model_part.NodesBegin();
								inode != mr_model_part.NodesEnd();
								inode++)	
				{
					double dist = inode->FastGetSolutionStepValue(DISTANCE);

					//calculated smoothed density and viscosity distribution
					double H;
					if(dist < -eps)	H = 0.0;
					else if(dist > eps)	H = 1.0;
					else	H = (dist+eps)/(2.0*eps) + sin(3.141592*dist/eps)/(2.0*3.141592);

					double rho_node = rho_dense + (rho_light-rho_dense)*H;
					inode->FastGetSolutionStepValue(DENSITY) = rho_node;

					double nu_node = nu_dense + (nu_light-nu_dense)*H;
					inode->FastGetSolutionStepValue(VISCOSITY) = nu_node;

					//reset variables outside of the fluid domain
					if( dist < 0 )
						noalias(inode->FastGetSolutionStepValue(BODY_FORCE)) = body_force;
					else
					{
						inode->FastGetSolutionStepValue(PRESSURE) = 0.0;
						noalias(inode->FastGetSolutionStepValue(BODY_FORCE)) = body_force;
						noalias(inode->FastGetSolutionStepValue(VELOCITY)) = ZeroVector(3);
						noalias(inode->FastGetSolutionStepValue(VELOCITY,1)) = ZeroVector(3);
					}
				}

				
				KRATOS_CATCH("")
			}



			//********************************
			//function to compute coefficients
			void CalculateCoefficients(ModelPart::NodesContainerType& rNodes)
			{
			KRATOS_TRY	

				unsigned int n_nodes = mPn1.size();

				//pressure coefficient
				#pragma omp parallel for 
				for ( int i_node = 0; i_node < n_nodes; i_node++)
					mCp[i_node] = (mPn1[i_node] - mPinf) / mQinf;
				mr_matrix_container.WriteScalarToDatabase(PRESSURE_COEFFICIENT, mCp, rNodes);

				//Mach number
				#pragma omp parallel for 
				for ( int i_node = 0; i_node < n_nodes; i_node++)
					mMach[i_node] = norm_2(mvel_n1[i_node]) * sqrt(mC2inv[i_node]);
				mr_matrix_container.WriteScalarToDatabase(MACH_NUMBER, mMach, rNodes);

			KRATOS_CATCH("")
			}

			//**************************************
			//function to calculate the area normals
			void CalculateNormals(ModelPart::ConditionsContainerType& rConditions)
			//void CalculateNormals(ModelPart::NodesContainerType& rNodes, MatrixContainer& matrix_container)
			{
			KRATOS_TRY

				//calculate area normals face-by-face
				array_1d<double,3> area_normal;
				//2D case
				if(TDim == 2)
				{
					for(ModelPart::ConditionsContainerType::iterator cond_it=rConditions.begin(); cond_it!=rConditions.end(); cond_it++)
						CalculateNormal2D(cond_it,area_normal);
				}
				//3D case
				else if(TDim == 3)
				{
					//help vectors for cross product
					array_1d<double,3> v1;
					array_1d<double,3> v2;
					for(ModelPart::ConditionsContainerType::iterator cond_it=rConditions.begin(); cond_it!=rConditions.end(); cond_it++)
						CalculateNormal3D(cond_it,area_normal,v1,v2);
				}

				//(re)initialize normals
				unsigned int n_nodes = mNodalFlag.size();
				mSlipNormal.resize(n_nodes);
				mPressureNormal.resize(n_nodes);
				for (unsigned int i_node = 0; i_node < n_nodes; i_node++) 
				{
					noalias(mSlipNormal[i_node]) = ZeroVector(TDim);
					noalias(mPressureNormal[i_node]) = ZeroVector(TDim);
				}

				//loop over all faces
				for(ModelPart::ConditionsContainerType::iterator cond_it=rConditions.begin(); cond_it!=rConditions.end(); cond_it++)
				{
					//get geometry data of the face
					Geometry<Node<3> >& face_geometry = cond_it->GetGeometry();

					//boolean variables to characterize faces
					bool is_slip_condition = true;
					bool is_pressure_face = true;
					bool is_velocity_inlet = true;
					for (unsigned int if_node = 0; if_node < TDim; if_node++)
					{
						unsigned int i_node = static_cast<unsigned int>(face_geometry[if_node].FastGetSolutionStepValue(AUX_INDEX));
						
						//if the face contains at least 1 node that is not of slip or mixed 
						//then it is not a slip face
						if (	static_cast<unsigned int>(mNodalFlag[i_node]) != 3 &&
							static_cast<unsigned int>(mNodalFlag[i_node]) != 4)
							is_slip_condition = false;
						
						//if the face contains at least one node of pressure it is  a pressure face
						if (	static_cast<unsigned int>(mNodalFlag[i_node]) != 5 && 
							static_cast<unsigned int>(mNodalFlag[i_node]) != 4)
							is_pressure_face = false;
						
						if (static_cast<unsigned int>(mNodalFlag[i_node]) != 1)
							is_velocity_inlet = false;
					}

					//reference for area normal of the face
					array_1d<double,3>& face_normal = cond_it->GetValue(NORMAL);
					
					double node_factor = 1.0/TDim;

					//slip condition
					if (is_slip_condition == true)
						for (unsigned int if_node = 0; if_node < TDim; if_node++)
						{
							unsigned int i_node = static_cast<unsigned int>(face_geometry[if_node].FastGetSolutionStepValue(AUX_INDEX));
							array_1d<double,TDim>& slip_normal = mSlipNormal[i_node];
							for (unsigned int comp = 0; comp < TDim; comp++)
								slip_normal[comp] += node_factor * face_normal[comp];
						}

					//pressure face
					if (is_pressure_face == true || is_velocity_inlet == true)
						for (unsigned int if_node = 0; if_node < TDim; if_node++)
						{
							unsigned int i_node = static_cast<unsigned int>(face_geometry[if_node].FastGetSolutionStepValue(AUX_INDEX));
							array_1d<double,TDim>& pressure_normal = mPressureNormal[i_node];
							for (unsigned int comp = 0; comp < TDim; comp++)
								pressure_normal[comp] += node_factor * face_normal[comp];
						}
						
					//remaining case ... add pressure to pressure nodes and slip to the others
					if(is_pressure_face == false && is_slip_condition == false && is_velocity_inlet == false)
						for (unsigned int if_node = 0; if_node < TDim; if_node++)
						{
							unsigned int i_node = static_cast<unsigned int>(face_geometry[if_node].FastGetSolutionStepValue(AUX_INDEX));
							
							if (	static_cast<unsigned int>(mNodalFlag[i_node]) == 5) //pressure node
							{
								array_1d<double,TDim>& pressure_normal = mPressureNormal[i_node];
								for (unsigned int comp = 0; comp < TDim; comp++)
									pressure_normal[comp] += node_factor * face_normal[comp];
							}
							else if (	static_cast<unsigned int>(mNodalFlag[i_node]) == 3) //slip node
							{
								array_1d<double,TDim>& slip_normal = mPressureNormal[i_node];
								for (unsigned int comp = 0; comp < TDim; comp++)
									slip_normal[comp] += node_factor * face_normal[comp];
							}
						}
				}

			KRATOS_CATCH("")
			}


			
			void SetSpeedOfSound(double c, ModelPart::NodesContainerType& rNodes)
			{
			KRATOS_TRY

				unsigned int n_nodes = mC2inv.size();
				double temp = 1.0 / (c * c);
				for (unsigned int i_node = 0; i_node < n_nodes; i_node++)
					mC2inv[i_node] = temp;
				//WriteScalarToDatabase(LIFT_COEFFICIENT, mC2inv, rNodes);

			KRATOS_CATCH("")
			}


			void SetFreeFlowConditions(array_1d<double, 3> velocity, double pressure, double density, double gamma)
			{
			KRATOS_TRY

				mUinf = velocity;
				mPinf = pressure;
				mRhoinf = density;

				mGamma = gamma;

				mQinf = 0.5 * mRhoinf * norm_2(mUinf) * norm_2(mUinf);
				mMachinf = norm_2(mUinf) / (sqrt(mGamma*mPinf/mRhoinf));

				unsigned int n_nodes = mPn1.size();
				for (unsigned int i_node = 0; i_node < n_nodes; i_node++)
					mC2inv[i_node] = mRho[i_node] / (mGamma * mPn1[i_node]);

				for (unsigned int i_velocity = 0; i_velocity < mVelocityInletList.size(); i_velocity++)
					noalias(mVelocityInlet[i_velocity]) = velocity;

			KRATOS_CATCH("")
			}
			

						
			//**********************************************************************
			void CalculateVelocity(	CalcVectorType& velocity,
							const CalcVectorType& momentum,
    							const ValuesVectorType& rho)
			{
				int loop_size = velocity.size();
				#pragma omp parallel for 
				for (int i_node = 0; i_node < loop_size; i_node++)
				{
					double inv_rho = 1.0/mRho[i_node];
					array_1d<double,TDim>& vel = velocity[i_node];
					const array_1d<double,TDim>& mom = momentum[i_node];
					for (unsigned int comp = 0; comp < TDim; comp++)
						vel[comp] = mom[comp] * inv_rho;
				}
			}


			void SetDissipationLength(double h)
			{
			KRATOS_TRY

				mDissipationLength = h;

			KRATOS_CATCH("")
			}

			void CalculateDrag (CalcVectorType& convective_velocity,
							CalcVectorType& momentum,
							CalcVectorType& drag,
							const ValuesVectorType& viscosity)
			{	
				mViscosity = viscosity;
				int n_nodes = mViscosity.size();
				for ( int i_node = 0; i_node < n_nodes; i_node++)
				{
					double dist = mdistances[i_node];
					if (dist < 0.0) //node is inside domain ---- if outside do nothing
					{
						const array_1d<double, TDim>& a_i = convective_velocity[i_node]; 
						const array_1d<double, TDim>& U_i = momentum[i_node];
						array_1d<double, TDim>& Drag_i = drag[i_node];
						const double& nu_i = viscosity[i_node];
						//porous contribution
						double eps = mEps[i_node];
						double d = mD[i_node]; //diameter of the particle	
						double kinv = 150.0*(1.0-eps)*(1.0-eps)/(eps*eps*eps*d*d);
						double norm_u_2 = 0.0;
						for (unsigned int comp = 0; comp < TDim; comp++)
							norm_u_2 = a_i[comp]*a_i[comp];

						//CORRECTED Term
						double nonlin_term = kinv * nu_i * eps + 1.75 * norm_u_2  * sqrt(kinv / ( eps * 150.0));
						
						for (unsigned int comp = 0; comp < TDim; comp++)
							Drag_i[comp] = nonlin_term * U_i[comp];
					}
				}
			}	
			//*******************************
			//function to free dynamic memory
			void Clear()
			{
			KRATOS_TRY
				mWork.clear();
				mvel_n.clear();
				mvel_n1.clear();
				mA.clear();
				mPn.clear();
				mPn1.clear();
				mHmin.clear();
				mHavg.clear();
				//mAreaNormal.clear();
				//mvel_nitNormal.clear();
				mPressureNormal.clear();
				mSlipNormal.clear();
				mNodalFlag.clear();
				mVelocityInletList.clear();
				mVelocityInlet.clear();
				mPressureOutletList.clear();
				mPressureOutlet.clear();
				mSlipBoundaryList.clear();
				mNoSlipBoundaryList.clear();
				mL.clear();
				mTauPressure.clear();
				mTauConvection.clear();
				mViscosity.clear();
				mEps.clear();
				mEpsOld.clear();

			KRATOS_CATCH("")
			}
			
			//******************************************
			void CalculateForces()
			{
				KRATOS_TRY
						
				//variables for node based data handling
				ModelPart::NodesContainerType& rNodes = mr_model_part.Nodes();
				int n_nodes = rNodes.size();

				//storage of nodal values in local variables
				CalcVectorType rhs;
				rhs.resize(n_nodes);

				//read velocity and pressure data from Kratos
				mr_matrix_container.FillVectorFromDatabase(VELOCITY, mvel_n1, rNodes);
				mr_matrix_container.FillOldVectorFromDatabase(VELOCITY, mvel_n, rNodes);
				
				mr_matrix_container.FillScalarFromDatabase(PRESSURE, mPn1, rNodes);
				mr_matrix_container.FillOldScalarFromDatabase(PRESSURE, mPn, rNodes);
				
				mr_matrix_container.FillScalarFromDatabase(DENSITY, mRho, rNodes);
				mr_matrix_container.FillOldScalarFromDatabase(DENSITY, mRhoOld, rNodes);
				
				mr_matrix_container.FillVectorFromDatabase(BODY_FORCE, mBodyForce, rNodes);
				mr_matrix_container.FillScalarFromDatabase(VISCOSITY, mViscosity, rNodes);

				//read time step size from Kratos
				ProcessInfo& CurrentProcessInfo = mr_model_part.GetProcessInfo();
				double delta_t = CurrentProcessInfo[DELTA_TIME];
				
#pragma omp parallel for 	
				for ( int i_node = 0; i_node < n_nodes; i_node++)
				{
					// -> mCurrMom
					//compute the momentum at the current step -> mCurrMom
					double& rho_i = mRho[i_node];
					const array_1d<double, TDim>& u_i = mvel_n1[i_node];
					array_1d<double, TDim>& U_i = mCurrMom[i_node];
					for (unsigned int comp = 0; comp < TDim; comp++)
						U_i[comp] = rho_i * u_i[comp];
					
					// -> mInitMom
					double& rho_i_old = mRhoOld[i_node];
					//compute the momentum at the beginning of the tep
					const array_1d<double, TDim>& u_i_old = mvel_n[i_node];
					array_1d<double, TDim>& U_i_old = mInitMom[i_node];
					for (unsigned int comp = 0; comp < TDim; comp++)
						U_i_old[comp] = rho_i_old * u_i_old[comp];
					
					//compute volumetric body force
					array_1d<double, TDim>& f_i = mBodyForce[i_node];
					for (unsigned int comp = 0; comp < TDim; comp++)
						f_i[comp] *= rho_i;
				}

				//compute advective velocity - area average of the current velocity
				CalculateAdvectiveVelocity(mvel_n1, mA, msmooth_convective_velocity);

				//compute intrinsic time
				double time_inv = 1.0/delta_t;
				
#pragma omp parallel for firstprivate(time_inv)	
				for (int i_node = 0; i_node < n_nodes; i_node++)
				{
// 					double& h_i = mHavg[i_node];
					double& h_i = mHmin[i_node];
					array_1d<double, TDim>& a_i = mA[i_node];
					const double nu_i = mViscosity[i_node];

					double vel_norm = norm_2(a_i);
					mTauPressure[i_node] = 1.0 / (2.0 * vel_norm/h_i + 0.01*time_inv + nu_i /(h_i*h_i) );
					mTauConvection[i_node] = 1.0 / (2.0 * vel_norm/h_i + 0.01*time_inv + nu_i /(h_i*h_i) );
					
					if (mTauPressure[i_node] < delta_t)
						mTauPressure[i_node] =  delta_t;
					else if(mTauPressure[i_node] > 100.0*delta_t)
						mTauPressure[i_node] = 100.0*delta_t;
				}

				//compute pressure switch
				if (mFirstStep == false)
					if(minclude_shock_capturing == true)
						ComputeMonotonicityPreserving();
				
				mr_matrix_container.SetToZero(rhs);
				CalculateRHS( mCurrMom, mPn1, mA, mBodyForce, mViscosity, rhs);
				
				ValuesVectorType& lumped_mass = mr_matrix_container.GetLumpedMass();
				//add inertia term
				#pragma omp parallel for firstprivate(time_inv)	
				for (int i_node = 0; i_node < n_nodes; i_node++)
				{
					array_1d<double, TDim>& rhs_i = rhs[i_node];
					const array_1d<double, TDim>& curr_mom_i = mCurrMom[i_node];
					const array_1d<double, TDim>& old_mom_i = mInitMom[i_node];
					
					for (unsigned int comp = 0; comp < TDim; comp++)
						rhs_i[comp]-=time_inv*lumped_mass[i_node]*(curr_mom_i[comp]-old_mom_i[comp]);
					
					//change of sign
/*					for (unsigned int comp = 0; comp < TDim; comp++)
						rhs_i[comp] = -rhs_i[comp];*/
					
				}
				
				mr_matrix_container.WriteVectorToDatabase(FORCE, rhs, mr_model_part.Nodes());
						
				KRATOS_CATCH("")
			}

		private:
			MatrixContainer& mr_matrix_container;
			ModelPart& mr_model_part;
			
   			bool msmooth_convective_velocity;
			bool minclude_shock_capturing;
			
			//nodal values
			//velocity vector U at time steps n and n+1
			CalcVectorType mWork, mvel_n, mvel_n1, mInitMom, mCurrMom, mFracMom, mx;
			//pressure vector p at time steps n and n+1
			ValuesVectorType mPn, mPn1, mViscosity;
			//monotony preserving term
			ValuesVectorType mBeta;
			//density
			ValuesVectorType mRho, mRhoOld;
			//compressibility parameter
			ValuesVectorType mC2inv;
			double mGamma;
			double mQinf;
			array_1d<double, TDim> mUinf;
			double mPinf;
			double mRhoinf;
			double mMachinf;
			//coefficients
			ValuesVectorType mCp, mMach, mdistances;
			//advective velocity vector
			CalcVectorType mA;
			//minimum length of the edges surrounding edges surrounding each nodal point
			ValuesVectorType mHmin;
			ValuesVectorType mHavg;
			ValuesVectorType mEps;
			ValuesVectorType mEpsOld;
			ValuesVectorType mD;
			CalcVectorType mEdgeDimensions;
			double mDissipationLength;
			//area normal
			//CalcVectorType mAreaNormal, mvel_nitNormal;
			CalcVectorType mPressureNormal, mSlipNormal;
			//projection terms
			CalcVectorType mPi, mXi;

			CalcVectorType mBodyForce, mDrag;

			//flag for first time step
			bool mFirstStep;

			//flag to differentiate interior and boundary nodes
			ValuesVectorType mNodalFlag;
			//lists of nodes with different types of boundary conditions
			IndicesVectorType mSlipBoundaryList, mNoSlipBoundaryList, mPressureOutletList, mVelocityInletList;
			IndicesVectorType mDissipationList;
			CalcVectorType mVelocityInlet;
			ValuesVectorType mPressureOutlet, mDensityInlet;
			//list for pressure boundary faces
			ModelPart::ConditionsContainerType mPressureFaces;

			//intrinsic time step size
			ValuesVectorType mTauPressure;
			ValuesVectorType mTauConvection;

			//variables for resolving pressure equation
			//laplacian matrix
			TSystemMatrixType mL;


			//***********************************************************
			//functions to calculate area normals for boundary conditions
			void CalculateNormal2D(ModelPart::ConditionsContainerType::iterator cond_it, array_1d<double,3>& area_normal)
			{
					Geometry<Node<3> >& face_geometry = (cond_it)->GetGeometry();
	
					area_normal[0] =    face_geometry[1].Y() - face_geometry[0].Y();
					area_normal[1] = - (face_geometry[1].X() - face_geometry[0].X());
					area_normal[2] =    0.00;
	
					noalias((cond_it)->GetValue(NORMAL)) = area_normal;
			}	
	
			void CalculateNormal3D(ModelPart::ConditionsContainerType::iterator cond_it, array_1d<double,3>& area_normal, array_1d<double,3>& v1,array_1d<double,3>& v2 )
			{
					Geometry<Node<3> >& face_geometry = (cond_it)->GetGeometry();
	
					v1[0] = face_geometry[1].X() - face_geometry[0].X();
					v1[1] = face_geometry[1].Y() - face_geometry[0].Y();
					v1[2] = face_geometry[1].Z() - face_geometry[0].Z();

					v2[0] = face_geometry[2].X() - face_geometry[0].X();
					v2[1] = face_geometry[2].Y() - face_geometry[0].Y();
					v2[2] = face_geometry[2].Z() - face_geometry[0].Z();
	
					MathUtils<double>::CrossProduct(area_normal,v1,v2);
					area_normal *= -0.5;
	
					noalias((cond_it)->GetValue(NORMAL)) = area_normal;
			}


			//******************************************
			//function to calculate advective velocities
			void CalculateAdvectiveVelocity(const CalcVectorType& rVelocity, CalcVectorType& rAdvectiveVelocity, bool smooth_convective_velocity)
			{
			KRATOS_TRY

				if(smooth_convective_velocity == true)
				{	
					//get number of nodes
					int n_nodes = rVelocity.size();
					//initialize advective velocities
/*					#pragma omp parallel for 
					for (int i_node = 0; i_node < n_nodes; i_node++)
						noalias(rAdvectiveVelocity[i_node]) = ZeroVector(TDim);*/
	
					//loop over all nodes
					#pragma omp parallel for 
					for (int i_node = 0; i_node < n_nodes; i_node++)
					{
						//reference for advective velocity of node i
						array_1d<double, TDim>&  a_i = rAdvectiveVelocity[i_node];
						noalias(a_i) = ZeroVector(TDim);

						//setting weighting mass to zero
						double mass_sum = 0.0;
	
						//loop over all neighbours
						for (unsigned int csr_index=mr_matrix_container.GetRowStartIndex()[i_node]; csr_index!=mr_matrix_container.GetRowStartIndex()[i_node+1]; csr_index++)
						{
							//add consistent mass of edge ij to denominator
							double& m_ij = mr_matrix_container.GetEdgeValues()[csr_index].Mass;
							mass_sum += m_ij;
							//reference for velocity of neighbouring node j
							unsigned int j_neighbour = mr_matrix_container.GetColumnIndex()[csr_index];
							const array_1d<double, TDim>& u_j = rVelocity[j_neighbour];
							//add contributions of numerator componentwisely
							for (unsigned int comp = 0; comp < TDim; comp++)
								a_i[comp] += m_ij * u_j[comp];
						}
	
						//for Dirichlet boundary nodes lumped values have to be included
						//attention: nodes with Neumann pressure condition are treated as interior points!
						if ((static_cast<unsigned int>(mNodalFlag[i_node]) != 0) && (static_cast<unsigned int>(mNodalFlag[i_node]) != 5) && (static_cast<unsigned int>(mNodalFlag[i_node]) != 4))
						{
							//taking into account diagonal matrix elements
							double m_ii = mr_matrix_container.GetLumpedMass()[i_node] - mass_sum;
							const array_1d<double, TDim>& u_i = rVelocity[i_node];
							//add contribution to advective velocity
							for (unsigned int comp = 0; comp < TDim; comp++)
								a_i[comp] +=  m_ii * u_i[comp];
							//add contribution to mass sum
							mass_sum += m_ii;
						}
	
						//weighting contributions by the mass sum of all (surrounding) edges
						for (unsigned int comp = 0; comp < TDim; comp++)
							a_i[comp] /= mass_sum;
					}
				}
				else
				{
					//get number of nodes
					int n_nodes = rVelocity.size();
					#pragma omp parallel for 
					for (int i_node = 0; i_node < n_nodes; i_node++)
					{
						 array_1d<double, TDim>& aaa = rAdvectiveVelocity[i_node];
						const array_1d<double, TDim>& u_i = rVelocity[i_node];
						for (unsigned int comp = 0; comp < TDim; comp++)
							aaa[comp] = u_i[comp];
					}
// 	noalias(rAdvectiveVelocity[i_node]) = mvel_n1[i_node];
				}
					
				
// for (unsigned int i_node = 0; i_node < n_nodes; i_node++)
// 	noalias(rAdvectiveVelocity[i_node]) = mvel_n1[i_node];


			KRATOS_CATCH("")
			}

			//*********************************************************
			//function to calculate minimum length of surrounding edges
			void CalculateEdgeLengths(ModelPart::NodesContainerType& rNodes)
			{
			KRATOS_TRY

				//get number of nodes
				unsigned int n_nodes = rNodes.size();
				//reserve memory for storage of nodal coordinates
				std::vector< array_1d<double, 3> > position;
				position.resize(n_nodes);

				//get position of all nodes
				for (typename ModelPart::NodesContainerType::iterator node_it=rNodes.begin(); node_it!=rNodes.end(); node_it++)
				{
					//get the global index of the node
					unsigned int i_node = static_cast<unsigned int>(node_it->FastGetSolutionStepValue(AUX_INDEX));
					//save its coordinates locally
					noalias(position[i_node]) = node_it->Coordinates();

					//initialize minimum edge length with relatively big values
					mHmin[i_node] = 1e10;
				}

				ValuesVectorType& aaa = mr_matrix_container.GetHmin();
				for (unsigned int i_node = 0; i_node < n_nodes; i_node++)
				{
					mHmin[i_node] = aaa[i_node];
				}

				//take unstructured meshes into account
				if(TDim == 2)
				{
					for (unsigned int i_node = 0; i_node < n_nodes; i_node++)
					{
						double& h_i = mHavg[i_node];
						double& m_i = mr_matrix_container.GetLumpedMass()[i_node];
// 						double& rho_i = mRho[i_node];

						h_i = sqrt(2.0*m_i);
					}
				}
				else if(TDim == 3)
				{
					for (unsigned int i_node = 0; i_node < n_nodes; i_node++)
					{
						double& h_i = mHavg[i_node];
						double& m_i = mr_matrix_container.GetLumpedMass()[i_node];
// 						double& rho_i = mRho[i_node];

						h_i = pow (6.0*m_i, 1.0/3.0);
					}
				}

				//compute edge coordinates
				for (unsigned int i_node = 0; i_node < n_nodes; i_node++)
				{
					array_1d<double, 3>& pos_i = position[i_node];

					for (unsigned int csr_index=mr_matrix_container.GetRowStartIndex()[i_node]; csr_index!=mr_matrix_container.GetRowStartIndex()[i_node+1]; csr_index++)
					{
						unsigned int j_neighbour = mr_matrix_container.GetColumnIndex()[csr_index];
						array_1d<double, 3>& pos_j = position[j_neighbour];

						array_1d<double, TDim>& l_k = mEdgeDimensions[csr_index];
						for (unsigned int comp = 0; comp < TDim; comp++)
							l_k[comp] = pos_i[comp] - pos_j[comp];
					}
				}

			KRATOS_CATCH("")
			}


			//*******************************************************
			//function to calculate monotonicity preserving term beta
			void ComputeMonotonicityPreserving()
			{
			KRATOS_TRY

				unsigned int n_nodes = mPn1.size();

				for (unsigned int i_node = 0; i_node < n_nodes; i_node++)
				{
					double& p_i = mPn1[i_node];
					array_1d<double, TDim>& xi_i = mXi[i_node];

					for (unsigned int csr_index=mr_matrix_container.GetRowStartIndex()[i_node]; csr_index!=mr_matrix_container.GetRowStartIndex()[i_node+1]; csr_index++)
					{
						unsigned int j_neighbour = mr_matrix_container.GetColumnIndex()[csr_index];
						double& p_j = mPn1[j_neighbour];
						array_1d<double, TDim>& l_k = mEdgeDimensions[csr_index];
						array_1d<double, TDim>& xi_j = mXi[j_neighbour];

						double press_diff = p_i - p_j;
						double proj_sum = 0.0;
						for (unsigned int comp = 0; comp < TDim; comp++)
							proj_sum += l_k[comp] * (xi_i[comp] + xi_j[comp]);
						proj_sum *= 0.5;

						double temp = fabs(press_diff) + fabs(proj_sum);
						if (temp <= 1e-10)
							mBeta[csr_index] = 1.0;
						else
// 							mBeta[csr_index] = 1.0 - fabs(fabs(press_diff) - fabs(proj_sum)) / temp;
 						mBeta[csr_index] = 1.0 - fabs(press_diff + proj_sum) / temp;
  /*mBeta[csr_index]=1.0;*/
/*						if (mNodalFlag[i_node] == 1.0 || mNodalFlag[i_node] == 4.0 || mNodalFlag[i_node] == 5.0 || mNodalFlag[j_neighbour] == 1.0 || mNodalFlag[j_neighbour] == 4.0 || mNodalFlag[j_neighbour] == 5.0)
							mBeta[csr_index] = 0.0;*/
/*if (mBeta[csr_index]<0.0 && mBeta[csr_index]>1.0)
	KRATOS_WATCH(mBeta[csr_index]);*/
					}
				}

			KRATOS_CATCH("")
			}




			
			inline double CalculateEdgeTau( const double time_inv, const double h_i, 
							const array_1d<double,TDim>& v_i, 
       							const double h_j,
       							const array_1d<double,TDim>& v_j)
			{
				double h_avg = 0.5 * (h_i+h_j);
				
				//calculating norm o
				double norm_avg = 0.0;
				for(unsigned int k=0; k<TDim; k++)
					norm_avg += pow(v_i[k] + v_j[k],2);
				norm_avg *= 0.25;
				norm_avg = sqrt(norm_avg);
				
				return 1.0 / (2.0 * norm_avg/h_avg + time_inv + 1e-6 /(h_avg*h_avg) );
			}

			void DivideByPorosity(CalcVectorType& r_destination,const CalcVectorType& r_origin, const ValuesVectorType& porosity)
			{
					int n_nodes = r_origin.size();
					#pragma omp parallel for 
					for (int i_node = 0; i_node < n_nodes; i_node++)
					{
						array_1d<double, TDim>& dest = r_destination[i_node];
						const array_1d<double, TDim>& orig = r_origin[i_node];
						double factor = 1.0/porosity[i_node];

						for (unsigned int comp = 0; comp < TDim; comp++)
							dest[comp] = factor * orig[comp];
					}
			}
			
			void MultiplyByPorosity(CalcVectorType& r_destination, const CalcVectorType& r_origin, const ValuesVectorType& porosity)
			{
					int n_nodes = r_origin.size();
					#pragma omp parallel for 
					for (int i_node = 0; i_node < n_nodes; i_node++)
					{
						array_1d<double, TDim>& dest = r_destination[i_node];
						const array_1d<double, TDim>& orig = r_origin[i_node];
						double factor = porosity[i_node];

						for (unsigned int comp = 0; comp < TDim; comp++)
							dest[comp] = factor * orig[comp];
					}
			}			
			
			
	};
} //namespace Kratos

#endif //KRATOS_LEVELSET_FLUID_SOLVER_H_INCLUDED defined


