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
//   Date:                $Date: 2011-02-27 13:42:51 $
//   Revision:            $Revision: 1.00 $
//
//

//
// opencl_fluid_solver.cl
//
// OpenCL kernels and functions used in opencl_fluid_solver.h


#include "opencl_edge_data_common.cl"


//
// SolveStep1_1
//
// Part of SolveStep1

__kernel void SolveStep1_1(__global ValueType *mHavg, __global VectorType *mvel_n1, __global VectorType *mTauPressure, __global VectorType *mTauConvection, __global VectorType *mTau2, double mViscosity, double time_inv_avg, double mstabdt_pressure_factor, double mstabdt_convection_factor, double mtau2_factor, const IndexType n_nodes)
{
	// Get work item index
	const size_t i_node = get_global_id(0);

	// Check if we are in the range
	if (i_node < n_nodes)
	{
		double h_avg_i = mHavg[i_node];
		double vel_norm = KRATOS_OCL_LENGTH3(mvel_n1[i_node]);

		mTauPressure[i_node] = 1.00 / (2.00 * vel_norm / h_avg_i + mstabdt_pressure_factor * time_inv_avg + (4.00 * mViscosity) / (h_avg_i * h_avg_i));
		mTauConvection[i_node] = 1.00 / (2.00 * vel_norm / h_avg_i + mstabdt_convection_factor * time_inv_avg + (4.00 * mViscosity) / (h_avg_i * h_avg_i));
		mTau2[i_node] = (mViscosity + h_avg_i * vel_norm * 0.5) * mtau2_factor;
	}
}

//
// SolveStep1_2
//
// Part of SolveStep1

__kernel void SolveStep1_2(__global VectorType *mPi, __global VectorType *mvel_n1, __global IndexType *RowStartIndex, __global IndexType *ColumnIndex, __read_only image2d_t EdgeValues, __global ValueType *InvertedMass, const IndexType n_nodes, __local IndexType *Bounds)
{
	// Get work item index
	const size_t i_node = get_global_id(0);
	const size_t i_thread = get_local_id(0);

	// Reading for loop bounds

	if (i_thread == 0)
	{
		Bounds[0] = RowStartIndex[i_node];
	}

	if (i_node < n_nodes)
	{
		Bounds[i_thread + 1] = RowStartIndex[i_node + 1];
	}

	barrier(CLK_LOCAL_MEM_FENCE);

	// Check if we are in the range
	if (i_node < n_nodes)
	{
		VectorType pi_i = 0.00;

		VectorType a_i = mvel_n1[i_node];

		for (IndexType csr_index = Bounds[i_thread]; csr_index != Bounds[i_thread + 1]; csr_index++)
		{
			IndexType j_neighbour = ColumnIndex[csr_index];
		    VectorType a_j = mvel_n1[j_neighbour];

			EdgeType CurrentEdge = ReadDouble16FromDouble16Image(EdgeValues, csr_index);
			VectorType Ni_DNj = KRATOS_OCL_VECTOR3(KRATOS_OCL_NI_DNJ_0(CurrentEdge), KRATOS_OCL_NI_DNJ_1(CurrentEdge), KRATOS_OCL_NI_DNJ_2(CurrentEdge));

		    Add_ConvectiveContribution(Ni_DNj, &pi_i, a_i, a_i, a_j);  // U_i = a_i, U_j = a_j
		}

		mPi[i_node] = InvertedMass[i_node] * pi_i;
	}
}

	//*********************************************************************
	//function to calculate right-hand side of fractional momentum equation
/*
	void CalculateRHS(
		const CalcVectorType& vel,
		const ValuesVectorType& pressure,
		const CalcVectorType& convective_velocity,
		CalcVectorType& rhs)
	{
	    KRATOS_TRY

	    int n_nodes = vel.size();

	    //calculating the RHS
	    array_1d<double, TDim> stab_low;
	    array_1d<double, TDim> stab_high;
	    const double nu_i = mViscosity;
	    const double nu_j = mViscosity;
	    double inverse_rho = 1.0 / mRho;
	    #pragma omp parallel for private(stab_low,stab_high)
	    for (int i_node = 0; i_node < n_nodes; i_node++) {
		    array_1d<double, TDim>& rhs_i = rhs[i_node];
		    const array_1d<double, TDim>& f_i = mBodyForce;
		    array_1d<double, TDim> a_i = convective_velocity[i_node];

		    const array_1d<double, TDim>& U_i = vel[i_node];
		    const array_1d<double, TDim>& pi_i = mPi[i_node];
		    const double& p_i = pressure[i_node];

		    double edge_tau = mTauConvection[i_node];

		    //initializing with the external forces (e.g. gravity)
		    double& m_i = mr_matrix_container.GetLumpedMass()[i_node];
		    for (unsigned int comp = 0; comp < TDim; comp++)
			rhs_i[comp] = m_i  * f_i[comp] ;

		    //convective term
		    for (unsigned int csr_index = mr_matrix_container.GetRowStartIndex()[i_node]; csr_index != mr_matrix_container.GetRowStartIndex()[i_node + 1]; csr_index++) {
			unsigned int j_neighbour = mr_matrix_container.GetColumnIndex()[csr_index];
			    array_1d<double, TDim> a_j = convective_velocity[j_neighbour];
			    const array_1d<double, TDim>& U_j = vel[j_neighbour];
			    const array_1d<double, TDim>& pi_j = mPi[j_neighbour];
			    const double& p_j = pressure[j_neighbour];

			    CSR_Tuple& edge_ij = mr_matrix_container.GetEdgeValues()[csr_index];

			    edge_ij.Sub_ConvectiveContribution(rhs_i, a_i, U_i, a_j, U_j);
			    edge_ij.Sub_grad_p(rhs_i, p_i*inverse_rho, p_j * inverse_rho);
			    edge_ij.Sub_ViscousContribution(rhs_i, U_i, nu_i, U_j, nu_j);

			    //add stabilization
			    edge_ij.CalculateConvectionStabilization_LOW(stab_low, a_i, U_i, a_j, U_j);
			    edge_ij.CalculateConvectionStabilization_HIGH(stab_high, a_i, pi_i, a_j, pi_j);

//                            double beta = 1.0;
//                             double beta = beta_i;
//                             if(beta_j > beta)
//                                 beta = beta_j;
//                            beta = 1.0;

//                            edge_ij.Sub_StabContribution(rhs_i, edge_tau*beta, 1.0, stab_low, stab_high);
//                             edge_ij.Sub_StabContribution(rhs_i, edge_tau, (1.0-beta), stab_low, stab_high);
			    edge_ij.Sub_StabContribution(rhs_i, edge_tau, 1.0, stab_low, stab_high);


			    //add tau2 term
//                             boost::numeric::ublas::bounded_matrix<double,TDim,TDim>& LL = edge_ij.LaplacianIJ;
//                             for (unsigned int k_comp = 0; k_comp < TDim; k_comp++)
//                             {
//                                 double aaa = 0.0;
//                                 for (unsigned int m_comp = 0; m_comp < TDim; m_comp++)
//                                     aaa +=  LL(k_comp,m_comp) * (U_j[m_comp] - U_i[m_comp]);
//                                 rhs_i[k_comp] -= tau2_i*aaa;
//                             }


		    }



	    }

	    //apply wall resistance
	  if(mWallLawIsActive == true)
		  ComputeWallResistance(vel,rhs);

	    ModelPart::NodesContainerType& rNodes = mr_model_part.Nodes();
	    mr_matrix_container.WriteVectorToDatabase(VELOCITY, mvel_n1, rNodes);
	    KRATOS_CATCH("")
	}

	//*************************************************************************
	//function to solve fluid equations - fractional step 2: calculate pressure

	void SolveStep2(typename TLinearSolver::Pointer pLinearSolver) {
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

	    //read the pressure projection from the database
#endif
	    mr_matrix_container.FillScalarFromDatabase(PRESSURE, mPn1, mr_model_part.Nodes());
	    mr_matrix_container.FillVectorFromDatabase(PRESS_PROJ, mXi, rNodes);
	    mr_matrix_container.FillVectorFromDatabase(VELOCITY, mvel_n1, rNodes);

	    #pragma omp parallel for firstprivate(time_inv)
	    for (int i_node = 0; i_node < n_nodes; i_node++) {

		double& rhs_i = rhs[i_node];
		rhs_i = 0.0;
		const double& p_i = mPn1[i_node];
		const double& p_old_i = mPn[i_node];
		const array_1d<double, TDim>& U_i_curr = mvel_n1[i_node];

		array_1d<double, TDim>& xi_i = mXi[i_node];

		double l_ii = 0.0;

		//loop over all neighbours
		for (unsigned int csr_index = mr_matrix_container.GetRowStartIndex()[i_node]; csr_index != mr_matrix_container.GetRowStartIndex()[i_node + 1]; csr_index++) {
		    unsigned int j_neighbour = mr_matrix_container.GetColumnIndex()[csr_index];
		    const double& p_j = mPn1[j_neighbour];
		    const double& p_old_j = mPn[j_neighbour];
		    const array_1d<double, TDim>& U_j_curr = mvel_n1[j_neighbour];
		    const array_1d<double, TDim>& xi_j = mXi[j_neighbour];

		    CSR_Tuple& edge_ij = mr_matrix_container.GetEdgeValues()[csr_index];

#ifdef SYMM_PRESS
		    double edge_tau = 0.5 * (mTauPressure[i_node] + mTauPressure[j_neighbour]);
#else
		    double edge_tau = mTauPressure[i_node];
#endif



		    //compute laplacian operator
		    double sum_l_ikjk;
		    edge_ij.CalculateScalarLaplacian(sum_l_ikjk);
		    double sum_l_ikjk_onlydt = sum_l_ikjk * (delta_t);
		    sum_l_ikjk *= (delta_t + edge_tau);

		    //assemble right-hand side
		    //pressure contribution
		    rhs_i -= sum_l_ikjk * (p_j - p_i);
		    rhs_i += sum_l_ikjk_onlydt * (p_old_j - p_old_i);


		    //calculating the divergence of the fract vel
		    edge_ij.Sub_D_v(rhs_i, U_i_curr*mRho, U_j_curr * mRho);

		    //high order stabilizing term
		    double temp = 0.0;
		    edge_ij.Add_div_v(temp, xi_i, xi_j);
		    rhs_i += edge_tau * temp;

		    //assemble laplacian matrix
		    mL(i_node, j_neighbour) = sum_l_ikjk;
		    l_ii -= sum_l_ikjk;
		}

		mL(i_node, i_node) = l_ii;
	    }

	    if(muse_mass_correction == true)
	    {
		#pragma omp parallel for
		for (int i_node = 0; i_node < n_nodes; i_node++)
		{
		    double& rhs_i = rhs[i_node];
		    rhs_i -= mdiv_error[i_node];
		}
	    }

	    //find the max diagonal term
	    double max_diag = 0.0;
	    for (int i_node = 0; i_node < n_nodes; i_node++) {
		double L_diag = mL(i_node, i_node);
		if (fabs(L_diag) > fabs(max_diag)) max_diag = L_diag;
	    }




	    //respect pressure boundary conditions by penalization
//            double huge = max_diag * 1e6;
//            for (unsigned int i_pressure = 0; i_pressure < mPressureOutletList.size(); i_pressure++) {
//                unsigned int i_node = mPressureOutletList[i_pressure];
//                mL(i_node, i_node) = huge;
//                rhs[i_node] = 0.0;
//            }
	    for (unsigned int i_pressure = 0; i_pressure < mPressureOutletList.size(); i_pressure++) {
		unsigned int i_node = mPressureOutletList[i_pressure];
		mL(i_node, i_node) = max_diag;
		rhs[i_node] = 0.0;
		for (unsigned int csr_index = mr_matrix_container.GetRowStartIndex()[i_node]; csr_index != mr_matrix_container.GetRowStartIndex()[i_node + 1]; csr_index++)
		{
		    unsigned int j_neighbour = mr_matrix_container.GetColumnIndex()[csr_index];
		    mL(i_node, j_neighbour) = 0.0;
		}
	    }

	    //set starting vector for iterative solvers
	    for (int i_node = 0; i_node < n_nodes; i_node++)
		dp[i_node] = 0.0;

	    pLinearSolver->Solve(mL, dp, rhs);
	    KRATOS_WATCH(*pLinearSolver)


	    //update pressure
	    for (int i_node = 0; i_node < n_nodes; i_node++)
		mPn1[i_node] += dp[i_node];

	    //write pressure and density to Kratos
	    mr_matrix_container.WriteScalarToDatabase(PRESSURE, mPn1, rNodes);


	    //compute pressure proj for the next step

	    #pragma omp parallel for firstprivate(time_inv), private(work_array)
	    for (int i_node = 0; i_node < n_nodes; i_node++) {
		array_1d<double, TDim>& xi_i = mXi[i_node];
		for (unsigned int comp = 0; comp < TDim; comp++)
		    xi_i[comp] = 0.0;


		    const double& p_i = mPn1[i_node];

		    for (unsigned int csr_index = mr_matrix_container.GetRowStartIndex()[i_node]; csr_index != mr_matrix_container.GetRowStartIndex()[i_node + 1]; csr_index++) {
			//get global index of neighbouring node j
			unsigned int j_neighbour = mr_matrix_container.GetColumnIndex()[csr_index];

			    const double& p_j = mPn1[j_neighbour];

			    //projection of pressure gradients
			    CSR_Tuple& edge_ij = mr_matrix_container.GetEdgeValues()[csr_index];

			    edge_ij.Add_grad_p(xi_i, p_i, p_j);
		    }

		    const double& m_inv = mr_matrix_container.GetInvertedMass()[i_node];
		    for (unsigned int l_comp = 0; l_comp < TDim; l_comp++)
			xi_i[l_comp] *= m_inv;


	    }

	    mr_matrix_container.WriteVectorToDatabase(PRESS_PROJ, mXi, rNodes);

	    KRATOS_CATCH("")
	}

	//**********************************************************************************
	//function to solve fluid equations - fractional step 3: correct fractional momentum

	void SolveStep3() {
	    KRATOS_TRY
	    //get number of nodes
	    ModelPart::NodesContainerType& rNodes = mr_model_part.Nodes();

	    int n_nodes = rNodes.size();

	    //define work array
	    array_1d<double, TDim> correction;
	    //read time step size from Kratos
	    ProcessInfo& CurrentProcessInfo = mr_model_part.GetProcessInfo();
	    double delta_t = CurrentProcessInfo[DELTA_TIME];

	    double factor = 0.5;
	    if(massume_constant_dp == true)
		factor = 1.0;

	    //compute end of step momentum
	    double rho_inv = 1.0 / mRho;
	    #pragma omp parallel for private(correction) firstprivate(delta_t,rho_inv,factor)
	    for (int i_node = 0; i_node < n_nodes; i_node++) {

		    array_1d<double, TDim>& U_i_curr = mvel_n1[i_node];
		    double delta_p_i = (mPn1[i_node] - mPn[i_node]) * rho_inv*factor;
		    const double m_inv = mr_matrix_container.GetInvertedMass()[i_node];

		    //setting to zero
		    for (unsigned int l_comp = 0; l_comp < TDim; l_comp++)
			correction[l_comp] = 0.0;

		    //compute edge contributions dt*M^(-1)Gp
		    for (unsigned int csr_index = mr_matrix_container.GetRowStartIndex()[i_node]; csr_index != mr_matrix_container.GetRowStartIndex()[i_node + 1]; csr_index++) {
			unsigned int j_neighbour = mr_matrix_container.GetColumnIndex()[csr_index];
			double delta_p_j = (mPn1[j_neighbour] - mPn[j_neighbour]) * rho_inv*factor;

			CSR_Tuple& edge_ij = mr_matrix_container.GetEdgeValues()[csr_index];

			edge_ij.Sub_grad_p(correction, delta_p_i, delta_p_j);
		    }
		    //compute prefactor
		    double coefficient = delta_t * m_inv;

		    //correct fractional momentum
		    for (unsigned int comp = 0; comp < TDim; comp++)
			U_i_curr[comp] += coefficient * correction[comp];

	    }

	    ApplyVelocityBC(mvel_n1);

	    //write velocity of time step n+1 to Kratos
	    mr_matrix_container.WriteVectorToDatabase(VELOCITY, mvel_n1, rNodes);



	    //calculate the error on the divergence
	    if(muse_mass_correction == true)
	    {
		#pragma omp parallel for private(correction) firstprivate(delta_t,rho_inv)
		for (int i_node = 0; i_node < n_nodes; i_node++)
		{
		    double& div_i_err = mdiv_error[i_node];
		    div_i_err = 0.0;

			const array_1d<double, TDim>& U_i_curr = mvel_n1[i_node];

			//compute edge contributions dt*M^(-1)Gp
			for (unsigned int csr_index = mr_matrix_container.GetRowStartIndex()[i_node]; csr_index != mr_matrix_container.GetRowStartIndex()[i_node + 1]; csr_index++)
			{
			    unsigned int j_neighbour = mr_matrix_container.GetColumnIndex()[csr_index];
			    array_1d<double, TDim>& U_j_curr = mvel_n1[j_neighbour];

			    CSR_Tuple& edge_ij = mr_matrix_container.GetEdgeValues()[csr_index];

			    edge_ij.Add_D_v(div_i_err, U_i_curr*mRho, U_j_curr * mRho);
			}

		}
	    }

	    KRATOS_CATCH("")
	}


	//************************************

	void ApplyVelocityBC(CalcVectorType& VelArray) {
	    KRATOS_TRY


	    if(mWallLawIsActive == false)
	    {
		//apply conditions on corner edges
		int edge_size = medge_nodes_direction.size();
		#pragma omp parallel for firstprivate(edge_size)
		for (int i = 0; i < edge_size; i++)
		{
		    int i_node = medge_nodes[i];
		    const array_1d<double, TDim>& direction = medge_nodes_direction[i];
			array_1d<double, TDim>& U_i = VelArray[i_node];
			double temp=0.0;
			for (unsigned int comp = 0; comp < TDim; comp++)
			    temp += U_i[comp] * direction[comp];

			for (unsigned int comp = 0; comp < TDim; comp++)
			    U_i[comp] = direction[comp]*temp;

		}

		//apply conditions on corners
		int corner_size = mcorner_nodes.size();
		for (int i = 0; i < corner_size; i++)
		{
		    int i_node = mcorner_nodes[i];

		    array_1d<double, TDim>& U_i = VelArray[i_node];
			for (unsigned int comp = 0; comp < TDim; comp++)
			    U_i[comp] = 0.0;
		}
	    }


	    //slip condition
	    int slip_size = mSlipBoundaryList.size();
	    #pragma omp parallel for firstprivate(slip_size)
	    for (int i_slip = 0; i_slip < slip_size; i_slip++)
	    {
		unsigned int i_node = mSlipBoundaryList[i_slip];
		    array_1d<double, TDim>& U_i = VelArray[i_node];
		    array_1d<double, TDim>& an_i = mSlipNormal[i_node];
		    double projection_length = 0.0;
		    double normalization = 0.0;
		    for (unsigned int comp = 0; comp < TDim; comp++) {
			projection_length += U_i[comp] * an_i[comp];
			normalization += an_i[comp] * an_i[comp];
		    }
		    projection_length /= normalization;
		    //tangential momentum as difference between original and normal momentum
		    for (unsigned int comp = 0; comp < TDim; comp++)
			U_i[comp] -= projection_length * an_i[comp];

	    }

	    //fixed condition
	    int fixed_size = mFixedVelocities.size();
	    #pragma omp parallel for firstprivate(fixed_size)
	    for (int i_velocity = 0; i_velocity < fixed_size; i_velocity++)
	    {
		unsigned int i_node = mFixedVelocities[i_velocity];

		    const array_1d<double, TDim>& u_i_fix = mFixedVelocitiesValues[i_velocity];
		    array_1d<double, TDim>& u_i = VelArray[i_node];

		    for (unsigned int comp = 0; comp < TDim; comp++)
			u_i[comp] = u_i_fix[comp];

	    }



	    KRATOS_CATCH("")
	}

void ComputeWallResistance(
		const CalcVectorType& vel,
		CalcVectorType& rhs
		)
	{
	    //parameters:
	    double k = 0.41;
	    double B = 5.1;
	    double density = mRho;
	    double mu = mViscosity;
	    double toll = 1e-6;
	    double ym = mY_wall; //0.0825877; //0.0093823
	    double y_plus_incercept = 10.9931899;
	    unsigned int itmax = 100;

	    if (mu == 0)
		KRATOS_ERROR(std::logic_error, "it is not possible to use the wall law with 0 viscosity", "");

	    //slip condition
	    int slip_size = mSlipBoundaryList.size();
#pragma omp parallel for firstprivate(slip_size,B,density,mu,toll,ym,y_plus_incercept,itmax)
	    for (int i_slip = 0; i_slip < slip_size; i_slip++)
	    {
		unsigned int i_node = mSlipBoundaryList[i_slip];

		    array_1d<double, TDim>& rhs_i = rhs[i_node];
		    const array_1d<double, TDim>& U_i = vel[i_node];
		    const array_1d<double, TDim>& an_i = mSlipNormal[i_node];

		    //compute the modulus of the velocity
		    double mod_vel = 0.0;
		    double area = 0.0;
		    for (unsigned int comp = 0; comp < TDim; comp++)
		    {
			mod_vel += U_i[comp] * U_i[comp];
			area += an_i[comp] * an_i[comp];
		    }
		    mod_vel = sqrt(mod_vel);
		    area = sqrt(area);

		    //now compute the skin friction
		    double mod_uthaw = sqrt(mod_vel * mu / ym);
		    const double y_plus = ym * mod_uthaw / mu;

		    if (y_plus > y_plus_incercept)
		    {
			//begin cicle to calculate the real u_thaw's module:
			unsigned int it = 0;
			double dx = 1e10;
			//                        KRATOS_WATCH(fabs(dx));
			while (fabs(dx) > toll * mod_uthaw && it < itmax)
			{
			    double a = 1.0 / k;
			    double temp = a * log(ym * mod_uthaw / mu) + B;
			    double y = mod_uthaw * (temp) - mod_vel;
			    double y1 = temp + a;
			    dx = y / y1;
			    mod_uthaw -= dx;
			    it = it + 1;
			}

			//                         KRATOS_WATCH(toll*mod_uthaw);
			//                         KRATOS_WATCH(area);
			//                        KRATOS_WATCH(it);
			if (it == itmax)
			    std::cout << "attention max number of iterations exceeded in wall law computation" << std::endl;


		    }
		    //                    else
		    //                    {
		    //                        for (unsigned int comp = 0; comp < TDim; comp++)
		    //                            rhs_i[comp] -= U_i[comp] * area * mu  / (density*ym) ;
		    //                    }

		    if (mod_vel > 1e-12)
			for (unsigned int comp = 0; comp < TDim; comp++)
			    rhs_i[comp] -= U_i[comp] * area * mod_uthaw * mod_uthaw * density / (mod_vel);




	    }
	}
*/
