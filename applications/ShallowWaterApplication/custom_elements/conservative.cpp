//
//   Project Name:        Kratos
//   Last modified by:    $Author:  Miguel Masó Sotomayor $
//   Date:                $Date:              july 3 2017 $
//   Revision:            $Revision:                  1.1 $
//
//

// Project includes
#include "includes/define.h"
#include "custom_elements/conservative.h"
#include "shallow_water_application.h"
#include "utilities/math_utils.h"
#include "utilities/geometry_utilities.h"

namespace Kratos
{


	//************************************************************************************
	//************************************************************************************
	Conservative::Conservative(IndexType NewId, GeometryType::Pointer pGeometry)
		: Element(NewId, pGeometry)
	{
		// DO NOT ADD DOFS HERE
	}

	//************************************************************************************
	//************************************************************************************
	Conservative::Conservative(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
		: Element(NewId, pGeometry, pProperties)
	{
	}

	Element::Pointer Conservative::Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const
	{
		return Element::Pointer(new Conservative(NewId, GetGeometry().Create(ThisNodes), pProperties));
	}

	Conservative::~Conservative()
	{
	}

	//************************************************************************************
	//************************************************************************************
	void Conservative::CalculateConsistentMassMatrix(boost::numeric::ublas::bounded_matrix<double,9,9>& rMassMatrix)
	{
		const unsigned int number_of_nodes = 3;
		//~ const unsigned int number_of_dof = 3;
		rMassMatrix = IdentityMatrix(number_of_nodes*3, number_of_nodes*3);
		for(unsigned int i = 0; i<number_of_nodes; i++){
			for(unsigned int j = 0; j<number_of_nodes; j++){
				for(unsigned int k = 0; k<3; k++)
					rMassMatrix(k+3*i, k+3*j) += 1.0;
			}
		}
		rMassMatrix *= 1.0/(12.0);
	}

	void Conservative::CalculateLumpedMassMatrix(boost::numeric::ublas::bounded_matrix<double,9,9>& rMassMatrix)
	{
		const unsigned int number_of_nodes = 3;
		rMassMatrix = IdentityMatrix(number_of_nodes*3, number_of_nodes*3);
		rMassMatrix *= 1.0/(3.0);
	}

	//************************************************************************************
	//************************************************************************************
	void Conservative::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY

		// Getting the BDF2 coefficients (not fixed to allow variable time step)
		// The coefficients INCLUDE the time step
		const double delta_t = rCurrentProcessInfo[DELTA_TIME];
		//~ const Vector& BDFcoeffs = rCurrentProcessInfo[BDF_COEFFICIENTS];
		//~ array_1d<double,2> BDFcoeffs = {1.0, 1.0};
		double BDFcoeffs[2] = {1.0/delta_t, 1.0/delta_t};

		boost::numeric::ublas::bounded_matrix<double,9,9> msMass   = ZeroMatrix(9,9);     // Mass matrix
		boost::numeric::ublas::bounded_matrix<double,3,2> msDN_DX  = ZeroMatrix(3,2);     // Shape functions gradients
		boost::numeric::ublas::bounded_matrix<double,3,1> ms_hU    = ZeroMatrix(3,1);     // Iteration matrix: momentum unknown
		boost::numeric::ublas::bounded_matrix<double,9,9> msC      = ZeroMatrix(9,9);     // Nt*A*B (LHS)
		boost::numeric::ublas::bounded_matrix<double,3,9> msN      = ZeroMatrix(3,9);     // Shape functions type
		//
		boost::numeric::ublas::bounded_matrix<double,2,9> msN_mom        = ZeroMatrix(2,9);   // Shape functions type matrix (for momentum unknown)
		boost::numeric::ublas::bounded_matrix<double,1,9> msN_height     = ZeroMatrix(1,9);   // Shape functions type matrix (for height unknown)
		boost::numeric::ublas::bounded_matrix<double,2,9> msDN_DX_height = ZeroMatrix(2,9);   // Shape functions for scalar gradients (for height unknown)
		boost::numeric::ublas::bounded_matrix<double,1,9> msDN_DX_mom    = ZeroMatrix(1,9);   // Shape functions for vectior divergence (for momentum unknown)
		boost::numeric::ublas::bounded_matrix<double,2,9> msGrad_mom     = ZeroMatrix(2,9);   // Shape functions for vector gradient (for momentum unknown)
		//
		array_1d<double,3> msNGauss;                                    // Dimension = number of nodes . Position of the gauss point
		array_1d<double,9> ms_depth;
		array_1d<double,9> ms_rain;
		array_1d<double,9> ms_unknown;
		array_1d<double,9> ms_proj_unknown;
		array_1d<double,9> ms_inv_unknown;
		array_1d<double,2> momentum;
		array_1d<double,2> inv_h_grad;
		//~ boost::numeric::ublas::bounded_matrix<double,2,1> inv_h_grad;
		double height;
		// double inv_height;
		double divU;
		//
		double Ctau;        // Stabilization parameter >0.005 (R.Codina, CMAME 197, 2008, 1305-1322)
		double depth;
		double tau_h;
		boost::numeric::ublas::bounded_matrix<double,2,2> tau_u = ZeroMatrix(2,2);
		double gradient_threshold;    // Shock capturing parameters
		double m_residual;
		double m_height_grad_norm;
		double k_dc;

		const unsigned int number_of_points = GetGeometry().size();
		if(rLeftHandSideMatrix.size1() != number_of_points*3)
			rLeftHandSideMatrix.resize(number_of_points*3,number_of_points*3,false); // Resizing the system in case it does not have the right size
		if(rRightHandSideVector.size() != number_of_points*3)
			rRightHandSideVector.resize(number_of_points*3,false);

		// Getting data for the given geometry
		double Area;
		GeometryUtils::CalculateGeometryData(GetGeometry(), msDN_DX, msNGauss, Area); // Asking for gradients and other info
		double elem_size = pow(Area,0.5);

		// Reading properties and conditions
		//~ const double gravity = rCurrentProcessInfo[GRAVITY_Z];
		const double gravity = 9.81;

		// Get current step and projected values
		int counter = 0;
		for(unsigned int iii = 0; iii<number_of_points; iii++){
			ms_depth[counter]        = 0;
			ms_rain[counter]         = 0;
			ms_unknown[counter]      = GetGeometry()[iii].FastGetSolutionStepValue(MOMENTUM_X);
			ms_proj_unknown[counter] = GetGeometry()[iii].FastGetSolutionStepValue(PROJECTED_MOMENTUM_X);
			counter++;

			ms_depth[counter]        = 0;
			ms_rain[counter]         = 0;
			ms_unknown[counter]      = GetGeometry()[iii].FastGetSolutionStepValue(MOMENTUM_Y);
			ms_proj_unknown[counter] = GetGeometry()[iii].FastGetSolutionStepValue(PROJECTED_MOMENTUM_Y);
			counter++;

			ms_depth[counter]        = GetGeometry()[iii].FastGetSolutionStepValue(BATHYMETRY);
			ms_rain[counter]         = GetGeometry()[iii].FastGetSolutionStepValue(RAIN);
			ms_unknown[counter]      = GetGeometry()[iii].FastGetSolutionStepValue(HEIGHT);
			ms_proj_unknown[counter] = GetGeometry()[iii].FastGetSolutionStepValue(PROJECTED_HEIGHT);
			counter++;
		}

		// Compute parameters and derivatives matrices on Gauss points

		// Height gradient
		msDN_DX_height(0,2) = msDN_DX(0,0);
		msDN_DX_height(0,5) = msDN_DX(1,0);
		msDN_DX_height(0,8) = msDN_DX(2,0);
		msDN_DX_height(1,2) = msDN_DX(0,1);
		msDN_DX_height(1,5) = msDN_DX(1,1);
		msDN_DX_height(1,8) = msDN_DX(2,1);
		// Momentum divergence
		msDN_DX_mom(0,0) = msDN_DX(0,0);
		msDN_DX_mom(0,1) = msDN_DX(0,1);
		msDN_DX_mom(0,3) = msDN_DX(1,0);
		msDN_DX_mom(0,4) = msDN_DX(1,1);
		msDN_DX_mom(0,6) = msDN_DX(2,0);
		msDN_DX_mom(0,7) = msDN_DX(2,1);
		// Momentum gradient
		msGrad_mom(0,0) = msDN_DX(0,0);
		msGrad_mom(0,1) = msDN_DX(0,0);
		msGrad_mom(1,0) = msDN_DX(0,1);
		msGrad_mom(1,1) = msDN_DX(0,1);
		msGrad_mom(0,3) = msDN_DX(1,0);
		msGrad_mom(0,4) = msDN_DX(1,0);
		msGrad_mom(1,3) = msDN_DX(1,1);
		msGrad_mom(1,4) = msDN_DX(1,1);
		msGrad_mom(0,6) = msDN_DX(2,0);
		msGrad_mom(0,7) = msDN_DX(2,0);
		msGrad_mom(1,6) = msDN_DX(2,1);
		msGrad_mom(1,7) = msDN_DX(2,1);

		// N matrix: shape functions
		for (unsigned int jjj = 0; jjj < number_of_points; jjj++) {
			for (unsigned int iii = 0; iii < number_of_points; iii++)
				msN(iii,iii+3*jjj) = msNGauss[jjj];
		}
		//noalias(msNmass) = row(msN, 2);
		msN_height(0,2) = msNGauss[0];
		msN_height(0,5) = msNGauss[1];
		msN_height(0,8) = msNGauss[2];
		//noalias(msNmoment) = row(msN, 0,1);
		msN_mom(0,0) = msNGauss[0];
		msN_mom(0,3) = msNGauss[1];
		msN_mom(0,6) = msNGauss[2];
		msN_mom(1,1) = msNGauss[0];
		msN_mom(1,5) = msNGauss[1];
		msN_mom(1,7) = msNGauss[2];


		// Previous height iteration at current time step
		height = norm_1(prod(msN_height,ms_unknown));
		//~ height = inner_prod(msN_height,ms_unknown);
		// Previous momentum iteration at current time step
		momentum = prod(msN_mom,ms_unknown);
		// Previous inv height gradient at current time step
		inv_h_grad[0] = msDN_DX(0,0)/ms_unknown[2] + msDN_DX(1,0)/ms_unknown[5] + msDN_DX(2,0)/ms_unknown[8];
		inv_h_grad[1] = msDN_DX(0,1)/ms_unknown[2] + msDN_DX(1,1)/ms_unknown[5] + msDN_DX(2,1)/ms_unknown[8];
		// Previous div(U) iteration
		divU  = msDN_DX(0,0)*ms_unknown[0]/ms_unknown[2];
		divU += msDN_DX(0,1)*ms_unknown[1]/ms_unknown[2];
		divU += msDN_DX(1,0)*ms_unknown[3]/ms_unknown[5];
		divU += msDN_DX(1,1)*ms_unknown[4]/ms_unknown[5];
		divU += msDN_DX(2,0)*ms_unknown[6]/ms_unknown[8];
		divU += msDN_DX(2,1)*ms_unknown[7]/ms_unknown[8];

		// U matrix
		ms_hU(0,0) = momentum[0];
		ms_hU(1,0) = momentum[1];
		ms_hU(2,0) = 0;


		// Main loop
		// LHS
		// Cross terms
		noalias(rLeftHandSideMatrix)  = ZeroMatrix(9,9);
		noalias(msC)                  = gravity*height*prod(trans(msN_mom),msDN_DX_height);  // Add <w,g*h*grad(h)> to Momentum Eq.
		noalias(rLeftHandSideMatrix) += msC;


		// Inertia terms
		//~ CalculateConsistentMassMatrix(msMass);
		CalculateLumpedMassMatrix(msMass);
		noalias(rLeftHandSideMatrix) += BDFcoeffs[0] * msMass;          // LHS += bdf*M


		noalias(rLeftHandSideMatrix) += divU * msMass;                  // Add <q,div(u)*h> to Mass Eq. and <w,div(u)*hu> to Momentum Eq.


		// Stabilization parameters
		Ctau = 0.01;
		depth = norm_1(prod(msN_height,ms_depth));
		tau_h = Ctau/elem_size*pow(depth/gravity,0.5);
		tau_u(0,0) = Ctau/elem_size*pow(gravity/depth,0.5);
		tau_u(1,1) = Ctau/elem_size*pow(gravity/depth,0.5);
		// Stabilization term
		noalias(rLeftHandSideMatrix) += tau_h * prod(trans(msDN_DX_mom), msDN_DX_mom);                    // Add artifficial diffusion to Mass eq.
		noalias(rLeftHandSideMatrix) += prod(trans(msDN_DX_height), Matrix(prod(tau_u,msDN_DX_height)));  // Add artifficial diffusion to Momentum eq.

		// Add discontinuity capturing term via adding artificial diffusion to height
		gradient_threshold = 1e-6;
		m_residual = norm_1(prod(msN_height,ms_unknown)) * norm_1(prod(msDN_DX_mom,ms_unknown)) + BDFcoeffs[1]*norm_1(prod(msN_height, (ms_unknown - ms_proj_unknown)));
		m_height_grad_norm = norm_2(prod(msDN_DX_height,ms_unknown));
		if (m_height_grad_norm < gradient_threshold){
			k_dc = 0;
		}
		else{
			k_dc = 0.5*0.4*elem_size*m_residual;//m_height_grad_norm;  // Residual formulation
		}
		// Print values on Gauss points
		this->SetValue(MIU,m_height_grad_norm);
		this->SetValue(RESIDUAL_NORM,m_residual);
		this->SetValue(PR_ART_VISC,k_dc);
		// Add discontinuity capturing to LHS
		noalias(rLeftHandSideMatrix) += k_dc * prod(trans(msDN_DX_height), msDN_DX_height);


		// RHS
		// Source terms
		noalias(rRightHandSideVector)  = -prod(msC, ms_depth);          // Add <w,-g*h*grad(H)> to RHS (Momentum Eq.)
		noalias(rRightHandSideVector) +=  prod(msMass, ms_rain);        // Add <q,rain>         to RHS (Mass Eq.)

		// Inertia terms
		// RHS += M*vhistory
		noalias(rRightHandSideVector) += BDFcoeffs[1] * prod(msMass, ms_proj_unknown);

		// Subtracting the dirichlet term
		// RHS -= LHS*UNKNOWNs
		noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix, ms_unknown);



		//~ KRATOS_WATCH(ms_unknown)
		//~ KRATOS_WATCH(rRightHandSideVector)



		rRightHandSideVector *= Area;
		rLeftHandSideMatrix *= Area;


		KRATOS_CATCH("");
	}

	//************************************************************************************
	//************************************************************************************
	void Conservative::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_THROW_ERROR(std::logic_error,  "method not implemented" , "");
	}

	//************************************************************************************
	//************************************************************************************
	// This subroutine calculates the nodal contributions for the explicit steps of the
	// Fractional step procedure
	void Conservative::InitializeSolutionStep(ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY

		KRATOS_CATCH("");
	}

	//************************************************************************************
	//************************************************************************************
	void Conservative::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo)
	{
		unsigned int number_of_nodes = GetGeometry().PointsNumber();
		if(rResult.size() != number_of_nodes*3)
			rResult.resize(number_of_nodes*3,false);
		int counter=0;

		for (unsigned int i=0;i<number_of_nodes;i++){
			rResult[counter++] = GetGeometry()[i].GetDof(MOMENTUM_X).EquationId();
			rResult[counter++] = GetGeometry()[i].GetDof(MOMENTUM_Y).EquationId();
			rResult[counter++] = GetGeometry()[i].GetDof(HEIGHT).EquationId();
		}
	}

	//************************************************************************************
	//************************************************************************************
	void Conservative::GetDofList(DofsVectorType& rElementalDofList,ProcessInfo& rCurrentProcessInfo)
	{
		unsigned int number_of_nodes = GetGeometry().PointsNumber();
		if(rElementalDofList.size() != number_of_nodes*3)
			rElementalDofList.resize(number_of_nodes*3);

		int counter=0;

		for (unsigned int i=0;i<number_of_nodes;i++){
			rElementalDofList[counter++] = GetGeometry()[i].pGetDof(MOMENTUM_X);
			rElementalDofList[counter++] = GetGeometry()[i].pGetDof(MOMENTUM_Y);
			rElementalDofList[counter++] = GetGeometry()[i].pGetDof(HEIGHT);
		}
	}

	void Conservative::GetValueOnIntegrationPoints(const Variable<double>& rVariable, std::vector<double>& rValues, const ProcessInfo& rCurrentProcessInfo)
	{
		if (rVariable == VEL_ART_VISC){
			for (unsigned int PointNumber = 0; PointNumber < 1; PointNumber++)
				rValues[PointNumber] = double(this->GetValue(VEL_ART_VISC));
		}
		if (rVariable == PR_ART_VISC){
			for (unsigned int PointNumber = 0; PointNumber < 1; PointNumber++)
				rValues[PointNumber] = double(this->GetValue(PR_ART_VISC));
		}
		if (rVariable == RESIDUAL_NORM){
			for (unsigned int PointNumber = 0; PointNumber < 1; PointNumber++)
				rValues[PointNumber] = double(this->GetValue(RESIDUAL_NORM));
		}
		if (rVariable == MIU){
			for (unsigned int PointNumber = 0; PointNumber < 1; PointNumber++)
				rValues[PointNumber] = double(this->GetValue(MIU));
		}
    }



} // namespace Kratos
