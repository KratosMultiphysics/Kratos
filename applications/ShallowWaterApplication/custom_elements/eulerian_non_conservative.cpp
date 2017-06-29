//
//   Project Name:        Kratos
//   Last modified by:    $Author:  Miguel Masó Sotomayor $
//   Date:                $Date:             june 14 2017 $
//   Revision:            $Revision:                  1.1 $
//
//

// Project includes
#include "includes/define.h"
#include "custom_elements/non_conservative_dc.h"
#include "shallow_water_application.h"
#include "utilities/math_utils.h"
#include "utilities/geometry_utilities.h"

namespace Kratos
{


	//************************************************************************************
	//************************************************************************************
	EulerianNonConservative::EulerianNonConservative(IndexType NewId, GeometryType::Pointer pGeometry)
		: Element(NewId, pGeometry)
	{
		// DO NOT ADD DOFS HERE
	}

	//************************************************************************************
	//************************************************************************************
	EulerianNonConservative::EulerianNonConservative(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
		: Element(NewId, pGeometry, pProperties)
	{
	}

	Element::Pointer EulerianNonConservative::Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const
	{
		return Element::Pointer(new EulerianNonConservative(NewId, GetGeometry().Create(ThisNodes), pProperties));
	}

	EulerianNonConservative::~EulerianNonConservative()
	{
	}

	//************************************************************************************
	//************************************************************************************
	void EulerianNonConservative::CalculateConsistentMassMatrix(boost::numeric::ublas::bounded_matrix<double,9,9>& rMassMatrix)
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

	void EulerianNonConservative::CalculateLumpedMassMatrix(boost::numeric::ublas::bounded_matrix<double,9,9>& rMassMatrix)
	{
		const unsigned int number_of_nodes = 3;
		rMassMatrix = IdentityMatrix(number_of_nodes*3, number_of_nodes*3);
		rMassMatrix *= 1.0/(3.0);
	}

	//************************************************************************************
	//************************************************************************************
	void EulerianNonConservative::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY

		// Getting the BDF2 coefficients (not fixed to allow variable time step)
		// The coefficients INCLUDE the time step
		const double delta_t = rCurrentProcessInfo[DELTA_TIME];
		//~ const Vector& BDFcoeffs = rCurrentProcessInfo[BDF_COEFFICIENTS];
		//~ array_1d<double,2> BDFcoeffs = {1.0, 1.0};
		double BDFcoeffs[2] = {1.0/delta_t, 1.0/delta_t};

		boost::numeric::ublas::bounded_matrix<double,9,9> msM      = ZeroMatrix(9,9);     // Mass matrix
		boost::numeric::ublas::bounded_matrix<double,3,2> msDN_DX  = ZeroMatrix(3,2);     // Shape functions gradients
		boost::numeric::ublas::bounded_matrix<double,2,2> msG      = ZeroMatrix(2,2);     // Gravity matrix
		boost::numeric::ublas::bounded_matrix<double,3,1> msU      = ZeroMatrix(3,1);     // Iteration matrix: velocity unknown
		boost::numeric::ublas::bounded_matrix<double,9,9> msC      = ZeroMatrix(9,9);     // Nt*A*B (LHS)
		//
		boost::numeric::ublas::bounded_matrix<double,2,9> msN_vel        = ZeroMatrix(2,9);   // Shape functions type matrix (for velocity unknown)
		boost::numeric::ublas::bounded_matrix<double,1,9> msN_height     = ZeroMatrix(1,9);   // Shape functions type matrix (for height unknown)
		boost::numeric::ublas::bounded_matrix<double,2,9> msDN_DX_height = ZeroMatrix(2,9);   // Shape functions gradients matrix (for height unknown)
		boost::numeric::ublas::bounded_matrix<double,1,9> msDN_DX_vel    = ZeroMatrix(1,9);   // Shape functions gradients matrix (for velocity unknown)
		//
		array_1d<double,3> msNGauss;                                    // Dimension = number of nodes . Position of the gauss point
		array_1d<double,9> ms_depth;
		array_1d<double,9> ms_unknown;
		array_1d<double,9> ms_prev_unknown;
		array_1d<double,9> ms_proj_unknown;  // TODO: this variable and dependencies should be removed
		array_1d<double,2> ms_velocity;
		double ms_height;
		double k_dc;                                                        // Discontinuity capturing

		const unsigned int number_of_points = GetGeometry().size();
		if(rLeftHandSideMatrix.size1() != number_of_points*3)
			rLeftHandSideMatrix.resize(number_of_points*3,number_of_points*3,false); // Resizing the system in case it does not have the right size
		if(rRightHandSideVector.size() != number_of_points*3)
			rRightHandSideVector.resize(number_of_points*3,false);

		// Getting data for the given geometry
		double Area;
		GeometryUtils::CalculateGeometryData(GetGeometry(), msDN_DX, msNGauss, Area); // Asking for gradients and other info
		// double elem_size = pow(Area,0.5);

		// Reading properties and conditions
		//~ const double gravity = rCurrentProcessInfo[GRAVITY_Z];
		const double gravity = 9.81;

		// Get current step and projected values
		int counter = 0;
		for(unsigned int iii = 0; iii<number_of_points; iii++){
			ms_depth[counter] = 0;
			ms_unknown[counter] = GetGeometry()[iii].FastGetSolutionStepValue(VELOCITY_X);
			ms_prev_unknown[counter] = GetGeometry()[iii].FastGetSolutionStepValue(VELOCITY_X,1);
			ms_proj_unknown[counter] = GetGeometry()[iii].FastGetSolutionStepValue(PROJECTED_VELOCITY_X);
			counter++;

			ms_depth[counter] = 0;
			ms_unknown[counter] = GetGeometry()[iii].FastGetSolutionStepValue(VELOCITY_Y);
			ms_prev_unknown[counter] = GetGeometry()[iii].FastGetSolutionStepValue(VELOCITY_Y,1);
			ms_proj_unknown[counter] = GetGeometry()[iii].FastGetSolutionStepValue(PROJECTED_VELOCITY_Y);
			counter++;

			ms_depth[counter]   = GetGeometry()[iii].FastGetSolutionStepValue(BATHYMETRY);
			ms_unknown[counter] = GetGeometry()[iii].FastGetSolutionStepValue(HEIGHT);
			ms_prev_unknown[counter] = GetGeometry()[iii].FastGetSolutionStepValue(HEIGHT,1);
			ms_proj_unknown[counter] = GetGeometry()[iii].FastGetSolutionStepValue(PROJECTED_HEIGHT);
			counter++;
		}

		// Compute parameters and derivatives matrices
		// Loop on Gauss points: ONE GAUSS POINT

		// Height gradient
		msDN_DX_height(0,2) = msDN_DX(0,0);
		msDN_DX_height(0,5) = msDN_DX(1,0);
		msDN_DX_height(0,8) = msDN_DX(2,0);
		msDN_DX_height(1,2) = msDN_DX(0,1);
		msDN_DX_height(1,5) = msDN_DX(1,1);
		msDN_DX_height(1,8) = msDN_DX(2,1);
		// Velocity gradient
		msDN_DX_vel(0,0) = msDN_DX(0,0);
		msDN_DX_vel(0,1) = msDN_DX(0,1);
		msDN_DX_vel(0,3) = msDN_DX(1,0);
		msDN_DX_vel(0,4) = msDN_DX(1,1);
		msDN_DX_vel(0,6) = msDN_DX(2,0);
		msDN_DX_vel(0,7) = msDN_DX(2,1);

		// Shape functions for height (Mass)
		msN_height(0,2) = msNGauss[0];
		msN_height(0,5) = msNGauss[1];
		msN_height(0,8) = msNGauss[2];

		// Shape functions for velocity (Momentum)
		msN_vel(0,0) = msNGauss[0];
		msN_vel(0,3) = msNGauss[1];
		msN_vel(0,6) = msNGauss[2];
		msN_vel(1,1) = msNGauss[0];
		msN_vel(1,5) = msNGauss[1];
		msN_vel(1,7) = msNGauss[2];

		// A matrix: gravity and previous iteration height
		msG(0,0) = gravity;
		msG(1,1) = gravity;
		ms_height = norm_1(prod(msN_height,ms_unknown));

		// B matrix: previous iteration velocity at current step
		ms_velocity = prod(msN_vel,ms_unknown);
		msU(0,0) = ms_velocity[0];
		msU(1,0) = ms_velocity[1];
		msU(2,0) = 0;

		// End loop on Gauss point


		// Main loop
		// LHS
		// Cross terms
		noalias(rLeftHandSideMatrix)  = ms_height*prod(trans(msN_height),msDN_DX_vel);          // Mass: q*h*div(u)
		noalias(msC)                  = prod(trans(msN_vel),Matrix(prod(msG,msDN_DX_height)));  // Momentum: w*g*grad(h)
		noalias(rLeftHandSideMatrix) += msC;

		// Convective terms
		noalias(rLeftHandSideMatrix) += prod(trans(msN_vel),Matrix(prod(msU,msDN_DX_vel)));                // Momentum: w*U*div(u)
		noalias(rLeftHandSideMatrix) += prod(trans(msN_height),Matrix(prod(trans(msU),msDN_DX_height)));   // Mass: q*U*grad(h)

		// Inertia terms
		// Compute the mass matrix
		CalculateConsistentMassMatrix(msM);
		//~ CalculateLumpedMassMatrix(msM);
		// LHS += bdf*M
		noalias(rLeftHandSideMatrix) += BDFcoeffs[0] * msM;

		// Add artificial diffusion
		k_dc = 0.4;
		noalias(rLeftHandSideMatrix) += k_dc * prod(trans(msDN_DX_height), msDN_DX_height);

		// RHS
		// Body force
		noalias(rRightHandSideVector) = prod(msC, ms_depth);

		// Inertia terms
		// RHS += M*vhistory
		noalias(rRightHandSideVector) += BDFcoeffs[1] * prod(msM, ms_prev_unknown);

		// Substracting the dirichlet term
		// RHS -= LHS*UNKNOWNs
		noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix, ms_unknown);

		rRightHandSideVector *= Area;
		rLeftHandSideMatrix *= Area;


		KRATOS_CATCH("");
	}

	//************************************************************************************
	//************************************************************************************
	void EulerianNonConservative::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_THROW_ERROR(std::logic_error,  "method not implemented" , "");
	}

	//************************************************************************************
	//************************************************************************************
	// This subroutine calculates the nodal contributions for the explicit steps of the
	// Fractional step procedure
	void EulerianNonConservative::InitializeSolutionStep(ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY

		KRATOS_CATCH("");
	}

	//************************************************************************************
	//************************************************************************************
	void EulerianNonConservative::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo)
	{
		unsigned int number_of_nodes = GetGeometry().PointsNumber();
		if(rResult.size() != number_of_nodes*3)
			rResult.resize(number_of_nodes*3,false);
		int counter=0;

		for (unsigned int i=0;i<number_of_nodes;i++){
			rResult[counter++] = GetGeometry()[i].GetDof(VELOCITY_X).EquationId();
			rResult[counter++] = GetGeometry()[i].GetDof(VELOCITY_Y).EquationId();
			rResult[counter++] = GetGeometry()[i].GetDof(HEIGHT).EquationId();
		}
	}

	//************************************************************************************
	//************************************************************************************
	void EulerianNonConservative::GetDofList(DofsVectorType& rElementalDofList,ProcessInfo& rCurrentProcessInfo)
	{
		unsigned int number_of_nodes = GetGeometry().PointsNumber();
		if(rElementalDofList.size() != number_of_nodes*3)
			rElementalDofList.resize(number_of_nodes*3);

		int counter=0;

		for (unsigned int i=0;i<number_of_nodes;i++){
			rElementalDofList[counter++] = GetGeometry()[i].pGetDof(VELOCITY_X);
			rElementalDofList[counter++] = GetGeometry()[i].pGetDof(VELOCITY_Y);
			rElementalDofList[counter++] = GetGeometry()[i].pGetDof(HEIGHT);
		}
	}

	void EulerianNonConservative::GetValueOnIntegrationPoints(const Variable<double>& rVariable, std::vector<double>& rValues, const ProcessInfo& rCurrentProcessInfo)
	{
		if (rVariable == VEL_ART_VISC || PR_ART_VISC || RESIDUAL_NORM || MIU){
			for (unsigned int PointNumber = 0; PointNumber < 1; PointNumber++)
				rValues[PointNumber] = double(this->GetValue(rVariable));
		}
  }


} // namespace Kratos
