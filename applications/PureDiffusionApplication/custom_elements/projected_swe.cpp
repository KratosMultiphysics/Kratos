//
//   Project Name:        Kratos
//   Last modified by:    $Author: it's me! $
//   Date:                $Date: 2008-08-08 23:58:38 $
//   Revision:            $Revision: 1.0 $
//
//

// Project includes
#include "includes/define.h"
#include "custom_elements/projected_swe.h"
#include "pure_diffusion_application.h"
#include "utilities/math_utils.h"
#include "utilities/geometry_utilities.h"

namespace Kratos
{


	//************************************************************************************
	//************************************************************************************
	ProjectedSWE::ProjectedSWE(IndexType NewId, GeometryType::Pointer pGeometry)
		: Element(NewId, pGeometry)
	{
		// DO NOT ADD DOFS HERE
	}

	//************************************************************************************
	//************************************************************************************
	ProjectedSWE::ProjectedSWE(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
		: Element(NewId, pGeometry, pProperties)
	{
	}

	Element::Pointer ProjectedSWE::Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const
	{
		return Element::Pointer(new ProjectedSWE(NewId, GetGeometry().Create(ThisNodes), pProperties));
	}

	ProjectedSWE::~ProjectedSWE()
	{
	}

	//************************************************************************************
	//************************************************************************************
	void ProjectedSWE::CalculateConsistentMassMatrix(boost::numeric::ublas::bounded_matrix<double,9,9>& rMassMatrix) 
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
		rMassMatrix *= 1.0/(number_of_nodes*3.0*2.0);
	}

	void ProjectedSWE::CalculateLumpedMassMatrix(boost::numeric::ublas::bounded_matrix<double,9,9>& rMassMatrix) 
	{
		const unsigned int number_of_nodes = 3;
		rMassMatrix = IdentityMatrix(number_of_nodes*3, number_of_nodes*3);
		rMassMatrix *= 1.0/(number_of_nodes*3.0);
	}

	//************************************************************************************
	//************************************************************************************
	void ProjectedSWE::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY

		// Getting the BDF2 coefficients (not fixed to allow variable time step)
		// The coefficients INCLUDE the time step
		const double delta_t = rCurrentProcessInfo[DELTA_TIME];
		//~ const int Stationary = rCurrentProcessInfo[STATIONARY];
		//~ const Vector& BDFcoeffs = rCurrentProcessInfo[BDF_COEFFICIENTS];
		//~ array_1d<double,2> BDFcoeffs = {1.0, 1.0};  
		double BDFcoeffs[2] = {1.0/delta_t, 1.0/delta_t};

		boost::numeric::ublas::bounded_matrix<double,9,9> msM     = ZeroMatrix(9,9);     // Mass matrix
		boost::numeric::ublas::bounded_matrix<double,3,2> msDN_DX = ZeroMatrix(3,2);     // Shape functions gradients
		boost::numeric::ublas::bounded_matrix<double,3,3> msA     = ZeroMatrix(3,3);     // Parameters matrix
		boost::numeric::ublas::bounded_matrix<double,3,9> msB     = ZeroMatrix(3,9);     // Gradients and derivatives matrix
		boost::numeric::ublas::bounded_matrix<double,9,9> msC     = ZeroMatrix(9,9);     // Nt*A*B (LHS)
		boost::numeric::ublas::bounded_matrix<double,3,9> msN     = ZeroMatrix(3,9);     // Shape functions matrix
		array_1d<double,3> msNGauss;                                    // Dimension = number of nodes . Position of the gauss point
		array_1d<double,3> ms_prev_height;
		array_1d<double,9> ms_proj_unknown;
		array_1d<double,9> ms_depth;
		array_1d<double,9> ms_unknown;
		int counter, counter1, counter2;
		
		array_1d<double,9> h_valor_arb;

		const unsigned int number_of_points = GetGeometry().size();
		if(rLeftHandSideMatrix.size1() != number_of_points*3)
			rLeftHandSideMatrix.resize(number_of_points*3,number_of_points*3,false); // Resizing the system in case it does not have the right size
		if(rRightHandSideVector.size() != number_of_points*3)
			rRightHandSideVector.resize(number_of_points*3,false);

		// Getting data for the given geometry
		double Area;
		GeometryUtils::CalculateGeometryData(GetGeometry(), msDN_DX, msNGauss, Area); // Asking for gradients and other info

		// Reading properties and conditions
		//~ const double gravity = rCurrentProcessInfo[GRAVITY_Z];
		const double gravity = 9.81;

		// Get previous step and projected values
		counter1 = 0;
		counter2 = 0;
		for(unsigned int iii = 0; iii<number_of_points; iii++){
			ms_depth[counter1++] = 0;
			ms_depth[counter1++] = 0;
			ms_depth[counter1++]   = GetGeometry()[iii].FastGetSolutionStepValue(BATHYMETRY);

			h_valor_arb[iii*3+2]   = GetGeometry()[iii].FastGetSolutionStepValue(HEIGHT,1);

			ms_prev_height[iii] = GetGeometry()[iii].FastGetSolutionStepValue(HEIGHT,1);

			ms_proj_unknown[counter2++] = GetGeometry()[iii].FastGetSolutionStepValue(PROJECTED_VELOCITY_X);
			ms_proj_unknown[counter2++] = GetGeometry()[iii].FastGetSolutionStepValue(PROJECTED_VELOCITY_Y);
			ms_proj_unknown[counter2++] = GetGeometry()[iii].FastGetSolutionStepValue(PROJECTED_HEIGHT);
		}


		// Compute the mass matrix
		//~ CalculateConsistentMassMatrix(msM);
		CalculateLumpedMassMatrix(msM);

		// Compute parameters and derivatives matrices
		// Loop on Gauss points: ONE GAUSS POINT

		// A matrix: gravity and previous height
		msA(0,0) = gravity;
		msA(1,1) = gravity;
		for (unsigned int iii = 0; iii < number_of_points; iii++)  // One Gauss point
			msA(2,2) += ms_prev_height(iii) * msNGauss(iii);
		//msA(2,2) *= Area;

		// B matrix: shape functions derivatives
		msB(0,2) = msDN_DX(0,0);
		msB(0,5) = msDN_DX(1,0);
		msB(0,8) = msDN_DX(2,0);
		msB(1,2) = msDN_DX(0,1);
		msB(1,5) = msDN_DX(1,1);
		msB(1,8) = msDN_DX(2,1);
		msB(2,0) = msDN_DX(0,0);
		msB(2,1) = msDN_DX(0,1);
		msB(2,3) = msDN_DX(1,0);
		msB(2,4) = msDN_DX(1,1);
		msB(2,6) = msDN_DX(2,0);
		msB(2,7) = msDN_DX(2,1);

		// N matrix: shape functions
		for (unsigned int jjj = 0; jjj < number_of_points; jjj++) {
			for (unsigned int iii = 0; iii < number_of_points; iii++)
				msN(iii,iii+3*jjj) = msNGauss[jjj];
		}
		// End loop on Gauss point


		// Main loop
		// LHS
		// const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints();
		msC = prod(trans(msN),Matrix(prod(msA,msB)));  // Nt*A*B
		noalias(rLeftHandSideMatrix) = msC;
		//~ noalias(rLeftHandSideMatrix) = prod(trans(msN),Matrix(prod(msA,msB)));  // Nt*A*B
		//~ KRATOS_WATCH(msN)
		//~ KRATOS_WATCH(msA)
		//~ KRATOS_WATCH(msB)
		//~ KRATOS_WATCH(msC)
		//~ KRATOS_WATCH(h_valor_arb)
		//~ KRATOS_WATCH(prod(msC,h_valor_arb))

		// Inertia terms
		// LHS += bdf*M
		noalias(rLeftHandSideMatrix) += BDFcoeffs[0] * msM;

		// RHS
		// TODO: SOURCE TERM
		noalias(rRightHandSideVector) = prod(msC, ms_depth);

		// Inertia terms
		// RHS += M*vhistory
		noalias(rRightHandSideVector) += BDFcoeffs[1] * prod(msM, ms_proj_unknown);

		// Substracting the dirichlet term
		// RHS -= LHS*UNKNOWNs
		counter = 0;
		for(unsigned int iii = 0; iii<number_of_points; iii++){
			ms_unknown[counter++] = GetGeometry()[iii].FastGetSolutionStepValue(VELOCITY_X);
			ms_unknown[counter++] = GetGeometry()[iii].FastGetSolutionStepValue(VELOCITY_Y);
			ms_unknown[counter++] = GetGeometry()[iii].FastGetSolutionStepValue(HEIGHT);
			//ms_unknown[iii*3+2] = GetGeometry()[iii].FastGetSolutionStepValue(HEIGHT);
		}
		noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix,ms_unknown);

		rRightHandSideVector *= Area;
		rLeftHandSideMatrix *= Area;


		KRATOS_CATCH("");
	}

	//************************************************************************************
	//************************************************************************************
	void ProjectedSWE::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_THROW_ERROR(std::logic_error,  "method not implemented" , "");
	}

	//************************************************************************************
	//************************************************************************************
	// This subroutine calculates the nodal contributions for the explicit steps of the
	// Fractional step procedure
	void ProjectedSWE::InitializeSolutionStep(ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY

		KRATOS_CATCH("");
	}

	//************************************************************************************
	//************************************************************************************
	void ProjectedSWE::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo)
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
	  void ProjectedSWE::GetDofList(DofsVectorType& rElementalDofList,ProcessInfo& rCurrentProcessInfo)
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



} // namespace Kratos
