//
//   Project Name:        Kratos
//   Last modified by:    Miguel Masó Sotomayor
//   Date:                June 28th 2017
//   Revision:            1.1
//
//

// Project includes
#include "includes/define.h"
#include "custom_elements/primitive_var_element.hpp"
#include "shallow_water_application.h"
#include "utilities/math_utils.h"
#include "utilities/geometry_utilities.h"

namespace Kratos
{


	//******************************************************************
	//******************************************************************
	template< unsigned int TNumNodes >
	PrimitiveVarElement<TNumNodes>::PrimitiveVarElement(IndexType NewId, GeometryType::Pointer pGeometry)
		: Element(NewId, pGeometry)
	{
	}

	//******************************************************************
	//******************************************************************
	template< unsigned int TNumNodes >
	PrimitiveVarElement<TNumNodes>::PrimitiveVarElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
		: Element(NewId, pGeometry, pProperties)
	{
	}

	template< unsigned int TNumNodes >
	Element::Pointer PrimitiveVarElement<TNumNodes>::Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const
	{
		return Element::Pointer(new PrimitiveVarElement(NewId, GetGeometry().Create(ThisNodes), pProperties));
	}

	template< unsigned int TNumNodes >
	PrimitiveVarElement<TNumNodes>::~PrimitiveVarElement()
	{
	}

	//******************************************************************
	//******************************************************************
	template< unsigned int TNumNodes >
	void PrimitiveVarElement<TNumNodes>::CalculateConsistentMassMatrix(boost::numeric::ublas::bounded_matrix<double,9,9>& rMassMatrix) 
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

	template< unsigned int TNumNodes >
	void PrimitiveVarElement<TNumNodes>::CalculateLumpedMassMatrix(boost::numeric::ublas::bounded_matrix<double,9,9>& rMassMatrix) 
	{
		const unsigned int number_of_nodes = 3;
		rMassMatrix = IdentityMatrix(number_of_nodes*3, number_of_nodes*3);
		rMassMatrix *= 1.0/(3.0);
	}

	//******************************************************************
	//******************************************************************
	template< unsigned int TNumNodes >
	void PrimitiveVarElement<TNumNodes>::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY

		const unsigned int number_of_points = GetGeometry().size();
		if(rLeftHandSideMatrix.size1() != number_of_points*3)
			rLeftHandSideMatrix.resize(number_of_points*3,number_of_points*3,false); // Resizing the system in case it does not have the right size
		if(rRightHandSideVector.size() != number_of_points*3)
			rRightHandSideVector.resize(number_of_points*3,false);

		// Getting the BDF2 coefficients (not fixed to allow variable time step)
		// The coefficients INCLUDE the time step
		const double delta_t = rCurrentProcessInfo[DELTA_TIME];
		//~ const Vector& BDFcoeffs = rCurrentProcessInfo[BDF_COEFFICIENTS];
		double BDFcoeffs[2] = {1.0/delta_t, 1.0/delta_t};

		boost::numeric::ublas::bounded_matrix<double,9,9> msMass   = ZeroMatrix(9,9);     // Mass matrix
		boost::numeric::ublas::bounded_matrix<double,3,2> msDN_DX  = ZeroMatrix(3,2);     // Shape functions gradients
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
		array_1d<double,9> ms_proj_unknown;
		double height;
		//
		bool stabilization = true;           // Stabilization parameter >0.005 (R.Codina, CMAME 197, 2008, 1305-1322) En caso de multiplicar: 0.02
		double height_threshold = 1e-6;
		double Ctau;
		double tau_h;
		boost::numeric::ublas::bounded_matrix<double,2,2> tau_u = ZeroMatrix(2,2);
		//
		bool discontinuity_capturing = true; // Discontinuity capturing parameters
		double gradient_threshold = 1e-6;
		double residual;
		double height_grad_norm;
		double k_dc;


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
			ms_unknown[counter]      = GetGeometry()[iii].FastGetSolutionStepValue(VELOCITY_X);
			ms_proj_unknown[counter] = GetGeometry()[iii].FastGetSolutionStepValue(PROJECTED_VELOCITY_X);
			counter++;

			ms_depth[counter]        = 0;
			ms_unknown[counter]      = GetGeometry()[iii].FastGetSolutionStepValue(VELOCITY_Y);
			ms_proj_unknown[counter] = GetGeometry()[iii].FastGetSolutionStepValue(PROJECTED_VELOCITY_Y);
			counter++;

			ms_depth[counter]        = GetGeometry()[iii].FastGetSolutionStepValue(BATHYMETRY);
			ms_unknown[counter]      = GetGeometry()[iii].FastGetSolutionStepValue(HEIGHT);
			ms_proj_unknown[counter] = GetGeometry()[iii].FastGetSolutionStepValue(PROJECTED_HEIGHT);
			counter++;
		}


		// Compute parameters and derivatives matrices at Gauss points

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

		// N matrix: shape functions
		msN_height(0,2) = msNGauss[0];
		msN_height(0,5) = msNGauss[1];
		msN_height(0,8) = msNGauss[2];

		msN_vel(0,0) = msNGauss[0];
		msN_vel(0,3) = msNGauss[1];
		msN_vel(0,6) = msNGauss[2];
		msN_vel(1,1) = msNGauss[0];
		msN_vel(1,4) = msNGauss[1];
		msN_vel(1,7) = msNGauss[2];
		
		// Previous iteration height at current step
		height = msNGauss[0]*ms_unknown[2] + msNGauss[1]*ms_unknown[5] + msNGauss[2]*ms_unknown[8];


		// Main loop
		// LHS
		// Cross terms
		noalias(rLeftHandSideMatrix)  = height*prod(trans(msN_height),msDN_DX_vel);      // Add <q*h*div(u)> to Mass Eq.
		noalias(msC)                  = gravity*prod(trans(msN_vel),msDN_DX_height);     // Add <w*g*grad(h)> to Momentum Eq.
		noalias(rLeftHandSideMatrix) += msC;

		// Inertia terms
		// Compute the mass matrix
		//~ CalculateConsistentMassMatrix(msMass);
		CalculateLumpedMassMatrix(msMass);
		// LHS += bdf*M
		noalias(rLeftHandSideMatrix) += BDFcoeffs[0] * msMass;

		// Stabilization parameters
		if (stabilization){
			Ctau = 0.01;
			tau_h = Ctau/elem_size*pow(height/gravity,0.5);
			if (height > height_threshold){
				tau_u(0,0) = Ctau/elem_size*pow(gravity/height,0.5);
				tau_u(1,1) = Ctau/elem_size*pow(gravity/height,0.5);
			}
			// Stabilization term
			noalias(rLeftHandSideMatrix) += tau_h * prod(trans(msDN_DX_vel), msDN_DX_vel);                    // Artifficial diffusion to Mass eq.
			noalias(rLeftHandSideMatrix) += prod(trans(msDN_DX_height), Matrix(prod(tau_u,msDN_DX_height)));  // Artifficial diffusion to Momentum eq.
		}
		if (discontinuity_capturing){			
			residual = norm_1(prod(msN_height,ms_unknown)) * norm_1(prod(msDN_DX_vel,ms_unknown)) + BDFcoeffs[1]*norm_1(prod(msN_height, (ms_unknown - ms_proj_unknown)));
			height_grad_norm = norm_2(prod(msDN_DX_height,ms_unknown));
			if (height_grad_norm < gradient_threshold){
				k_dc = 0.0;
			}
			else{
				k_dc = 0.5*0.4*elem_size*residual;//height_grad_norm;  // Residual formulation
			}
			// Add discontinuity capturing to LHS
			noalias(rLeftHandSideMatrix) += k_dc * prod(trans(msDN_DX_vel), msDN_DX_vel); 
			// Print values on Gauss points
			this->SetValue(MIU,height_grad_norm);
			this->SetValue(RESIDUAL_NORM,residual);
			this->SetValue(PR_ART_VISC,k_dc);
		}


		// RHS
		// Source terms
		noalias(rRightHandSideVector)  = -prod(msC, ms_depth);          // Add <w,-g*h*grad(H)> to RHS (Momentum Eq.)

		// Inertia terms
		// RHS += M*vhistory
		noalias(rRightHandSideVector) += BDFcoeffs[1] * prod(msMass, ms_proj_unknown);

		// Substracting the dirichlet term
		// RHS -= LHS*UNKNOWNs
		noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix, ms_unknown);

		rRightHandSideVector *= Area;
		rLeftHandSideMatrix *= Area;

		KRATOS_CATCH("");
	}

	//******************************************************************
	//******************************************************************
	template< unsigned int TNumNodes >
	void PrimitiveVarElement<TNumNodes>::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_THROW_ERROR(std::logic_error,  "method not implemented" , "");
	}

	//******************************************************************
	//******************************************************************
	// This subroutine calculates the nodal contributions for the explicit steps of the
	// Fractional step procedure
	template< unsigned int TNumNodes >
	void PrimitiveVarElement<TNumNodes>::InitializeSolutionStep(ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY

		KRATOS_CATCH("");
	}

	//******************************************************************
	//******************************************************************
	template< unsigned int TNumNodes >
	void PrimitiveVarElement<TNumNodes>::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo)
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

	//******************************************************************
	//******************************************************************
	template< unsigned int TNumNodes >
	void PrimitiveVarElement<TNumNodes>::GetDofList(DofsVectorType& rElementalDofList,ProcessInfo& rCurrentProcessInfo)
	{
		GeometryType& rGeom = GetGeometry();
		const unsigned int number_of_nodes = rGeom.size();
		const unsigned int element_size = 3*number_of_nodes;
		if(rElementalDofList.size() != element_size)
			rElementalDofList.resize(element_size);
		
		int counter=0;
		
		for (unsigned int i=0;i<number_of_nodes;i++){
			rElementalDofList[counter++] = rGeom[i].pGetDof(VELOCITY_X);
			rElementalDofList[counter++] = rGeom[i].pGetDof(VELOCITY_Y);
			rElementalDofList[counter++] = rGeom[i].pGetDof(HEIGHT);
		}
	}

	//******************************************************************
	//******************************************************************
	template< unsigned int TNumNodes >
	void PrimitiveVarElement<TNumNodes>::GetValueOnIntegrationPoints(const Variable<double>& rVariable, std::vector<double>& rValues, const ProcessInfo& rCurrentProcessInfo)
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

template class PrimitiveVarElement<3>;

} // namespace Kratos
