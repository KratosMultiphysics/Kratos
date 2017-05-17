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
		return Element::Pointer(new Poisson2D(NewId, GetGeometry().Create(ThisNodes), pProperties));
	}

	ProjectedSWE::~ProjectedSWE()
	{
	}

	//************************************************************************************
	//************************************************************************************
	void ProjectedSWE::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY
		
		// Getting the BDF2 coefficients (not fixed to allow variable time step)
		// The coefficients INCLUDE the time step
		const int Stationary = rCurrentProcessInfo[STATIONARY];
		const Vector& BDFcoeffs = rCurrentProcessInfo[BDF_COEFFICIENTS];
		const double delta_t = rCurrentProcessInfo[DELTA_TIME];
		
		boost::numeric::ublas::bounded_matrix<double,3,3> msMassFactors = 1.0 / 3.0 * IdentityMatrix(3, 3);
		boost::numeric::ublas::bounded_matrix<double,3,2> msDN_DX;      // Shape functions gradients matrix 
		boost::numeric::ublas::bounded_matrix<double,2,2> msD;          // Conductivity matrix 
		msD = ZeroMatrix(2,2);                                          // Initializing the matrix as zero
  		array_1d<double,3> msN;                                         // Dimension = number of nodes . Position of the gauss point 
		array_1d<double,3> ms_temp;                                     // Dimension = number of nodes . . since we are using a residualbased approach 
		array_1d<double,3> ms_temp_vec_np;
		//~ boost::numeric::ublas::bounded_matrix<double, 2, 2 > Identity = IdentityMatrix(2, 2);
		//~ boost::numeric::ublas::bounded_matrix<double, 2, 2 > First;
		//~ boost::numeric::ublas::bounded_matrix<double, 2, 2 > Second;
		//~ boost::numeric::ublas::bounded_matrix<double, 2, 3 > Third;
		const unsigned int number_of_points = GetGeometry().size();
		if(rLeftHandSideMatrix.size1() != number_of_points)
			rLeftHandSideMatrix.resize(number_of_points,number_of_points,false); // Resizing the system in case it does not have the right size 

		if(rRightHandSideVector.size() != number_of_points)
			rRightHandSideVector.resize(number_of_points,false);

		// Getting data for the given geometry
		double Area;
		GeometryUtils::CalculateGeometryData(GetGeometry(), msDN_DX, msN, Area); // Asking for gradients and other info 

		// Reading properties and conditions
		double permittivity = GetProperties()[CONDUCTIVITY];
		msD(0,0)=permittivity;
		msD(1,1)=permittivity;
		double density = GetProperties()[DENSITY];
		double specific_heat = GetProperties()[SPECIFIC_HEAT];          // see commented lines below:
		//~ const Variable<double>& rDensityVar = my_settings->GetDensityVariable();
		//~ density += GetGeometry()[i].FastGetSolutionStepValue(rDensityVar);
		// End For; density *= lumping_factor;

		// Main loop
		// LHS
		// const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints();
		noalias(rLeftHandSideMatrix) = prod(msDN_DX,Matrix(prod(msD,trans(msDN_DX))));  // Bt*D*B
		
		// Inertia terms
		// LHS += bdf*M
		if(Stationary!=1){
			noalias(rLeftHandSideMatrix) += BDFcoeffs[0] * (density * specific_heat) * msMassFactors;
		}

		// RHS
		// No source terms
		noalias(rRightHandSideVector) = 0 * rRightHandSideVector;

		// Inertia terms
		// RHS += M*vhistory
		if(Stationary!=1){
			for (unsigned int iii = 0; iii < number_of_points; iii++)
				ms_temp_vec_np[iii] = BDFcoeffs[1] * GetGeometry()[iii].FastGetSolutionStepValue(TEMPERATURE, 1);
			for (unsigned int jjj = 2; jjj < BDFcoeffs.size(); jjj++){
				for (unsigned int iii = 0; iii < number_of_points; iii++)
					ms_temp_vec_np[iii] += BDFcoeffs[jjj] * GetGeometry()[iii].FastGetSolutionStepValue(TEMPERATURE, jjj);
			}
			noalias(rRightHandSideVector) += density * specific_heat * prod(msMassFactors, ms_temp_vec_np);
		}

		// Substracting the dirichlet term
		// RHS -= LHS*UNKNOWNs
		for(unsigned int iii = 0; iii<number_of_points; iii++)
			ms_temp[iii] = GetGeometry()[iii].FastGetSolutionStepValue(TEMPERATURE);
		noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix,ms_temp);
		
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
		if(rResult.size() != number_of_nodes)
			rResult.resize(number_of_nodes,false);	

		for (unsigned int i=0;i<number_of_nodes;i++)
				rResult[i] = GetGeometry()[i].GetDof(TEMPERATURE).EquationId();
	}

	//************************************************************************************
	//************************************************************************************
	  void ProjectedSWE::GetDofList(DofsVectorType& rElementalDofList,ProcessInfo& rCurrentProcessInfo)
	{
		unsigned int number_of_nodes = GetGeometry().PointsNumber();
		if(rElementalDofList.size() != number_of_nodes)
			rElementalDofList.resize(number_of_nodes);	

		for (unsigned int i=0;i<number_of_nodes;i++)
			rElementalDofList[i] = GetGeometry()[i].pGetDof(TEMPERATURE);

	}



} // namespace Kratos
