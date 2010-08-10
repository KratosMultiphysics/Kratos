/*
==============================================================================
KratosStructuralApplication 
A library based on:
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi, Janosch Stascheit, Felix Nagel 
pooyan@cimne.upc.edu 
rrossi@cimne.upc.edu
janosch.stascheit@rub.de
nagel@sd.rub.de
- CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain
- Ruhr-University Bochum, Institute for Structural Mechanics, Germany


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
 
/* **************************************************************************************
*          
*   Last Modified by:    $Author: rrossi $
*   Date:                $Date: 2008-10-13 07:00:53 $
*   Revision:            $Revision: 1.12 $
*
* ***************************************************************************************/


// System includes


// External includes


// Project includes
#include "includes/define.h"
#include "custom_elements/membrane_element.h"
#include "includes/constitutive_law.h"
#include "structural_application.h"


namespace Kratos
{
	namespace MembraneAuxiliaries
	{
	    Matrix  msB(0,0);
	    #pragma omp threadprivate(msB)
	    
	    boost::numeric::ublas::bounded_matrix<double,3,3>  msQ = ZeroMatrix(3,3);
	    #pragma omp threadprivate(msQ)
	    
	    Matrix	msD = ZeroMatrix(3,3);
	    #pragma omp threadprivate(msD)
	    
	    Vector	msStrainVector = ZeroVector(3);
	    #pragma omp threadprivate(msStrainVector)
	    
	    Vector	msStressVector = ZeroVector(3);
	    #pragma omp threadprivate(msStressVector)
	    
	    boost::numeric::ublas::bounded_matrix<double,2,2>  msC = ZeroMatrix(2,2);
	    #pragma omp threadprivate(msC)
	    
	    Matrix	msDN_DX(0,0);
	    #pragma omp threadprivate(msDN_DX)
	}
	using namespace MembraneAuxiliaries;

	//***********************************************************************************
	//***********************************************************************************
		  // -------- //
		 //  PUBLIC  //
		// -------- //

	// Constructor
	MembraneElement::MembraneElement()
	{
	}

	// Constructor
	MembraneElement::MembraneElement(IndexType NewId, GeometryType::Pointer pGeometry)
		: Element(NewId, pGeometry)
	{
	}

	// Constructor
	MembraneElement::MembraneElement(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
		: Element(NewId, pGeometry, pProperties)

	{
		//initializing static variables
		unsigned int number_of_nodes = GetGeometry().size();
		unsigned int dim = number_of_nodes*3;
		msB.resize(3,dim);
	        msDN_DX.resize(number_of_nodes,3);
	}

	//***********************************************************************************
	//***********************************************************************************

	Element::Pointer MembraneElement::Create(
		IndexType NewId,
		NodesArrayType const& ThisNodes,
		PropertiesType::Pointer pProperties) const

	{
		return Element::Pointer(new MembraneElement(NewId, GetGeometry().Create(ThisNodes), pProperties));
	}

	//***********************************************************************************
	//***********************************************************************************
	// Destructor
	MembraneElement::~MembraneElement()
	{
	}

	//***********************************************************************************
	//***********************************************************************************

	void MembraneElement::EquationIdVector(
		EquationIdVectorType& rResult,
		ProcessInfo& rCurrentProcessInfo)

	{
		KRATOS_TRY
		unsigned int number_of_nodes = GetGeometry().size();
		unsigned int dim = number_of_nodes*3;
		if(rResult.size() != dim)
			rResult.resize(dim);

		for (unsigned int i=0;i<number_of_nodes;i++)
		{
			int index = i*3;
			rResult[index]   = GetGeometry()[i].GetDof(DISPLACEMENT_X).EquationId();
			rResult[index+1] = GetGeometry()[i].GetDof(DISPLACEMENT_Y).EquationId();
			rResult[index+2] = GetGeometry()[i].GetDof(DISPLACEMENT_Z).EquationId();
		}
		KRATOS_CATCH("")
	}

	//***********************************************************************************
	//***********************************************************************************

	void MembraneElement::GetDofList(
		DofsVectorType& ElementalDofList,
		ProcessInfo& rCurrentProcessInfo)
	
	{
		ElementalDofList.resize(0);

		for (unsigned int i=0;i<GetGeometry().size();i++)
		{
			ElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_X));
			ElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_Y));
			ElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_Z));
		}
	}

	//***********************************************************************************
	//***********************************************************************************

	void MembraneElement::Initialize()

	{
		KRATOS_TRY

		//getting all "Actual" info from the geometry
		const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints();
	
		//resizing jacobian inverses containers
		mDetJ0.resize(integration_points.size());
		
		GeometryType::JacobiansType J0;
		J0 = GetGeometry().Jacobian(J0);  
		mTotalDomainInitialSize = 0.00;

		mStrainsVector.resize(integration_points.size());
		mStressesVector.resize(integration_points.size());

		mV1.resize(integration_points.size());
		mV2.resize(integration_points.size());
		mG_Vector.resize(integration_points.size(),ZeroMatrix(2,2));
		
		array_1d<double,3> g0e;
		array_1d<double,3> g0n;
		array_1d<double,3> V1;
		array_1d<double,3> V2;
		array_1d<double,3> V3;
		array_1d<double,3> N;
		
		// Initialize Variables
		mdensity = GetProperties()[DENSITY];
		mThickness0 = GetProperties()[THICKNESS];
		mThickness.resize(integration_points.size(),0.00);

		// Initialize Material
		if(mConstitutiveLawVector.size() == 0)
		{
			mConstitutiveLawVector.resize(integration_points.size());

			for (unsigned int i=0; i<mConstitutiveLawVector.size(); i++)
			{
/*				ConstitutiveLaw<Node<3> >::Pointer material = ConstitutiveLaw<Node<3> >::Pointer( new Isotropic2D() );
				mConstitutiveLawVector[i] = material;*/
				
				
				mConstitutiveLawVector[i] = GetProperties()[CONSTITUTIVE_LAW]->Clone();
				
				mConstitutiveLawVector[i]->InitializeMaterial( 	GetProperties(), GetGeometry(),	row(GetGeometry().ShapeFunctionsValues(), i) ); 
			}
		}

		//calculating the inverse J0
		for(unsigned int PointNumber = 0; PointNumber<integration_points.size(); PointNumber++)
		{
			//getting informations for integration
			double IntegrationWeight = integration_points[PointNumber].Weight();

			//calculating and storing inverse of the jacobian and the parameters needed
			g0e[0]=J0[PointNumber](0,0); g0n[0]=J0[PointNumber](0,1);
			g0e[1]=J0[PointNumber](1,0); g0n[1]=J0[PointNumber](1,1);
			g0e[2]=J0[PointNumber](2,0); g0n[2]=J0[PointNumber](2,1);

			//calculate base vectors
			CrossProduct(V3,g0e,g0n);
			N = V3;
			N /= norm_2(V3);
			V1 = g0e; 
			V1 /= norm_2(V1);
			CrossProduct(V2,N,V1);

			//saving the initial local base vectors
			mV1[PointNumber] = V1;
			mV2[PointNumber] = V2;

			// Calculation of Matrix G (a sort of inverse jacobian)
			double J11 = norm_2(g0e);
			double J12 = inner_prod(g0e,g0n)/norm_2(g0e);
			double J22 = norm_2(V3) / norm_2(g0e);

			Matrix G(2,2);
			G(0,0) = 1/J11;    G(0,1) = -J12/(J11*J22);
			G(1,0) = 0.00;     G(1,1) = 1/J22;

			//saving the G matrix for the point number
			noalias(mG_Vector[PointNumber]) = G;

			//Calculate the reduced mass matrix
			mDetJ0[PointNumber] = norm_2(V3);

			//calculating the total area
			mTotalDomainInitialSize += mDetJ0[PointNumber]*IntegrationWeight;
		}

		KRATOS_CATCH("")
	}

	//***********************************************************************************
	//***********************************************************************************
	
	void MembraneElement::CalculateRightHandSide(
		VectorType& rRightHandSideVector,
		ProcessInfo& rCurrentProcessInfo)

	{
		//calculation flags
		bool CalculateStiffnessMatrixFlag = false;
		bool CalculateResidualVectorFlag = true;
		MatrixType temp = Matrix();
		
		CalculateAll(temp, rRightHandSideVector, rCurrentProcessInfo, CalculateStiffnessMatrixFlag, CalculateResidualVectorFlag);
//std::cout << "RHS della membrana" << std::endl;
//KRATOS_WATCH(Id() )
//KRATOS_WATCH(rRightHandSideVector)
//std::cout << "*******************" << std::endl;
	}

	//***********************************************************************************
	//***********************************************************************************
	
	void MembraneElement::CalculateLocalSystem(
		MatrixType& rLeftHandSideMatrix,
		VectorType& rRightHandSideVector,
		ProcessInfo& rCurrentProcessInfo)

	{
		//calculation flags
		bool CalculateStiffnessMatrixFlag = true;
		bool CalculateResidualVectorFlag = true;
		
		CalculateAll(rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo, CalculateStiffnessMatrixFlag, CalculateResidualVectorFlag);
	}

	//***********************************************************************************
	//***********************************************************************************

	void MembraneElement::CalculateOnIntegrationPoints(
		const Variable<Matrix>& rVariable,
		std::vector<Matrix>& Output,
		const ProcessInfo& rCurrentProcessInfo)

	{
		
		//reading integration points and local gradients
		const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints();
		//const Matrix& Ncontainer = GetGeometry().ShapeFunctionsValues();

		//calculating actual jacobian
		GeometryType::JacobiansType J;
		J = GetGeometry().Jacobian(J);

		//auxiliary terms
		boost::numeric::ublas::bounded_matrix<double,2,2> j;
		boost::numeric::ublas::bounded_matrix<double,2,2> g;
		boost::numeric::ublas::bounded_matrix<double,2,2> C;
		array_1d<double,3> ge;
		array_1d<double,3> gn;
		array_1d<double,3> v3;


		if(Output.size() != integration_points.size())
			Output.resize(integration_points.size());

		for(unsigned int PointNumber=0;PointNumber<integration_points.size();PointNumber++)
		{
			//double IntegrationWeight = GetGeometry().IntegrationPoints()[PointNumber].Weight();

			ge[0]=J[PointNumber](0,0); gn[0]=J[PointNumber](0,1);
			ge[1]=J[PointNumber](1,0); gn[1]=J[PointNumber](1,1);
			ge[2]=J[PointNumber](2,0); gn[2]=J[PointNumber](2,1);

			CrossProduct(v3,ge,gn);
			CalculateJ(j,ge,gn,v3);

			// Calculation of matrix g = jtrans*j;
			noalias(g) = prod(trans(j),j);

			// calculation of the Right Cauchy-Green Tensor C = Gtrans*g*G
			boost::numeric::ublas::bounded_matrix<double,2,2> tmp;
			tmp = prod(g,mG_Vector[PointNumber]);
			noalias(C) = prod(trans(mG_Vector[PointNumber]),tmp);

			// Calculation of the StrainVector
			CalculateStrain(msStrainVector, C);
            
            Matrix dummy = ZeroMatrix(0,0);
//            KRATOS_WATCH(msStrainVector);
            mConstitutiveLawVector[PointNumber]->CalculateMaterialResponse(
                    msStrainVector,
                    dummy,
                    msStressVector,
                    dummy,
                    rCurrentProcessInfo,
                    GetProperties(),
                    GetGeometry(),
                    row(GetGeometry().ShapeFunctionsValues(),PointNumber),
                    true,
                    false,
                    false );
		    
//		    KRATOS_WATCH(msStrainVector);
//		    KRATOS_WATCH(msStressVector);

			// calculation of the StressVector (PK2)
// 			mConstitutiveLawVector[PointNumber]->UpdateMaterial( msStrainVector,
// 						GetProperties(),
// 						GetGeometry(),
// 						row(GetGeometry().ShapeFunctionsValues(),PointNumber),
// 						rCurrentProcessInfo );
			
// 			mConstitutiveLawVector[PointNumber]->CalculateStress(msStrainVector,msStressVector); 

			noalias(mStressesVector[PointNumber]) = ZeroVector(6);
			Calculate_GlobalStressVector(mStressesVector[PointNumber], msStressVector, mV1[PointNumber], mV2[PointNumber]);	//saving the stress vector


			if(rVariable==GREEN_LAGRANGE_STRAIN_TENSOR)
			{
				if(Output[PointNumber].size2() != msStrainVector.size())
					Output[PointNumber].resize(1,msStrainVector.size());
				for(unsigned int ii = 0; ii<msStrainVector.size(); ii++)
					Output[PointNumber](0,ii) = msStrainVector[ii];
			}
			else if(rVariable==PK2_STRESS_TENSOR)
			{
				if(Output[PointNumber].size2() != 6)
					Output[PointNumber].resize(1,6);
				
//			KRATOS_WATCH(mStressesVector[PointNumber]);
/*
				mConstitutiveLawVector[PointNumber]->UpdateMaterial( msStrainVector,
							GetProperties(),
							GetGeometry(),
							row(GetGeometry().ShapeFunctionsValues(),PointNumber),
							rCurrentProcessInfo );
				mConstitutiveLawVector[PointNumber]->CalculateStress(msStrainVector,msStressVector); 

				Calculate_GlobalStressVector(mStressesVector[PointNumber], msStressVector, mV1[PointNumber], mV2[PointNumber]);	//saving the stress vector*/

				
				for(unsigned int ii = 0; ii<6; ii++)
					Output[PointNumber](0,ii) = mStressesVector[PointNumber][ii];
			}
		}

	}

	//***********************************************************************************
	//***********************************************************************************

	void MembraneElement::MassMatrix(
		MatrixType& rMassMatrix,
		ProcessInfo& rCurrentProcessInfo)

	{
		KRATOS_TRY

		//rMassMatrix.resize(0,0);
		// LUMPED MASS MATRIX
		unsigned int number_of_nodes = GetGeometry().size();
		unsigned int MatSize = number_of_nodes * 3;
		if(rMassMatrix.size1() != MatSize)
			rMassMatrix.resize(MatSize,MatSize);
		rMassMatrix = ZeroMatrix(MatSize,MatSize);

		double TotalMass = mTotalDomainInitialSize * GetProperties()[THICKNESS] * GetProperties()[DENSITY];
		Vector LumpFact;
		LumpFact = GetGeometry().LumpingFactors(LumpFact);
		
		for(unsigned int i=0; i<number_of_nodes; i++)
		{
			double temp = LumpFact[i]*TotalMass;
			for(unsigned int j=0; j<3; j++)
			{
				unsigned int index = i*3 + j;
				rMassMatrix(index,index) = temp;
			}
		}

		KRATOS_CATCH("")
	}   

	//***********************************************************************************
	//***********************************************************************************

	void MembraneElement::DampMatrix(
		MatrixType& rDampMatrix,
		ProcessInfo& rCurrentProcessInfo)

	{
		KRATOS_TRY

		// LUMPED DAMPING MATRIX
		unsigned int number_of_nodes = GetGeometry().size();
		unsigned int MatSize = number_of_nodes * 3;
		if(rDampMatrix.size1() != MatSize)
			rDampMatrix.resize(MatSize,MatSize);
		rDampMatrix = ZeroMatrix(MatSize,MatSize);

		double TotalMass = mTotalDomainInitialSize * GetProperties()[THICKNESS] * GetProperties()[DENSITY];
		Vector LumpFact;
		LumpFact = GetGeometry().LumpingFactors(LumpFact);

		for(unsigned int i=0; i<number_of_nodes; i++)
		{
			double temp = LumpFact[i]*TotalMass;
			for(unsigned int j=0; j<3; j++)
			{
				unsigned int index = i*3 + j;
				rDampMatrix(index,index) = temp;
			}
		}

		KRATOS_CATCH("")
	}

	//***********************************************************************************
	//***********************************************************************************
	
	void MembraneElement::FinalizeSolutionStep(
		ProcessInfo& rCurrentProcessInfo)
	{
		for(unsigned int i = 0; i<mConstitutiveLawVector.size(); i++)
			mConstitutiveLawVector[i]->FinalizeSolutionStep( GetProperties(),
				GetGeometry(),
				row( GetGeometry().ShapeFunctionsValues(), i ),
				rCurrentProcessInfo);	
	}

	//***********************************************************************************
	//***********************************************************************************
	
	void MembraneElement::GetValuesVector(
		  Vector& values,
		  int Step)

	{
		const unsigned int number_of_nodes = GetGeometry().size();
		const unsigned int MatSize = number_of_nodes*3;
		if(values.size() != MatSize)
			values.resize(MatSize);
		for (unsigned int i=0;i<number_of_nodes;i++)
		{
			const array_1d<double,3>& disp = GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT,Step);
			unsigned int index = i*3;
			values[index]   = disp[0];
			values[index+1] = disp[1];
			values[index+2] = disp[2];
		}
	}

	//***********************************************************************************
	//***********************************************************************************

	void MembraneElement::GetFirstDerivativesVector(
		Vector& values,
		int Step)

	{
		const unsigned int number_of_nodes = GetGeometry().size();
		const unsigned int MatSize = number_of_nodes*3;
		if(values.size() != MatSize)
			values.resize(MatSize);
		for (unsigned int i=0;i<number_of_nodes;i++)
		{
			const array_1d<double,3>& vel = GetGeometry()[i].FastGetSolutionStepValue(VELOCITY,Step);
			unsigned int index = i*3;			
			values[index]   = vel[0];
			values[index+1] = vel[1];
			values[index+2] = vel[2];
		}

	}

	//***********************************************************************************
	//***********************************************************************************

	void MembraneElement::GetSecondDerivativesVector(
		Vector& values,
		int Step)
	
	{
		const unsigned int number_of_nodes = GetGeometry().size();
		const unsigned int MatSize = number_of_nodes*3;
		if(values.size() != MatSize)
			values.resize(MatSize);
		for (unsigned int i=0;i<number_of_nodes;i++)
		{
			const array_1d<double,3>& acc = GetGeometry()[i].FastGetSolutionStepValue(ACCELERATION,Step);
			unsigned int index = i*3;
			values[index]   = acc[0];
			values[index+1] = acc[1];
			values[index+2] = acc[2];
		}
	}
	
	//***********************************************************************************
	//***********************************************************************************
		  // --------- //
		 //  PRIVATE  //
		// --------- //


	//***********************************************************************************
	//This function reads the follower pressure applied to the element
	//if corrections have to be applied to this local pressure this is the best 
	//place to do it
	//***********************************************************************************
	//double MembraneElement::GetElementalPressure(
	//	const ProcessInfo& rCurrentProcessInfo)
	//
	//{
	//	KRATOS_TRY
	//		
	//	//first of all i read the esternal pressure
	//	ElementSources::SourcesArrayType sourceP;
	//	sourceP = (*mpSources)(ST_SURFACE_1);
	//	PropertiesType& properties = *mProperties[0];

	//	double pext = 0.00; 
	//	for (ElementSources::SourcesArrayType::iterator isource = sourceP.begin();
	//		isource != sourceP.end(); isource++)
	//	{
	//		if ((*isource)->Has(PRESSURE))
	//		{
	//			pext = (**isource)(PRESSURE);
	//		}
	//	}
	//	return pext;

	//	KRATOS_CATCH("")
	//}

	//***********************************************************************************
	//***********************************************************************************

	void MembraneElement::CalculateAndAddKm(
		Matrix& K,
		Matrix& msB,
		Matrix& msD,
		double weight)
	
	{
		KRATOS_TRY
		
		unsigned int dim = msB.size2();
		Matrix temp(3,dim);
		noalias(temp) = prod(msD,msB);
		temp *= weight;
		Matrix Km(dim,dim);
		noalias(Km) = prod(trans(msB),temp);
		noalias(K) += Km;

		KRATOS_CATCH("")
	}

	//***********************************************************************************
	//***********************************************************************************

	void MembraneElement::CalculateAndAddKg(
		Matrix& K,
		boost::numeric::ublas::bounded_matrix<double,3,3>& msQ,
		const Matrix& DN_De,
		Vector& msStressVector,
		double weight)

	{
		KRATOS_TRY
		
		unsigned int number_of_nodes = GetGeometry().size();
		Vector s(3);
		noalias(s) = prod(trans(msQ),msStressVector);

		double s11 = s[0];
		double s22 = s[1];
		double s12 = s[2];

		Matrix Kloc(number_of_nodes,number_of_nodes);
		Vector a = ZeroVector(number_of_nodes);
		Vector b = ZeroVector(number_of_nodes);

		for(unsigned int i=0;i<number_of_nodes;i++)
		{
			a[i] = DN_De(i,0);
			b[i] = DN_De(i,1);
		}

		for(unsigned int i=0;i<number_of_nodes;i++)
		{
			for(unsigned int j=0;j<number_of_nodes;j++) 
			{
				Kloc(i,j) = a[i]*a[j]*s11 + b[i]*b[j]*s22 + (b[i]*a[j]+a[i]*b[j])*s12;
			}
		}

		Kloc *= weight;
		ExpandReducedMatrix(K,Kloc);

		KRATOS_CATCH("")
	}

	//***********************************************************************************
	//***********************************************************************************

	void MembraneElement::CalculateAndSubKp(
		Matrix& K,
		array_1d<double,3>& ge,
		array_1d<double,3>& gn,
		const Matrix& DN_De,
		const Vector& N,
		double pressure,
		double weight)

	{
		KRATOS_TRY

		boost::numeric::ublas::bounded_matrix<double,3,3> Kij;
		boost::numeric::ublas::bounded_matrix<double,3,3> Cross_ge;
		boost::numeric::ublas::bounded_matrix<double,3,3> Cross_gn;
		double coeff;
		unsigned int number_of_nodes = GetGeometry().size();

		MakeCrossMatrix(Cross_ge,ge);
		MakeCrossMatrix(Cross_gn,gn);

		for(unsigned int i=0;i<number_of_nodes;i++)
		{
			int RowIndex = i*3;
			for(unsigned int j=0;j<number_of_nodes;j++)
			{
				int ColIndex = j*3;

				coeff = pressure*N[i]*DN_De(j,1)*weight;
				noalias(Kij)  = coeff * Cross_ge;

				coeff = pressure*N[i]*DN_De(j,0)*weight;

				noalias(Kij) -= coeff * Cross_gn;
//Kij *= -1;
				//TAKE CARE: the load correction matrix should be SUBTRACTED not added
				SubtractMatrix(K,Kij,RowIndex,ColIndex);	  
			}
		}

		KRATOS_CATCH("")
	}

	//***********************************************************************************
	//***********************************************************************************

	void MembraneElement::MakeCrossMatrix(
		boost::numeric::ublas::bounded_matrix<double,3,3>& M,
		array_1d<double,3>& U)

	{
		M(0,0) =  0.00;
		M(0,1) = -U[2];
		M(0,2) =  U[1];
		M(1,0) =  U[2];
		M(1,1) =  0.00;
		M(1,2) = -U[0];
		M(2,0) = -U[1];
		M(2,1) =  U[0];
		M(2,2) =  0.00;
	}

	//***********************************************************************************
	//***********************************************************************************

	void MembraneElement::CrossProduct(
		array_1d<double,3>& cross,
		array_1d<double,3>& a,
		array_1d<double,3>& b)

	{
		cross[0] = a[1]*b[2] - a[2]*b[1];
		cross[1] = a[2]*b[0] - a[0]*b[2];
		cross[2] = a[0]*b[1] - a[1]*b[0];
	}

	//***********************************************************************************
	//***********************************************************************************

	void  MembraneElement::ExpandReducedMatrix(
		Matrix& Destination,
		Matrix& ReducedMatrix)

	{
		KRATOS_TRY

		unsigned int size=ReducedMatrix.size2();
		
		for (unsigned int i=0;i<size;i++)
		{
			int rowindex = i*3;
			for (unsigned int j=0;j<size;j++)
			{
				unsigned int colindex = j*3;
				for(unsigned int ii=0;ii<3;ii++)
					Destination(rowindex+ii,colindex+ii)+=ReducedMatrix(i,j);
			}
		}

		KRATOS_CATCH("")
	}

	//***********************************************************************************
	//***********************************************************************************

	void  MembraneElement::SubtractMatrix(
		MatrixType& Destination,
		boost::numeric::ublas::bounded_matrix<double,3,3>& InputMatrix, 
		int InitialRow,
		int InitialCol)

	{
		KRATOS_TRY

		for(unsigned int i=0; i<3; i++)
			for(unsigned int j=0; j<3; j++)
				Destination(InitialRow+i, InitialCol+j) -= InputMatrix(i,j);

		KRATOS_CATCH("")
	}

	//***********************************************************************************
	//***********************************************************************************

	void MembraneElement::CalculateQ(
		boost::numeric::ublas::bounded_matrix<double,3,3>& msQ,
		Matrix& mG)

	{
		KRATOS_TRY

		msQ(0,0)=pow(mG(0,0),2);  
		msQ(1,0)=pow(mG(0,1),2);       msQ(1,1)=pow(mG(1,1),2);  msQ(1,2)=mG(0,1)*mG(1,1);
		msQ(2,0)=2.00*mG(0,0)*mG(0,1);					         msQ(2,2)=mG(0,0)*mG(1,1);

		KRATOS_CATCH("")
	}

	//***********************************************************************************
	//***********************************************************************************

	void MembraneElement::CalculateB(
		Matrix& msB,
		boost::numeric::ublas::bounded_matrix<double,3,3>& msQ,
		const Matrix& DN_De,
		array_1d<double,3>& ge,
		array_1d<double,3>& gn)

	{
		KRATOS_TRY
		
		const unsigned int number_of_nodes = GetGeometry().size();
		Matrix b(3,number_of_nodes*3);

		for (unsigned int i=0;i<number_of_nodes;i++)
		{
			unsigned int index = 3*i;

			//first line
			b(0,index)   = DN_De(i,0)*ge[0];
			b(0,index+1) = DN_De(i,0)*ge[1];
			b(0,index+2) = DN_De(i,0)*ge[2];

			//second line
			b(1,index)   = DN_De(i,1)*gn[0];
			b(1,index+1) = DN_De(i,1)*gn[1];
			b(1,index+2) = DN_De(i,1)*gn[2];

			//third line
			b(2,index)   = DN_De(i,1)*ge[0]+DN_De(i,0)*gn[0];
			b(2,index+1) = DN_De(i,1)*ge[1]+DN_De(i,0)*gn[1];
			b(2,index+2) = DN_De(i,1)*ge[2]+DN_De(i,0)*gn[2];
		}
		msB = prod(msQ,b);

		KRATOS_CATCH("")
	}

	//***********************************************************************************
	//***********************************************************************************

	void MembraneElement::CalculateJ(
		boost::numeric::ublas::bounded_matrix<double,2,2>& j,
		array_1d<double,3>& ge,
		array_1d<double,3>& gn,
		array_1d<double,3>& v3)
	
	{
		double Norm_v3 = norm_2(v3);
		double Norm_ge = norm_2(ge);

		double j11 = Norm_ge; 
		double j12 = inner_prod(ge,gn)/Norm_ge;
		double j22 = Norm_v3 / Norm_ge;

		j(0,0) = j11;    j(0,1) = j12;
		j(1,0) = 0.00;   j(1,1) = j22;
	}

	//***********************************************************************************
	//***********************************************************************************

	void MembraneElement::CalculateStrain(
		Vector& StrainVector,
		boost::numeric::ublas::bounded_matrix<double,2,2>& C)

	{
		KRATOS_TRY

		StrainVector[0] = 0.5*(C(0,0) - 1.00);
		StrainVector[1] = 0.5*(C(1,1) - 1.00);
		StrainVector[2] = C(0,1);

		KRATOS_CATCH("")
	}

	//***********************************************************************************
	//***********************************************************************************

	void MembraneElement::CalculateAndAdd_BodyForce(
		const Vector& N,
		const ProcessInfo& rCurrentProcessInfo,
		array_1d<double,3>& BodyForce,
		VectorType& rRightHandSideVector,
		double weight)

	{
		KRATOS_TRY

		unsigned int number_of_nodes = GetGeometry().size();

		for (unsigned int i=0;i<number_of_nodes;i++)
		{
			int index = 3*i;
			for (unsigned int j=0; j<3; j++)
				rRightHandSideVector[index+j] += weight*N[i]*BodyForce[j];		  
		} 

		KRATOS_CATCH("")
	}

	//***********************************************************************************
	//***********************************************************************************

	void MembraneElement::CalculateAndAdd_PressureForce(
		VectorType& residualvector,
		const Vector& N,
		const array_1d<double,3>& v3,
		double pressure,
		double weight,
		const ProcessInfo& rCurrentProcessInfo)

	{
		KRATOS_TRY

		unsigned int number_of_nodes = GetGeometry().size();

		for (unsigned int i=0;i<number_of_nodes;i++)
		{
			int index = 3*i;
			double coeff = pressure * N[i] * weight;
			residualvector[index]   += coeff*v3[0];
			residualvector[index+1] += coeff*v3[1];
			residualvector[index+2] += coeff*v3[2];
		}

		KRATOS_CATCH("")
	}

	//***********************************************************************************
	//***********************************************************************************

	void MembraneElement::CalculateAll(
		MatrixType& rLeftHandSideMatrix,
		VectorType& rRightHandSideVector,
		const ProcessInfo& rCurrentProcessInfo,
		bool CalculateStiffnessMatrixFlag,
		bool CalculateResidualVectorFlag)
		
	{
		KRATOS_TRY

		const unsigned int number_of_nodes = GetGeometry().size();
		unsigned int MatSize = number_of_nodes * 3;

		//resizing as needed the LHS
		if (CalculateStiffnessMatrixFlag == true) //calculation of the matrix is required
		{
			if(rLeftHandSideMatrix.size1() != MatSize)
				rLeftHandSideMatrix.resize(MatSize,MatSize);
			noalias(rLeftHandSideMatrix) = ZeroMatrix(MatSize,MatSize); //resetting LHS
		}

		//resizing as needed the RHS
		if (CalculateResidualVectorFlag == true) //calculation of the matrix is required
		{
			if(rRightHandSideVector.size() != MatSize)
				rRightHandSideVector.resize(MatSize);
			rRightHandSideVector = ZeroVector(MatSize); //resetting RHS
		}

		//reading integration points and local gradients
		const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints();
		const GeometryType::ShapeFunctionsGradientsType& DN_DeContainer = GetGeometry().ShapeFunctionsLocalGradients();
		const Matrix& Ncontainer = GetGeometry().ShapeFunctionsValues();

		//calculating actual jacobian
		GeometryType::JacobiansType J;
		J = GetGeometry().Jacobian(J);

		//auxiliary terms
		array_1d<double,3> BodyForce;
		double elem_positive_face_pressure = 0.0;
		double elem_negative_face_pressure = 0.0;
		bool zero_positive_nodal_pressure = false;
		bool zero_negative_nodal_pressure = false;

		for(unsigned int k = 0; k<GetGeometry().size();k++)
		{
			double temp = GetGeometry()[k].FastGetSolutionStepValue(NEGATIVE_FACE_PRESSURE);
			elem_negative_face_pressure += temp;
			if(temp == 0)
				zero_negative_nodal_pressure = true;

			temp = GetGeometry()[k].FastGetSolutionStepValue(POSITIVE_FACE_PRESSURE);
			elem_positive_face_pressure += temp;
			if(temp == 0)
				zero_positive_nodal_pressure = true;
		}
		if(zero_negative_nodal_pressure == true)
			elem_negative_face_pressure = 0.0;
		else
			elem_negative_face_pressure /= (double) GetGeometry().size();
		if(zero_positive_nodal_pressure == true)
			elem_positive_face_pressure = 0.0;
		else
			elem_positive_face_pressure /= (double) GetGeometry().size();

/*	
		for(unsigned int k = 0; k<PressureOnNodes.size();k++)
		{
			//double temp = GetGeometry()[i].FastGetSolutionStepValue(PRESSURE);
			double temp = GetGeometry()[k].FastGetSolutionStepValue(NEGATIVE_FACE_PRESSURE);
			temp -= GetGeometry()[k].FastGetSolutionStepValue(POSITIVE_FACE_PRESSURE);
			PressureOnNodes[k] += temp;
		}

		bool zero_nodal_pressure = false;
		for(unsigned int k = 0; k<PressureOnNodes.size();k++)
		{
 			if(GetGeometry()[k].FastGetSolutionStepValue(NEGATIVE_FACE_PRESSURE) == 0 ||
			   GetGeometry()[k].FastGetSolutionStepValue(POSITIVE_FACE_PRESSURE) == 0 )
				zero_nodal_pressure = true;
                }
*/
		array_1d<double,3> ge;
		array_1d<double,3> gn;
		array_1d<double,3> v3;
		for(unsigned int PointNumber=0;PointNumber<integration_points.size();PointNumber++)
		{
			double IntegrationWeight = GetGeometry().IntegrationPoints()[PointNumber].Weight();

			ge[0]=J[PointNumber](0,0); gn[0]=J[PointNumber](0,1);
			ge[1]=J[PointNumber](1,0); gn[1]=J[PointNumber](1,1);
			ge[2]=J[PointNumber](2,0); gn[2]=J[PointNumber](2,1);

			CrossProduct(v3,ge,gn);
			boost::numeric::ublas::bounded_matrix<double,2,2> j;
			CalculateJ(j,ge,gn,v3);

			// calculation of matrix g = jtrans*j;
			boost::numeric::ublas::bounded_matrix<double,2,2> g;
			noalias(g) = prod(trans(j),j);

			// calculation of the Right Cauchy-Green Tensor C = Gtrans*g*G
			boost::numeric::ublas::bounded_matrix<double,2,2> tmp;
			tmp = prod(g,mG_Vector[PointNumber]);
			noalias(msC) = prod(trans(mG_Vector[PointNumber]),tmp);

			// calculation of the StrainVector
			CalculateStrain(msStrainVector, msC);
			mStrainsVector[PointNumber] = msStrainVector;	//saving the strain vector
            
            mConstitutiveLawVector[PointNumber]->CalculateMaterialResponse(
                    msStrainVector,
                    ZeroMatrix(1),
                    msStressVector,
                    msD,
                    rCurrentProcessInfo,
                    GetProperties(),
                    GetGeometry(),
                    row(Ncontainer,PointNumber),
                    true,
                    (int)CalculateStiffnessMatrixFlag,
                    true );

			// material update (considering the level of strain achieved)
// 			mConstitutiveLawVector[PointNumber]->UpdateMaterial( msStrainVector,
// 						GetProperties(),
// 						GetGeometry(),
// 						row(Ncontainer,PointNumber),
// 						rCurrentProcessInfo );

			// calculation of the StressVector (PK2)			
// 			mConstitutiveLawVector[PointNumber]->CalculateStress(msStrainVector,msStressVector); 

			noalias(mStressesVector[PointNumber]) = ZeroVector(6);
			Calculate_GlobalStressVector(mStressesVector[PointNumber], msStressVector, mV1[PointNumber], mV2[PointNumber]);	//saving the stress vector

			// calculation of the actual thickness
//			mThickness[PointNumber]  = mConstitutiveLawVector[PointNumber]->CalculateThicknessRatio(msStrainVector);
//			mThickness[PointNumber] *= mThickness0;

/////////////////////////
//if (rCurrentProcessInfo[TIME] > 4.8e-006)
//{
//KRATOS_WATCH(Id() );
//KRATOS_WATCH(msStrainVector );
//KRATOS_WATCH(msStressVector );
//KRATOS_WATCH(mThickness[PointNumber] );
//}
/////////////////////////

			// calculating the pressure on the gauss point
			double pressure = elem_negative_face_pressure - 				  elem_positive_face_pressure;
//			double pressure = 0.00;
//			for (unsigned int ii=0; ii<number_of_nodes; ii++) 
//				pressure += Ncontainer(PointNumber,ii)*PressureOnNodes[ii];

//			if(zero_nodal_pressure == true) //case in which a node has no pressure applied
//				pressure = 0.0;

			CalculateQ(msQ, mG_Vector[PointNumber]);
			CalculateB(msB, msQ, DN_DeContainer[PointNumber], ge, gn);

			// integration on the REFERENCE CONFIGURATION
			double DetJ0 = mDetJ0[PointNumber];
			double IntToReferenceWeight = IntegrationWeight * DetJ0 * mThickness0;

			// LEFT HAND SIDE MATRIX
			if (CalculateStiffnessMatrixFlag == true)
			{
// 				mConstitutiveLawVector[PointNumber]->CalculateConstitutiveMatrix(msStrainVector,msD);

				//adding contributions to the stiffness matrix
				CalculateAndAddKm(rLeftHandSideMatrix, msB, msD, IntToReferenceWeight);
				CalculateAndAddKg(rLeftHandSideMatrix, msQ, DN_DeContainer[PointNumber], msStressVector, IntToReferenceWeight);
				if (pressure != 0.00) 
				{
					const Matrix& dN_dE = DN_DeContainer[PointNumber];
					CalculateAndSubKp(rLeftHandSideMatrix,ge,gn,dN_dE,row(Ncontainer,PointNumber),pressure,IntegrationWeight);
				}
			}
			
			// RIGHT HAND SIDE VECTOR
			if (CalculateResidualVectorFlag == true) //calculation of the matrix is required
			{
				//contribution to external forces
				BodyForce = GetProperties()[BODY_FORCE];
				// operation performed: rRightHandSideVector += ExtForce*IntToReferenceWeight
				CalculateAndAdd_BodyForce(row(Ncontainer,PointNumber),rCurrentProcessInfo,BodyForce,rRightHandSideVector,IntToReferenceWeight);
				// operation performed: rRightHandSideVector -= Weight*IntForce
				noalias(rRightHandSideVector) -= IntToReferenceWeight * prod(trans(msB),msStressVector);
				if (pressure!= 0.00)
					CalculateAndAdd_PressureForce(rRightHandSideVector,row(Ncontainer,PointNumber),v3,pressure,IntegrationWeight,rCurrentProcessInfo);
			}
		}
//KRATOS_WATCH(Id() )
//KRATOS_WATCH(rRightHandSideVector)

		KRATOS_CATCH("")
	}

	//***********************************************************************************
	//***********************************************************************************

	void MembraneElement::Calculate_GlobalStressVector(
		array_1d<double,6>& GlobalVector,
		Vector& LocalStressVector,
		array_1d<double,3>& v1,
		array_1d<double,3>& v2)

	{	
		KRATOS_TRY

		array_1d<double,6> temp;

		//adding the component S11
		noalias(temp)  = VoigtTensorComponents(v1,v1);
		temp *= LocalStressVector[0];
		noalias(GlobalVector) += temp;

		//adding the component S22
		noalias(temp)  = VoigtTensorComponents(v2,v2);
		temp *= LocalStressVector[1];
		noalias(GlobalVector) += temp;

		//adding the component S12 (& S21)
		noalias(temp)  = VoigtTensorComponents(v1,v2);
		noalias(temp) += VoigtTensorComponents(v2,v1);
		temp *= LocalStressVector[2];
		noalias(GlobalVector) += temp;

		KRATOS_CATCH("")
	}

	//***********************************************************************************
	//***********************************************************************************

	//auxiliary function needed in the calculation of output stresses
	inline array_1d<double,6> MembraneElement::VoigtTensorComponents(
		array_1d<double,3>& a,
		array_1d<double,3>& b)

	{
		array_1d<double,6> v;

		v[0] = a[0]*b[0];
		v[1] = a[1]*b[1];
		v[2] = a[2]*b[2];
		v[3] = a[0]*b[1];
		v[4] = a[1]*b[2];
		v[5] = a[0]*b[2];

		return v;
	}
	
        //************************************************************************************
        //************************************************************************************
	void MembraneElement::GetValueOnIntegrationPoints(const Variable<Matrix>& rVariable,
			std::vector<Matrix>& rValues, const ProcessInfo& rCurrentProcessInfo)
			{
				if(rVariable==GREEN_LAGRANGE_STRAIN_TENSOR)
				{
					CalculateOnIntegrationPoints( rVariable, rValues, rCurrentProcessInfo );
				}
				if(rVariable==PK2_STRESS_TENSOR)
				{
					CalculateOnIntegrationPoints( rVariable, rValues, rCurrentProcessInfo );
				}
			}

	//***********************************************************************************
	//***********************************************************************************

}	// Namespace Kratos.
