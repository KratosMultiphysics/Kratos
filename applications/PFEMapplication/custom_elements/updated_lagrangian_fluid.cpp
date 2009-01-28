//   
//   Project Name:        Kratos       
//   Last modified by:    $Author: rrossi $
//   Date:                $Date: 2007-01-23 15:38:00 $
//   Revision:            $Revision: 1.1 $
//
//
//this is an element for 3D

// System includes 


// External includes 


// Project includes 
#include "includes/define.h"
#include "custom_elements/updated_lagrangian_fluid.h"
#include "utilities/math_utils.h"

//#include "constitutive_laws/constitutive_laws.h"
//#include "constitutive_laws/isotropic_planestress_wrinkling_new.h"
#include "PFEM_application.h"

namespace Kratos
{
    Matrix UpdatedLagrangianFluid::msB(0,0);
	Matrix UpdatedLagrangianFluid::msF(0,0);
	Matrix UpdatedLagrangianFluid::msD(0,0);
	Matrix UpdatedLagrangianFluid::msC(0,0);
	Matrix UpdatedLagrangianFluid::msCapx(0,0);
	
	Vector UpdatedLagrangianFluid::msStrainVector(0);
	Vector UpdatedLagrangianFluid::msStressVector(0);
	
	Matrix UpdatedLagrangianFluid::msDN_DX(0,0);

	//************************************************************************************
	//************************************************************************************
	UpdatedLagrangianFluid::UpdatedLagrangianFluid(IndexType NewId, GeometryType::Pointer pGeometry)
		: Element(NewId, pGeometry)
	{		
		//DO NOT ADD DOFS HERE!!!
	}

	//************************************************************************************
	//************************************************************************************
	UpdatedLagrangianFluid::UpdatedLagrangianFluid(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
		: Element(NewId, pGeometry, pProperties)
	{
		unsigned int dim = GetGeometry().WorkingSpaceDimension();

		mpReferenceGeometry = GetGeometry().Clone();

		//setting up the nodal degrees of freedom
		for(unsigned int i = 0 ; i != GetGeometry().size() ; ++i)
		{
			(GetGeometry()[i].pAddDof(DISPLACEMENT_X,REACTION_X));
			(GetGeometry()[i].pAddDof(DISPLACEMENT_Y,REACTION_Y));
			if(dim == 3)
				(GetGeometry()[i].pAddDof(DISPLACEMENT_Z,REACTION_Z));
		}

		//initializing static variables
		unsigned int dimension = GetGeometry().WorkingSpaceDimension();
		unsigned int dim2 = GetGeometry().size()*dimension;
		unsigned int StrainSize;
		if (dimension==2) 
			StrainSize = 3;
		else 
			StrainSize = 6;
		msB.resize(StrainSize,dim2);
		msF.resize(dimension,dimension);
		msD.resize(StrainSize,StrainSize);
		msC.resize(dimension,dimension);
		msCapx.resize(StrainSize,StrainSize);
		msStrainVector.resize(StrainSize);
		msStressVector.resize(StrainSize);
		msDN_DX.resize(GetGeometry().size(),dimension);
	}

	Element::Pointer UpdatedLagrangianFluid::Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const
	{
		return Element::Pointer(new UpdatedLagrangianFluid(NewId, GetGeometry().Create(ThisNodes), pProperties));
	}

	UpdatedLagrangianFluid::~UpdatedLagrangianFluid()
	{
	}
	//************************************************************************************
	//************************************************************************************
	void UpdatedLagrangianFluid::InitializeSolutionStep(ProcessInfo &CurrentProcessInfo)
		{
		KRATOS_TRY

		//calculating the values of Jacobian at the integration points 
		//JACOBIAN IS CONST FOR A LINEAR ELEMENT
		const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints();

		//resizing jacobian inverses containers
		mInvJ0.resize(integration_points.size());
		mDetJ0.resize(integration_points.size());

		GeometryType::JacobiansType J0;
		J0 = GetGeometry().Jacobian(J0);  
		mTotalDomainInitialSize = 0.00;
	
		//calculating the inverse J0
		for(unsigned int PointNumber = 0; PointNumber<integration_points.size(); PointNumber++)
		{
			//getting informations for integration
			double IntegrationWeight = integration_points[PointNumber].Weight();
			//calculating and storing inverse of the jacobian and the parameters needed
			MathUtils<double>::InvertMatrix(J0[PointNumber],mInvJ0[PointNumber],mDetJ0[PointNumber]);				
			//calculating the total area
			mTotalDomainInitialSize += mDetJ0[PointNumber]*IntegrationWeight;
		}
		KRATOS_CATCH("")
	}
	
	//************************************************************************************
	//************************************************************************************
	void UpdatedLagrangianFluid::CalculateAll(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, 
										ProcessInfo& rCurrentProcessInfo,
										bool CalculateStiffnessMatrixFlag,
										bool CalculateResidualVectorFlag)
	{
		KRATOS_TRY

		unsigned int number_of_nodes = GetGeometry().size();
		unsigned int dim = GetGeometry().WorkingSpaceDimension();
		//resizing as needed the LHS
		//e.g. in 2D and using linear triangles => MatSize = 3*2=6
		unsigned int MatSize=number_of_nodes*dim;

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
		const GeometryType::ShapeFunctionsGradientsType& DN_De = GetGeometry().ShapeFunctionsLocalGradients();
		const Matrix& Ncontainer = GetGeometry().ShapeFunctionsValues();

		//calculating actual jacobian
		GeometryType::JacobiansType J;
		GetGeometry().Jacobian(J); 
		
		double dd=0.0;
		Matrix Jinvv (2,2);
		
		MathUtils<double>::InvertMatrix(J[0],Jinvv,dd);				
		
		if (dd<=0.0)
		{
			std::cout<<"ERROR........ negative det"<<std::endl;
			std::cout<<"ERROR........ negative det"<<std::endl;
			std::cout<<"ERROR........ negative det"<<std::endl;
			std::cout<<"ERROR........ negative det"<<std::endl;
			std::cout<<"ERROR........ negative det"<<std::endl;
			std::cout<<"ERROR........ negative det"<<std::endl;
			std::cout<<"ERROR........ negative det"<<std::endl;
		}
		//auxiliary terms
		Vector BodyForce;		
		//sizing work matrices
		for(unsigned int PointNumber = 0; PointNumber<integration_points.size(); PointNumber++)
		{
			//Calculating the cartesian derivatives (it is avoided storing them to minimize storage)
			msDN_DX = ZeroMatrix(msDN_DX.size1(),msDN_DX.size2());
			noalias(msDN_DX) = prod(DN_De[PointNumber],mInvJ0[PointNumber]);

			//deformation gradient
			noalias(msF) = prod(J[PointNumber],mInvJ0[PointNumber]);

			//strain calculation
			noalias(msC) = prod(trans(msF),msF);
			CalculateStrain(msC,msStrainVector);

		//	double E = GetProperties()[YOUNG_MODULUS];
			double NU = GetProperties()[VISCOSITY];
			double RO = GetProperties()[DENSITY];
			double K = GetProperties()[BULK_MODULUS];
			K = 100000.0;


			//Obtain PK2 stress vector - first parameter is the bulk modulus K
			CalculateIncPK2Stress(K, msStressVector, msF);
							
			if (dim==2)
			{
			msCapx(0,0) = 4.0/3.0;		msCapx(0,1)=-2.0/3.0;		msCapx(0,2)=0.0;
			msCapx(1,0) =-2.0/3.0;		msCapx(1,1)=4.0/3.0;		msCapx(1,2)=0.0;
			msCapx(2,0)=0.0;			msCapx(2,1)=0.0;			msCapx(2,2)=1.0;
			}
			if (dim==3)
			{
			//BELOW - BULLSHIT!
			msCapx(0,0)=	 4.0/3.0;	 msCapx(0,1)= -2.0/3.0;		msCapx(0,2)= -2.0/3.0;		msCapx(0,3)=0;	msCapx(0,4)=	0;	msCapx(0,5)=0;
			msCapx(1,0)=	-2.0/3.0;	 msCapx(1,1)= 4.0/3.0;		msCapx(1,2)= -2.0/3.0;		msCapx(1,3)=0;	msCapx(1,4)=	0;	msCapx(1,5)=0;
			msCapx(2,0)=	-2.0/3.0;	 msCapx(2,1)= -2.0/3.0;		msCapx(2,2)=  4.0/3.0;		msCapx(2,3)=0;	msCapx(2,4)=	0;	msCapx(2,5)=0;
			msCapx(3,0)=	0;			 msCapx(3,1)=0;				msCapx(3,2)=0;				msCapx(3,3)=1;	msCapx(3,4)=	0;	msCapx(3,5)=0;
			msCapx(4,0)=	0;			 msCapx(4,1)=0;				msCapx(4,2)=0;				msCapx(4,3)=0;	msCapx(4,4)=	1;	msCapx(4,5)=0;
			msCapx(5,0)=	0;			 msCapx(5,1)=0;				msCapx(5,2)=0;				msCapx(5,3)=0;	msCapx(5,4)=    0;	msCapx(5,5)=1;
			}
			msCapx*=RO*NU/(rCurrentProcessInfo[DELTA_TIME]);
			//and now adding the terms accounting for the volume change
			if (dim==2)
			{
			msCapx(0,0)+=K/2.0;	msCapx(0,1)+=K/2.0;
			msCapx(1,1)+=K/2.0;	msCapx(1,1)+=K/2.0;
			}
			if (dim==3)
			{
			msCapx(0,0)+=K/3.0;	msCapx(0,1)+=K/3.0;	msCapx(0,2)+=K/3.0;
			msCapx(1,0)+=K/3.0;	msCapx(1,1)+=K/3.0;	msCapx(1,2)+=K/3.0;
			msCapx(2,0)+=K/3.0;	msCapx(2,1)+=K/3.0;	msCapx(2,2)+=K/3.0;
			}
						
			//calculating operator B
			CalculateB(msB,msF,msDN_DX,msStrainVector.size());

			//calculating weights for integration on the reference configuration
			double IntToReferenceWeight = integration_points[PointNumber].Weight() * mDetJ0[PointNumber];
			//if (dim == 2) IntToReferenceWeight *= GetProperties()[THICKNESS];
			if (dim == 2) IntToReferenceWeight *= 1;
			
			//WHAT TO DO WITH 3D : IntToReferenceWeight??!?!?


			if (CalculateStiffnessMatrixFlag == true) //calculation of the matrix is required
			{			
				//contributions to stiffness matrix calculated on the reference config
				//Kmat=Inetgral(BT*Capx*B)
				noalias(rLeftHandSideMatrix) += prod(trans(msB),(IntToReferenceWeight)*Matrix(prod(msCapx,msB)) ); //to be optimized to remove the temporary
				//here we use the incremental stress (only the increment  - no histroy part)
				CalculateAndAddKg(rLeftHandSideMatrix,msDN_DX,msStressVector,IntToReferenceWeight);
			}

			if (CalculateResidualVectorFlag == true) //calculation of the matrix is required
			{
				//contribution to external forces 
				BodyForce = GetProperties()[BODY_FORCE];
				// operation performed: rRightHandSideVector += ExtForce*IntToReferenceWeight
				CalculateAndAdd_ExtForceContribution(row(Ncontainer,PointNumber),rCurrentProcessInfo,BodyForce,rRightHandSideVector,IntToReferenceWeight);
				noalias(rRightHandSideVector) -= IntToReferenceWeight * prod(trans(msB),msStressVector);
			/*	KRATOS_WATCH("INSIDE ELEMENT")
				KRATOS_WATCH(rRightHandSideVector)*/
			}	
			//CHAPUZA!!!!!!!!! - if inverted elements - dont add anything
			if (dd<=0)
			{
			std::cout<<"BAD ONE - not adding contribs"<<std::endl;
			std::cout<<"BAD ONE - not adding contribs"<<std::endl;
			rRightHandSideVector=ZeroVector(MatSize);
			rLeftHandSideMatrix=ZeroMatrix(MatSize,MatSize);
			KRATOS_WATCH("LEFT AND RIGHT HAND SIDE VECTOR");
			KRATOS_WATCH(rRightHandSideVector);
			KRATOS_WATCH(rLeftHandSideMatrix);
			}
			if (dd==0.0 || dd==0)
				{
					std::cout<<"THREE NODES OF ELEMENT ON ONE LINE"<<std::endl;
					std::cout<<"THREE NODES OF ELEMENT ON ONE LINE"<<std::endl;
					std::cout<<"THREE NODES OF ELEMENT ON ONE LINE"<<std::endl;
					std::cout<<"THREE NODES OF ELEMENT ON ONE LINE"<<std::endl;
					std::cout<<"THREE NODES OF ELEMENT ON ONE LINE"<<std::endl;
					std::cout<<"TIME IS"<<rCurrentProcessInfo[TIME]<<std::endl;
				}

		}
		
		KRATOS_CATCH("")
	}

	//************************************************************************************
	//************************************************************************************
	void UpdatedLagrangianFluid::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
	{
//KRATOS_WATCH("UpdatedLagrangianFluid::CalculateRightHandSide");
		//calculation flags
		bool CalculateStiffnessMatrixFlag = false;
		bool CalculateResidualVectorFlag = true;
		MatrixType temp = Matrix();
		
		CalculateAll(temp, rRightHandSideVector, rCurrentProcessInfo, CalculateStiffnessMatrixFlag, CalculateResidualVectorFlag);
//KRATOS_WATCH("UpdatedLagrangianFluid::CalculateRightHandSide finished");
	}

	//************************************************************************************
	//************************************************************************************
	void UpdatedLagrangianFluid::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
	{
		//calculation flags
		bool CalculateStiffnessMatrixFlag = true;
		bool CalculateResidualVectorFlag = true;
		
		CalculateAll(rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo, CalculateStiffnessMatrixFlag, CalculateResidualVectorFlag);
	}

	//************************************************************************************
	//************************************************************************************
	double UpdatedLagrangianFluid::CalculateIntegrationWeight(double GaussPointWeight, double DetJ0)
	{
		//to permorm the integration over the reference domain we need to include 
		// the thickness in 2D
		unsigned int dimension = mpReferenceGeometry->WorkingSpaceDimension();
		double weight = GaussPointWeight;
	  
		weight *= DetJ0;
		//if (dimension == 2) weight *= GetProperties()[THICKNESS];
		if (dimension == 2) weight *= 1;
		return weight;
	}

	////************************************************************************************
	////************************************************************************************
	void UpdatedLagrangianFluid::FinalizeSolutionStep(ProcessInfo& CurrentProcessInfo)
	{
	/*	for(unsigned int i = 0; i < mConstitutiveLawVector.size(); i++)
			mConstitutiveLawVector[i]->FinalizeSolutionStep(CurrentProcessInfo);*/
	}
	//************************************************************************************
	//************************************************************************************
	inline void UpdatedLagrangianFluid::CalculateAndAdd_ExtForceContribution(
		const Vector& N,
		const ProcessInfo& CurrentProcessInfo,
		Vector& BodyForce,
		VectorType& rRightHandSideVector,
		double weight
		)
	{
		KRATOS_TRY
		unsigned int number_of_nodes = mpReferenceGeometry->PointsNumber();
		unsigned int dimension = mpReferenceGeometry->WorkingSpaceDimension();

		for (unsigned int i=0;i<number_of_nodes;i++)
		{
			int index = dimension*i;
			for (unsigned int j=0; j<dimension; j++)  rRightHandSideVector[index+j] += weight*N[i]*BodyForce[j];		  
		}
		//calculating actual jacobian
		
		GeometryType::JacobiansType J;
		GetGeometry().Jacobian(J); 
		
		double dd=0.0;
		Matrix Jinvv (2,2);
		
		MathUtils<double>::InvertMatrix(J[0],Jinvv,dd);				
		
		if (dd<=0.0)
		{
			for (unsigned int i=0;i<number_of_nodes;i++)
			{
			int index = dimension*i;
			for (unsigned int j=0; j<dimension; j++)  rRightHandSideVector[index+j] = 0.0;		  
			}
		}
		/*KRATOS_WATCH(BodyForce)
		KRATOS_WATCH(rRightHandSideVector)*/
		KRATOS_CATCH("")
	}



	//************************************************************************************
	//************************************************************************************
	//!!! THIS FUNCTION calculates the INCREMENT of PK2-stress from t=n to n+1
	void UpdatedLagrangianFluid::CalculateIncPK2Stress(double& BulkModulus, Vector& StressVector, Matrix& F)
	{
		KRATOS_TRY
			unsigned int dim = mpReferenceGeometry->WorkingSpaceDimension();
			Matrix StressTensor (dim,dim);
			Matrix PK2 (dim, dim);
			Matrix VelGrad (dim, dim);
			Matrix invF(dim, dim);
			Matrix invTransF(dim, dim);
			
			Matrix aux(dim,dim+1);
			//calculating viscosity and density
			double nu = GetProperties()[VISCOSITY]; double ro = GetProperties()[DENSITY];
			//calculating velocity gradient
			for(unsigned int iii = 0; iii < GetGeometry().size(); iii++)
				{
				const array_1d<double,3>& v = GetGeometry()[iii].FastGetSolutionStepValue(VELOCITY);
				for(unsigned int j=0; j <dim; j++)
				aux(j,iii) = v[j];
				//aux contains nodal velocities
				}
			//std::cout<<"AUX from PK2 stress FCT"<<aux<<std::endl;
			
			noalias(VelGrad) = prod(aux,msDN_DX);
			//this is the deviatoric part of Cauchy stress
			noalias(StressTensor) = 0.5 * nu * ro * (VelGrad + trans(VelGrad));
					
			if (dim==2)
			{
			//double av_pressure = 0.333333333333333*(GetGeometry()[0].FastGetSolutionStepValue(PRESSURE,1)+GetGeometry()[1].FastGetSolutionStepValue(PRESSURE,1)+GetGeometry()[2].FastGetSolutionStepValue(PRESSURE,1));
			//HISTORY (of the volumetric stress) IS NOT ADDED - we are calculating increment
			const double volumetric_stress = BulkModulus*0.5*(msStrainVector[0]+msStrainVector[1]);//+ av_pressure;
			//	 and will transform it: S=J*F-1*sigma*F-T
			double detF;
			MathUtils<double>::InvertMatrix(F,invF,detF);
			
			Matrix temp = prod (   StressTensor, trans(invF ));
			noalias(PK2) = detF * prod(  invF    ,       temp  );
			//now we apply Voigt rule:
			//and add the volumetric part - resulting from the trace of G-L strain tensor
			//std::cout<<" STRESS VEC " <<PK2<<std::endl;
			StressVector[0] = PK2(0,0) + volumetric_stress;
			StressVector[1] = PK2(1,1) + volumetric_stress;
			StressVector[2] = PK2(0,1);
			}
			if (dim == 3)
			{
			//adding the volumetric part to the tensor (k*EPSILONv*Identity) 
			//volumetric strain = trace of the Gr-Lagr strain tensor
			const double volumetric_stress = BulkModulus*(1.0/3.0)*(msStrainVector[0]+msStrainVector[1]+msStrainVector[2]);
			StressTensor(0,0) += volumetric_stress;
			StressTensor(1,1) += volumetric_stress;
			StressTensor(2,2) += volumetric_stress;
			//now we have a Cauchy stress tensor and will transform it: S=J*F-1*sigma*F-T
			
			//double detF;
			/*MathUtils<double>::InvertMatrix(F,invF,detF);
			MathUtils<double>::InvertMatrix(trans(F),invTransF,detTransF);

			PK2 = prod(  invF    ,       prod (   StressTensor, invTransF )  );
			PK2 *= detF;*/
						
			//now we apply Voigt rule:
			StressVector[0] = PK2(0,0);
			StressVector[1] = PK2(1,1);
			StressVector[2] = PK2(2,2);
			StressVector[3] = PK2(1,2);
			StressVector[4] = PK2(0,2);
			StressVector[5] = PK2(0,1);
			}
		KRATOS_CATCH("")		
	}

	
	void UpdatedLagrangianFluid::CalculateAndAddKg(
		MatrixType& K,
		Matrix& DN_DX,
		Vector& StressVector,
		double weight)
	{
		KRATOS_TRY
		unsigned int dimension = GetGeometry().WorkingSpaceDimension(); 
		Matrix StressTensor = MathUtils<double>::StressVectorToTensor(StressVector);
		Matrix ReducedKg = prod(DN_DX, weight * Matrix(prod(StressTensor,trans(DN_DX)) ) ); //to be optimized
		MathUtils<double>::ExpandAndAddReducedMatrix(K,ReducedKg,dimension);
		//std::cout<<"Kg "<<ReducedKg<<std::endl;
		KRATOS_CATCH("")
	}  

	//************************************************************************************
	//************************************************************************************
	void UpdatedLagrangianFluid::CalculateStrain(
		const Matrix& C,
		Vector& StrainVector)
	{
		KRATOS_TRY
			unsigned int dimension = mpReferenceGeometry->WorkingSpaceDimension();

		if (dimension==2) 
		{
			if(StrainVector.size() != 3) StrainVector.resize(3);
			StrainVector[0] = 0.5*(C(0,0) - 1.00);
			StrainVector[1] = 0.5*(C(1,1) - 1.00);
			StrainVector[2] = C(0,1);
		}
		if (dimension==3) 
		{
			if(StrainVector.size() != 6) StrainVector.resize(6);
			StrainVector[0] = 0.5*(C(0,0) - 1.00);
			StrainVector[1] = 0.5*(C(1,1) - 1.00);
			StrainVector[2] = 0.5*(C(2,2) - 1.00);
			StrainVector[3] = C(0,1);
			StrainVector[4] = C(1,2);
			StrainVector[5] = C(0,2);
		}
		KRATOS_CATCH("")
	}
	//************************************************************************************
	//************************************************************************************
	void UpdatedLagrangianFluid::CalculateB(
		Matrix& B,
		Matrix& F,
		Matrix& DN_DX,
		unsigned int StrainSize)
	{
		KRATOS_TRY
			const unsigned int number_of_nodes = mpReferenceGeometry->PointsNumber();
		unsigned int dimension = mpReferenceGeometry->WorkingSpaceDimension();
		//
		//unsigned int dim2 = number_of_nodes*dimension;
		//if(B.size1() != StrainSize || B.size2()!=dim2)
		//	B.resize(StrainSize,dim2);
		//Matrix Bi;
		for (unsigned int i=0;i<number_of_nodes;i++)
		{
			unsigned int index = dimension*i;

			if (dimension == 2)
			{
				B(0,index+0)=F(0,0)*DN_DX(i,0);						B(0,index+1)=F(1,0)*DN_DX(i,0);
				B(1,index+0)=F(0,1)*DN_DX(i,1);						B(1,index+1)=F(1,1)*DN_DX(i,1);
				B(2,index+0)=F(0,0)*DN_DX(i,1)+F(0,1)*DN_DX(i,0);	B(2,index+1)=F(1,0)*DN_DX(i,1)+F(1,1)*DN_DX(i,0); 
			}
			else
			{
				B(0,index+0)=F(0,0)*DN_DX(i,0);		B(0,index+1)=F(1,0)*DN_DX(i,0);		B(0,index+2)=F(2,0)*DN_DX(i,0);
				B(1,index+0)=F(0,1)*DN_DX(i,1);		B(1,index+1)=F(1,1)*DN_DX(i,1);		B(1,index+2)=F(2,1)*DN_DX(i,1);
				B(2,index+0)=F(0,2)*DN_DX(i,2);		B(2,index+1)=F(1,2)*DN_DX(i,2);		B(2,index+2)=F(2,2)*DN_DX(i,2);
				B(3,index+0)=F(0,0)*DN_DX(i,1)+F(0,1)*DN_DX(i,0);	B(3,index+1)=F(1,0)*DN_DX(i,1)+F(1,1)*DN_DX(i,0);	B(3,index+2)=F(2,0)*DN_DX(i,1)+F(2,1)*DN_DX(i,0);
				B(4,index+0)=F(0,1)*DN_DX(i,2)+F(0,2)*DN_DX(i,1);	B(4,index+1)=F(1,1)*DN_DX(i,2)+F(1,2)*DN_DX(i,1);	B(4,index+2)=F(2,1)*DN_DX(i,2)+F(2,2)*DN_DX(i,1);
				B(5,index+0)=F(0,2)*DN_DX(i,0)+F(0,0)*DN_DX(i,2);	B(5,index+1)=F(1,2)*DN_DX(i,0)+F(1,0)*DN_DX(i,2);	B(5,index+2)=F(2,2)*DN_DX(i,0)+F(2,0)*DN_DX(i,2);
			}

			//CalculateBi(Bi,F,DN_DX,i);
			//MathUtils<double>::WriteMatrix(B,Bi,0,index);
		}
		KRATOS_CATCH("")
	}

	//************************************************************************************
	//************************************************************************************
	void UpdatedLagrangianFluid::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo)
	{
		int number_of_nodes = GetGeometry().size();
		int dim = GetGeometry().WorkingSpaceDimension();
		unsigned int dim2 = number_of_nodes*dim;
		if(rResult.size() != dim2)
			rResult.resize(dim2);

		for (int i=0;i<number_of_nodes;i++)
		{
			int index = i*dim;
			rResult[index] = GetGeometry()[i].GetDof(DISPLACEMENT_X).EquationId();
			rResult[index + 1] = GetGeometry()[i].GetDof(DISPLACEMENT_Y).EquationId();
			if(dim == 3)
				rResult[index + 2] = GetGeometry()[i].GetDof(DISPLACEMENT_Z).EquationId();
		}

	}

	//************************************************************************************
	//************************************************************************************
	void UpdatedLagrangianFluid::GetDofList(DofsVectorType& ElementalDofList,ProcessInfo& CurrentProcessInfo)
	{
		ElementalDofList.resize(0);

		for (unsigned int i=0;i<GetGeometry().size();i++)
		{
			ElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_X));
			ElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_Y));
			if(GetGeometry().WorkingSpaceDimension() == 3)
				ElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_Z));
		}
//std::cout << "el = " << Id() << std::endl;
//for(int j=0; j<ElementalDofList.size(); j++)
//std::cout << "dof id = " << ElementalDofList[j]->Id() << std::endl;
//std::cout << "ripetizione = " << std::endl;
//for(int j=0; j<ElementalDofList.size(); j++)
//std::cout << "dof id = " << ElementalDofList[j]->Id() << std::endl;
	}

	//************************************************************************************
	//************************************************************************************
	void UpdatedLagrangianFluid::MassMatrix(MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY

		//lumped
		unsigned int dimension = GetGeometry().WorkingSpaceDimension();
		unsigned int NumberOfNodes = GetGeometry().size();
		unsigned int MatSize = dimension * NumberOfNodes;
		if(rMassMatrix.size1() != MatSize)
			rMassMatrix.resize(MatSize,MatSize);

		rMassMatrix = ZeroMatrix(MatSize,MatSize);

		double TotalMass = mTotalDomainInitialSize * GetProperties()[DENSITY];;
		if (TotalMass<=0)
		{
			std::cout<<"??? INVERTED ELEMENT AT THE BEGINNING OF T-STEP!!"<<std::endl;
			std::cout<<"??? INVERTED ELEMENT AT THE BEGINNING OF T-STEP!!"<<std::endl;
			TotalMass=0;
		}

		//if(dimension == 2) TotalMass *= GetProperties()[THICKNESS];
		if(dimension == 2) TotalMass *= 1;
		//KRATOS_WATCH(mTotalDomainInitialSize)
		//KRATOS_WATCH(GetProperties()[DENSITY])
		//KRATOS_WATCH(GetProperties()[THICKNESS])
		Vector LumpFact;
		LumpFact = GetGeometry().LumpingFactors(LumpFact);

		for(unsigned int i=0; i<NumberOfNodes; i++)
		{
			double temp = LumpFact[i]*TotalMass;
			for(unsigned int j=0; j<dimension; j++)
			{
				unsigned int index = i*dimension + j;
				rMassMatrix(index,index) = temp;
			}
		}
		//std::cout<<" MASS MAT "<<rMassMatrix;
//KRATOS_WATCH(rMassMatrix);
		KRATOS_CATCH("")
	}

	//************************************************************************************
	//************************************************************************************
	void UpdatedLagrangianFluid::DampMatrix(MatrixType& rDampMatrix, ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY
			rDampMatrix.resize(0,0);
		KRATOS_CATCH("")
	}

	//************************************************************************************
	//************************************************************************************
	void UpdatedLagrangianFluid::CalculateOnIntegrationPoints(const Variable<Matrix >& rVariable, std::vector< Matrix >& Output, const ProcessInfo& rCurrentProcessInfo)
	{
		//unsigned int number_of_nodes = GetGeometry().size();
		unsigned int dim = GetGeometry().WorkingSpaceDimension();

		//resizing as needed the LHS
		//unsigned int MatSize=number_of_nodes*dim;

		//reading integration points and local gradients
		const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints();
		//const Matrix& Ncontainer = mpReferenceGeometry->ShapeFunctionsValues();

		//calculating actual jacobian
		GeometryType::JacobiansType J;
		J = GetGeometry().Jacobian(J);  

		//auxiliary terms
		Matrix F(dim,dim);
		Matrix C(dim,dim);
		Vector StrainVector;
		Vector StressVector;

		if(Output.size() != integration_points.size())
			Output.resize(integration_points.size());

		for(unsigned int PointNumber = 0; PointNumber<integration_points.size(); PointNumber++)
		{
			//deformation gradient
			noalias(F) = prod(J[PointNumber],mInvJ0[PointNumber]);

			//strain calculation
			noalias(C) = prod(trans(F),F);
			CalculateStrain(C,StrainVector);

			if(rVariable==GREEN_LAGRANGE_STRAIN_TENSOR)
			{
				if(Output[PointNumber].size2() != StrainVector.size())
					Output[PointNumber].resize(1,StrainVector.size());
				for(unsigned int ii = 0; ii<StrainVector.size(); ii++)
					Output[PointNumber](0,ii) = StrainVector[ii];
			}
			else if(rVariable==PK2_STRESS_TENSOR)
			{
				if(Output[PointNumber].size2() != StrainVector.size())
					Output[PointNumber].resize(1,StrainVector.size());
			//	mConstitutiveLawVector(PointNumber)->UpdateMaterial(StrainVector,rCurrentProcessInfo);
			//	StressVector = mConstitutiveLawVector[PointNumber]->CalculateStress(StrainVector); //saving the stress vector

				for(unsigned int ii = 0; ii<StrainVector.size(); ii++)
					Output[PointNumber](0,ii) = StressVector[ii];
			}
		}

	}
	//************************************************************************************
	//************************************************************************************

void UpdatedLagrangianFluid::Calculate(const Variable< array_1d<double,3> >& rVariable, array_1d<double,3>& Output, const ProcessInfo& rCurrentProcessInfo)
{
	KRATOS_TRY

		unsigned int number_of_nodes = GetGeometry().size();
	unsigned int dim = GetGeometry().WorkingSpaceDimension();

	//resizing as needed the LHS
	//e.g. in 2D and using linear triangles => MatSize = 3*2=6
	unsigned int MatSize=number_of_nodes*dim;

	Vector aaa(MatSize);


	//reading integration points and local gradients
	const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints();
	const GeometryType::ShapeFunctionsGradientsType& DN_De = GetGeometry().ShapeFunctionsLocalGradients();
	//const Matrix& Ncontainer = GetGeometry().ShapeFunctionsValues();

	//calculating actual jacobian
	GeometryType::JacobiansType J;
	GetGeometry().Jacobian(J);  

	for(unsigned int PointNumber = 0; PointNumber<integration_points.size(); PointNumber++)
	{

	//Calculating the cartesian derivatives (it is avoided storing them to minimize storage)
		msDN_DX = ZeroMatrix(msDN_DX.size1(),msDN_DX.size2());

		noalias(msDN_DX) = prod(DN_De[PointNumber],mInvJ0[PointNumber]);


		//deformation gradient
		noalias(msF) = prod(J[PointNumber],mInvJ0[PointNumber]);

		//strain calculation
		noalias(msC) = prod(trans(msF),msF);
		CalculateStrain(msC,msStrainVector);

		//PAVEL
		//Obtain PK2 stress vector - first parameter is the bulk modulus K
		double K = GetProperties()[BULK_MODULUS];
		K = 100000.0;

		noalias(msStressVector) = ZeroVector(msStrainVector.size());
		double vol_var  =0.00;
		for(unsigned int i = 0; i< dim; i++)
		{
			vol_var += msStrainVector[i];
		}

		vol_var /= double(dim);
		for(unsigned int i = 0; i< dim; i++)
		{
			msStressVector[i] = vol_var * K;
		}

		//calculating operator B
		CalculateB(msB,msF,msDN_DX,msStrainVector.size());
		//calculating weights for integration on the reference configuration
		double IntToReferenceWeight = integration_points[PointNumber].Weight() * mDetJ0[PointNumber];
		if (dim == 2) IntToReferenceWeight *= 1;

		//contribution to external forces 
		// operation performed: rRightHandSideVector += ExtForce*IntToReferenceWeight
		noalias(aaa) = IntToReferenceWeight * prod(trans(msB),msStressVector);

		//CHAPUZA: checking if the element was inverted within the nonlin iter
		double dd=0.0;
		Matrix Jinvv (2,2);
		
		MathUtils<double>::InvertMatrix(J[0],Jinvv,dd);				
		
		//and now we store the pressure forces on the nodes
		for(unsigned int inode = 0; inode<GetGeometry().size(); inode++)
		{
			array_1d<double,3>& tmp = GetGeometry()[inode].FastGetSolutionStepValue(PRESSURE_FORCE);
			for(unsigned k=0; k<dim; k++)
			{
				tmp[k] += aaa[inode*dim + k];
				if (dd<=0.0)
				{
					std::cout<<"inside calculate - inv element"<<std::endl;
					tmp[k]=0.0;
				}
				if (dd==0.0 || dd==0)
				{
					std::cout<<"THREE NODES OF ELEMENT ON ONE LINE"<<std::endl;
					std::cout<<"THREE NODES OF ELEMENT ON ONE LINE"<<std::endl;
					std::cout<<"TIME IS"<<rCurrentProcessInfo[TIME]<<std::endl;
				}
			}
		}

	}
	KRATOS_CATCH("")
}

	//************************************************************************************
	//************************************************************************************
	  void UpdatedLagrangianFluid::GetValuesVector(Vector& values, int Step)

	{
		unsigned int number_of_nodes = GetGeometry().size();
		unsigned int dim = GetGeometry().WorkingSpaceDimension();
		unsigned int MatSize = number_of_nodes*dim;
		if(values.size() != MatSize)	values.resize(MatSize);
		for (unsigned int i=0;i<number_of_nodes;i++)
		{
			unsigned int index = i*dim;
			values[index] = GetGeometry()[i].GetSolutionStepValue(DISPLACEMENT_X,Step);
			values[index + 1] = GetGeometry()[i].GetSolutionStepValue(DISPLACEMENT_Y,Step);
			if(dim == 3)
				values[index + 2] = GetGeometry()[i].GetSolutionStepValue(DISPLACEMENT_Z,Step);
		}
//std::cout << "GetValuesVector " << Id() << " step = " << Step << std::endl;
//KRATOS_WATCH(Step)
//KRATOS_WATCH(values)
	}
	//************************************************************************************
	//************************************************************************************
	  void UpdatedLagrangianFluid::GetFirstDerivativesVector(Vector& values, int Step)
	{
		unsigned int number_of_nodes = GetGeometry().size();
		unsigned int dim = GetGeometry().WorkingSpaceDimension();
		unsigned int MatSize = number_of_nodes*dim;
		if(values.size() != MatSize)   values.resize(MatSize);
		for (unsigned int i=0;i<number_of_nodes;i++)
		{
			unsigned int index = i*dim;
			values[index] = GetGeometry()[i].GetSolutionStepValue(VELOCITY_X,Step);
			values[index + 1] = GetGeometry()[i].GetSolutionStepValue(VELOCITY_Y,Step);
			if(dim == 3)
				values[index + 2] = GetGeometry()[i].GetSolutionStepValue(VELOCITY_Z,Step);
		}
//std::cout << "GetFirstDerivativesVector " << Id() << " step = " << Step << std::endl;
//KRATOS_WATCH(Step)
//KRATOS_WATCH(values)

	}
	//************************************************************************************
	//************************************************************************************

	
	
	
	
	void UpdatedLagrangianFluid::GetSecondDerivativesVector(Vector& values, int Step)
	{
		unsigned int number_of_nodes = GetGeometry().size();
		unsigned int dim = GetGeometry().WorkingSpaceDimension();
		unsigned int MatSize = number_of_nodes*dim;
		if(values.size() != MatSize) values.resize(MatSize);
		for (unsigned int i=0;i<number_of_nodes;i++)
		{
			unsigned int index = i*dim;
			values[index] = GetGeometry()[i].GetSolutionStepValue(ACCELERATION_X,Step);
			values[index + 1] = GetGeometry()[i].GetSolutionStepValue(ACCELERATION_Y,Step);
			if(dim == 3)
				values[index + 2] = GetGeometry()[i].GetSolutionStepValue(ACCELERATION_Z,Step);
		}
//std::cout << "GetSecondDerivativesVector " << Id() << " step = " << Step << std::endl;
//KRATOS_WATCH(Step)
//KRATOS_WATCH(values)
	}
} // Namespace Kratos


