//
//   Project Name:        Kratos       
//   Last modified by:    $Author: kazem $
//   Date:                $Date: 2009-01-15 18:49:25 $
//   Revision:            $Revision: 1.3 $



// System includes 
// External includes 
// Project includes 

#include "includes/define.h"
#include "custom_elements/Hypoelastic_element.h"
#include "utilities/math_utils.h"
//#include "constitutive_laws/isotropic_planestress_wrinkling_new.h"
#include "constitutive_laws/isotropic_2d.h"
#include "constitutive_laws/isotropic_3d.h"
#include "structural_application.h"



namespace Kratos

{

    Matrix HypoelasticElement::msB(0,0);
	Matrix HypoelasticElement::msD(0,0);
	Matrix HypoelasticElement::msDN_DX(0,0);
	Vector HypoelasticElement::msStrainVector(0);
	Vector HypoelasticElement::msStressVector(0);
	Matrix HypoelasticElement::msF(0,0);

	//************************************************************************************

	//************************************************************************************

	HypoelasticElement::HypoelasticElement(IndexType NewId, GeometryType::Pointer pGeometry)
		: Element(NewId, pGeometry)
	{		
		//DO NOT ADD DOFS HERE!!!
	}

	//************************************************************************************
	//************************************************************************************

	HypoelasticElement::HypoelasticElement(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
		: Element(NewId, pGeometry, pProperties)
	{

		//const unsigned int dim = GetGeometry().WorkingSpaceDimension();

		//setting up the nodal degrees of freedom
/*		for(unsigned int i = 0 ; i != GetGeometry().size() ; ++i)
		{
			(GetGeometry()[i].pAddDof(DISPLACEMENT_X,REACTION_X));
			(GetGeometry()[i].pAddDof(DISPLACEMENT_Y,REACTION_Y));
			if(dim == 3)
				(GetGeometry()[i].pAddDof(DISPLACEMENT_Z,REACTION_Z));
		}*/
		//initializing static variables
		unsigned int dimension = GetGeometry().WorkingSpaceDimension();
		unsigned int dim2 = GetGeometry().size()*dimension;
		unsigned int StrainSize;
		if (dimension==2) 
			StrainSize = 3;
		else 
			StrainSize = 6;
		msB.resize(StrainSize,dim2);
		msD.resize(StrainSize,StrainSize);
		msDN_DX.resize(GetGeometry().size(),dimension);
		msStrainVector.resize(StrainSize);
		msStressVector.resize(StrainSize);
		msF.resize(dimension,dimension);
		((*this).GetValue(PK2_STRESS_TENSOR)).resize(1,3);
		((*this).GetValue(CAUCHY_STRESS_TENSOR)).resize(1,3);

		noalias( (*this).GetValue(PK2_STRESS_TENSOR)) = ZeroMatrix(1,3);
		noalias( (*this).GetValue(CAUCHY_STRESS_TENSOR)) = ZeroMatrix(1,3);

	}
	Element::Pointer HypoelasticElement::Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const
	{

		return Element::Pointer(new HypoelasticElement(NewId, GetGeometry().Create(ThisNodes), pProperties));

	}
	HypoelasticElement::~HypoelasticElement()
	{
	}

	//************************************************************************************
	//************************************************************************************

	void HypoelasticElement::Initialize()
	{
		KRATOS_TRY
			const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints();
		//resizing jacobian inverses containers
		mInvJ0.resize(integration_points.size());
		mDetJ0.resize(integration_points.size());
		mInvJ.resize(integration_points.size());
		mDetJ.resize(integration_points.size());

		GeometryType::JacobiansType J0;
		J0 = GetGeometry().Jacobian(J0);  
		mTotalDomainInitialSize = 0.00;

		//Constitutive Law initialisation
		if(mConstitutiveLawVector.size() == 0)
		{
			mConstitutiveLawVector.resize(integration_points.size());
			InitializeMaterial();
		}

		//calculating the inverse J0
		for(unsigned int PointNumber = 0; PointNumber<integration_points.size(); PointNumber++)
		{
			//getting informations for integration
			double IntegrationWeight = integration_points[PointNumber].Weight();
			//calculating and storing inverse of the jacobian and the parameters needed

			MathUtils<double>::InvertMatrix(J0[PointNumber],mInvJ0[PointNumber],mDetJ0[PointNumber]);				
			//KRATOS_WATCH("inside element initialize");
			//KRATOS_WATCH(mInvJ0[0]);
	
			//KRATOS_WATCH(J0[PointNumber]);
			//KRATOS_WATCH(mInvJ0[0]);

			//calculating the total area

			mTotalDomainInitialSize += mDetJ0[PointNumber]*IntegrationWeight;

		}



		KRATOS_CATCH("")

	}



	//************************************************************************************
	//************************************************************************************
	void HypoelasticElement::InitializeSolutionStep(ProcessInfo& CurrentProcessInfo)
		{
		
		for(unsigned int i=0; i< mConstitutiveLawVector.size(); ++i)
			mConstitutiveLawVector[i]->InitializeSolutionStep(GetProperties(), GetGeometry(),	row(GetGeometry().ShapeFunctionsValues(), i),CurrentProcessInfo);


		}

	//************************************************************************************
	//************************************************************************************


	void HypoelasticElement::CalculateAll(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, 
										ProcessInfo& rCurrentProcessInfo,
										bool CalculateStiffnessMatrixFlag,
										bool CalculateResidualVectorFlag)
	{

		KRATOS_TRY

		const unsigned int number_of_nodes = GetGeometry().size();
		const unsigned int dim = GetGeometry().WorkingSpaceDimension();

		

		//resizing as needed the LHS

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
		
		//stress calls
		Matrix vectorcauchyN;
		vectorcauchyN = (*this).GetValue(PK2_STRESS_TENSOR);
		Vector GLstrain;
		//auxiliary terms

		Vector BodyForce;		

		for(unsigned int PointNumber = 0; PointNumber<integration_points.size(); PointNumber++)
		{
			//Calculating the cartesian derivatives (it is avoided storing them to minimize storage)
			MathUtils<double>::InvertMatrix(J[PointNumber],mInvJ[PointNumber],mDetJ[PointNumber]);	
			msDN_DX = ZeroMatrix(msDN_DX.size1(),msDN_DX.size2());
			noalias(msDN_DX) = prod(DN_De[PointNumber],mInvJ[PointNumber]);
			

			
			//material update (considering the level of strain achieved)

			//Calculation of constitutive
			//Hypoconstitutive(msD);

			//calculating operator B
			CalculateB(msB,msDN_DX);
			
			//calculate strain using B
			//CalculateStrain(msB,msStrainVector);
			//
				//CALCULATE J
			
			noalias(msF) = prod(J[PointNumber],mInvJ0[PointNumber]);

			//KRATOS_WATCH("FOR CALCULATEALL");
			//KRATOS_WATCH(msF);
			//CalculateStress(msStrainVector,msStressVector,msF,msD, rCurrentProcessInfo);
// 			mConstitutiveLawVector[PointNumber]->CalculateCauchyStresses(msStressVector, msF, row(vectorcauchyN,PointNumber), GLstrain);
            Matrix dummy(0,0);
            Vector cauchy_stress = row(vectorcauchyN,PointNumber);
            mConstitutiveLawVector[PointNumber]->CalculateMaterialResponse(
                    GLstrain,
                    msF,
                    cauchy_stress,
                    dummy,
                    rCurrentProcessInfo,
                    GetProperties(),
                    GetGeometry(),
                    row(Ncontainer,PointNumber),
                    true,
                    0,//no material matrix is computed
                    true );
            noalias(row(vectorcauchyN,PointNumber)) = cauchy_stress;

			//CalculateRotatedStress(msStressVector,msF, rCurrentProcessInfo);
			
			//calculating weights for integration on the reference configuration
			
			double IntToReferenceWeight = integration_points[PointNumber].Weight() * mDetJ[PointNumber];
			

			
			if (CalculateResidualVectorFlag == true) //calculation of the matrix is required

			{
				//contribution to external forces 
				
				BodyForce=ZeroVector(2);
				BodyForce[0]=0.0;
				BodyForce[1]= -10.0*GetProperties()[DENSITY];
			
				
				
				
				// operation performed: rRightHandSideVector += ExtForce*IntToReferenceWeight
				CalculateAndAdd_ExtForceContribution(row(Ncontainer,PointNumber),rCurrentProcessInfo,BodyForce,rRightHandSideVector,IntToReferenceWeight);
				// operation performed: rRightHandSideVector -= IntForce*IntToReferenceWeight
				//msStressVector = prod(msD, msStrainVector);
				
				noalias(rRightHandSideVector) -= IntToReferenceWeight * prod(trans(msB),msStressVector);
				//KRATOS_WATCH(prod(trans(msB),msStressVector));
			}	

		}
		//KRATOS_WATCH(rRightHandSideVector);

		KRATOS_CATCH("")

	}

	//************************************************************************************

	//************************************************************************************

	void HypoelasticElement::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)

	{

		//calculation flags

		bool CalculateStiffnessMatrixFlag = false;
		bool CalculateResidualVectorFlag = true;
		MatrixType temp = Matrix();
		CalculateAll(temp, rRightHandSideVector, rCurrentProcessInfo, CalculateStiffnessMatrixFlag, CalculateResidualVectorFlag);

	}



	//************************************************************************************

	//************************************************************************************

	void HypoelasticElement::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)

	{
		//calculation flags
		bool CalculateStiffnessMatrixFlag = true;
		bool CalculateResidualVectorFlag = true;
		
		CalculateAll(rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo, CalculateStiffnessMatrixFlag, CalculateResidualVectorFlag);

	}
	//**************************************************************************
	//**************************************************************************

	double HypoelasticElement::CalculateIntegrationWeight(double GaussPointWeight, double DetJ0)
	{
		//to permorm the integration over the reference domain we need to include 
		// the thickness in 2D
		unsigned int dimension = GetGeometry().WorkingSpaceDimension();
		double weight = GaussPointWeight;
	  
		weight *= DetJ0;
		if (dimension == 2) weight *= GetProperties()[THICKNESS];
		return weight;
	}



	////**********************************************************************************
	//************************************************************************************

	void HypoelasticElement::InitializeMaterial()

	{
		KRATOS_TRY

		for (unsigned int i=0; i<mConstitutiveLawVector.size(); i++)

		{
			//temporary

			mConstitutiveLawVector[i] = GetProperties()[CONSTITUTIVE_LAW];
			if( mConstitutiveLawVector[i]->GetStressMeasure() != ConstitutiveLaw::StressMeasure_Cauchy )
                KRATOS_ERROR( std::logic_error, "The specified constitutive law does not support Cauchy stresses, please choose a suitable constitutive model", "" );
            mConstitutiveLawVector[i]->InitializeMaterial( 	GetProperties(), GetGeometry(),	row(GetGeometry().ShapeFunctionsValues(), i) ); 

		}

		KRATOS_CATCH("")

	} 



	//************************************************************************************

	//************************************************************************************

	inline void HypoelasticElement::CalculateAndAdd_ExtForceContribution(

		const Vector& N,
		const ProcessInfo& CurrentProcessInfo,
		Vector& BodyForce,
		VectorType& rRightHandSideVector,
		double weight
		)

	{

		KRATOS_TRY
		
			unsigned int number_of_nodes = GetGeometry().PointsNumber();

		unsigned int dimension = GetGeometry().WorkingSpaceDimension();
		for (unsigned int i=0;i<number_of_nodes;i++)
		{
			int index = dimension*i;
			for (unsigned int j=0; j<dimension; j++)  rRightHandSideVector[index+j] += weight*N[i]*BodyForce[j];		  
		}

		KRATOS_CATCH("")

	}



	//************************************************************************************

	//************************************************************************************

	void HypoelasticElement::CalculateStrain(
		const Matrix& B,
		Vector& StrainVector)

	{

		KRATOS_TRY

			int dimension = GetGeometry().WorkingSpaceDimension();
			int number_of_nodes = GetGeometry().size();
			Vector Velocity(dimension*number_of_nodes);
			Velocity=ZeroVector(dimension*number_of_nodes);
	
			for (int i=0; i<number_of_nodes; i++)
				{
					int index = dimension*i;
					array_1d<double,3>& vel = GetGeometry()[i].FastGetSolutionStepValue(VELOCITY);
			//KRATOS_WATCH(disp);
					Velocity[index+0] = vel[0];
					Velocity[index+1] = vel[1];
				//KRATOS_WATCH(GetGeometry()[i].GetDof(DISPLACEMENT_X).EquationId());
				//KRATOS_WATCH(GetGeometry()[i].Id());											
				if (dimension==3)
					Velocity[index+2] = vel[2];
				
				}	
				
				//KRATOS_WATCH(Velocity);
				//KRATOS_WATCH(B);
		noalias(StrainVector) = prod(B,Velocity);
			
		//KRATOS_WATCH(prod(B,Velocity));
		//KRATOS_WATCH(StrainVector);
		KRATOS_CATCH("")
		

	}

	//************************************************************************************

	//************************************************************************************

	void HypoelasticElement::CalculateB(

		Matrix& B,
		Matrix& DN_DX)

	{

		KRATOS_TRY

			const unsigned int number_of_nodes = GetGeometry().PointsNumber();
		unsigned int dimension = GetGeometry().WorkingSpaceDimension();
		
		for (unsigned int i=0;i<number_of_nodes;i++)
		{
			unsigned int index = dimension*i;
			if (dimension == 2)
			{
				B(0,index+0)=DN_DX(i,0);				B(0,index+1)=0;
				B(1,index+0)=0;						B(1,index+1)=DN_DX(i,1);
				B(2,index+0)=DN_DX(i,1);				B(2,index+1)=DN_DX(i,0); 
			}
			else
			{
		           B(0,index+0)=DN_DX(i,0);                B(0,index+1)=0;		B(0,index+2)=0;
			   B(1,index+0)=0;			   B(1,index+1)=DN_DX(i,1);	B(1,index+2)=0;
			   B(2,index+0)=0;			   B(2,index+1)=0;		B(2,index+2)=DN_DX(i,2);
			B(3,index+0)=DN_DX(i,1);		B(3,index+1)=DN_DX(i,0);	B(3,index+2)=0;
			B(4,index+0)=0;				B(4,index+1)=DN_DX(i,2);	B(4,index+2)=DN_DX(i,1);
			B(5,index+0)=DN_DX(i,2);		B(5,index+1)=0;			B(5,index+2)=DN_DX(i,0);

			}

		}

		KRATOS_CATCH("")

	}
	//********************************************************************************
	//********************************************************************************
// 	void HypoelasticElement::CalculateBPV(
// 
// 		Matrix& BPV,
// 		const Vector& N,
// 		Matrix& DN_DX)
// 
// 	{
// 
// 		KRATOS_TRY
// 
// 		const unsigned int number_of_nodes = GetGeometry().PointsNumber();
// 		unsigned int dimension = GetGeometry().WorkingSpaceDimension();
// 		
// 		for (unsigned int i=0;i<number_of_nodes;i++)
// 		{
// 			unsigned int index = dimension*i;
// 			if (dimension == 2)
// 			{
// 		BPV(index+0,0)=N[0]*DN_DX(i,0);	BPV(index+0,1)=N[1]*DN_DX(i,0); BPV(index+0,2)=N[2]*DN_DX(i,0);
// 		BPV(index+1,0)=N[0]*DN_DX(i,1);	BPV(index+1,1)=N[1]*DN_DX(i,1); BPV(index+1,2)=N[2]*DN_DX(i,1);					
// 			}
// 			else
// 			{
// 		BPV(index+0,0)=N[0]*DN_DX(i,0);	BPV(index+0,1)=N[1]*DN_DX(i,0); BPV(index+0,2)=N[2]*DN_DX(i,0);BPV(index+0,3)=N[3]*DN_DX(i,0);
// 
// 		BPV(index+1,0)=N[0]*DN_DX(i,1);	BPV(index+1,1)=N[1]*DN_DX(i,1); BPV(index+1,2)=N[2]*DN_DX(i,1);BPV(index+1,3)=N[3]*DN_DX(i,1);
// 
// 			}
// 
// 		}
// 
// 		KRATOS_CATCH("")
// 
// 	}
	//*********************************************************************************
	//*********************************************************************************

	void HypoelasticElement::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo)

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
		//KRATOS_WATCH("FIRST MEMBER OF EQUATIONID");
		//KRATOS_WATCH(rResult[0]);
		
	}

	//************************************************************************************
	//************************************************************************************

	void HypoelasticElement::GetDofList(DofsVectorType& ElementalDofList,ProcessInfo& CurrentProcessInfo)

	{
		ElementalDofList.resize(0);
		for (unsigned int i=0;i<GetGeometry().size();i++)
		{
			
			ElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_X));
			ElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_Y));
			if(GetGeometry().WorkingSpaceDimension() == 3)
				ElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_Z));

		}

	}



	//************************************************************************************
	//************************************************************************************

	void HypoelasticElement::Calculate(const Variable<Matrix >& rVariable, Matrix& Output, const ProcessInfo& rCurrentProcessInfo)
	{

		KRATOS_TRY
		//KRATOS_WATCH(*this);
		//KRATOS_WATCH(Output);
		//KRATOS_WATCH("inside calculate");
		//KRATOS_WATCH(msStressVector);
		
		//ResizeAndInitializeAuxiliaries();
		const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints();
		const GeometryType::ShapeFunctionsGradientsType& DN_De = GetGeometry().ShapeFunctionsLocalGradients();
		
		//KRATOS_WATCH(integration_points.size());
		//calculating actual jacobian
		GeometryType::JacobiansType J;
		J = GetGeometry().Jacobian(J);  

		mInvJ.resize(integration_points.size());
		mDetJ.resize(integration_points.size());

		
		
		//Output.resize(integration_points.size());
		
		//KRATOS_WATCH(J[0]);
		
		MathUtils<double>::InvertMatrix(J[0],mInvJ[0],mDetJ[0]);
		//KRATOS_WATCH(mInvJ[0]);
		//KRATOS_WATCH(mDetJ[0]);
		noalias(msF) = prod(J[0],mInvJ0[0]);

		msDN_DX = ZeroMatrix(msDN_DX.size1(),msDN_DX.size2());
		noalias(msDN_DX) = prod(DN_De[0],mInvJ[0]);
		
// 		Matrix StressTensor;
// 		Matrix StrainTensor;
// 		mConstitutiveLawVector[0]->CalculateStressAndTangentMatrix(StressTensor, msF, StrainTensor, msD);
        mConstitutiveLawVector[0]->CalculateMaterialResponse(
                    msStrainVector,
                    msF,
                    msStressVector,
                    msD,
                    rCurrentProcessInfo,
                    GetProperties(),
                    GetGeometry(),
                    row( GetGeometry().ShapeFunctionsValues(), 0 ),
                    true,
                    1,
                    true );
        
		//Hypoconstitutive(msD,msF,rCurrentProcessInfo);
		CalculateB(msB,msDN_DX);
		
		CalculateStrain(msB,msStrainVector);
		
		
		//KRATOS_WATCH(mInvJ0[0]);
		
		
		
		//KRATOS_WATCH(msF);
		
		CalculateStress(msStrainVector,msStressVector,msF,msD, rCurrentProcessInfo);
		
		//KRATOS_WATCH("inside caalculate gg");
		//KRATOS_WATCH(Output);
		//KRATOS_WATCH(msStressVector);
		if(rVariable == CAUCHY_STRESS_TENSOR)
		 {
			if (Output.size2()!=msStressVector.size())
				Output.resize(1,msStressVector.size(),false);
			//KRATOS_WATCH("inside CAUCHY loop");
			//KRATOS_WATCH(msStressVector.size());
			for(unsigned int ii=0; ii<msStressVector.size(); ii++)
				Output(0,ii) = msStressVector[ii];
				
		 }
		//KRATOS_WATCH(Output);
		//KRATOS_WATCH("end of Calculate");
		KRATOS_CATCH("")



	}

	//************************************************************************************
	//************************************************************************************
	void HypoelasticElement::CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo)
		{
		
		bool CalculateStiffnessMatrixFlag = true;
		bool CalculateResidualVectorFlag = false;
		Vector temp = Vector();
		CalculateAll(rLeftHandSideMatrix, temp, rCurrentProcessInfo, CalculateStiffnessMatrixFlag, CalculateResidualVectorFlag);
		
		}
	//************************************************************************************
	//************************************************************************************
	void HypoelasticElement::CalculateRotatedStress(Vector& StressVector,Matrix& F,const ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY
		
		const unsigned int dim = GetGeometry().WorkingSpaceDimension();
		double detF = MathUtils<double>::Det(F);
		
		Matrix cauchyN;
		Matrix vectorcauchyN;
		Matrix transformed_cauchy;
		Matrix tempstr;
		
		StressVector=ZeroVector(3);
		
		cauchyN.resize(dim,dim,false);
		//KRATOS_WATCH(StrainVector);
		//KRATOS_WATCH(D);
		//KRATOS_WATCH(Hypofac);
			
		
		//KRATOS_WATCH(StressVector);
		
		vectorcauchyN = (*this).GetValue(PK2_STRESS_TENSOR);
		cauchyN(0,0)=vectorcauchyN(0,0);
		cauchyN(1,1)=vectorcauchyN(0,1);
		cauchyN(0,1)=vectorcauchyN(0,2);
		cauchyN(1,0)=vectorcauchyN(0,2);

		
		//KRATOS_WATCH(vectorcauchyN);
		//KRATOS_WATCH(cauchyN);
		tempstr = prod(cauchyN,trans(F));

		
		//KRATOS_WATCH(cauchyN);
		transformed_cauchy = prod(F,tempstr);
		

		StressVector[0]= transformed_cauchy(0,0)/detF;
		StressVector[1]= transformed_cauchy(1,1)/detF;
		StressVector[2]= transformed_cauchy(0,1)/detF;	
		//KRATOS_WATCH("IN ROTATED STRESS");
		//KRATOS_WATCH(StressVector);
		KRATOS_CATCH("")
	}
	//************************************************************************************
	//************************************************************************************
	void HypoelasticElement::CalculateStress(Vector& StrainVector, Vector& StressVector, Matrix& F, Matrix& D, const ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY
		//double deltaT= rCurrentProcessInfo[DELTA_TIME];
		const unsigned int dim = GetGeometry().WorkingSpaceDimension();
		//double detF = MathUtils<double>::Det(F);
		
		Matrix cauchyN;
		Vector vectorcauchyN;
		Vector transformed_cauchy;
		Vector GLstrain;
		Matrix tempstr;
		double Hypofac = 1.00;
		
		cauchyN.resize(dim,dim,false);
		StressVector=transformed_cauchy=ZeroVector(3);
		//KRATOS_WATCH(D);
		//KRATOS_WATCH(Hypofac);
		//KRATOS_WATCH(StressVector);	
		StressVector[0] = Hypofac*(D(0,0)*StrainVector[0] + D(0,1) * (StrainVector[1]))	;
		StressVector[1] = Hypofac*(D(1,0)*StrainVector[0] + D(1,1) * (StrainVector[1]))	;
		StressVector[2] = Hypofac*D(2,2)*StrainVector[2];
		//KRATOS_WATCH("IT IS SUPPOSED TO BE ZERO");
		//KRATOS_WATCH(StrainVector);
		//KRATOS_WATCH(StressVector);
		
		vectorcauchyN = row((*this).GetValue(PK2_STRESS_TENSOR),0);
		//KRATOS_WATCH("IT IS SUPPOSED TO BE (1000,0,0)");
		//KRATOS_WATCH(vectorcauchyN);
        mConstitutiveLawVector[0]->CalculateMaterialResponse(
                    GLstrain,
                    F,
                    vectorcauchyN,
                    tempstr,
                    rCurrentProcessInfo,
                    GetProperties(),
                    GetGeometry(),
                    row(GetGeometry().ShapeFunctionsValues(),0),
                    true,
                    0,
                    true );

// 		mConstitutiveLawVector[0]->CalculateCauchyStresses(transformed_cauchy, F, row(vectorcauchyN,0), GLstrain);

/*
		cauchyN(0,0)=vectorcauchyN(0,0);
		cauchyN(1,1)=vectorcauchyN(0,1);
		cauchyN(0,1)=vectorcauchyN(0,2);
		cauchyN(1,0)=vectorcauchyN(0,2);
		
		//KRATOS_WATCH(cauchyN);
		tempstr = prod(cauchyN,trans(F));

		//KRATOS_WATCH(F);
		//KRATOS_WATCH(cauchyN);
		transformed_cauchy = prod(F,tempstr);
		noalias(transformed_cauchy)/=detF;*/
		//KRATOS_WATCH("IT IS SUPPOSED TO BE (0,1000,0)");
		//KRATOS_WATCH(transformed_cauchy);
		noalias(StressVector)+=  transformed_cauchy;
// 		StressVector[0]=StressVector[0] + transformed_cauchy(0,0);
// 		StressVector[1]=StressVector[1] + transformed_cauchy(1,1);
// 		StressVector[2]=StressVector[2] + transformed_cauchy(0,1);
		//KRATOS_WATCH("IT IS SUPPOSED TO REMAIN (0,1000,0)");
		//KRATOS_WATCH(detF);
		KRATOS_CATCH("")
	}
	//************************************************************************************
	//************************************************************************************	
	void HypoelasticElement::Hypoconstitutive(Matrix& D, Matrix& F,const ProcessInfo& rCurrentProcessInfo)
	{
		
	double miu = GetProperties()[MIU];
	double lambda = GetProperties()[LAMBDA];
	double deltaT= rCurrentProcessInfo[DELTA_TIME];
	double detF = MathUtils<double>::Det(F);
	//KRATOS_WATCH(deltaT);
	miu = miu*deltaT/detF;
	lambda= lambda*deltaT/detF;
	
/*	double miu = 8.1e9;
	double lambda = 7.1e9;*/
	//KRATOS_WATCH(detF);
	double Ehypo = miu*(3*lambda+2*miu)/(lambda+miu);
	double NUhypo = lambda/(2*(lambda+miu));
	double c1 = Ehypo*(1.00-NUhypo) / ((1.00 + NUhypo)*(1.00 - 2*NUhypo));
	double c2 = c1 * NUhypo/(1.00-NUhypo);
	double c3 = 0.5*Ehypo / (1 + NUhypo);
	
	D(0,0) = c1;	D(0,1) = c2;	D(0,2) = 0.0;
	D(1,0) = c2;	D(1,1) = c1;	D(1,2) = 0.0;
	D(2,0) = 0.0;	D(2,1) = 0.0;	D(2,2) = c3;
	//KRATOS_WATCH(D);
	}
	//************************************************************************************
	//************************************************************************************
	void HypoelasticElement::FinalizeSolutionStep(ProcessInfo& CurrentProcessInfo)
	{
		KRATOS_TRY
//         std::cout << "in TL: calling FinalizeSolutionStep" << std::endl;
		for(unsigned int i = 0; i < mConstitutiveLawVector.size(); i++)
			mConstitutiveLawVector[i]->FinalizeSolutionStep( GetProperties(),
				GetGeometry(),
                            row( GetGeometry().ShapeFunctionsValues(), i ),
				CurrentProcessInfo);	
		
		//Updating stress
		Matrix temp;	
		int dim=GetGeometry().WorkingSpaceDimension();
		temp.resize(1,(dim*dim+dim)/2,false);
		temp=ZeroMatrix(1,(dim*dim+dim)/2);
				
		Calculate(CAUCHY_STRESS_TENSOR,temp,CurrentProcessInfo);
		//KRATOS_WATCH("IN FINALIZE AFTER CALCULATE TEMP SUPPOSED TO BE (0,1000,0)");
		//KRATOS_WATCH(temp);
		SetValue(PK2_STRESS_TENSOR,temp);
		//KRATOS_WATCH((*this).GetValue(PK2_STRESS_TENSOR));
	
		const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints();
		//resizing jacobian inverses containers
		mInvJ0.resize(integration_points.size());
		mDetJ0.resize(integration_points.size());

		GeometryType::JacobiansType J0;
		J0 = GetGeometry().Jacobian(J0);  
		//KRATOS_WATCH("IN FINALIZE AFTER CALCULATE J0[0] SUPPOSED TO BE ((0,-2),(2,2))");
		//KRATOS_WATCH(J0[0]);

		//mTotalDomainInitialSize = 0.00;


		//calculating the inverse J0
		for(unsigned int PointNumber = 0; PointNumber<integration_points.size(); PointNumber++)
		{
			//getting informations for integration
			//double IntegrationWeight = integration_points[PointNumber].Weight();
			//calculating and storing inverse of the jacobian and the parameters needed

			MathUtils<double>::InvertMatrix(J0[PointNumber],mInvJ0[PointNumber],mDetJ0[PointNumber]);				



			//calculating the total area

			//mTotalDomainInitialSize += mDetJ0[PointNumber]*IntegrationWeight;

		}





		KRATOS_CATCH("")

	}

	//*****************************************************************************
	//*****************************************************************************
	  void HypoelasticElement::GetSecondDerivativesVector(Vector& values, int Step)
	{
		const unsigned int number_of_nodes = GetGeometry().size();
		const unsigned int dim = GetGeometry().WorkingSpaceDimension();
		unsigned int MatSize = number_of_nodes*dim;
		if(values.size() != MatSize) values.resize(MatSize,false);
		for (unsigned int i=0;i<number_of_nodes;i++)
		{
			unsigned int index = i*dim;
			values[index] = GetGeometry()[i].GetSolutionStepValue(ACCELERATION_X,Step);
			values[index + 1] = GetGeometry()[i].GetSolutionStepValue(ACCELERATION_Y,Step);
			if(dim == 3)
				values[index + 2] = GetGeometry()[i].GetSolutionStepValue(ACCELERATION_Z,Step);
		}
	}
	//************************************************************************************
	//************************************************************************************
	  void HypoelasticElement::GetFirstDerivativesVector(Vector& values, int Step)
	{
		const unsigned int number_of_nodes = GetGeometry().size();
		const unsigned int dim = GetGeometry().WorkingSpaceDimension();
		unsigned int MatSize = number_of_nodes*dim;
		if(values.size() != MatSize)   values.resize(MatSize,false);
		for (unsigned int i=0;i<number_of_nodes;i++)
		{
			unsigned int index = i*dim;
			values[index] = GetGeometry()[i].GetSolutionStepValue(VELOCITY_X,Step);
			values[index + 1] = GetGeometry()[i].GetSolutionStepValue(VELOCITY_Y,Step);
			if(dim == 3)
				values[index + 2] = GetGeometry()[i].GetSolutionStepValue(VELOCITY_Z,Step);
		//KRATOS_WATCH("IN GETFIRSTDERIVETIVE VELOCITY VALUE");
		//KRATOS_WATCH(values);
		}
	}
	//************************************************************************************
	//************************************************************************************

	void HypoelasticElement::ResizeAndInitializeAuxiliaries()
	{
		//initializing static variables
		unsigned int dimension = GetGeometry().WorkingSpaceDimension();
		unsigned int dim2 = GetGeometry().size()*dimension;
		unsigned int StrainSize;
		if (dimension==2) 
			StrainSize = 3;
		else 
			StrainSize = 6;

		if(msB.size2() != dim2)
		{
			msB.resize(StrainSize,dim2,false);
			msF.resize(dimension,dimension,false);
			msD.resize(StrainSize,StrainSize,false);
			
			msStrainVector.resize(StrainSize,false);
			msStressVector.resize(StrainSize,false);
			msDN_DX.resize(GetGeometry().size(),dimension,false);
		}
		
	}

	//************************************************************************************
	//************************************************************************************
	void HypoelasticElement::DampMatrix(MatrixType& rDampMatrix, ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY
		unsigned int number_of_nodes = GetGeometry().size();
		unsigned int dim = GetGeometry().WorkingSpaceDimension();
		
		//resizing as needed the LHS
		unsigned int MatSize=number_of_nodes*dim;

		if(rDampMatrix.size1() != MatSize)
			rDampMatrix.resize(MatSize,MatSize,false);

		noalias(rDampMatrix)= ZeroMatrix(MatSize,MatSize);
	
		
		const GeometryType::IntegrationPointsArrayType& intg_points = GetGeometry().IntegrationPoints();
		const GeometryType::ShapeFunctionsGradientsType& DNN = GetGeometry().ShapeFunctionsLocalGradients();
		

		GeometryType::JacobiansType nowJ;
		GetGeometry().Jacobian(nowJ);
		noalias(msF) = prod(nowJ[0],mInvJ0[0]);
		//KRATOS_WATCH("IN dampmatrix nowJ[0] SUPPOSED TO BE ((0,-2),(2,2))");
		//KRATOS_WATCH(nowJ);

		for(unsigned int intpn = 0; intpn < intg_points.size(); intpn++)
        {
            MathUtils<double>::InvertMatrix(nowJ[intpn],mInvJ[intpn],mDetJ[intpn]);
			msDN_DX = ZeroMatrix(msDN_DX.size1(),msDN_DX.size2());
			
			noalias(msDN_DX) = prod(DNN[intpn],mInvJ[intpn]);

			//Hypoconstitutive(msD,msF,rCurrentProcessInfo);;
			Vector StressTensor;
			Vector StrainTensor;
// 			mConstitutiveLawVector[0]->CalculateStressAndTangentMatrix(StressTensor, msF, StrainTensor, msD);
            mConstitutiveLawVector[intpn]->CalculateMaterialResponse(
                    StrainTensor,
                    msF,
                    StressTensor,
                    msD,
                    rCurrentProcessInfo,
                    GetProperties(),
                    GetGeometry(),
                    row(GetGeometry().ShapeFunctionsValues(),intpn),
                    true,
                    1,
                    true );
			CalculateB(msB,msDN_DX);
			//KRATOS_WATCH(msB);
			double weightdv = intg_points[intpn].Weight()*mDetJ[intpn];
			noalias(rDampMatrix) += prod(trans(msB),(weightdv)*Matrix(prod(msD,msB)));
			//KRATOS_WATCH(rDampMatrix);
		     }
		//KRATOS_WATCH(rDampMatrix);
		KRATOS_CATCH("")
	}
        
	//*************************************************************************************
	//*************************************************************************************
	void HypoelasticElement::MassMatrix(MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY
		unsigned int number_of_nodes = GetGeometry().size();
		unsigned int dim = GetGeometry().WorkingSpaceDimension();
		
		//resizing as needed the LHS
		unsigned int MatSize=number_of_nodes*dim;

		if(rMassMatrix.size1() != MatSize)
			rMassMatrix.resize(MatSize,MatSize,false);

		noalias(rMassMatrix)= ZeroMatrix(MatSize,MatSize);
		double massfac = mTotalDomainInitialSize * GetProperties()[DENSITY];
		if(dim == 2) massfac *= GetProperties()[THICKNESS];

		for(unsigned int ii=0;ii<= dim;ii++)
			{
			  int rowind = dim*(ii+1)-1;
			    for(unsigned int jj=0; jj<=dim; jj++)
				{
				  int colind = dim*(jj+1)-1;
					if(rowind == colind)
						rMassMatrix(rowind,colind) = 2.0*(massfac/12.0);
					else
						rMassMatrix(rowind,colind) = 1.0*(massfac/12.0);
				}
			}
				
		KRATOS_CATCH("")
	}

	//************************************************************************************
	//************************************************************************************
	void HypoelasticElement::CalculateOnIntegrationPoints(const Variable<Matrix >& rVariable, std::vector< Matrix >& Output, const ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY

		
		const GeometryType::IntegrationPointsArrayType& intg_points = GetGeometry().IntegrationPoints();
		const unsigned int dim = GetGeometry().WorkingSpaceDimension();
		
		Matrix vectorcauchyN;
		Matrix transformed_cauchy;
		Matrix tempstr;
		
		
			
		vectorcauchyN=ZeroMatrix(1,msStrainVector.size());
		vectorcauchyN = (*this).GetValue(PK2_STRESS_TENSOR);
		

		if(Output.size() !=intg_points .size())
			Output.resize(intg_points .size());
	
		if(rVariable==PK2_STRESS_TENSOR)
			{
				if (Output[0].size2()!=msStressVector.size())
				Output[0].resize(1,msStressVector.size(),false);
				
				for(unsigned int ii = 0; ii<msStrainVector.size(); ii++)
					Output[0](0,ii) = vectorcauchyN(0,ii);
			
			}
		if(rVariable==AUXILIARY_MATRIX_1)
			{
				
				Output[0].resize(dim,dim,false);	
				Output[0] = mInvJ0[0];
			
			}		

			
		
		KRATOS_CATCH("")

	}
	//***********************************************************************************
	//***********************************************************************************
// 	 void HypoelasticElement::FinalizeNonLinearIteration(ProcessInfo& CurrentProcessInfo)
// 	     {
// 		const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints();
// 		const GeometryType::ShapeFunctionsGradientsType& DN_De = GetGeometry().ShapeFunctionsLocalGradients();
// 		const Matrix& Ncontainer = GetGeometry().ShapeFunctionsValues();
// 		unsigned int number_of_nodes = GetGeometry().size();
// 		unsigned int dim = GetGeometry().WorkingSpaceDimension();
// 		int matsize = number_of_nodes*dim;
// 		
// 		//calculating actual jacobian
// 		GeometryType::JacobiansType J;
// 		GetGeometry().Jacobian(J);  
// 		
// 		
// 		Vector Prforce=ZeroVector(matsize);
// 		Vector pressure=ZeroVector(number_of_nodes);
// 		//Get pressure
// 		for(int ii=0; ii<number_of_nodes; ++ii)
// 			pressure[ii]=GetGeometry()[ii].FastGetSolutionStepValue(PRESSURE);
// 		
// 		for(unsigned int PointNumber = 0; PointNumber<integration_points.size(); PointNumber++)
// 		{
// 			//Calculating the cartesian derivatives (it is avoided storing them to minimize storage)
// 			MathUtils<double>::InvertMatrix(J[PointNumber],mInvJ[PointNumber],mDetJ[PointNumber]);	
// 			msDN_DX = ZeroMatrix(msDN_DX.size1(),msDN_DX.size2());
// 			noalias(msDN_DX) = prod(DN_De[PointNumber],mInvJ[PointNumber]);
// 			
// 
// 			
// 			//material update (considering the level of strain achieved)
// 
// 			//Calculation of constitutive
// 			//Hypoconstitutive(msD);
// 			
// 			Matrix	BPV = ZeroMatrix(matsize,dim);
// 			CalculateBPV(BPV,row(Ncontainer,PointNumber),msDN_DX);
// 			double AreaWeight = integration_points[PointNumber].Weight() * mDetJ[PointNumber];
// 			noalias(Prforce)+=(AreaWeight*prod(BPV,pressure));
// 		}
// 		for(int ii=0; ii<number_of_nodes; ++ii)
// 		     {
// 			Vector nodeforce = ZeroVector(dim);
// 			int index=ii*dim;
// 
// 			for(int jj=0; jj<dim; ++jj)
// 				nodeforce[jj]+=Prforce[index+jj];
// 
// 			(GetGeometry()[ii].FastGetSolutionStepValue(FORCE)) += nodeforce;
// 		     }
// 	     }
	//************************************************************************************
	//************************************************************************************
	/*  void HypoelasticElement::GetValuesVector(Vector& values, int Step)



	{

		const unsigned int number_of_nodes = GetGeometry().size();
		const unsigned int dim = GetGeometry().WorkingSpaceDimension();
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

	}*/

	//************************************************************************************

	//************************************************************************************
		
	
} // Namespace Kratos




