//
//   Project Name:        Kratos       
//   Last modified by:    $Author: kazem $
//   Date:                $Date: 2009-01-15 18:49:33 $
//   Revision:            $Revision: 1.6 $



// System includes 
// External includes 
// Project includes 

#include "includes/define.h"
#include "custom_elements/linear_element.h"
#include "utilities/math_utils.h"
//#include "constitutive_laws/isotropic_planestress_wrinkling_new.h"
#include "constitutive_laws/plane_strain.h"
#include "constitutive_laws/isotropic_3d.h"
#include "structural_application.h"



namespace Kratos

{

    Matrix LinearElement::msB(0,0);
	Matrix LinearElement::msD(0,0);
	Matrix LinearElement::msDN_DX(0,0);
	Vector LinearElement::msStrainVector(0);
	Vector LinearElement::msStressVector(0);

	//************************************************************************************

	//************************************************************************************

	LinearElement::LinearElement(IndexType NewId, GeometryType::Pointer pGeometry)
		: Element(NewId, pGeometry)
	{		
		//DO NOT ADD DOFS HERE!!!
	}

	//************************************************************************************
	//************************************************************************************

	LinearElement::LinearElement(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
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
	}
	Element::Pointer LinearElement::Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const
	{

		return Element::Pointer(new LinearElement(NewId, GetGeometry().Create(ThisNodes), pProperties));

	}
	LinearElement::~LinearElement()
	{
	}

	//************************************************************************************
	//************************************************************************************

	void LinearElement::Initialize()
	{
		KRATOS_TRY
			const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints();
		//resizing jacobian inverses containers
		mInvJ0.resize(integration_points.size());
		mDetJ0.resize(integration_points.size());

		GeometryType::JacobiansType J0;
		J0 = GetGeometry().Jacobian(J0);  
		mTotalDomainInitialSize = 0.00;

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



			//calculating the total area

			mTotalDomainInitialSize += mDetJ0[PointNumber]*IntegrationWeight;

		}





		KRATOS_CATCH("")

	}



	//************************************************************************************

	//************************************************************************************

	void LinearElement::CalculateAll(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, 
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

		//auxiliary terms

		Vector BodyForce;		

		for(unsigned int PointNumber = 0; PointNumber<integration_points.size(); PointNumber++)
		{
			//Calculating the cartesian derivatives (it is avoided storing them to minimize storage)
			msDN_DX = ZeroMatrix(msDN_DX.size1(),msDN_DX.size2());
			noalias(msDN_DX) = prod(DN_De[PointNumber],mInvJ0[PointNumber]);
			
			
			//material update (considering the level of strain achieved)

			//Calculation of stress
			

			//calculating operator B
			CalculateB(msB,msDN_DX);
			
			//calculate strain using B
			CalculateStrain(msB,msStrainVector);

			mConstitutiveLawVector[PointNumber]->CalculateStress(msStrainVector,msStressVector);

			//KRATOS_WATCH(msStrainVector);
			//calculating weights for integration on the reference configuration

			double IntToReferenceWeight = integration_points[PointNumber].Weight() * mDetJ0[PointNumber];
			if (CalculateStiffnessMatrixFlag == true) //calculation of the matrix is required

			{
				mConstitutiveLawVector[PointNumber]->CalculateConstitutiveMatrix(msStrainVector,msD);
				//mConstitutiveLawVector[PointNumber]->PlaneStrainConstitutiveMatrix(msStrainVector,msD);
				
				//KRATOS_WATCH(msD);
				//contributions to stiffness matrix calculated on the reference config

				noalias(rLeftHandSideMatrix) += prod(trans(msB),(IntToReferenceWeight)*Matrix(prod(msD,msB)) ); //to be optimized to remove the temporary			
				//KRATOS_WATCH(msD);
				//KRATOS_WATCH(rLeftHandSideMatrix)

			}



			if (CalculateResidualVectorFlag == true) //calculation of the matrix is required

			{
				//contribution to external forces 
				BodyForce = GetProperties()[BODY_FORCE];
				// operation performed: rRightHandSideVector += ExtForce*IntToReferenceWeight
				CalculateAndAdd_ExtForceContribution(row(Ncontainer,PointNumber),rCurrentProcessInfo,BodyForce,rRightHandSideVector,IntToReferenceWeight);
				// operation performed: rRightHandSideVector -= IntForce*IntToReferenceWeight
				msStressVector = prod(msD, msStrainVector);
				noalias(rRightHandSideVector) -= IntToReferenceWeight * prod(trans(msB),msStressVector);
				
			}	

		}
		
	
			

		KRATOS_CATCH("")

	}

	//************************************************************************************

	//************************************************************************************

	void LinearElement::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)

	{

		//calculation flags

		bool CalculateStiffnessMatrixFlag = false;
		bool CalculateResidualVectorFlag = true;
		MatrixType temp = Matrix();
		CalculateAll(temp, rRightHandSideVector, rCurrentProcessInfo, CalculateStiffnessMatrixFlag, CalculateResidualVectorFlag);

	}



	//************************************************************************************

	//************************************************************************************

	void LinearElement::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)

	{
		//calculation flags
		bool CalculateStiffnessMatrixFlag = true;
		bool CalculateResidualVectorFlag = true;
		
		CalculateAll(rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo, CalculateStiffnessMatrixFlag, CalculateResidualVectorFlag);

	}
	//**************************************************************************
	//**************************************************************************

	double LinearElement::CalculateIntegrationWeight(double GaussPointWeight, double DetJ0)
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

	void LinearElement::InitializeMaterial()

	{
		KRATOS_TRY
			unsigned int dimension = GetGeometry().WorkingSpaceDimension();

		//temporary - should be set from outside

		if(dimension == 2)

		{
			ConstitutiveLaw<Node<3> >::Pointer material = ConstitutiveLaw<Node<3> >::Pointer( new PlaneStrain() );GetProperties()[CONSTITUTIVE_LAW] = material;

		}

		else

		{

			ConstitutiveLaw<Node<3> >::Pointer material = ConstitutiveLaw<Node<3> >::Pointer( new Isotropic3D() );GetProperties()[CONSTITUTIVE_LAW] = material;

		}

		for (unsigned int i=0; i<mConstitutiveLawVector.size(); i++)

		{
			//temporary

			mConstitutiveLawVector[i] = GetProperties()[CONSTITUTIVE_LAW];

			mConstitutiveLawVector[i]->InitializeMaterial( 	GetProperties(), GetGeometry(),	row(GetGeometry().ShapeFunctionsValues(), i) ); 

		}

		KRATOS_CATCH("")

	} 



	//************************************************************************************

	//************************************************************************************

	inline void LinearElement::CalculateAndAdd_ExtForceContribution(

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

	void LinearElement::CalculateStrain(
		const Matrix& B,
		Vector& StrainVector)

	{

		KRATOS_TRY

			unsigned int dimension = GetGeometry().WorkingSpaceDimension();
			unsigned int number_of_nodes = GetGeometry().size();
			Vector Displacement(dimension*number_of_nodes);
	
			for (unsigned int i=0; i<number_of_nodes; i++)
				{
					int index = dimension*i;
					array_1d<double,3>& disp = GetGeometry()[i].FastGetSolutionStepValue(DISPLACEMENT);
			//KRATOS_WATCH(disp);
					Displacement[index+0] = disp[0];
					Displacement[index+1] = disp[1];
																
				if (dimension==3)
					Displacement[index+2] = disp[2];
				}	
				
		noalias(StrainVector) = prod(B,Displacement);
			
		//KRATOS_WATCH(Displacement);
		//KRATOS_WATCH(StrainVector);
		KRATOS_CATCH("")
		

	}

	//************************************************************************************

	//************************************************************************************

	void LinearElement::CalculateB(

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

	//*********************************************************************************
	//*********************************************************************************

	void LinearElement::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo)

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

	void LinearElement::GetDofList(DofsVectorType& ElementalDofList,ProcessInfo& CurrentProcessInfo)

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

	void LinearElement::CalculateOnIntegrationPoints(const Variable<Matrix >& rVariable, std::vector< Matrix >& Output, const ProcessInfo& rCurrentProcessInfo)

	{

		KRATOS_TRY
			//const unsigned int number_of_nodes = GetGeometry().size();
		//const unsigned int dim = GetGeometry().WorkingSpaceDimension();
		//resizing as needed the LHS
		//unsigned int MatSize=number_of_nodes*dim;
		//reading integration points and local gradients

		const GeometryType::IntegrationPointsArrayType& integration_points = GetGeometry().IntegrationPoints();
		const Matrix& Ncontainer = GetGeometry().ShapeFunctionsValues();

		//calculating actual jacobian
		GeometryType::JacobiansType J;
		J = GetGeometry().Jacobian(J);  

		if(Output.size() != integration_points.size())
			Output.resize(integration_points.size());

		for(unsigned int PointNumber = 0; PointNumber<integration_points.size(); PointNumber++)

		{
			CalculateStrain(msB,msStrainVector);
			if(rVariable==GREEN_LAGRANGE_STRAIN_TENSOR)
			{
				if(Output[PointNumber].size2() != msStrainVector.size())
					Output[PointNumber].resize(1,msStrainVector.size());
				for(unsigned int ii = 0; ii<msStrainVector.size(); ii++)
					Output[PointNumber](0,ii) = msStrainVector[ii];
			}

			else if(rVariable==PK2_STRESS_TENSOR)

			{

				if(Output[PointNumber].size2() != msStrainVector.size())
					Output[PointNumber].resize(1,msStrainVector.size());


				mConstitutiveLawVector[PointNumber]->UpdateMaterial( msStrainVector,
						GetProperties(),
						GetGeometry(),
						row(Ncontainer,PointNumber),
						rCurrentProcessInfo );

				mConstitutiveLawVector[PointNumber]->CalculateStress(msStrainVector,msStressVector); //saving 

				for(unsigned int ii = 0; ii<msStrainVector.size(); ii++)
					Output[PointNumber](0,ii) = msStressVector[ii];
			}

		}

		KRATOS_CATCH("")



	}



	//************************************************************************************
	void LinearElement::CalculateLeftHandSide(MatrixType& rLeftHandSideMatrix, ProcessInfo& rCurrentProcessInfo)
		{
		
		bool CalculateStiffnessMatrixFlag = true;
		bool CalculateResidualVectorFlag = false;
		Vector temp = Vector();
		CalculateAll(rLeftHandSideMatrix, temp, rCurrentProcessInfo, CalculateStiffnessMatrixFlag, CalculateResidualVectorFlag);
		
		}
	//************************************************************************************

	/*  void LinearElement::GetValuesVector(Vector& values, int Step)



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




