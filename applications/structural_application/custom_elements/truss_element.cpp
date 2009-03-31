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
/* *********************************************************   
*          
*   Last Modified by:    $Author: hurga $
*   Date:                $Date: 2008-02-15 10:37:19 $
*   Revision:            $Revision: 1.2 $
*
* ***********************************************************/

// System includes 


// External includes 

// Project includes 
// #include "includes/define.h"
#include "custom_elements/truss_element.h"
#include "utilities/math_utils.h"
#include "custom_utilities/sd_math_utils.h"

namespace Kratos
{

    TrussElement::TrussElement(IndexType NewId, 
            GeometryType::Pointer pGeometry)
    : Element(NewId, pGeometry)
    {		
		//DO NOT ADD DOFS HERE!!!
		//THIS IS THE DEFAULT CONSTRUCTOR
    }

	//************************************************************************************
	//************************************************************************************
	//THIS IS A TRUSS ELEMENT FOR FINITE STRAINS (ALSO CALLED CRISFIELD-ELEMENT)
	// it consists of a 2Node Element with neo-Hook Material behaviour
	// all operations are done inside the elements	
	// it does not contain works only with line2d and does not contain a material
	//************************************************************************************
	//************************************************************************************
    TrussElement::TrussElement(IndexType NewId, 
            GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
    : Element(NewId, pGeometry, pProperties)
    {
		//DEFINE HERE THE INTEGRATION METHOD (the number of integration points to be used)
		//see /kratos/source/integration_rules.cpp for further details
		if(GetGeometry().size() != 2)
		{
			std::cout<<"this element works only with a 2 node line"<<std::endl;
		}

		GetGeometry()[0].pAddDof(DISPLACEMENT_X, REACTION_X);
    	GetGeometry()[0].pAddDof(DISPLACEMENT_Y, REACTION_Y);
   		GetGeometry()[0].pAddDof(DISPLACEMENT_Z, REACTION_Z);
   		GetGeometry()[1].pAddDof(DISPLACEMENT_X, REACTION_X);
 		GetGeometry()[1].pAddDof(DISPLACEMENT_Y, REACTION_Y);
        GetGeometry()[1].pAddDof(DISPLACEMENT_Z, REACTION_Z);

		//set up the Matrix A
		if(mA_Matrix.size1()!=6 || mA_Matrix.size2()!=6)
			mA_Matrix.resize(6,6,false);
		noalias(mA_Matrix)= ZeroMatrix(6,6);
		mA_Matrix(0,0)=0.25;mA_Matrix(1,1)=0.25;mA_Matrix(2,2)=0.25;		mA_Matrix(0,3)=-0.25;mA_Matrix(1,4)=-0.25;mA_Matrix(2,5)=-0.25;
		mA_Matrix(3,0)=-0.25;mA_Matrix(4,1)=-0.25;mA_Matrix(5,2)=-0.25;		mA_Matrix(3,3)=0.25;mA_Matrix(4,4)=0.25;mA_Matrix(5,5)=0.25;
	
		if(mCurrentDisplacement.size() != 6)
			mCurrentDisplacement.resize(6,false);
		if(mAtimesU.size() != 6)
			mAtimesU.resize(6,false);
    }

    Element::Pointer TrussElement::Create(IndexType NewId, 
            NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const
    {
        return Element::Pointer(new TrussElement(NewId, GetGeometry().Create(ThisNodes), 
                                pProperties));
    }

	//************************************************************************************
	//************************************************************************************
	//THIS IS THE DESTRUCTOR, as we use UBLAS classes, we don't have to care about garbage collecting
	//************************************************************************************
	//************************************************************************************
    TrussElement::~TrussElement()
    {
    }
	//************************************************************************************
	//THIS IS THE INITIALIZATION OF THE ELEMENT (CALLED AT THE BEGIN OF EACH CALCULATION)
	//************************************************************************************
    void TrussElement::Initialize()
    {
       	KRATOS_TRY
		//Material Porperties
		mArea= 0.01;
	    mYoungs= GetProperties()[YOUNG_MODULUS];
// mPrescribedStrain= 0.1;
// 		initializing mb
		Vector x_zero(6);

		x_zero(0)= GetGeometry()[0].X0();x_zero(1)= GetGeometry()[0].Y0();x_zero(2)= GetGeometry()[0].Z0();
		x_zero(3)= GetGeometry()[1].X0();x_zero(4)= GetGeometry()[1].Y0();x_zero(5)= GetGeometry()[1].Z0();
// 		calculate reference length
		mReference_length= 2*sqrt(VectorMatrixVector(x_zero));
// 		initializing mb
		//set up the Matrix b
		if(mb_Vector.size()!=6)
			mb_Vector.resize(6,false);

		noalias(mb_Vector)= 4/(mReference_length*mReference_length)*MatrixVector(x_zero);

        KRATOS_CATCH("")
    }	

	//************************************************************************************
	//************************************************************************************
	//this is to finalize the solution step ("time step")
	//************************************************************************************
	//************************************************************************************
 
    void TrussElement::InitializeSolutionStep(ProcessInfo& CurrentProcessInfo)
    {
// 		mPrescribedStrain= GetProperties()[POISSON_RATIO];
    }
	//************************************************************************************
	//************************************************************************************
	//this is to finalize the solution step ("time step")
	//************************************************************************************
	//************************************************************************************
 
    void TrussElement::FinalizeSolutionStep(ProcessInfo& CurrentProcessInfo)
    {

    }

	//************************************************************************************
	//THIS is the main method here the integration in space (loop over the integration points) is done
	//************************************************************************************
    void TrussElement::CalculateAll(MatrixType& rLeftHandSideMatrix, 
            VectorType& rRightHandSideVector,ProcessInfo& rCurrentProcessInfo,
            bool CalculateStiffnessMatrixFlag,bool CalculateResidualVectorFlag)
    {
        KRATOS_TRY
		//update the current displacement
		mCurrentDisplacement(0)= GetGeometry()[0].GetSolutionStepValue(DISPLACEMENT_X);	
		mCurrentDisplacement(1)= GetGeometry()[0].GetSolutionStepValue(DISPLACEMENT_Y);
		mCurrentDisplacement(2)= GetGeometry()[0].GetSolutionStepValue(DISPLACEMENT_Z);
		mCurrentDisplacement(3)= GetGeometry()[1].GetSolutionStepValue(DISPLACEMENT_X);
		mCurrentDisplacement(4)= GetGeometry()[1].GetSolutionStepValue(DISPLACEMENT_Y);
		mCurrentDisplacement(5)= GetGeometry()[1].GetSolutionStepValue(DISPLACEMENT_Z);
		//calculate AtimesU
		CalcAtimesU();
		//update the current strain
		GreenStrain();

		//resizing and calculate the LHS contribution if needed
		if (CalculateStiffnessMatrixFlag == true) 
        {
                 if(rLeftHandSideMatrix.size1() != 6)
                 	rLeftHandSideMatrix.resize(6,6,false);
                 noalias(rLeftHandSideMatrix) = ZeroMatrix(6,6); 
				 CalculateLHS(rLeftHandSideMatrix);
         }
		 //resizing and calculate the RHS contribution if needed
         if (CalculateResidualVectorFlag == true) 
         {
                 if(rRightHandSideVector.size() != 6)
                       rRightHandSideVector.resize(6,false);
                 noalias(rRightHandSideVector) = ZeroVector(6); 
				 CalculateRHS(rRightHandSideVector);
         }

        KRATOS_CATCH("")
    }

	//************************************************************************************
	//************************************************************************************
	//This method is called from outside the element
	//************************************************************************************
 	//************************************************************************************
     void TrussElement::CalculateRightHandSide(VectorType& rRightHandSideVector, 
             ProcessInfo& rCurrentProcessInfo)
      {
	   //calculation flags
           bool CalculateStiffnessMatrixFlag = false;
           bool CalculateResidualVectorFlag = true;
           MatrixType temp = Matrix();
		
           CalculateAll(temp, rRightHandSideVector, rCurrentProcessInfo, CalculateStiffnessMatrixFlag,  CalculateResidualVectorFlag);
       }	

	//************************************************************************************
	//************************************************************************************
	//This method is called from outside the element
	//************************************************************************************
 	//************************************************************************************
 
        void TrussElement::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, 
                VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
        {
	   //calculation flags
            bool CalculateStiffnessMatrixFlag = true;
            bool CalculateResidualVectorFlag = true;
            CalculateAll(rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo, 
                     CalculateStiffnessMatrixFlag,CalculateResidualVectorFlag);
        }

    //************************************************************************************
	//************************************************************************************
	// returns the used integration method 
    //************************************************************************************
	//************************************************************************************
	TrussElement::IntegrationMethod TrussElement::GetIntegrationMethod()
	{
			return GeometryData::GI_GAUSS_1;
	}

    //************************************************************************************
	//Informations to assemble the global vectors and matrices
	//************************************************************************************

        void TrussElement::EquationIdVector(EquationIdVectorType& rResult, 
                ProcessInfo& CurrentProcessInfo)
        {
            if(rResult.size() != 6)
                    rResult.resize(6);

 			rResult[0] = GetGeometry()[0].GetDof(DISPLACEMENT_X).EquationId();
   			rResult[1] = GetGeometry()[0].GetDof(DISPLACEMENT_Y).EquationId();
			rResult[2] = GetGeometry()[0].GetDof(DISPLACEMENT_Z).EquationId();
 			rResult[3] = GetGeometry()[1].GetDof(DISPLACEMENT_X).EquationId();
   			rResult[4] = GetGeometry()[1].GetDof(DISPLACEMENT_Y).EquationId();
			rResult[5] = GetGeometry()[1].GetDof(DISPLACEMENT_Z).EquationId();
        }

	//************************************************************************************
	//************************************************************************************
        
        void TrussElement::GetDofList(DofsVectorType& ElementalDofList,ProcessInfo& 
                CurrentProcessInfo)
        {
                ElementalDofList.resize(0);

 				ElementalDofList.push_back(GetGeometry()[0].pGetDof(DISPLACEMENT_X));
          		ElementalDofList.push_back(GetGeometry()[0].pGetDof(DISPLACEMENT_Y));
				ElementalDofList.push_back(GetGeometry()[0].pGetDof(DISPLACEMENT_Z));
 				ElementalDofList.push_back(GetGeometry()[1].pGetDof(DISPLACEMENT_X));
          		ElementalDofList.push_back(GetGeometry()[1].pGetDof(DISPLACEMENT_Y));
				ElementalDofList.push_back(GetGeometry()[1].pGetDof(DISPLACEMENT_Z));  
         }

	//************************************************************************************
	//************************************************************************************
        void TrussElement::GetValuesVector(Vector& values, int Step)
        {
                if(values.size() != 6)
                    values.resize(6);

				values(0) = GetGeometry()[0].GetSolutionStepValue(DISPLACEMENT_X,Step);
                values(1) = GetGeometry()[0].GetSolutionStepValue(DISPLACEMENT_Y,Step); 
                values(2) = GetGeometry()[0].GetSolutionStepValue(DISPLACEMENT_Z,Step); 
				values(3) = GetGeometry()[1].GetSolutionStepValue(DISPLACEMENT_X,Step);
                values(4) = GetGeometry()[1].GetSolutionStepValue(DISPLACEMENT_Y,Step); 
                values(5) = GetGeometry()[1].GetSolutionStepValue(DISPLACEMENT_Z,Step); 
          }

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////
               
		double TrussElement::VectorMatrixVector(Vector& vector)
		{
			if(vector.size()!= 6)
				std::cout<<"this works only for the displacementVectors together with the A-Matrix"<<std::endl;

			double result= 0.25*(
				vector(0)*(vector(0)-vector(3))+
				vector(1)*(vector(1)-vector(4))+
				vector(2)*(vector(2)-vector(5))+
				vector(3)*(-vector(0)+vector(3))+
				vector(4)*(-vector(1)+vector(4))+
				vector(5)*(-vector(2)+vector(5)));

			return result;
		}
//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////
               
		Vector TrussElement::MatrixVector(Vector& vector)
		{
			Vector result(6);

			if(vector.size()!= 6)
				std::cout<<"this works only for the displacementVectors together with the A-Matrix"<<std::endl;
		
			result(0) = 0.25*(vector(0)-vector(3));
			result(1) = 0.25*(vector(1)-vector(4));
			result(2) = 0.25*(vector(2)-vector(5));
			result(3) = 0.25*(-vector(0)+vector(3));
			result(4) = 0.25*(-vector(1)+vector(4));
			result(5) = 0.25*(-vector(2)+vector(5));

			return result;
		}
//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////
  
		void TrussElement::GreenStrain()
		{
			mCurrentStrain= mb_Vector(0)*mCurrentDisplacement(0)+ mb_Vector(1)*mCurrentDisplacement(1) +mb_Vector(2)*mCurrentDisplacement(2)+mb_Vector(3)*mCurrentDisplacement(3)+ mb_Vector(4)*mCurrentDisplacement(4) +mb_Vector(5)*mCurrentDisplacement(5);

			mCurrentStrain+= 2/(mReference_length*mReference_length)*VectorMatrixVector(mCurrentDisplacement);

			mCurrentStrain-=mPrescribedStrain;
		}
//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////
		void TrussElement::CalcAtimesU()
		{
			mAtimesU(0) = 0.25*(mCurrentDisplacement(0)-mCurrentDisplacement(3));
			mAtimesU(1) = 0.25*(mCurrentDisplacement(1)-mCurrentDisplacement(4));
			mAtimesU(2) = 0.25*(mCurrentDisplacement(2)-mCurrentDisplacement(5));
			mAtimesU(3) = 0.25*(-mCurrentDisplacement(0)+mCurrentDisplacement(3));
			mAtimesU(4) = 0.25*(-mCurrentDisplacement(1)+mCurrentDisplacement(4));
			mAtimesU(5) = 0.25*(-mCurrentDisplacement(2)+mCurrentDisplacement(5));
		}
//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////
		void TrussElement::CalculateRHS(Vector& rRightHandSideVector)
		{
			for(unsigned int prim=0; prim<6; prim++)
			{
				rRightHandSideVector(prim)+= mArea*mReference_length*mYoungs*(mb_Vector(prim)+4/(mReference_length*mReference_length)*mAtimesU(prim))*mCurrentStrain;
			}
			return;
		}
//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////
		void TrussElement::CalculateLHS(Matrix& rLeftHandSideMatrix)
		{
			for(unsigned int prim=0; prim<6; prim++)
			{
				for(unsigned int sec=0; sec<6; sec++)
				{
					rLeftHandSideMatrix(prim,sec)+= (-1)*mArea*mYoungs*4/mReference_length*mCurrentStrain*mA_Matrix(prim,sec);
					rLeftHandSideMatrix(prim,sec)+= (-1)*mArea*mYoungs*mReference_length*(mb_Vector(prim)+4/(mReference_length*mReference_length)*mAtimesU(prim))*(mb_Vector(sec)+4/(mReference_length*mReference_length)*mAtimesU(sec));
				}
			}

			return;
		}

	  void TrussElement::SetValueOnIntegrationPoints(const Variable<double>& rVariable, std::vector<double>& rValues, const ProcessInfo& rCurrentProcessInfo)
	  	{
std::cout<<"hier isser"<<std::endl;
			if( rVariable== THICKNESS)
			{
std::cout<<"hier isser"<<std::endl;
				mPrescribedStrain= rValues[0];
KRATOS_WATCH(mPrescribedStrain);
			}
		}
} // Namespace Kratos


