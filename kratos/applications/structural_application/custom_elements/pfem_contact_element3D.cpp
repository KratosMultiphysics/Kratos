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
 
//   
//   Project Name:        Kratos       
//   Last modified by:    $Author: virginia $
//   Date:                $Date: 2009-01-23 14:39:59 $
//   Revision:            $Revision: 1.27 $
//
//


// System includes 


// External includes 


// Project includes 
#include "includes/define.h"
#include "custom_elements/pfem_contact_element.h"
#include "utilities/math_utils.h"


namespace Kratos
{

	//************************************************************************************
	//************************************************************************************
	PfemContactElement3D::PfemContactElement3D(IndexType NewId, GeometryType::Pointer pGeometry)
		: Element(NewId, pGeometry)
	{		
		//DO NOT ADD DOFS HERE!!!
	}

	//************************************************************************************
	//************************************************************************************
	PfemContactElement3D::PfemContactElement3D(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
		: Element(NewId, pGeometry, pProperties)
	{
	}

	Element::Pointer PfemContactElement3D::Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const
	{
		return Element::Pointer(new PfemContactElement3D(NewId, GetGeometry().Create(ThisNodes), pProperties));
	}

	PfemContactElement3D::~PfemContactElement3D()
	{
	}

	
	//************************************************************************************
	//************************************************************************************
	void PfemContactElement3D::Initialize()
	{
		KRATOS_TRY

		KRATOS_CATCH("")
	}

	//************************************************************************************
	//************************************************************************************
	void PfemContactElement3D::CalculateAll(MatrixType& rLeftHandSideMatrix, 
                                           VectorType& rRightHandSideVector, 
                                           ProcessInfo& rCurrentProcessInfo,
                                           bool CalculateStiffnessMatrixFlag,
                                           bool CalculateResidualVectorFlag)
        {
		KRATOS_TRY

		
		double l0 = 2*GetProperties()[THICKNESS];
		int number_of_nodes = 4;
		if(CalculateResidualVectorFlag==false) rRightHandSideVector.resize(12,false);
		if(CalculateStiffnessMatrixFlag==false) rLeftHandSideMatrix.resize(12,12,false);

		//find reference face. If a reference face is not found, then do nothing
		int reference_face = -1;
		WeakPointerVector< Element > neigb_el = GetValue(NEIGHBOUR_ELEMENTS);
		for(unsigned int i=0;i<4; i++)
		{
			if(neigb[i].Id() == Id()) //missing neighbour
			{
				reference_face=i;
				break;
			}
		}
		
		if(reference_face != -1) //note that if a reference face is not encountered nothing is done
		{
			//calculate h, n
			Geometry< Node<3> >& geom = i->GetGeometry();
			double Volume;				
			GeometryUtils::CalculateGeometryData(geom,DN_DX,N,Volume);

			array_1d<double,3> n;
			noalias(n) = row(DN_DX,reference_face);
			norm_n = norm_2(n);
			double h=1.0/norm_n;
			n/=norm_n;
			
			if(h>l0) //in this case it is not in contact
			{
				if(CalculateResidualVectorFlag==false) noalias(rRightHandSideVector) = ZeroVector(12);
				if(CalculateStiffnessMatrixFlag==false) noalias(rLeftHandSideMatrix) = ZeroMatrix(12,12);
			}
			else
			{
				//calculate 1d strain
				double eps= 0.5*(h*h-l0*l0)/(l0*l0);

				//calculate effective strain 
				array_1d<double,6> strain(6);
				noalias(strain)  = VoigtTensorComponents(n,n);
				strain *= eps;

				//compute elastic matrix
				//setting up material matrix
				const double E = GetProperties()[YOUNG_MODULUS];
				const double NU = GetProperties()[POISSON_RATIO];
				const double c1 = E / ((1.00+NU)*(1-2*NU));
				const double c2 = c1 * (1-NU);
				const double c3 = c1 * NU;
				const double c4 = c1 * 0.5 * (1 - 2*NU);
				//filling material matrix
				C(0,0) = c2;    C(0,1) = c3;    C(0,2) = c3;    C(0,3) = 0.0;   C(0,4) = 0.0;   C(0,5) = 0.0;
				C(1,0) = c3;    C(1,1) = c2;    C(1,2) = c3;    C(1,3) = 0.0;   C(1,4) = 0.0;   C(1,5) = 0.0;
				C(2,0) = c3;    C(2,1) = c3;    C(2,2) = c2;    C(2,3) = 0.0;   C(2,4) = 0.0;   C(2,5) = 0.0;
				C(3,0) = 0.0;   C(3,1) = 0.0;   C(3,2) = 0.0;   C(3,3) = c4;    C(3,4) = 0.0;   C(3,5) = 0.0;
				C(4,0) = 0.0;   C(4,1) = 0.0;   C(4,2) = 0.0;   C(4,3) = 0.0;   C(4,4) = c4;    C(4,5) = 0.0;
				C(5,0) = 0.0;   C(5,1) = 0.0;   C(5,2) = 0.0;   C(5,3) = 0.0;   C(5,4) = 0.0;   C(5,5) = c4;

				//computer stresses
				array_1d<double,6> stress;
				noalias(stress) = prod(C,strain);
				
				//calculate B
				for (unsigned int i=0;i<number_of_nodes;i++)
				{
					unsigned int index = dimension*i;
					B(0,index+0)=DN_DX(i,0);	B(0,index+1)=0.0;		B(0,index+2)=0.0;
					B(1,index+0)=0.0;		B(1,index+1)=DN_DX(i,1);	B(1,index+2)=0.0;
					B(2,index+0)=0.0;		B(2,index+1)=0.0;		B(2,index+2)=DN_DX(i,2);
					B(3,index+0)=DN_DX(i,1);	B(3,index+1)=DN_DX(i,0);	B(3,index+2)=0.0;
					B(4,index+0)=0.0;		B(4,index+1)=DN_DX(i,2);	B(4,index+2)=DN_DX(i,1);
					B(5,index+0)=DN_DX(i,2);	B(5,index+1)=0.0;		B(5,index+2)=DN_DX(i,0);
				}
				
				//calculate RHS and LHS;
				if(CalculateResidualVectorFlag==false) 
					noalias(rRightHandSideVector) = Volume*prod(trans(B),stress);
				
				if(CalculateStiffnessMatrixFlag==false)
				{
					noalias(aux) = prod(C,B);
					noalias(rLeftHandSideMatrix) = prod(trans(B),aux);
					noalias(rLeftHandSideMatrix) *= Volume;
				}
			}

			

			

			

			
		}

		

		KRATOS_CATCH("")
	}

	//************************************************************************************
	//************************************************************************************
	void PfemContactElement3D::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
	{
		//calculation flags
		bool CalculateStiffnessMatrixFlag = false;
		bool CalculateResidualVectorFlag = true;
		MatrixType temp = Matrix();
		
		CalculateAll(temp, rRightHandSideVector, rCurrentProcessInfo, CalculateStiffnessMatrixFlag, CalculateResidualVectorFlag);
	}

	//************************************************************************************
	//************************************************************************************
	void PfemContactElement3D::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
	{
		//calculation flags
		bool CalculateStiffnessMatrixFlag = true;
		bool CalculateResidualVectorFlag = true;
		
		CalculateAll(rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo, CalculateStiffnessMatrixFlag, CalculateResidualVectorFlag);
                

	}

	

        ////************************************************************************************
	////************************************************************************************
	
        void PfemContactElement3D::InitializeSolutionStep(ProcessInfo& CurrentProcessInfo)	
 	    {
 	    }

	////************************************************************************************
	////************************************************************************************
	void PfemContactElement3D::FinalizeSolutionStep(ProcessInfo& CurrentProcessInfo)
	{
	}

	//************************************************************************************
	//************************************************************************************
	void PfemContactElement3D::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo)
	{
		int number_of_nodes = GetGeometry().size();
		int dim = GetGeometry().WorkingSpaceDimension();
		unsigned int dim2 = number_of_nodes*dim;
		if(rResult.size() != dim2)
			rResult.resize(dim2,false);

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
	void PfemContactElement3D::GetDofList(DofsVectorType& ElementalDofList,ProcessInfo& CurrentProcessInfo)
	{
		ElementalDofList.resize(0);

		for (unsigned int i=0;i<GetGeometry().size();i++)
		{
			ElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_X));
			ElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_Y));
			if(GetGeometry().WorkingSpaceDimension() == 3)
                        {
				ElementalDofList.push_back(GetGeometry()[i].pGetDof(DISPLACEMENT_Z));
                        }
		}
	}

	//************************************************************************************
	//************************************************************************************
	void PfemContactElement3D::MassMatrix(MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY

			//lumped
			unsigned int dimension = GetGeometry().WorkingSpaceDimension();
		unsigned int NumberOfNodes = GetGeometry().size();
		unsigned int MatSize = dimension * NumberOfNodes;
		if(rMassMatrix.size1() != MatSize)
			rMassMatrix.resize(MatSize,MatSize,false);

		rMassMatrix = ZeroMatrix(MatSize,MatSize);


		KRATOS_CATCH("")
	}

	//************************************************************************************
	//************************************************************************************
	void PfemContactElement3D::DampMatrix(MatrixType& rDampMatrix, ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY
		unsigned int number_of_nodes = GetGeometry().size();
		unsigned int dim = GetGeometry().WorkingSpaceDimension();
		
		//resizing as needed the LHS
		unsigned int MatSize=number_of_nodes*dim;

		if(rDampMatrix.size1() != MatSize)
			rDampMatrix.resize(MatSize,MatSize,false);

		noalias(rDampMatrix)= ZeroMatrix(MatSize,MatSize);

		KRATOS_CATCH("")
	}
        
        




	//************************************************************************************
	//************************************************************************************
	void PfemContactElement3D::GetValuesVector(Vector& values, int Step)
	{
		const unsigned int number_of_nodes = GetGeometry().size();
		const unsigned int dim = GetGeometry().WorkingSpaceDimension();
		unsigned int MatSize = number_of_nodes*dim;
		if(values.size() != MatSize)	values.resize(MatSize,false);
		for (unsigned int i=0;i<number_of_nodes;i++)
		{
			unsigned int index = i*dim;
			values[index] = GetGeometry()[i].GetSolutionStepValue(DISPLACEMENT_X,Step);
			values[index + 1] = GetGeometry()[i].GetSolutionStepValue(DISPLACEMENT_Y,Step);
			if(dim == 3)
				values[index + 2] = GetGeometry()[i].GetSolutionStepValue(DISPLACEMENT_Z,Step);
		}
	}
	
	
	//************************************************************************************
	//************************************************************************************
	  void PfemContactElement3D::GetFirstDerivativesVector(Vector& values, int Step)
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
		}
	}

	//************************************************************************************
	//************************************************************************************
	  void PfemContactElement3D::GetSecondDerivativesVector(Vector& values, int Step)
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



} // Namespace Kratos


