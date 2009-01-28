//   
//   Project Name:        Kratos       
//   Last modified by:    $Author: pooyan $
//   Date:                $Date: 2006-11-27 16:07:42 $
//   Revision:            $Revision: 1.1.1.1 $
//
// 


// System includes 


// External includes 


// Project includes 
#include "includes/define.h"
#include "custom_elements/coupling_face2D.h"
#include "utilities/math_utils.h"
#include "fsi_application.h"

namespace Kratos
{
	//************************************************************************************
	//************************************************************************************
	CouplingFace2D::CouplingFace2D(IndexType NewId, GeometryType::Pointer pGeometry)
		: Element(NewId, pGeometry)
	{		
		//DO NOT ADD DOFS HERE!!!
	}

	//************************************************************************************
	//************************************************************************************
	CouplingFace2D::CouplingFace2D(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
		: Element(NewId, pGeometry, pProperties)
	{
	}

	Element::Pointer CouplingFace2D::Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const
	{
		return Element::Pointer(new CouplingFace2D(NewId, GetGeometry().Create(ThisNodes), pProperties));
	}

	CouplingFace2D::~CouplingFace2D()
	{
	}


	//************************************************************************************
	//************************************************************************************
	void CouplingFace2D::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
	{
		//calculation flags
		bool CalculateStiffnessMatrixFlag = false;
		bool CalculateResidualVectorFlag = true;
		MatrixType temp = Matrix();

		CalculateAll(temp, rRightHandSideVector, rCurrentProcessInfo, CalculateStiffnessMatrixFlag, CalculateResidualVectorFlag);

		//KRATOS_WATCH(rRightHandSideVector);
	}

	//************************************************************************************
	//************************************************************************************
	void CouplingFace2D::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
	{
		//calculation flags
		bool CalculateStiffnessMatrixFlag = true;
		bool CalculateResidualVectorFlag = true;
		
		CalculateAll(rLeftHandSideMatrix, rRightHandSideVector, rCurrentProcessInfo, CalculateStiffnessMatrixFlag, CalculateResidualVectorFlag);

		//KRATOS_WATCH(rLeftHandSideMatrix);
		//KRATOS_WATCH(rRightHandSideVector);
	}

	//************************************************************************************
	//************************************************************************************
	void CouplingFace2D::CalculateAll(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, 
										ProcessInfo& rCurrentProcessInfo,
										bool CalculateStiffnessMatrixFlag,
										bool CalculateResidualVectorFlag)
	{
		KRATOS_TRY

			unsigned int number_of_nodes = GetGeometry().size();
		unsigned int dim = GetGeometry().WorkingSpaceDimension();

		//resizing as needed the LHS
		unsigned int MatSize=number_of_nodes*dim;

		//calculating the lenght
		double lx = GetGeometry()[1].X() - GetGeometry()[0].X();
		double ly = GetGeometry()[1].Y() - GetGeometry()[0].Y();
		double node_area = 0.5*sqrt(lx*lx + ly*ly);

		double coupling_density  = GetProperties()[FICTITIOUS_FLUID_DENSITY];


		if (CalculateStiffnessMatrixFlag == true) //calculation of the matrix is required
		{
			if(rLeftHandSideMatrix.size1() != MatSize)
				rLeftHandSideMatrix.resize(MatSize,MatSize,false);
			noalias(rLeftHandSideMatrix) = ZeroMatrix(MatSize,MatSize); //resetting LHS
		}

		//resizing as needed the RHS - this assumes that the acceleration was saved at the previouws iteration
		if (CalculateResidualVectorFlag == true) //calculation of the matrix is required
		{
			if(rRightHandSideVector.size() != MatSize)
				rRightHandSideVector.resize(MatSize,false);
			for (unsigned int i=0;i<number_of_nodes;i++)
			{
				unsigned int index = i*dim;
				const array_1d<double,3>& acc_old = GetGeometry()[i].GetValue(ACCELERATION);
				rRightHandSideVector[index] = coupling_density * node_area * acc[0];
				rRightHandSideVector[index+1] = coupling_density * node_area * acc[1];
			}

		}

		KRATOS_CATCH("")
	}

	//************************************************************************************
	//************************************************************************************
	void CouplingFace2D::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo)
	{
		int number_of_nodes = GetGeometry().PointsNumber();
		unsigned int index;
		unsigned int dim = 2;
		if(rResult.size() != number_of_nodes*dim)
			rResult.resize(number_of_nodes*dim);
		for (int i=0;i<number_of_nodes;i++)
		{
			index = i*dim;
			rResult[index] = (GetGeometry()[i].GetDof(DISPLACEMENT_X)).EquationId();
			rResult[index+1] = (GetGeometry()[i].GetDof(DISPLACEMENT_Y).EquationId());
		}
	}
	//***********************************************************************************
	//***********************************************************************************

	void CouplingFace2D::MassMatrix(
		MatrixType& rMassMatrix,
		ProcessInfo& rCurrentProcessInfo)

	{
		KRATOS_TRY

		if(rMassMatrix.size1() != 4)
			rMassMatrix.resize(4);

		//calculating the lenght
		double lx = GetGeometry()[1].X() - GetGeometry()[0].X();
		double ly = GetGeometry()[1].Y() - GetGeometry()[0].Y();
		double node_area = 0.5*sqrt(lx*lx + ly*ly);

		double coupling_density  = GetProperties()[FICTITIOUS_FLUID_DENSITY];





		//calculating the lenght
		double lx = GetGeometry()[1].X() - GetGeometry()[0].X();
		double ly = GetGeometry()[1].Y() - GetGeometry()[0].Y();
		double node_mass = 0.5*sqrt(lx*lx + ly*ly);

double coupling_density  = 1000;

		double TotalMass = mTotalDomainInitialSize * GetProperties()[THICKNESS] * GetProperties()[DENSITY];
		Vector LumpFact;
		LumpFact = GetGeometry().LumpingFactors(LumpFact);
		
		for(unsigned int i=0; i<number_of_nodes; i++)
		{
			for(unsigned int j=0; j<2; j++)
			{
				unsigned int index = i*2 + j;
				rMassMatrix(index,index) = coupling_density * node_mass;
			}
		}

		KRATOS_CATCH("")
	}
	*/

	//************************************************************************************
	//************************************************************************************
	  void CouplingFace2D::GetDofList(DofsVectorType& ElementalDofList,ProcessInfo& CurrentProcessInfo)
	{
		unsigned int dim = 2;
		ElementalDofList.resize(GetGeometry().size()*dim);
		unsigned int index;
		for (unsigned int i=0;i<GetGeometry().size();i++)
		{
			index = i*dim;
			ElementalDofList[index] = (GetGeometry()[i].pGetDof(DISPLACEMENT_X));
			ElementalDofList[index+1] = (GetGeometry()[i].pGetDof(DISPLACEMENT_Y));
		}
	}
} // Namespace Kratos


