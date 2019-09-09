// Project includes 
#include "includes/define.h"
#include "custom_conditions/fixed_pressure_2d.h"
#include "pfem_2_application_variables.h"
#include "utilities/math_utils.h"

namespace Kratos
{
	//************************************************************************************
	//************************************************************************************
	FixedPressure2D::FixedPressure2D(IndexType NewId, GeometryType::Pointer pGeometry)
		: Condition(NewId, pGeometry)
	{		
		//DO NOT ADD DOFS HERE!!!
	}

	//************************************************************************************
	//************************************************************************************
	FixedPressure2D::FixedPressure2D(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
		: Condition(NewId, pGeometry, pProperties)
	{
	}
	Condition::Pointer FixedPressure2D::Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const
	{
		return Condition::Pointer(new FixedPressure2D(NewId, GetGeometry().Create(ThisNodes), pProperties));
	}

	FixedPressure2D::~FixedPressure2D()
	{
	}

	//************************************************************************************
	//************************************************************************************
	void FixedPressure2D::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY
		{
			KRATOS_THROW_ERROR(std::logic_error,  "method not implemented" , "");
		}
		KRATOS_CATCH("")
	}

	//************************************************************************************
	//************************************************************************************
	void FixedPressure2D::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY
		{
			const unsigned int TDim = 2;
			
			if(rLeftHandSideMatrix.size1() != TDim)
				rLeftHandSideMatrix.resize(TDim,TDim,false);
			if(rRightHandSideVector.size() != TDim)
				rRightHandSideVector.resize(TDim,false);

				
			const array_1d<double,3> & normal = GetGeometry()[0].FastGetSolutionStepValue(NORMAL);
			const double pressure = GetGeometry()[0].FastGetSolutionStepValue(PRESSURE);

			rLeftHandSideMatrix =ZeroMatrix(TDim,TDim);
			
			for (unsigned int i=0; i!=TDim ; i++)
				rRightHandSideVector[i] = - pressure*normal(i);
			
		}
		KRATOS_CATCH("")
	}


	//************************************************************************************
	//************************************************************************************
	void FixedPressure2D::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo)
	{
		{
			int number_of_nodes = GetGeometry().PointsNumber();
			const unsigned int TDim = 2;
			rResult.resize(number_of_nodes*TDim);
			rResult[0] = (GetGeometry()[0].GetDof(VELOCITY_X).EquationId());			
			rResult[1] = (GetGeometry()[0].GetDof(VELOCITY_Y).EquationId());			
			
		}
	}

	//************************************************************************************
	//************************************************************************************
	  void FixedPressure2D::GetDofList(DofsVectorType& ConditionalDofList,ProcessInfo& rCurrentProcessInfo)
	{
		{
			const unsigned int TDim = 2;
			ConditionalDofList.resize(GetGeometry().size()*TDim);
			ConditionalDofList[0] = (GetGeometry()[0].pGetDof(VELOCITY_X));
			ConditionalDofList[1] = (GetGeometry()[0].pGetDof(VELOCITY_Y));

		}
	}
} // Namespace Kratos
