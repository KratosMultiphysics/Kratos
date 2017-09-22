// Project includes 
#include "includes/define.h"
#include "custom_conditions/lagrange_multiplier_velocity_gradients_2d.h"
#include "fluid_rve_lagrange_multipliers_application.h"
#include "utilities/math_utils.h"
#include "utilities/geometry_utilities.h" 


namespace Kratos
{
	//************************************************************************************
	//************************************************************************************
	VelocityGradientsLagrangeMultiplierCondition2D::VelocityGradientsLagrangeMultiplierCondition2D(IndexType NewId, GeometryType::Pointer pGeometry)
		: Condition(NewId, pGeometry)
	{		
		//DO NOT ADD DOFS HERE!!!
	}

	//************************************************************************************
	//************************************************************************************
	VelocityGradientsLagrangeMultiplierCondition2D::VelocityGradientsLagrangeMultiplierCondition2D(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
		: Condition(NewId, pGeometry, pProperties)
	{
	}
	Condition::Pointer VelocityGradientsLagrangeMultiplierCondition2D::Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const
	{
		return Condition::Pointer(new VelocityGradientsLagrangeMultiplierCondition2D(NewId, GetGeometry().Create(ThisNodes), pProperties));
	}

	VelocityGradientsLagrangeMultiplierCondition2D::~VelocityGradientsLagrangeMultiplierCondition2D()
	{
	}

	//************************************************************************************
	//************************************************************************************
	void VelocityGradientsLagrangeMultiplierCondition2D::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY
		/*
		if(rCurrentProcessInfo[FRACTIONAL_STEP]==2)
		{
			KRATOS_THROW_ERROR(std::logic_error,  "method not implemented" , "");
		}
		*/
		KRATOS_THROW_ERROR(std::logic_error,  "method not implemented" , "");
		
		KRATOS_CATCH("")
	}

	//************************************************************************************
	//************************************************************************************
	void VelocityGradientsLagrangeMultiplierCondition2D::CalculateLocalVelocityContribution(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)

	{
		KRATOS_TRY

		//we only have to calculate if 
		const unsigned int dim = 2; // 2 = only 1 crossed gradients + du/dx . 3 = 2 crossed gradients (du/dy and dv/dx)  + du/dx

		if(dim==2) //(rCurrentProcessInfo[FRACTIONAL_STEP]==2)
		{
			const unsigned int dim = 2;
			if(rLeftHandSideMatrix.size1() != dim*4)
				rLeftHandSideMatrix.resize(dim*4,dim*4,false);
			noalias(rLeftHandSideMatrix) = ZeroMatrix(dim*4,dim*4);
			if(rRightHandSideVector.size() != dim*4)
				rRightHandSideVector.resize(dim*4,false);
			rRightHandSideVector = ZeroVector(dim*4);
		
		
			const array_1d<double,3> & target_gradients = rCurrentProcessInfo[LAGRANGE_MULTIPLIER_VELOCITY_GRADIENTS];

			double Area;
			array_1d<double, dim+1> N;
			boost::numeric::ublas::bounded_matrix<double, dim+1, dim> DN_DX;
			boost::numeric::ublas::bounded_matrix<double, 2, (dim+1)*dim > LM_matrix = ZeroMatrix(2, (dim+1)*dim) ; //(gradient)
			Geometry<Node<3> >& geom = this->GetGeometry();
			GeometryUtils::CalculateGeometryData(geom, DN_DX, N, Area);
			
			for (unsigned int j = 0; j < dim+1; j++) //3 nodes
			{
				LM_matrix(0, (j*dim+0) ) =  DN_DX(j,1)*Area; // du/dy : we use DNi/DNy
				LM_matrix(1, (j*dim+1) ) =  DN_DX(j,0)*Area; // dv/dx : we use DNi/DNx
			}
			
			
			for (unsigned int i = 0; i < (dim+1)*2; i++) //velocity dofs
			{
				for (unsigned int j = 0; j < 2 ; j++) // LM	
				{
					rLeftHandSideMatrix(i,j+6) = LM_matrix(j,i);
					rLeftHandSideMatrix(j+6,i) = LM_matrix(j,i);
				}
			}
			
			rRightHandSideVector(6) = target_gradients[0]*Area;
			rRightHandSideVector(7) = target_gradients[1]*Area;

			
			//substracting:
			array_1d<double,dim*4>  temp_vector=ZeroVector(dim*4);
			for (unsigned int j = 0; j < dim+1; j++) //3 nodes
			{
				const array_1d<double,3> & velocity = GetGeometry()[j].GetSolutionStepValue(VELOCITY,0);
				temp_vector[j*2]=velocity[0];
				temp_vector[j*2+1]=velocity[1];
			}
			const array_1d<double,3> & velocity_aux_node = (this->GetValue(NEIGHBOUR_NODES)[0]).GetSolutionStepValue(LAGRANGE_MULTIPLIER_VELOCITY_GRADIENTS,0);
			temp_vector[6]=velocity_aux_node[0];
			temp_vector[7]=velocity_aux_node[1];
 			noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix,temp_vector);
			
			//KRATOS_WATCH(rLeftHandSideMatrix);
			//KRATOS_WATCH(rRightHandSideVector);

			
		}

		else
		{
			if(rLeftHandSideMatrix.size1() != 2*3+3) // two dofs per node + 3 multipliers
				rLeftHandSideMatrix.resize(2*3+3,2*3+3,false);
			noalias(rLeftHandSideMatrix) = ZeroMatrix(2*3+3,2*3+3);
			if(rRightHandSideVector.size() != 2*3+3)
				rRightHandSideVector.resize(2*3+3,false);
			rRightHandSideVector = ZeroVector(2*3+3);
		
		
			const array_1d<double,3> & target_gradients = rCurrentProcessInfo[LAGRANGE_MULTIPLIER_VELOCITY_GRADIENTS];

			double Area;
			array_1d<double, 3> N;
			boost::numeric::ublas::bounded_matrix<double, 3, 2> DN_DX;
			boost::numeric::ublas::bounded_matrix<double, 3, 3*2 > LM_matrix = ZeroMatrix(3, 3*2) ; //(gradient)
			Geometry<Node<3> >& geom = this->GetGeometry();
			GeometryUtils::CalculateGeometryData(geom, DN_DX, N, Area);
			
			for (unsigned int j = 0; j < 3 ; j++) //3 nodes
			{
				LM_matrix(0, (j*2+0) ) =  DN_DX(j,1)*Area; // du/dy : we use DNi/DNy
				LM_matrix(1, (j*2+1) ) =  DN_DX(j,0)*Area; // dv/dx : we use DNi/DNx
				LM_matrix(2, (j*2) ) =  DN_DX(j,0)*Area; // dv/dx : we use DNi/DNx

				
			}
			
			
			for (unsigned int i = 0; i < 3*2; i++) //velocity dofs
			{
				for (unsigned int j = 0; j < 3 ; j++) // LM	
				{
					rLeftHandSideMatrix(i,j+6) = LM_matrix(j,i);
					rLeftHandSideMatrix(j+6,i) = LM_matrix(j,i);
				}
			}
			
			rRightHandSideVector(6) = target_gradients[0]*Area;
			rRightHandSideVector(7) = target_gradients[1]*Area;
			rRightHandSideVector(8) = target_gradients[2]*Area;


			
			//substracting:
			array_1d<double,dim*4>  temp_vector=ZeroVector(2*3+3);
			for (unsigned int j = 0; j < 3 ; j++) //3 nodes
			{
				const array_1d<double,3> & velocity = GetGeometry()[j].GetSolutionStepValue(VELOCITY,0);
				temp_vector[j*2]=velocity[0];
				temp_vector[j*2+1]=velocity[1];
			}
			const array_1d<double,3> & velocity_aux_node = (this->GetValue(NEIGHBOUR_NODES)[0]).GetSolutionStepValue(LAGRANGE_MULTIPLIER_VELOCITY_GRADIENTS,0);
			temp_vector[6]=velocity_aux_node[0];
			temp_vector[7]=velocity_aux_node[1];
			temp_vector[8]=velocity_aux_node[2];

 			noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix,temp_vector);
			
			//KRATOS_WATCH(rLeftHandSideMatrix);
			//KRATOS_WATCH(rRightHandSideVector);

			
		}
	
		KRATOS_CATCH("")
	}


	void VelocityGradientsLagrangeMultiplierCondition2D::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)

	{
		KRATOS_TRY

		//we only have to calculate if 
		
		const unsigned int dim = 2;

		if(dim == 2) //(rCurrentProcessInfo[FRACTIONAL_STEP]==2)
		{
			if(rLeftHandSideMatrix.size1() != dim*4)
				rLeftHandSideMatrix.resize(dim*4,dim*4,false);
			noalias(rLeftHandSideMatrix) = ZeroMatrix(dim*4,dim*4);
			if(rRightHandSideVector.size() != dim*4)
				rRightHandSideVector.resize(dim*4,false);
			rRightHandSideVector = ZeroVector(dim*4);
		}
		else
		{
			if(rLeftHandSideMatrix.size1() != 2*3+3)
				rLeftHandSideMatrix.resize(2*3+3,2*3+3,false);
			noalias(rLeftHandSideMatrix) = ZeroMatrix(2*3+3,2*3+3);
			if(rRightHandSideVector.size() != 2*3+3)
				rRightHandSideVector.resize(2*3+3,false);
			rRightHandSideVector = ZeroVector(2*3+3);
		}
		KRATOS_CATCH("")
	}


	//************************************************************************************
	//************************************************************************************
	void VelocityGradientsLagrangeMultiplierCondition2D::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo)
	{
		unsigned int dim = 2;
		if (dim==2)
		{
            if (rResult.size()!=(4*dim))
			    rResult.resize(4*dim);
			    
			unsigned int LocalIndex = 0;    
			for (unsigned int iNode = 0; iNode < 3; ++iNode)
			{
				rResult[LocalIndex++] = (GetGeometry()[iNode].GetDof(VELOCITY_X).EquationId());
				rResult[LocalIndex++] = (GetGeometry()[iNode].GetDof(VELOCITY_Y).EquationId());
			}
			rResult[LocalIndex++] = (this->GetValue(NEIGHBOUR_NODES)[0]).GetDof(LAGRANGE_MULTIPLIER_VELOCITY_GRADIENTS_X).EquationId();			
			rResult[LocalIndex++] = (this->GetValue(NEIGHBOUR_NODES)[0]).GetDof(LAGRANGE_MULTIPLIER_VELOCITY_GRADIENTS_Y).EquationId();			
		}
		else
		{
            if (rResult.size()!=(3*2+3))
			    rResult.resize(3*2+3);
			    
			unsigned int LocalIndex = 0;    
			for (unsigned int iNode = 0; iNode < 3; ++iNode)
			{
				rResult[LocalIndex++] = (GetGeometry()[iNode].GetDof(VELOCITY_X).EquationId());
				rResult[LocalIndex++] = (GetGeometry()[iNode].GetDof(VELOCITY_Y).EquationId());
			}
			rResult[LocalIndex++] = (this->GetValue(NEIGHBOUR_NODES)[0]).GetDof(LAGRANGE_MULTIPLIER_VELOCITY_GRADIENTS_X).EquationId();			
			rResult[LocalIndex++] = (this->GetValue(NEIGHBOUR_NODES)[0]).GetDof(LAGRANGE_MULTIPLIER_VELOCITY_GRADIENTS_Y).EquationId();			
			rResult[LocalIndex++] = (this->GetValue(NEIGHBOUR_NODES)[0]).GetDof(LAGRANGE_MULTIPLIER_VELOCITY_GRADIENTS_Z).EquationId();			

		}
	}

	//************************************************************************************
	//************************************************************************************
	  void VelocityGradientsLagrangeMultiplierCondition2D::GetDofList(DofsVectorType& ConditionalDofList,ProcessInfo& rCurrentProcessInfo)
	{
				unsigned int dim = 3;
		if (dim==2)
		{
            if (ConditionalDofList.size()!=4*dim)
			     ConditionalDofList.resize(4*dim);
			unsigned int LocalIndex = 0;    
			for (unsigned int iNode = 0; iNode < 3; ++iNode)
			{
				ConditionalDofList[LocalIndex++] = GetGeometry()[iNode].pGetDof(VELOCITY_X);
				ConditionalDofList[LocalIndex++] = GetGeometry()[iNode].pGetDof(VELOCITY_Y);
			}
			ConditionalDofList[LocalIndex++] = (this->GetValue(NEIGHBOUR_NODES)[0]).pGetDof(LAGRANGE_MULTIPLIER_VELOCITY_GRADIENTS_X);
			ConditionalDofList[LocalIndex++] = (this->GetValue(NEIGHBOUR_NODES)[0]).pGetDof(LAGRANGE_MULTIPLIER_VELOCITY_GRADIENTS_Y);
		}
		else
		{
            if (ConditionalDofList.size()!=3*2+3)
			     ConditionalDofList.resize(3*2+3);
			unsigned int LocalIndex = 0;    
			for (unsigned int iNode = 0; iNode < 3; ++iNode)
			{
				ConditionalDofList[LocalIndex++] = GetGeometry()[iNode].pGetDof(VELOCITY_X);
				ConditionalDofList[LocalIndex++] = GetGeometry()[iNode].pGetDof(VELOCITY_Y);
			}
			ConditionalDofList[LocalIndex++] = (this->GetValue(NEIGHBOUR_NODES)[0]).pGetDof(LAGRANGE_MULTIPLIER_VELOCITY_GRADIENTS_X);
			ConditionalDofList[LocalIndex++] = (this->GetValue(NEIGHBOUR_NODES)[0]).pGetDof(LAGRANGE_MULTIPLIER_VELOCITY_GRADIENTS_Y);
			ConditionalDofList[LocalIndex++] = (this->GetValue(NEIGHBOUR_NODES)[0]).pGetDof(LAGRANGE_MULTIPLIER_VELOCITY_GRADIENTS_Z);

		}
	}
} // Namespace Kratos
