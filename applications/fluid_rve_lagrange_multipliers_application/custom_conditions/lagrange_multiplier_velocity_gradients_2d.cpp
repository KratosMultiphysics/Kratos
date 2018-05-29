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
		const unsigned int dim = 3; // 2 = only 1 crossed gradients + du/dx . 3 = 2 crossed gradients (du/dy and dv/dx)  + du/dx
		if(dim==1) //(rCurrentProcessInfo[FRACTIONAL_STEP]==2)
		{
			if(rLeftHandSideMatrix.size1() != 2*3+1) // two dofs per node + 3 multipliers
				rLeftHandSideMatrix.resize(2*3+1,2*3+1,false);
			noalias(rLeftHandSideMatrix) = ZeroMatrix(2*3+1,2*3+1);
			if(rRightHandSideVector.size() != 2*3+1)
				rRightHandSideVector.resize(2*3+1,false);
			rRightHandSideVector = ZeroVector(2*3+1);
		
		
			const array_1d<double,3> & target_gradients = (this->GetValue(NEIGHBOUR_NODES)[0]).GetValue(LAGRANGE_MULTIPLIER_VELOCITY_GRADIENTS);

			double Area;
			array_1d<double, 3 > N;
			boost::numeric::ublas::bounded_matrix<double, 3, 2> DN_DX;
			boost::numeric::ublas::bounded_matrix<double, 1, 3*2+1 > LM_matrix = ZeroMatrix(1, 3*2+1) ; //(gradient)
			Geometry<Node<3> >& geom = this->GetGeometry();
			GeometryUtils::CalculateGeometryData(geom, DN_DX, N, Area);
			
			for (unsigned int j = 0; j < 3 ; j++) //3 nodes
			{
				LM_matrix(0, (j*2+0) ) =  DN_DX(j,1)*Area; // du/dy : we use DNi/DNy
				LM_matrix(0, (j*2+1) ) =  DN_DX(j,0)*Area; // dv/dx : we use DNi/DNx
				//LM_matrix(1, (j*2) ) =  DN_DX(j,0)*Area; // du/dx : we use DNi/DNx
				//LM_matrix(2, (j*2)+1 ) =  DN_DX(j,1)*Area; // dv/dy : we use DNi/DNy
			}
			
			
			KRATOS_WATCH(LM_matrix);
			
			for (unsigned int i = 0; i < (3)*2; i++) //velocity dofs
			{
				for (unsigned int j = 0; j < 1 ; j++) // LM	
				{
					rLeftHandSideMatrix(i,j+6) = LM_matrix(j,i);
					rLeftHandSideMatrix(j+6,i) = LM_matrix(j,i);
				}
			}
			
			rRightHandSideVector(6) = target_gradients[0]*Area;
			//rRightHandSideVector(7) = target_gradients[1]*Area;
			//rRightHandSideVector(8) = target_gradients[2]*Area;

			
			//substracting:
			array_1d<double,3*2+1>  temp_vector=ZeroVector(3*2+1);
			for (unsigned int j = 0; j < 3; j++) //3 nodes
			{
				const array_1d<double,3> & velocity = GetGeometry()[j].GetSolutionStepValue(VELOCITY,0);
				temp_vector[j*2]=velocity[0];
				temp_vector[j*2+1]=velocity[1];
			}
			const array_1d<double,3> & velocity_aux_node = (this->GetValue(NEIGHBOUR_NODES)[0]).GetSolutionStepValue(LAGRANGE_MULTIPLIER_VELOCITY_GRADIENTS,0);
			temp_vector[6]=velocity_aux_node[0];
			//temp_vector[7]=velocity_aux_node[1];
			//temp_vector[8]=velocity_aux_node[2];
 			noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix,temp_vector);
			
			//KRATOS_WATCH(rLeftHandSideMatrix);
			//KRATOS_WATCH(rRightHandSideVector);
			//KRATOS_WATCH( rLeftHandSideMatrix )

			
		}



		else if(dim==2) //(rCurrentProcessInfo[FRACTIONAL_STEP]==2)
		{
			if(rLeftHandSideMatrix.size1() != 2*3+2) // two dofs per node + 3 multipliers
				rLeftHandSideMatrix.resize(2*3+2,2*3+2,false);
			noalias(rLeftHandSideMatrix) = ZeroMatrix(2*3+2,2*3+2);
			if(rRightHandSideVector.size() != 2*3+2)
				rRightHandSideVector.resize(2*3+2,false);
			rRightHandSideVector = ZeroVector(2*3+2);
		
		
			const array_1d<double,3> & target_gradients = (this->GetValue(NEIGHBOUR_NODES)[0]).GetValue(LAGRANGE_MULTIPLIER_VELOCITY_GRADIENTS); //rCurrentProcessInfo[LAGRANGE_MULTIPLIER_VELOCITY_GRADIENTS];

			double Area;
			array_1d<double, 3 > N;
			boost::numeric::ublas::bounded_matrix<double, 3, 2> DN_DX;
			boost::numeric::ublas::bounded_matrix<double, 2, 3*2+2 > LM_matrix = ZeroMatrix(2, 3*2+2) ; //(gradient)
			Geometry<Node<3> >& geom = this->GetGeometry();
			GeometryUtils::CalculateGeometryData(geom, DN_DX, N, Area);
			
			for (unsigned int j = 0; j < 3 ; j++) //3 nodes
			{
				LM_matrix(0, (j*2+0) ) =  DN_DX(j,1)*Area; // du/dy : we use DNi/DNy
				LM_matrix(0, (j*2+1) ) =  DN_DX(j,0)*Area; // dv/dx : we use DNi/DNx
				LM_matrix(1, (j*2) ) =  DN_DX(j,0)*Area; // du/dx : we use DNi/DNx
				LM_matrix(1, (j*2)+1 ) =  DN_DX(j,1)*Area; // du/dx : we use DNi/DNx

				//LM_matrix(2, (j*2)+1 ) =  DN_DX(j,1)*Area; // dv/dy : we use DNi/DNy
			}
			
			
			//KRATOS_WATCH(LM_matrix);
			
			for (unsigned int i = 0; i < (3)*2; i++) //velocity dofs
			{
				for (unsigned int j = 0; j < 2 ; j++) // LM	
				{
					rLeftHandSideMatrix(i,j+6) = LM_matrix(j,i);
					rLeftHandSideMatrix(j+6,i) = LM_matrix(j,i);
				}
			}
			
			rRightHandSideVector(6) = target_gradients[0]*Area;
			rRightHandSideVector(7) = target_gradients[1]*Area;
			//rRightHandSideVector(8) = target_gradients[2]*Area;

			
			//substracting:
			array_1d<double,3*2+2>  temp_vector=ZeroVector(3*2+2);
			for (unsigned int j = 0; j < 3; j++) //3 nodes
			{
				const array_1d<double,3> & velocity = GetGeometry()[j].GetSolutionStepValue(VELOCITY,0);
				temp_vector[j*2]=velocity[0];
				temp_vector[j*2+1]=velocity[1];
			}
			const array_1d<double,3> & velocity_aux_node = (this->GetValue(NEIGHBOUR_NODES)[0]).GetSolutionStepValue(LAGRANGE_MULTIPLIER_VELOCITY_GRADIENTS,0);
			temp_vector[6]=velocity_aux_node[0];
			temp_vector[7]=velocity_aux_node[1];
			//temp_vector[8]=velocity_aux_node[2];
 			noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix,temp_vector);
			
			//KRATOS_WATCH(rLeftHandSideMatrix);
			//KRATOS_WATCH(rRightHandSideVector);
			//KRATOS_WATCH( rLeftHandSideMatrix )

			
		}
		else if(dim==3)
		{
			if(rLeftHandSideMatrix.size1() != 2*3+3) // two dofs per node + 3 multipliers
				rLeftHandSideMatrix.resize(2*3+3,2*3+3,false);
			noalias(rLeftHandSideMatrix) = ZeroMatrix(2*3+3,2*3+3);
			if(rRightHandSideVector.size() != 2*3+3)
				rRightHandSideVector.resize(2*3+3,false);
			rRightHandSideVector = ZeroVector(2*3+3);
		
		
			const array_1d<double,3> & target_gradients = (this->GetValue(NEIGHBOUR_NODES)[0]).GetValue(LAGRANGE_MULTIPLIER_VELOCITY_GRADIENTS);//rCurrentProcessInfo[LAGRANGE_MULTIPLIER_VELOCITY_GRADIENTS];

			double Area;
			array_1d<double, 3 > N;
			boost::numeric::ublas::bounded_matrix<double, 3, 2> DN_DX;
			boost::numeric::ublas::bounded_matrix<double, 3, 3*2+3 > LM_matrix = ZeroMatrix(3, 3*2+3) ; //(gradient)
			Geometry<Node<3> >& geom = this->GetGeometry();
			GeometryUtils::CalculateGeometryData(geom, DN_DX, N, Area);
			
			for (unsigned int j = 0; j < 3 ; j++) //3 nodes
			{
				LM_matrix(0, (j*2+0) ) =  DN_DX(j,1)*Area; //
				LM_matrix(0, (j*2+1) ) =  DN_DX(j,0)*Area; // //first row is the shear (Id)
				//LM_matrix(1, (j*2+0) ) =  DN_DX(j,1)*Area; // 
				LM_matrix(1, (j*2+0) ) = DN_DX(j,1)*Area; // //second column is the rotation
				LM_matrix(1, (j*2+1) ) = - DN_DX(j,0)*Area; // //second column is the rotation

				LM_matrix(2, (j*2+0) ) =  DN_DX(j,0)*Area; // dU/dx
				//LM_matrix(2, (j*2)+1 ) =  DN_DX(j,1)*Area; // //third column is the global incompressibility
			}
			
			
			//KRATOS_WATCH(LM_matrix);
			
			for (unsigned int i = 0; i < (3)*2; i++) //velocity dofs
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
			array_1d<double,3*2+3>  temp_vector=ZeroVector(3*2+3);
			for (unsigned int j = 0; j < 3; j++) //3 nodes
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
			//KRATOS_WATCH( rLeftHandSideMatrix )

			
		}
		else if(dim==4) //sending the strains to the RHS as a force : Sigma = div ( mu * symm_strains) 
		{
			if(rLeftHandSideMatrix.size1() != 2*3+3) // two dofs per node + 3 multipliers
				rLeftHandSideMatrix.resize(2*3+3,2*3+3,false);
			noalias(rLeftHandSideMatrix) = ZeroMatrix(2*3+3,2*3+3);
			if(rRightHandSideVector.size() != 2*3+3)
				rRightHandSideVector.resize(2*3+3,false);
			rRightHandSideVector = ZeroVector(2*3+3);
		
			Geometry<Node<3> >& geom = this->GetGeometry();
		
			const array_1d<double,3> & target_gradients = (this->GetValue(NEIGHBOUR_NODES)[0]).GetValue(LAGRANGE_MULTIPLIER_VELOCITY_GRADIENTS);//rCurrentProcessInfo[LAGRANGE_MULTIPLIER_VELOCITY_GRADIENTS];
			
			const double viscosity = 1.0/3.0* (  geom[0].FastGetSolutionStepValue(VISCOSITY,0) + geom[1].FastGetSolutionStepValue(VISCOSITY,0) + geom[2].FastGetSolutionStepValue(VISCOSITY,0) ); 
			array_1d<double,3> input_stresses = ZeroVector(3);
			input_stresses[0] = viscosity*(rCurrentProcessInfo[STRAIN])[0];
			input_stresses[1] = viscosity*(rCurrentProcessInfo[STRAIN])[1];
			input_stresses[2] = viscosity*(rCurrentProcessInfo[STRAIN])[2];
			

			double Area;
			array_1d<double, 3 > N;
			boost::numeric::ublas::bounded_matrix<double, 3, 2> DN_DX;
			boost::numeric::ublas::bounded_matrix<double, 3, 3*2+3 > LM_matrix = ZeroMatrix(3, 3*2+3) ; //(gradient)
			
			GeometryUtils::CalculateGeometryData(geom, DN_DX, N, Area);
			
			for (unsigned int j = 0; j < 3 ; j++) //3 nodes
			{
				LM_matrix(0, (j*2+0) ) =  DN_DX(j,1)*Area; //
				//LM_matrix(0, (j*2+1) ) =  DN_DX(j,0)*Area; // //first row is the shear (Id)
				//LM_matrix(1, (j*2+0) ) =  DN_DX(j,1)*Area; // 
				LM_matrix(1, (j*2+1) ) = DN_DX(j,0)*Area; // //second column is the rotation
				LM_matrix(2, (j*2+0) ) =  DN_DX(j,0)*Area; // dU/dx
				//LM_matrix(2, (j*2)+1 ) =  DN_DX(j,1)*Area; // //third column is the global incompressibility
			}
			
			
			//KRATOS_WATCH(LM_matrix);
			
			for (unsigned int i = 0; i < (3)*2; i++) //velocity dofs
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


			//adding the stresses on the RHS.
			for (unsigned int i = 0; i < 3; i++) //3 nodes
			{
				rRightHandSideVector(i*2+0) -= (DN_DX(i,0)*input_stresses(0)+DN_DX(i,1)*input_stresses(2))*Area;
				rRightHandSideVector(i*2+1) -= (DN_DX(i,1)*input_stresses(1)+DN_DX(i,0)*input_stresses(2))*Area;
			}
			
			//KRATOS_WATCH(input_stresses);
			
			//substracting:
			array_1d<double,3*2+3>  temp_vector=ZeroVector(3*2+3);
			for (unsigned int j = 0; j < 3; j++) //3 nodes
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
			//KRATOS_WATCH( rLeftHandSideMatrix )

			
		}
		else if(dim==5) // just saving the value of the gradients, but without touching the system at all.
		{
			if(rLeftHandSideMatrix.size1() != 2*3+3) // two dofs per node + 3 multipliers
				rLeftHandSideMatrix.resize(2*3+3,2*3+3,false);
			noalias(rLeftHandSideMatrix) = ZeroMatrix(2*3+3,2*3+3);
			if(rRightHandSideVector.size() != 2*3+3)
				rRightHandSideVector.resize(2*3+3,false);
			rRightHandSideVector = ZeroVector(2*3+3);
		
			double Area;
			array_1d<double, 3 > N;
			boost::numeric::ublas::bounded_matrix<double, 3, 2> DN_DX;
			boost::numeric::ublas::bounded_matrix<double, 3, 3*2+3 > LM_matrix = ZeroMatrix(3, 3*2+3) ; //(gradient)
			Geometry<Node<3> >& geom = this->GetGeometry();
			GeometryUtils::CalculateGeometryData(geom, DN_DX, N, Area);
			
			for (unsigned int j = 0; j < 3 ; j++) //3 nodes
			{
				LM_matrix(0, (j*2+0) ) =  DN_DX(j,1)*Area; //
				//LM_matrix(0, (j*2+1) ) =  DN_DX(j,0)*Area; // //first row is the shear (Id)
				//LM_matrix(1, (j*2+0) ) =  DN_DX(j,1)*Area; // 
				LM_matrix(1, (j*2+1) ) = DN_DX(j,0)*Area; // //second column is the rotation
				LM_matrix(2, (j*2+0) ) =  DN_DX(j,0)*Area; // dU/dx
				//LM_matrix(2, (j*2)+1 ) =  DN_DX(j,1)*Area; // //third column is the global incompressibility
			}
			
			
			//KRATOS_WATCH(LM_matrix);
			for (unsigned int j = 0; j < 3 ; j++) // LM	
			{
				for (unsigned int i = 0; i < (3)*2; i++) //velocity dofs
				{
					//rLeftHandSideMatrix(i,j+6) = LM_matrix(j,i); //we do not modify the original equations, we dont want to change the system, just read the value!
					rLeftHandSideMatrix(j+6,i) = LM_matrix(j,i);
				}
				rLeftHandSideMatrix(j+6,j+6) = Area;
			}
			
			
			rRightHandSideVector(6) = 0.0;
			rRightHandSideVector(7) = 0.0;
			rRightHandSideVector(8) = 0.0;
			
			//substracting:
			array_1d<double,3*2+3>  temp_vector=ZeroVector(3*2+3);
			for (unsigned int j = 0; j < 3; j++) //3 nodes
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
			//KRATOS_WATCH( rLeftHandSideMatrix )

			
		}
		KRATOS_CATCH("")
	}


	void VelocityGradientsLagrangeMultiplierCondition2D::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)

	{
		KRATOS_TRY

		//we only have to calculate if 
		
		const unsigned int dim = 3;
		if(dim == 1) //(rCurrentProcessInfo[FRACTIONAL_STEP]==2)
		{
			if(rLeftHandSideMatrix.size1() != 3*2+1)
				rLeftHandSideMatrix.resize(3*2+1,3*2+1,false);
			noalias(rLeftHandSideMatrix) = ZeroMatrix(3*2+1,3*2+1);
			if(rRightHandSideVector.size() != 3*2+1)
				rRightHandSideVector.resize(3*2+1,false);
			rRightHandSideVector = ZeroVector(3*2+1);
		}
		else if(dim == 2) //(rCurrentProcessInfo[FRACTIONAL_STEP]==2)
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
		unsigned int dim = 3;
		if (dim==1)
		{
            if (rResult.size()!=(3*2+1))
			    rResult.resize(3*2+1);
			    
			unsigned int LocalIndex = 0;    
			for (unsigned int iNode = 0; iNode < 3; ++iNode)
			{
				rResult[LocalIndex++] = (GetGeometry()[iNode].GetDof(VELOCITY_X).EquationId());
				rResult[LocalIndex++] = (GetGeometry()[iNode].GetDof(VELOCITY_Y).EquationId());
			}
			rResult[LocalIndex++] = (this->GetValue(NEIGHBOUR_NODES)[0]).GetDof(LAGRANGE_MULTIPLIER_VELOCITY_GRADIENTS_X).EquationId();			
			//rResult[LocalIndex++] = (this->GetValue(NEIGHBOUR_NODES)[0]).GetDof(LAGRANGE_MULTIPLIER_VELOCITY_GRADIENTS_Y).EquationId();			
		}
		else if (dim==2)
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
		if (dim==1)
		{
            if (ConditionalDofList.size()!=3*2+1)
			     ConditionalDofList.resize(3*2+1);
			unsigned int LocalIndex = 0;    
			for (unsigned int iNode = 0; iNode < 3; ++iNode)
			{
				ConditionalDofList[LocalIndex++] = GetGeometry()[iNode].pGetDof(VELOCITY_X);
				ConditionalDofList[LocalIndex++] = GetGeometry()[iNode].pGetDof(VELOCITY_Y);
			}
			ConditionalDofList[LocalIndex++] = (this->GetValue(NEIGHBOUR_NODES)[0]).pGetDof(LAGRANGE_MULTIPLIER_VELOCITY_GRADIENTS_X);
			//ConditionalDofList[LocalIndex++] = (this->GetValue(NEIGHBOUR_NODES)[0]).pGetDof(LAGRANGE_MULTIPLIER_VELOCITY_GRADIENTS_Y);
		}
		else if (dim==2)
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
