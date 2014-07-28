//   
//   Project Name:        Kratos       
//   Last modified by:    $Author: it's me! $
//   Date:                $Date: 2008-08-08 23:58:38 $
//   Revision:            $Revision: 1.0 $
//
//

// Project includes 
#include "includes/define.h"
#include "custom_elements/no_particles_solid_only_2d.h"
#include "custom_utilities/pfem_particle.h"
#include "pfem_2_application.h"
#include "utilities/math_utils.h"
#include "utilities/geometry_utilities.h" 
#include "utilities/enrichment_utilities.h"
#include "utilities/discont_utils.h"

#include "utilities/enrich_2d_2dofs.h"

#include "processes/process.h"
#include "includes/node.h"
//#include "includes/element.h"
#include "includes/model_part.h"



namespace Kratos
{
	


	//************************************************************************************
	//************************************************************************************
	NoParticlesSolidOnlyPFEM22D::NoParticlesSolidOnlyPFEM22D(IndexType NewId, GeometryType::Pointer pGeometry)
		: Element(NewId, pGeometry)
	{		
		//DO NOT ADD DOFS HERE!!!
	}

	//************************************************************************************
	//************************************************************************************
	NoParticlesSolidOnlyPFEM22D::NoParticlesSolidOnlyPFEM22D(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
		: Element(NewId, pGeometry, pProperties)
	{

	}

	Element::Pointer NoParticlesSolidOnlyPFEM22D::Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const
	{
		return Element::Pointer(new NoParticlesSolidOnlyPFEM22D(NewId, GetGeometry().Create(ThisNodes), pProperties));
	}

	NoParticlesSolidOnlyPFEM22D::~NoParticlesSolidOnlyPFEM22D()
	{
	}

	//************************************************************************************
	//************************************************************************************
	
	
	void NoParticlesSolidOnlyPFEM22D::AddExplicitContribution(ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY;

		switch ( rCurrentProcessInfo[FRACTIONAL_STEP] )
		{
			case 3: //calculating pressure projection. notthing returned. saving data in PRESS_PROJ, PRESS_PROJ_NO_RO , NODAL_MASS and NODAL_AREA
			{
				this->CalculatePressureProjection(rCurrentProcessInfo);
				break;
			}
			
			default:
			{
				KRATOS_ERROR(std::logic_error,"Unexpected value for FRACTIONAL_STEP index: ",rCurrentProcessInfo[FRACTIONAL_STEP]);
			}
		}

		KRATOS_CATCH("");
	}

	
	//************************************************************************************
	//************************************************************************************
	
	void NoParticlesSolidOnlyPFEM22D::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY
		array_1d<double,3> gravity= rCurrentProcessInfo[GRAVITY];

		const double pressure_increasing_factor=1.0e0;

		//WE MUST CALCULATE Mass and G(or D) matrixes

	    double delta_t = rCurrentProcessInfo[DELTA_TIME];	
		//const bool split_element = (this->GetValue(SPLIT_ELEMENT));

		double theta=1.0;
		
		//element mean stress, not used for calculations but rather only to reseed particles if needed.
		Vector & stresses = this->GetValue(ELEMENT_MEAN_STRESS);
		if (stresses.size()!=3)
		{
			stresses.resize(3);
			stresses=ZeroVector(3);
		}

		const unsigned int LocalSize = GetGeometry().size()*3; //2 vel + PRESSURE!
		unsigned int TNumNodes = GetGeometry().size();
		
        // Check sizes and initialize
        if (rRightHandSideVector.size() != LocalSize)
            rRightHandSideVector.resize(LocalSize, false);
       
        if(rLeftHandSideMatrix.size1() != LocalSize)
			rLeftHandSideMatrix.resize(LocalSize,LocalSize,false);
	    // Calculate this element's geometric parameters
		double Area;
		array_1d<double, 3> N; //dimension = number of nodes
		boost::numeric::ublas::bounded_matrix<double, 3, 2> DN_DX;
		GeometryUtils::CalculateGeometryData(this->GetGeometry(), DN_DX, N, Area);

		noalias(rRightHandSideVector) = ZeroVector(LocalSize);	
		noalias(rLeftHandSideMatrix) = ZeroMatrix(LocalSize,LocalSize);	
		
		
		//we start by restting the LHS and RHS

		
		array_1d<double,6>  previous_vel_in_mesh = ZeroVector(6); //to calculate the deformation Gradient F. Dimension = velocity dofs
		//array_1d<double,6>  previous_accel_in_mesh= ZeroVector(6); //to improve the computation of the displacements calculation, (used in the computation of the deformation gradient F too) Dimension = velocity dofs

		array_1d<double,9>  previous_vel_and_press = ZeroVector(9); //to add to the righthandside the unkwnonws from the previous time step. Dimension = total dofs
		//array_1d<double,9>  previous_accel= ZeroVector(9);          //same as above, but pressure spaces are left blank (unused). Dimension = total dofs
		
		for(unsigned int iii = 0; iii<3; iii++)
		{
			array_1d<double,3>& velocity =  GetGeometry()[iii].FastGetSolutionStepValue(VELOCITY);  //reference
			//array_1d<double,3>& accel =  GetGeometry()[iii].FastGetSolutionStepValue(ACCELERATION); //reference
			//array_1d<double,3> accel = (GetGeometry()[iii].FastGetSolutionStepValue(VELOCITY,1)- GetGeometry()[iii].GetSolutionStepValue(VELOCITY,2) )/delta_t;
			//velocity -= gravity*delta_t; //it was added in the particles, but since we need it here on the rhs instead of inside the velocity, we subtract it in order to add it correctly in the rhs.
			
			//saving everything
			previous_vel_in_mesh(iii*2) = velocity[0];
			previous_vel_in_mesh(iii*2+1) = velocity[1];
			//previous_accel_in_mesh(iii*2) = accel[0];
			//previous_accel_in_mesh(iii*2+1) = accel[1];
			//
			previous_vel_and_press(iii*3) = velocity[0];
			previous_vel_and_press(iii*3+1) = velocity[1];
			previous_vel_and_press(iii*3+2) = GetGeometry()[iii].FastGetSolutionStepValue(PRESSURE);
			//
			//previous_accel(iii*3+0) = accel[0];
			//previous_accel(iii*3+1) = accel[1]; 			

		}
		
		const double one_third = 1.0/3.0; //in 3d we would need one_quarter
	
		boost::numeric::ublas::bounded_matrix<double, 3, 3> Laplacian_matrix; //pressure dof^2
		noalias(Laplacian_matrix) = ZeroMatrix(3,3);	
		boost::numeric::ublas::bounded_matrix<double, 3, 6 > D_matrix; //(divergence)
		noalias(D_matrix) = ZeroMatrix(3,6);	
		boost::numeric::ublas::bounded_matrix<double, 3, 6 > G_matrix; //(gradient)
		noalias(G_matrix) = ZeroMatrix(3,6);	
		boost::numeric::ublas::bounded_matrix<double, 9, 9 > Mass_matrix; //2 vel + 1 pressure per node
		noalias(Mass_matrix) = ZeroMatrix(9,9);	
		//boost::numeric::ublas::bounded_matrix<double, 9, 9 > Pressure_Mass_matrix; //2 vel + 1 pressure per node
		//noalias(Pressure_Mass_matrix) = ZeroMatrix(9,9);	
		boost::numeric::ublas::bounded_matrix<double, 9, 9 > Lumped_Pressure_Mass_matrix; //2 vel + 1 pressure per node
		noalias(Lumped_Pressure_Mass_matrix) = ZeroMatrix(9,9);	
		
		//we start by calculating the non-enriched functions of divergence and gradient. as i write this, gradient integrated by parts.
		for (unsigned int i = 0; i < 3; i++)
		{
			for (unsigned int j = 0; j < 3; j++)
			{
				D_matrix(i, (j*2) ) =  DN_DX(j,0)*one_third;     
				D_matrix(i, (j*2+1) ) = DN_DX(j,1)*one_third;
				G_matrix(i, (j*2) ) =  -DN_DX(i,0)*one_third;     
				G_matrix(i, (j*2+1) ) = - DN_DX(i,1)*one_third;
			}
		}
		
		
		//area/volume integrals using the particles  in the elements:
		
		array_1d<double,3> total_integral_stresses = stresses*Area;
		double solid_area=Area;
		double density_solid_integral=this->GetValue(DENSITY)*Area;
		double density=this->GetValue(DENSITY);
		const double lambda = this->GetValue(YOUNG_MODULUS) * this->GetValue(POISSON_RATIO) / ( (1.0+this->GetValue(POISSON_RATIO))*(1.0-2.0*this->GetValue(POISSON_RATIO)) );	 
		const double mu = this->GetValue(YOUNG_MODULUS)/ (2.0*(1.0+this->GetValue(POISSON_RATIO))); 
		const double mu_integral = mu*solid_area;
		const double Bulk_modulus = (2.0/3.0 * mu + lambda);	
		
		//this->UpdateStressesToNewConfiguration(particle_stresses,particle_pressure,DN_DX,previous_vel_in_mesh,previous_accel_in_mesh, particle_mu , delta_t); //updating stresses
		//this->GetValue(ELEMENT_MEAN_STRESS)=stresses;

		double add_accel_stab_term=0.0; //(true or false)
		
		this->AddElasticityTerm(rLeftHandSideMatrix,DN_DX, (mu_integral*delta_t));

		
		for (unsigned int i=0; i<TNumNodes ; i++)
		{

			Mass_matrix (i*3,i*3) = density*Area*one_third; //nodal_mass_fluid(i)+nodal_mass_solid(i);
			Mass_matrix (i*3+1,i*3+1) = density*Area*one_third; //nodal_mass_fluid(i)+nodal_mass_solid(i);

			Lumped_Pressure_Mass_matrix(i*3+2,i*3+2) = one_third * solid_area /Bulk_modulus;
			//for (unsigned int j=0; j<TNumNodes ; j++)
			//	Pressure_Mass_matrix(i*3+2,j*3+2) = one_third * 0.25 * (Area /Bulk_modulus );
			//Pressure_Mass_matrix(i*3+2,i*3+2) = one_third * 0.5 * (Area /Bulk_modulus );
		}
			
		/*
		/////////////////// calculating tau
		double AdvVelNorm = sqrt(pow(this->GetGeometry()[0].FastGetSolutionStepValue(VELOCITY_X),2)+pow(this->GetGeometry()[0].FastGetSolutionStepValue(VELOCITY_Y),2));
		//const double TimeFactor = rCurrentProcessInfo.GetValue(DYNAMIC_TAU);
		double mElemSize;
		array_1d<double,3> Edge(3,0.0);
		Edge = this->GetGeometry()[1].Coordinates() - this->GetGeometry()[0].Coordinates();
		mElemSize = Edge[0]*Edge[0];
		for (SizeType d = 1; d < 2; d++)
			mElemSize += Edge[d]*Edge[d];
		for (SizeType i = 2; i < 3; i++)
			for(SizeType j = 0; j < i; j++)
			{
				Edge = this->GetGeometry()[i].Coordinates() - this->GetGeometry()[j].Coordinates();
				double Length = Edge[0]*Edge[0];
				for (SizeType d = 1; d < 2; d++)
					Length += Edge[d]*Edge[d];
				if (Length < mElemSize) mElemSize = Length;
			}
			
		mElemSize = sqrt(mElemSize);

		double apparent_viscosity = (viscosity * fluid_area +  0.0001*mu*solid_area*delta_t)/Area;
		if (viscosity>apparent_viscosity)
			apparent_viscosity=viscosity;
		const double TauOneFluid = (1.0 / ( ( 1.0 / delta_t + 4.0 * apparent_viscosity / (density * mElemSize * mElemSize) + 2.0 * AdvVelNorm / mElemSize) ) );
		
		double TauOneSolid = 0.0 /  (1.0 / delta_t +  10.0 * mu * delta_t / (density * mElemSize * mElemSize) + 2.0 * AdvVelNorm / mElemSize );
		double TauOne = (solid_area*TauOneSolid+fluid_area*TauOneFluid)/Area;
		*/
		const double TauOne=0.0;
		const double TauOneFluid=0.0;
		const double TauOneSolid=0.0;

			
		//this->GetValue(TAU)=TauOne;
		//this->GetValue(DENSITY)=density;
		//this->GetValue(VISCOSITY)=viscosity;
		
		Laplacian_matrix = prod(DN_DX,trans(DN_DX))/density;  // B B^T  (standard laplacian, now we must add the contribution of the condensed new nodes.
		for (unsigned int i = 0; i < 3; i++)
		{
			for (unsigned int j = 0; j < 3; j++)
			{
				
				rLeftHandSideMatrix(i*3+2, j*3+0 ) -= (D_matrix(i,j*2) + add_accel_stab_term *TauOne * DN_DX(i,0)*one_third/delta_t)*Area;     
				rLeftHandSideMatrix(i*3+2, j*3+1 ) -= (D_matrix(i,j*2+1) + add_accel_stab_term *TauOne * DN_DX(i,1)*one_third*Area/delta_t)*Area;    
				
				rLeftHandSideMatrix(j*3+0, i*3+2 ) -=  D_matrix(i,j*2)*Area;     
				rLeftHandSideMatrix(j*3+1, i*3+2 ) -= D_matrix(i,j*2+1)*Area; 
				//rLeftHandSideMatrix(i*3+2, j*3+2 ) -= Laplacian_matrix(i,j)*(fluid_area*TauOneFluid+solid_area*TauOneSolid);

				}
			}

		double mean_pressure= 0.0;
		for (unsigned int i = 0; i < 3; i++)
		{
			//mean_velocity += GetGeometry()[i].FastGetSolutionStepValue(VELOCITY)*one_third;
			//divergence_n += one_third*Area*(DN_DX(i,0)*GetGeometry()[i].FastGetSolutionStepValue(VELOCITY_X)+DN_DX(i,1)*GetGeometry()[i].FastGetSolutionStepValue(VELOCITY_Y));
			mean_pressure +=one_third*GetGeometry()[i].FastGetSolutionStepValue(PRESSURE);
		}	

		array_1d<double,2> extra_solid_forces =ZeroVector(2);
		for (unsigned int i = 0; i < 3; i++)
		{
			rRightHandSideVector(i*3+0) += gravity(0)*Mass_matrix (i*3,i*3);
			rRightHandSideVector(i*3+1) += gravity(1)*Mass_matrix (i*3+1,i*3+1);
				
			rRightHandSideVector(i*3+0) -= (DN_DX(i,0)*total_integral_stresses(0)+DN_DX(i,1)*total_integral_stresses(2));
			rRightHandSideVector(i*3+1) -= (DN_DX(i,1)*total_integral_stresses(1)+DN_DX(i,0)*total_integral_stresses(2));
				
			extra_solid_forces(0) += - (DN_DX(i,0)*total_integral_stresses(0)+DN_DX(i,1)*total_integral_stresses(2)) +  DN_DX(i,0) * (mean_pressure);
			extra_solid_forces(1) += - (DN_DX(i,1)*total_integral_stresses(1)+DN_DX(i,0)*total_integral_stresses(2)) +  DN_DX(i,1) * (mean_pressure);
		}
		//having all the forces, we calculate the mean x and y forces:
		/*
		array_1d<double, 2 > vel_gauss = ZeroVector(2);
		Vector pressure_rhs_stab(3);
		for (unsigned int i = 0; i < 3; i++)
		{
			array_1d<double,3>& node_press_proj = GetGeometry()[i].FastGetSolutionStepValue(PRESS_PROJ);
			vel_gauss[0] += node_press_proj(0)*one_third;
			vel_gauss[1] += node_press_proj(1)*one_third;
		}
		noalias(pressure_rhs_stab) = prod(DN_DX, vel_gauss);
		
		for (unsigned int i = 0; i < 3; i++) 
		{
			
			rRightHandSideVector(i*3+2) -= TauOneFluid*(DN_DX(i,0)*gravity(0)*fluid_area+DN_DX(i,1)*gravity(1)*fluid_area);
			rRightHandSideVector(i*3+2) -= TauOneSolid*(DN_DX(i,0)*(gravity(0)*solid_area)+DN_DX(i,1)*(gravity(1)*solid_area));
		}
		*/


		noalias(rRightHandSideVector) += prod((Mass_matrix),((1.0)*previous_vel_and_press/delta_t));
		noalias(rLeftHandSideMatrix) += (1.0)*Mass_matrix/delta_t;   
			
		Lumped_Pressure_Mass_matrix *= 1.0+TauOneSolid;	
	
		noalias(rRightHandSideVector) -= prod((Lumped_Pressure_Mass_matrix),(previous_vel_and_press/delta_t));
		noalias(rLeftHandSideMatrix) -= Lumped_Pressure_Mass_matrix/delta_t;  

		

		noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix,previous_vel_and_press);

		KRATOS_CATCH("");
	}
	
	
	
	//************************************************************************************
	//************************************************************************************
	// this subroutine calculates the nodal contributions for the explicit steps of the 
	// fractional step procedure
	void NoParticlesSolidOnlyPFEM22D::InitializeSolutionStep(ProcessInfo& CurrentProcessInfo)
	{
		KRATOS_TRY
		KRATOS_CATCH("");
	}
	
	//************************************************************************************
	//************************************************************************************


	//************************************************************************************
	//************************************************************************************
	void NoParticlesSolidOnlyPFEM22D::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo)
	{
		const SizeType NumNodes = 3;
		//const SizeType LocalSize = 6;
		const SizeType LocalSize = 9;
		GeometryType& rGeom = this->GetGeometry();

		SizeType LocalIndex = 0;

		if (rResult.size() != LocalSize)
			rResult.resize(LocalSize, false);

		for (SizeType i = 0; i < NumNodes; ++i)
		{
			rResult[LocalIndex++] = rGeom[i].GetDof(VELOCITY_X).EquationId();
			rResult[LocalIndex++] = rGeom[i].GetDof(VELOCITY_Y).EquationId();
			rResult[LocalIndex++] = rGeom[i].GetDof(PRESSURE).EquationId();
			//rResult[LocalIndex++] = rGeom[i].GetDof(SOLID_PRESSURE).EquationId();
		}
	}

	//************************************************************************************
	//************************************************************************************
	void NoParticlesSolidOnlyPFEM22D::GetDofList(DofsVectorType& ElementalDofList,ProcessInfo& CurrentProcessInfo)
	{
		const SizeType NumNodes = 3;
		//const SizeType LocalSize = 6;
		const SizeType LocalSize = 9;
		GeometryType& rGeom = this->GetGeometry();

		if (ElementalDofList.size() != LocalSize)
			ElementalDofList.resize(LocalSize);

		SizeType LocalIndex = 0;

		for (SizeType i = 0; i < NumNodes; ++i)
		{
			ElementalDofList[LocalIndex++] = rGeom[i].pGetDof(VELOCITY_X);
			ElementalDofList[LocalIndex++] = rGeom[i].pGetDof(VELOCITY_Y);
			ElementalDofList[LocalIndex++] = rGeom[i].pGetDof(PRESSURE);
			//ElementalDofList[LocalIndex++] = rGeom[i].pGetDof(SOLID_PRESSURE);
		}
	}
	
	//************************************************************************************
	
	//************************************************************************************
	
	//non partitioned elements using MatrixType
	void NoParticlesSolidOnlyPFEM22D::AddElasticityTerm(MatrixType& OutputMatrix,
                                       const boost::numeric::ublas::bounded_matrix<double, 3, 2>& rShapeDeriv,
                                       const double Weight)
	{

		
		//boost::numeric::ublas::bounded_matrix<double, 3, 3 > Laplacian_matrix =    prod(rShapeDeriv,trans(rShapeDeriv));
		boost::numeric::ublas::bounded_matrix<double, 6, 3 > B_matrix = ZeroMatrix(6,3);
		for (unsigned int i=0; i!=3; i++) //i node
		{
				B_matrix(i*2,0)=rShapeDeriv(i,0);
				B_matrix(i*2+1,1)=rShapeDeriv(i,1);
				
				B_matrix(i*2,2)=rShapeDeriv(i,1);
				B_matrix(i*2+1,2)=rShapeDeriv(i,0);
		}
		
		boost::numeric::ublas::bounded_matrix<double, 3, 3 > C_matrix = ZeroMatrix(3,3);
		
		C_matrix(0,0)=4.0/3.0;
		C_matrix(1,1)=4.0/3.0;
		C_matrix(2,2)=1.0;
		
		C_matrix(1,0)=-2.0/3.0;
		C_matrix(0,1)=-2.0/3.0;
		
		C_matrix *= Weight;
		
		boost::numeric::ublas::bounded_matrix<double, 3, 6 > temp_matrix = prod(C_matrix,trans(B_matrix));
		boost::numeric::ublas::bounded_matrix<double, 6, 6 > viscosity_matrix = prod(B_matrix, temp_matrix );

		
		for (unsigned int i=0; i!=3; i++) //i node
		{
			for (unsigned int j=0; j!=3; j++) //j neighbour
			{
					OutputMatrix(i*3+0,j*3+0)+=viscosity_matrix(i*2+0,j*2+0);   //0,0
					OutputMatrix(i*3+0,j*3+1)+=viscosity_matrix(i*2+0,j*2+1);   //0,1
					OutputMatrix(i*3+1,j*3+0)+=viscosity_matrix(i*2+1,j*2+0);  //1,0
					OutputMatrix(i*3+1,j*3+1)+=viscosity_matrix(i*2+1,j*2+1);      //1,1
			}
		}
		
		//OutputMatrix += viscosity_matrix;
	}
	

	
		//non partitioned elements using MatrixType
	void NoParticlesSolidOnlyPFEM22D::AddViscousTerm(MatrixType& OutputMatrix,
                                       const boost::numeric::ublas::bounded_matrix<double, 3, 2>& rShapeDeriv,
                                       const double Weight)
	{
		boost::numeric::ublas::bounded_matrix<double, 6, 3 > B_matrix = ZeroMatrix(6,3);
		for (unsigned int i=0; i!=3; i++) //i node
		{
				B_matrix(i*2,0)=rShapeDeriv(i,0);
				B_matrix(i*2+1,1)=rShapeDeriv(i,1);
				
				B_matrix(i*2,2)=rShapeDeriv(i,1);
				B_matrix(i*2+1,2)=rShapeDeriv(i,0);
		}
		
		boost::numeric::ublas::bounded_matrix<double, 3, 3 > C_matrix = ZeroMatrix(3,3);
		
		C_matrix(0,0)=2.0;
		C_matrix(1,1)=2.0;
		C_matrix(2,2)=1.0;

		C_matrix *= Weight;
		
		boost::numeric::ublas::bounded_matrix<double, 3, 6 > temp_matrix = prod(C_matrix,trans(B_matrix));
		boost::numeric::ublas::bounded_matrix<double, 6, 6 > viscosity_matrix = prod(B_matrix, temp_matrix );

		
		for (unsigned int i=0; i!=3; i++) //i node
		{
			for (unsigned int j=0; j!=3; j++) //j neighbour
			{
					//OutputMatrix(i*4+0,j*4+0)+=viscosity_matrix(i*2+0,j*2+0);   //0,0
					//OutputMatrix(i*4+0,j*4+1)+=viscosity_matrix(i*2+0,j*2+1);   //0,1
					//OutputMatrix(i*4+1,j*4+0)+=viscosity_matrix(i*2+1,j*2+0);  //1,0
					//OutputMatrix(i*4+1,j*4+1)+=viscosity_matrix(i*2+1,j*2+1);      //1,1
					OutputMatrix(i*3+0,j*3+0)+=viscosity_matrix(i*2+0,j*2+0);   //0,0
					OutputMatrix(i*3+0,j*3+1)+=viscosity_matrix(i*2+0,j*2+1);   //0,1
					OutputMatrix(i*3+1,j*3+0)+=viscosity_matrix(i*2+1,j*2+0);  //1,0
					OutputMatrix(i*3+1,j*3+1)+=viscosity_matrix(i*2+1,j*2+1);      //1,1
			}
		}
		
		//OutputMatrix += viscosity_matrix;
	}
	
		//non partitioned elements using MatrixType
	void NoParticlesSolidOnlyPFEM22D::UpdateStressesToNewConfiguration(array_1d<double,3>& OutputVector,
									   double& pressure,
                                       const boost::numeric::ublas::bounded_matrix<double, 3, 2>& rShapeDeriv,
                                       const array_1d<double,6>&  previous_step_vel,
                                       const array_1d<double,6>&  previous_step_accel,
                                       const double ShearModulus,
                                       const double delta_t)
	{
		boost::numeric::ublas::bounded_matrix<double, 6, 1 > velocities;
		for (unsigned int i=0;i!=6;i++)
			velocities(i,0)=previous_step_vel(i)+previous_step_accel(i)*delta_t*0.0;	 //hopefully this is something like V(n+1)
		
	
		
		/*	
		for (unsigned int i=0;i!=2;i++)
		{
			OutputVector(i) -=pressure;
		}	
		*/
		
		//having added the new component, we do sigma_new_config = F.S.F where F is the DeformationGradient
		boost::numeric::ublas::bounded_matrix<double, 3 , 2 > displacements;
		for (unsigned int i=0;i!=3;i++)
		{
			displacements(i,0)=velocities(i*2+0,0)*delta_t;
			displacements(i,1)=velocities(i*2+1,0)*delta_t;
		}	
		
		boost::numeric::ublas::bounded_matrix<double, 2 , 2 > DeformationGradient = identity_matrix<double> (2) ; //ZeroMatrix(2,2);// = prod(trans(displacements),rShapeDeriv);
		
		DeformationGradient+=prod(trans(displacements),rShapeDeriv);

        double J = DeformationGradient ( 0 , 0 ) *  DeformationGradient ( 1 , 1 ) - DeformationGradient ( 0 , 1 ) * DeformationGradient ( 1 , 0 ); //J is the determinant of the deformation gradient;
		
		if (J<=0.8 )
		{
			double amplifying_factor=1.0;
			double reducing_factor;
			while (J<=0.8 )
			{
				reducing_factor=1.0/(1.0+0.1*amplifying_factor);
				DeformationGradient = identity_matrix<double> (2) + reducing_factor*prod(trans(displacements),rShapeDeriv);
				J = DeformationGradient ( 0 , 0 ) *  DeformationGradient ( 1 , 1 ) - DeformationGradient ( 0 , 1 ) * DeformationGradient ( 1 , 0 );
				amplifying_factor++;
			}
			std::cout << "elem " << this->Id() << ", reduced timestep by " << reducing_factor << std::endl;
			
		}
		boost::numeric::ublas::bounded_matrix<double, 2 , 2 > tensorial_stresses;
		tensorial_stresses(0,0)=OutputVector(0);
		tensorial_stresses(1,1)=OutputVector(1);
		tensorial_stresses(1,0)=OutputVector(2);
		tensorial_stresses(0,1)=OutputVector(2);
		
		boost::numeric::ublas::bounded_matrix<double, 2 , 2 > temp_matrix = prod(tensorial_stresses, trans(DeformationGradient));
		tensorial_stresses = prod((DeformationGradient),temp_matrix);
		tensorial_stresses /=J;
		OutputVector(0)=tensorial_stresses(0,0);
		OutputVector(1)=tensorial_stresses(1,1);
		OutputVector(2)=tensorial_stresses(1,0);
		//pressure /= J ;
		
		/*
		const double residual_pressure =  - (tensorial_stresses(0,0)+tensorial_stresses(1,1))*0.5;
		pressure =  - (tensorial_stresses(0,0)+tensorial_stresses(1,1))*0.5;

		for (unsigned int i=0;i!=2;i++)
		{
			OutputVector(i) +=residual_pressure;
		}
		 */
		
		//KRATOS_WATCH(OutputVector)
		
		//KRATOS_WATCH(delta_strains)
		//KRATOS_WATCH(delta_stresses)
		/*
		for (unsigned int i=0;i!=3;i++)
			OutputVector(i)+=delta_strains(i,0);
		*/

	}
	
	
	//using MatriType (for the implicit step)
	void NoParticlesSolidOnlyPFEM22D::AddViscousTerm(MatrixType& rDampMatrix,
                                       const boost::numeric::ublas::bounded_matrix<double, 3, 2>& rShapeDeriv,
                                       const double viscosity_air,
                                       const double viscosity_water,
                                       array_1d<double,3>  volumes,
                                       array_1d<double,3>  signs,
                                       Matrix Ngauss,
                                       const double Area)
	{
		
		double Weight=0.0;
		for (unsigned int i=0;i<3;i++) //we go over the three partitions;
		{
			if (signs(i)>0.0)
				Weight += volumes(i)*viscosity_air;
			else
				Weight += volumes(i)*viscosity_water;
		}
		
		boost::numeric::ublas::bounded_matrix<double, 6, 3 > B_matrix = ZeroMatrix(6,3);
		for (unsigned int i=0; i!=3; i++) //i node
		{
				B_matrix(i*2,0)=rShapeDeriv(i,0);
				B_matrix(i*2+1,1)=rShapeDeriv(i,1);
				
				B_matrix(i*2,2)=rShapeDeriv(i,1);
				B_matrix(i*2+1,2)=rShapeDeriv(i,0);
		}
		
		boost::numeric::ublas::bounded_matrix<double, 3, 3 > C_matrix = ZeroMatrix(3,3);
		C_matrix(0,0)=2.0;
		C_matrix(1,1)=2.0;
		C_matrix(2,2)=1.0;
		
		C_matrix*=Weight;
		
		boost::numeric::ublas::bounded_matrix<double, 3, 6 > temp_matrix = prod(C_matrix,trans(B_matrix));
		rDampMatrix = prod(B_matrix, temp_matrix );
	}
	
	
	inline bool NoParticlesSolidOnlyPFEM22D::CalculatePosition(const bounded_matrix<double, 3, 3 > & coordinates,
                                         const double xc, const double yc, const double zc,
                                         array_1d<double, 3 > & N
                                        )
    {
        double x0 = coordinates(0,0);
        double y0 = coordinates(0,1);
        double x1 = coordinates(1,0);
        double y1 = coordinates(1,1);
        double x2 = coordinates(2,0);
        double y2 = coordinates(2,1);

        double area = CalculateVol(x0, y0, x1, y1, x2, y2);
        double inv_area = 0.0;
        if (area == 0.0)
        {
            KRATOS_ERROR(std::logic_error, "element with zero area found", "");
        }
        else
        {
            inv_area = 1.0 / area;
        }


        N[0] = CalculateVol(x1, y1, x2, y2, xc, yc) * inv_area;
        N[1] = CalculateVol(x2, y2, x0, y0, xc, yc) * inv_area;
        N[2] = CalculateVol(x0, y0, x1, y1, xc, yc) * inv_area;
        //KRATOS_WATCH(N);
        
        if (N[0] >= 0.0 && N[1] >= 0.0 && N[2] >= 0.0 && N[0] <= 1.0 && N[1] <= 1.0 && N[2] <= 1.0) //if the xc yc is inside the triangle return true
            return true;

        return false;
        
    }

    inline bool NoParticlesSolidOnlyPFEM22D::CalculatePosition(Geometry<Node < 3 > >&geom,
                const double xc, const double yc, const double zc,
                array_1d<double, 3 > & N
                )
        {
            double x0 = geom[0].X();
            double y0 = geom[0].Y();
            double x1 = geom[1].X();
            double y1 = geom[1].Y();
            double x2 = geom[2].X();
            double y2 = geom[2].Y();

            double area = CalculateVol(x0, y0, x1, y1, x2, y2);
            double inv_area = 0.0;
            if (area == 0.0)
            {
                KRATOS_ERROR(std::logic_error, "element with zero area found", "");
            } else
            {
                inv_area = 1.0 / area;
            }


            N[0] = CalculateVol(x1, y1, x2, y2, xc, yc) * inv_area;
            N[1] = CalculateVol(x2, y2, x0, y0, xc, yc) * inv_area;
            N[2] = CalculateVol(x0, y0, x1, y1, xc, yc) * inv_area;
			//KRATOS_WATCH(N);

            if (N[0] >= 0.0 && N[1] >= 0.0 && N[2] >= 0.0 && N[0] <= 1.0 && N[1] <= 1.0 && N[2] <= 1.0) //if the xc yc is inside the triangle return true
                return true;

            return false;
        }


    inline double NoParticlesSolidOnlyPFEM22D::CalculateVol(const double x0, const double y0,
                                      const double x1, const double y1,
                                      const double x2, const double y2
                                     )
    {
        return 0.5 * ((x1 - x0)*(y2 - y0)- (y1 - y0)*(x2 - x0));
    }
    
    void NoParticlesSolidOnlyPFEM22D::GetValueOnIntegrationPoints(const Variable<double>& rVariable,
            std::vector<double>& rValues,
            const ProcessInfo& rCurrentProcessInfo)
	{
		if (rVariable == YP)
		{
			// Set output vector (for a single integration point)
			rValues.resize(1, false);
			rValues[0]=double( this->GetValue(NUMBER_OF_PARTICLES) );
		}
		else if (rVariable == TEMPERATURE)
		{
			// Set output vector (for a single integration point)
			rValues.resize(1, false);
			rValues[0]=this->GetValue(TEMPERATURE);
		}
		else if (rVariable == DENSITY)
		{
			// Set output vector (for a single integration point)
			rValues.resize(1, false);
			rValues[0]=this->GetValue(DENSITY);
		}
		else if (rVariable == VISCOSITY)
		{
			// Set output vector (for a single integration point)
			rValues.resize(1, false);
			rValues[0]=this->GetValue(VISCOSITY);
		}
		else if (rVariable == TAU)
		{
			// Set output vector (for a single integration point)
			rValues.resize(1, false);
			rValues[0]=this->GetValue(TAU);
		}
		
		else // Default behaviour (returns elemental data)
		{
			rValues.resize(1, false);
			//const VMS<Dim,NumNodes>* const_this = static_cast< const VMS<Dim,NumNodes>* >(this);
			//rOutput[0] = const_this->GetValue(rVariable);
		}
		
	}
	
	   void NoParticlesSolidOnlyPFEM22D::GetValueOnIntegrationPoints(const Variable<array_1d<double, 3 > >& rVariable,
            std::vector<array_1d<double, 3 > >& rValues,
            const ProcessInfo& rCurrentProcessInfo)
	{
		if (rVariable == VELOCITIES)
		{
			// Set output vector (for a single integration point)
			rValues.resize(1);
			rValues[0]=this->GetValue(ELEMENT_MEAN_STRESS);
		}
		else // Default behaviour (returns elemental data)
		{
			rValues.resize(1);
			//const VMS<Dim,NumNodes>* const_this = static_cast< const VMS<Dim,NumNodes>* >(this);
			//rOutput[0] = const_this->GetValue(rVariable);
		}
		
	}

	template<class T>
	bool NoParticlesSolidOnlyPFEM22D::InvertMatrix(const T& input, T& inverse)
	{
		typedef permutation_matrix<std::size_t> pmatrix;

		// create a working copy of the input
		T A(input);

		// create a permutation matrix for the LU-factorization
		pmatrix pm(A.size1());

		// perform LU-factorization
		int res = lu_factorize(A, pm);
		if (res != 0)
			return false;

		// create identity matrix of "inverse"
		inverse.assign(identity_matrix<double> (A.size1()));

		// backsubstitute to get the inverse
		lu_substitute(A, pm, inverse);

		return true;
	}
	
	
	void NoParticlesSolidOnlyPFEM22D::CalculateInterfaceNormal(boost::numeric::ublas::bounded_matrix<double, 3, 2 >& rPoints, array_1d<double,3>&  rDistances, array_1d<double,2>&  normal, double & interface_area, array_1d<double,3>&  Ninterface)
	{
		double sign_correction=1.0;
		
		boost::numeric::ublas::bounded_matrix<double, 2, 2 > InterfacePoints;
		array_1d<bool,3>  cut_edges;
		array_1d<double,2>  interface_segment=ZeroVector(2);
		if ((rDistances(0)*rDistances(1))<0.0) cut_edges[0]=true;//edge 12 is cut	
		else 	cut_edges[0]=false;
			
		if ((rDistances(1)*rDistances(2))<0.0) cut_edges[1]=true;//edge 23 is cut. 
		else 	cut_edges[1]=false;
			
		if ((rDistances(2)*rDistances(0))<0.0) cut_edges[2]=true;//edge 13 is cut. 
		else 	cut_edges[2]=false;
					
					
		if (cut_edges[0])
		{
			if (rDistances(0)>0.0) sign_correction=1.0;
			else sign_correction=-1.0;
			
			const double relative_position = fabs(rDistances(1)/(rDistances(1)-rDistances(0) ) ); 
			InterfacePoints(0,0) = relative_position*rPoints(0,0) +  (1.0-relative_position)*rPoints(1,0);
			InterfacePoints(0,1) = relative_position*rPoints(0,1) +  (1.0-relative_position)*rPoints(1,1);
			
			if (cut_edges[1])
			{
				const double relative_position2 = fabs(rDistances(2)/(rDistances(1)-rDistances(2) ) ); 
				InterfacePoints(1,0) = relative_position2*rPoints(1,0) +  (1.0-relative_position2)*rPoints(2,0);
				InterfacePoints(1,1) = relative_position2*rPoints(1,1) +  (1.0-relative_position2)*rPoints(2,1);
			}
			else
			{
				const double relative_position2 = fabs(rDistances(0)/(rDistances(2)-rDistances(0) ) ); 
				InterfacePoints(1,0) = relative_position2*rPoints(2,0) +  (1.0-relative_position2)*rPoints(0,0);
				InterfacePoints(1,1) = relative_position2*rPoints(2,1) +  (1.0-relative_position2)*rPoints(0,1);
			}
		}
		else
		{
			if (rDistances(1)>0.0) sign_correction=1.0;
			else sign_correction=-1.0;
			
			const double relative_position = fabs(rDistances(2)/(rDistances(2)-rDistances(1) ) ); 
			InterfacePoints(0,0) = relative_position*rPoints(1,0) +  (1.0-relative_position)*rPoints(2,0);
			InterfacePoints(0,1) = relative_position*rPoints(1,1) +  (1.0-relative_position)*rPoints(2,1);
			
			const double relative_position2 = fabs(rDistances(0)/(rDistances(2)-rDistances(0) ) ); 
			InterfacePoints(1,0) = relative_position2*rPoints(2,0) +  (1.0-relative_position2)*rPoints(0,0);
			InterfacePoints(1,1) = relative_position2*rPoints(2,1) +  (1.0-relative_position2)*rPoints(0,1);
		}	
		
		interface_segment[0] = (InterfacePoints(1,0)-InterfacePoints(0,0));
		interface_segment[1] = (InterfacePoints(1,1)-InterfacePoints(0,1));	
		
		const double norm = sqrt(  pow((interface_segment[0]),2) + pow((interface_segment[1]),2));
		
		normal(0)= -interface_segment[1]*sign_correction/norm;
		normal(1)= interface_segment[0]*sign_correction/norm;
		//KRATOS_WATCH(interface_segment)
		//KRATOS_WATCH(InterfacePoints)
		interface_area=norm;
		
		const double x_interface = 0.5*(InterfacePoints(0,0)+InterfacePoints(1,0));
		const double y_interface = 0.5*(InterfacePoints(0,1)+InterfacePoints(1,1));
		
		CalculatePosition(rPoints, x_interface, y_interface, 0.0, Ninterface );
		
	}
	
	
	void NoParticlesSolidOnlyPFEM22D::CalculatePressureProjection(ProcessInfo& CurrentProcessInfo)
	{
		KRATOS_TRY
		
		const int number_of_particles_in_elem = this->GetValue(NUMBER_OF_PARTICLES);
		
		if( (number_of_particles_in_elem>0))
		//if ((this->GetValue(IS_INACTIVE))==false) //elements can be inactive to add temporary walls. fractional velocity is integrated by parts so walls are seen as having zero velocity without doing anything 
				{
					//const double delta_t = CurrentProcessInfo[DELTA_TIME];
					const double mDENSITY_AIR = CurrentProcessInfo[DENSITY_AIR]; // * (1.0-theta) ;
					const double mDENSITY_WATER = CurrentProcessInfo[DENSITY_WATER]; // * (1.0-theta) ;
					const double mINV_DENSITY_AIR = 1.0/mDENSITY_AIR;
					const double mINV_DENSITY_WATER = 1.0/mDENSITY_WATER;
					
					const double delta_t = CurrentProcessInfo[DELTA_TIME];
					//array_1d<double,3> & gravity= CurrentProcessInfo[GRAVITY];
					//const double x_force  = gravity(0);
					//const double y_force  = gravity(1);	
					const array_1d<double,3> zero3 = ZeroVector(3);
					const double mass_factor = 1.0/ (1.0 + double (2) );
					
					boost::numeric::ublas::bounded_matrix<double, (2+1), 2 > DN_DX;
					array_1d<double, (2+1) > N;
					array_1d<double, 2 > vel_gauss;
					
					double Area;
					Geometry<Node<3> >& geom = this->GetGeometry();
					//const unsigned int LocalSize = geom.size()*2;
					unsigned int TNumNodes = geom.size();
					//boost::numeric::ublas::bounded_matrix<double, (2+1),1 > enrich_rhs;
					GeometryUtils::CalculateGeometryData(geom, DN_DX, N, Area);
					
					array_1d<double,(2+1)> nodal_masses = ZeroVector(2+1);

					array_1d<double,(2+1)> pressure;
					for (unsigned int i=0; i<(2+1);i++)
						pressure(i) = geom[i].FastGetSolutionStepValue(PRESSURE);
						
					
					boost::numeric::ublas::bounded_matrix<double, (2+1), (2-1)*6 > G_matrix; //(gradient)
					noalias(G_matrix) = ZeroMatrix((2+1), (2-1)*6);	
					//new (upper) part of D, corresponding the first enrichment function (gradient discontinuity)
					boost::numeric::ublas::bounded_matrix<double, 1, (2-1)*6 > G_matrix_mixed;
					noalias(G_matrix_mixed) = ZeroMatrix(1,(2-1)*6);	
					//lower part of D, the one corresponding to the second enrichmend function (jump)
					boost::numeric::ublas::bounded_matrix<double, 1, (2-1)*6 > G_matrix_mixed_jump;
					noalias(G_matrix_mixed_jump) = ZeroMatrix(1,(2-1)*6);	
					
					boost::numeric::ublas::bounded_matrix<double, (2+1), (2-1)*6 > G_matrix_no_ro; //(gradient)
					noalias(G_matrix_no_ro) = ZeroMatrix((2+1),(2-1)*6);	
					//new (upper) part of D, corresponding the first enrichment function (gradient discontinuity
					boost::numeric::ublas::bounded_matrix<double, 1, (2-1)*6 > G_matrix_mixed_no_ro;
					noalias(G_matrix_mixed_no_ro) = ZeroMatrix(1,(2-1)*6);	
					//lower part of D, the one corresponding to the second enrichmend function (jump)
					boost::numeric::ublas::bounded_matrix<double, 1, (2-1)*6> G_matrix_mixed_jump_no_ro;
					noalias(G_matrix_mixed_jump_no_ro) = ZeroMatrix(1,(2-1)*6);		
					
					
					//for the enrichment:
					array_1d<double,(2+1)> distances;
					boost::numeric::ublas::bounded_matrix<double,3*(2-1), 2> Nenriched;
					array_1d<double,(3*(2-1))> volumes;
					array_1d<double,(3*(2-1))> densities;
					array_1d<double,(3*(2-1))> inv_densities;
					boost::numeric::ublas::bounded_matrix<double,(2+1), 2 > coords;
					boost::numeric::ublas::bounded_matrix<double, 3*(2-1), (2+1) > Ngauss;
					array_1d<double,(3*(2-1))> signs;
					std::vector< Matrix > gauss_gradients(((2-1)*3));
					//fill coordinates
				   
					unsigned int single_triangle_node;

					
			
						
						
					
					
					//if ((this->GetValue(SPLIT_ELEMENT))==true)
					if (false==true)
					{
						//to begin with we calculate the gradient discontinuity:
						const int iteration_number =  CurrentProcessInfo[NL_ITERATION_NUMBER];
						if (true) 
						{
							
							double gradient_discontinuity = this->GetValue(ENRICH_RHS) ;
							double gradient_discontinuity_complete_rhs = this->GetValue(GRADIENT_DISCONTINUITY)/((this->GetValue(INV_LAPLACIAN_ENRICH)));
							
							array_1d<double,(2+1)> enrich_lhs;
							enrich_lhs = this->GetValue(ENRICH_LHS_ROW);

							array_1d<double,(2+1)> pressure_array;
							array_1d<double,(2+1)> previous_iteration_pressure_array;
							for (unsigned int i=0; i<(2+1);i++)
							{
								pressure_array(i) = geom[i].FastGetSolutionStepValue(PRESSURE);//-geom[i].GetSolutionStepValue(PRESSURE,1);
								previous_iteration_pressure_array(i) = geom[i].FastGetSolutionStepValue(PREVIOUS_ITERATION_PRESSURE);
							}
							for (unsigned int i = 0; i < (2+1); i++)  // node i
							{ 
								gradient_discontinuity -= pressure_array(i)*enrich_lhs(i);
								if (iteration_number>1)
									gradient_discontinuity_complete_rhs += previous_iteration_pressure_array(i)*enrich_lhs(i);
							}
							//gradient_discontinuity_complete_rhs *= delta_t;
							gradient_discontinuity += gradient_discontinuity_complete_rhs;
							gradient_discontinuity *= this->GetValue(INV_LAPLACIAN_ENRICH);
							this->GetValue(GRADIENT_DISCONTINUITY) = gradient_discontinuity;
						}
						
						
						
						double input_jump_value = this->GetValue(PRESS_DISCONTINUITY);
						double jump_value=-input_jump_value;
						//if (iteration_number>1) jump_value *=double(iteration_number);
						double gradient_discontinuity; // = this->GetValue(ENRICH_RHS);
						//array_1d<double,(2+1)> & enrich_lhs = this->GetValue(ENRICH_LHS_ROW);
						
						/*
						noalias(G_matrix) = ZeroMatrix((2+1), (2-1)*6);	
						noalias(G_matrix_mixed) = ZeroMatrix(1,(2-1)*6);	
						noalias(G_matrix_mixed_jump) = ZeroMatrix(1,(2-1)*6);	
						*/
						noalias(G_matrix_no_ro) = ZeroMatrix((2+1),6*(2-1)*6);	
						noalias(G_matrix_mixed_no_ro) = ZeroMatrix(1,(2-1)*6);	
						noalias(G_matrix_mixed_jump_no_ro) = ZeroMatrix(1,(2-1)*6);	
						
						
						
						//TOOLS NEEDED TO FIND THE ENRICHMENT SHAPE FUNCTIONS
						//get position of the cut surface

						for (unsigned int i = 0; i < TNumNodes; i++)
						{
							const array_1d<double, 3 > & xyz = geom[i].Coordinates();
							distances[i] = geom[i].FastGetSolutionStepValue(DISTANCE);
							//KRATOS_WATCH(distances(i));
							for (unsigned int j = 0; j < (2); j++)
								coords(i, j) = xyz[j];
						}

						for (unsigned int i = 0; i < ((2-1)*3) ; i++)
							gauss_gradients[i].resize(2, 2, false);  //2 values of the 2 shape functions, and derivates in (xy) direction).
						unsigned int ndivisions;
						ndivisions = EnrichmentUtilities::CalculateEnrichedShapeFuncions(coords, DN_DX, distances, volumes, Ngauss, signs, gauss_gradients, Nenriched); //, face_gauss_N, face_gauss_Nenriched, face_Area, face_n  ,single_triangle_node);


						//in case we've changed the distance function:
						for (unsigned int i = 0; i < ndivisions; i++)
						{
							//geom[i].FastGetSolutionStepValue(DISTANCE)=distances[i];
							
							if (signs[i]>0)
							{
								densities(i) = mDENSITY_AIR;
								inv_densities(i) = mINV_DENSITY_AIR;
							}
							else
							{
								densities(i) = mDENSITY_WATER;
								inv_densities(i) = mINV_DENSITY_WATER;
							}
						}

						gradient_discontinuity = this->GetValue(GRADIENT_DISCONTINUITY);

						///*****************************************************************************************
						for (unsigned int i = 0; i < ndivisions; i++)  // partiton i
						{ 
							for (unsigned int j = 0; j < (2+1); j++) //shape function (Nj)
							{
								for (unsigned int k=0;k!=2;k++) //x,y,(z)
								{
									//G_matrix_mixed(0,j*(2)+k) += gauss_gradients[i](0,k) *volumes(i)*Ngauss(i,j)*inv_densities(i);
									G_matrix_mixed_no_ro(0,j*(2)+k) += gauss_gradients[i](0,k) *volumes(i)*Ngauss(i,j);
									

									//G_matrix_mixed_jump(0,j*(2)+k) += gauss_gradients[i](1,k) *volumes(i)*Ngauss(i,j)*inv_densities(i);
									G_matrix_mixed_jump_no_ro(0,j*(2)+k) += gauss_gradients[i](1,k) *volumes(i)*Ngauss(i,j);

									
									//for (unsigned int l = 0; l < (2+1); l++) //shape function derivative (dNl/dk)
									//{
									//	G_matrix(l, (j*(2)+k) ) += DN_DX(l,k)*volumes(i)*Ngauss(i,j)*inv_densities(i);
									//}
									
								}
								
								nodal_masses(j) += volumes(i)*Ngauss(i,j)*densities(i);
							}
						}
						//for the massless g matrix (to calculate the press proj for the stabilization) we do not need to loop the partitions, so we create a new loop
						for (unsigned int j = 0; j < (2+1); j++) //shape function (Nj)
						{
							for (unsigned int k=0;k!=2;k++) //x,y,(z)
							{
								for (unsigned int l = 0; l < (2+1); l++) //shape function derivative (dNl/dk)
								{
									G_matrix_no_ro(l, (j*(2)+k) ) = DN_DX(l,k)*Area*mass_factor;
								}
							}
						}
								
						//gradient_discontinuity = this->GetValue(GRADIENT_DISCONTINUITY);
						
						const double reduction_factor=0.00000001; //the splitted elements should not add anything to the stabilization terms to avoid problems. we keep a small value to avoid problems due to zero nodal area.
						//now we save the data:
						for (unsigned int i=0; i!=(2+1); i++) //loop around the nodes of the element to add contribution to node i
						{
							
							array_1d<double, 3 > & current_press_proj_no_ro = geom[i].FastGetSolutionStepValue(PRESS_PROJ_NO_RO);
							//array_1d<double, 3 > & current_press_proj = geom[i].FastGetSolutionStepValue(PRESS_PROJ);
								
							for (unsigned int j = 0; j < (2+1) ; j++) //loop around the nodes of the element to add contribution of node j to node i
							{
								
								for (unsigned int k = 0; k < (2) ; k++) //x,y,(z)
								{
									//current_press_proj[k] += reduction_factor*G_matrix(j,i*(2)+k)*(pressure(j));///-old_pressures(j)); //gamma=0!
									current_press_proj_no_ro[k] += G_matrix_no_ro(j,i*(2)+k)*(pressure(j));///-old_pressures(j)); //gamma=0!
								}
							}
							
							for (unsigned int k = 0; k < (2) ; k++) //x,y,(z)
							{
								//current_press_proj[k] += reduction_factor* (G_matrix_mixed(0,2*i+k)*gradient_discontinuity + G_matrix_mixed_jump(0,2*i+k)*jump_value);
								current_press_proj_no_ro[k] += G_matrix_mixed_no_ro(0,2*i+k)*gradient_discontinuity + G_matrix_mixed_jump_no_ro(0,2*i+k)*jump_value;
							}
							
							geom[i].FastGetSolutionStepValue(NODAL_AREA) += Area*mass_factor*reduction_factor;
							geom[i].FastGetSolutionStepValue(NODAL_MASS) += nodal_masses(i)/delta_t;
							
						}
						
					}
					else //normal (uncut) element
					{
						noalias(G_matrix_no_ro) = ZeroMatrix((2+1),(2-1)*6);
						//noalias(G_matrix) = ZeroMatrix((2+1),(2-1)*6);
						
						double density;
						if (geom[0].FastGetSolutionStepValue(DISTANCE)<0.0)
							density = mDENSITY_WATER;
						else
							density = mDENSITY_AIR;
							

						
						for (unsigned int i = 0; i < (2+1); i++) //loop in shape functions (row)
						{
							for (unsigned int j = 0; j < (2+1) ; j++) //loop in shape functions (column)
							{
								for (unsigned int k = 0; k < (2) ; k++) //x,y,(z)
								{
								G_matrix(i, (j*2)+k ) = DN_DX(i,k)*Area*mass_factor;
								//G_matrix_no_ro(i, (j*2)+k ) = DN_DX(i,k)*Area*mass_factor;
								}
							}
						}
						
						G_matrix /=density;

						//now we save the data:
						for (unsigned int i=0; i!=(2+1); i++) //loop around the nodes of the element to add contribution to node i
						{
							
							//array_1d<double, 3 > & current_press_proj_no_ro = geom[i].FastGetSolutionStepValue(PRESS_PROJ_NO_RO);
							array_1d<double, 3 > & current_press_proj = geom[i].FastGetSolutionStepValue(PRESS_PROJ);
								
							for (unsigned int j = 0; j < (2+1) ; j++) //loop around the nodes of the element to add contribution of node j to node i
							{
								
								for (unsigned int k = 0; k < (2) ; k++) //x,y,(z)
								{
									current_press_proj[k] += G_matrix(j,i*(2)+k)*(pressure(j));///-old_pressures(j)); //gamma=0!
									//current_press_proj_no_ro[k] += G_matrix_no_ro(j,i*(2)+k)*(pressure(j));///-old_pressures(j)); //gamma=0!
								}
							}
							
							geom[i].FastGetSolutionStepValue(NODAL_AREA) += Area*mass_factor;
							//geom[i].FastGetSolutionStepValue(NODAL_MASS) += Area*mass_factor*density/delta_t;
							
						}

					} //closing the normal (uncut) element
					
				} //closing the if(is_inactive==false)

		
		KRATOS_CATCH("");
	}
	
	
	void NoParticlesSolidOnlyPFEM22D::Calculate(const Variable<Vector> &rVariable,
                                     Vector &rOutput,
                                     const ProcessInfo &rCurrentProcessInfo)
	{
		
		
		if (rVariable == ELEMENT_MEAN_STRESS)
		{
				const bool use_failure_criteria=true;
			
				double delta_t = rCurrentProcessInfo[DELTA_TIME];	
	
			//we will add the delta stresses to each of the particles inside the elements.
				Geometry<Node<3> >& geom = this->GetGeometry(); 
				
				
				//std::cout << "elem " << this->Id() << " with " << (unsigned int)number_of_particles_in_elem << " particles" << std::endl;
				////////////////////////
				boost::numeric::ublas::bounded_matrix<double, 6, 3 > B_matrix = ZeroMatrix(6,3);
				boost::numeric::ublas::bounded_matrix<double, 6, 1 >velocities = ZeroMatrix(6,1);
				double Area;
				array_1d<double, 3> N;
				array_1d<double, 3> distances;
				boost::numeric::ublas::bounded_matrix<double, 3, 2> DN_DX;
				GeometryUtils::CalculateGeometryData(geom, DN_DX, N, Area);
				
				const double mu = this->GetValue(YOUNG_MODULUS)/ (2.0*(1.0+this->GetValue(POISSON_RATIO))); 

				for (unsigned int i=0; i!=3; i++) //i node
				{
						B_matrix(i*2,0)=DN_DX(i,0);
						B_matrix(i*2+1,1)=DN_DX(i,1);
						
						B_matrix(i*2,2)=DN_DX(i,1);
						B_matrix(i*2+1,2)=DN_DX(i,0);
						
						array_1d<double, 3 >velocity =  ( geom[i].FastGetSolutionStepValue(VELOCITY));//+geom[i].FastGetSolutionStepValue(MESH_VELOCITY) ) * 0.5;
						velocities(i*2+0,0)=velocity(0);
						velocities(i*2+1,0)=velocity(1);
						
						distances(i)=geom[i].FastGetSolutionStepValue(DISTANCE);
				}
				
				double reduced_delta_t=delta_t;
				CalculateReducedTimeStepForPositiveJacobian(DN_DX,velocities,delta_t,reduced_delta_t);
				
				
				double number_of_solid_particles=1.e-10; //almost zero
				Vector & element_stresses= this->GetValue(ELEMENT_MEAN_STRESS);
				//array_1d<double,3> sum_particle_stresses=ZeroVector(3);
				
				UpdateElementStresses(element_stresses,mu,geom,velocities,B_matrix,reduced_delta_t);

				if(use_failure_criteria)
					TestElementWithDruckerPrager(geom);
		}
	}
	
	

	void NoParticlesSolidOnlyPFEM22D::UpdateElementStresses(
						 Vector & element_stresses,
						 const double mu,
						 Geometry< Node<3> >& geom,
						 const boost::numeric::ublas::bounded_matrix<double, 6, 1 >& velocities,
						 const boost::numeric::ublas::bounded_matrix<double, 6, 3 >& B_matrix,
						 const double delta_t
						 )
	{
		//if(mesh_distance<0.0)	
		const int TDim=2;
		if(true)
		{
			boost::numeric::ublas::bounded_matrix<double, (TDim-1)*3, (TDim-1)*3 > C_matrix = ZeroMatrix((TDim-1)*3,(TDim-1)*3);
			
			C_matrix(0,0)=4.0/3.0;
			C_matrix(1,1)=4.0/3.0;
			C_matrix(2,2)=1.0;
			
			C_matrix(1,0)=-2.0/3.0;
			C_matrix(0,1)=-2.0/3.0;
			
			C_matrix *= mu;
		
			boost::numeric::ublas::bounded_matrix<double, 3, 1 > delta_strains = prod(trans(B_matrix),velocities)*delta_t;
			
			boost::numeric::ublas::bounded_matrix<double, 3, 1 > delta_stresses = prod(C_matrix,delta_strains); //for the moment is DELTA_stresses
			
			for (unsigned int i=0;i<3;i++)
				//particle_stress(i)=(pelement->GetValue(ELEMENT_MEAN_STRESS))(i)+delta_stresses(i,0);
				element_stresses(i) += delta_stresses(i,0); 
				
		}
	}
	
	
	void NoParticlesSolidOnlyPFEM22D::TestElementWithDruckerPrager( Geometry< Node<3> >& geom)
	{
			Vector& stress = this->GetValue(ELEMENT_MEAN_STRESS);
			const double mean_pressure = 1.0/3.0*(geom[0].FastGetSolutionStepValue(PRESSURE)+geom[1].FastGetSolutionStepValue(PRESSURE)+geom[2].FastGetSolutionStepValue(PRESSURE));
			const double von_misses_trial_stress =  sqrt(1.0/3.0* (0.25 * pow((stress(0)-stress(1)),2) + pow((stress(2)),2) ) ) ;

		    const double  theta =  this->GetValue(INTERNAL_FRICTION_ANGLE);
		    const double cohesion =  this->GetValue(YIELD_STRESS);
		    //if (pparticle.GetDistance()> -0.9);
			//	cohesion = 0.0;
			double A = 6.0 * cohesion *sin(theta) / ( 3.0 - sin(theta));
			double B = 2.0 *sin(theta) / ( 3.0 - sin(theta));
			
			double failure = von_misses_trial_stress - A - B * mean_pressure;
			if (failure>0.0)
			{				
				//cohesion*=0.8;
				
				double pressure=(von_misses_trial_stress+(1.0/B)*mean_pressure - A)/(B+1.0/B);
				double von_misses_stress = A + B*pressure;
				double reducing_factor= von_misses_stress/von_misses_trial_stress ;
				if (reducing_factor>1.0)
					KRATOS_WATCH("ARTT");
				if (reducing_factor>0.0)
				{
					if (pressure<mean_pressure)
						KRATOS_WATCH("molt malament");
					stress=stress*reducing_factor;
					//particle_pressure=pressure;					
				}
				else
				{
					stress = stress*0.0;
					pressure = - A/B;
				}
			}
	}
	
	
	void NoParticlesSolidOnlyPFEM22D::CalculateReducedTimeStepForPositiveJacobian(const boost::numeric::ublas::bounded_matrix<double, 3, 2>& rShapeDeriv,const boost::numeric::ublas::bounded_matrix<double, 6, 1 >& velocities, const double& delta_t, double& reduced_delta_t)
    {
		boost::numeric::ublas::bounded_matrix<double, 3 , 2 > displacements;
		for (unsigned int i=0;i!=3;i++)
		{
			displacements(i,0)=velocities(i*2+0,0)*delta_t;
			displacements(i,1)=velocities(i*2+1,0)*delta_t;
		}	
		
		boost::numeric::ublas::bounded_matrix<double, 2 , 2 > DeformationGradient = identity_matrix<double> (2) ; //ZeroMatrix(2,2);// = prod(trans(displacements),rShapeDeriv);
		
		DeformationGradient+=prod(trans(displacements),rShapeDeriv);

        double J = DeformationGradient ( 0 , 0 ) *  DeformationGradient ( 1 , 1 ) - DeformationGradient ( 0 , 1 ) * DeformationGradient ( 1 , 0 ); //J is the determinant of the deformation gradient;
		
		double reducing_factor=1.0;
		if (J<=0.8 )
		{
			double amplifying_factor=1.0;
			
			while (J<=0.8 )
			{
				reducing_factor=1.0/(1.0+0.1*amplifying_factor);
				DeformationGradient = identity_matrix<double> (2) + reducing_factor*prod(trans(displacements),rShapeDeriv);
				J = DeformationGradient ( 0 , 0 ) *  DeformationGradient ( 1 , 1 ) - DeformationGradient ( 0 , 1 ) * DeformationGradient ( 1 , 0 );
				amplifying_factor++;
			}
			//std::cout << "elem " << this->Id() << ", reduced timestep by " << reducing_factor << std::endl;
		}
		reduced_delta_t = reducing_factor*delta_t;
	}
	

} // Namespace Kratos
