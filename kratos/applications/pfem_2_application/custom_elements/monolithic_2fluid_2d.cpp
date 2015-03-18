//   
//   Project Name:        Kratos       
//   Last modified by:    $Author: it's me! $
//   Date:                $Date: 2008-08-08 23:58:38 $
//   Revision:            $Revision: 1.0 $
//
//

// Project includes 
#include "includes/define.h"
#include "custom_elements/monolithic_2fluid_2d.h"
#include "pfem_2_application.h"
#include "utilities/math_utils.h"
#include "utilities/geometry_utilities.h" 
#include "utilities/enrichment_utilities.h"
#include "utilities/enrich_2d_2dofs.h"

#include "processes/process.h"
#include "includes/node.h"
//#include "includes/element.h"
#include "includes/model_part.h"



namespace Kratos
{
	


	//************************************************************************************
	//************************************************************************************
	MonolithicPFEM22D::MonolithicPFEM22D(IndexType NewId, GeometryType::Pointer pGeometry)
		: Element(NewId, pGeometry)
	{		
		//DO NOT ADD DOFS HERE!!!
	}

	//************************************************************************************
	//************************************************************************************
	MonolithicPFEM22D::MonolithicPFEM22D(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
		: Element(NewId, pGeometry, pProperties)
	{

	}

	Element::Pointer MonolithicPFEM22D::Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const
	{
		return Element::Pointer(new MonolithicPFEM22D(NewId, GetGeometry().Create(ThisNodes), pProperties));
	}

	MonolithicPFEM22D::~MonolithicPFEM22D()
	{
	}

	//************************************************************************************
	//************************************************************************************
	
	
	void MonolithicPFEM22D::AddExplicitContribution(ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY;

		switch ( rCurrentProcessInfo[FRACTIONAL_STEP] )
		{
		case 0: //calculating viscous term and the lumped mass matrix. // nothing returned. saving data in RHS and NODAL_MASS
		{
			this->CalculateViscousRHS(rCurrentProcessInfo);
			break;
		}	
			
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
	
	void MonolithicPFEM22D::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY
		array_1d<double,3> gravity= rCurrentProcessInfo[GRAVITY];

		//WE MUST CALCULATE Mass and G(or D) matrixes

	    double delta_t = rCurrentProcessInfo[DELTA_TIME];	
		//KRATOS_WATCH(gradient_discontinuity);
		const bool & split_element = (this->GetValue(SPLIT_ELEMENT));

		double theta=1.0;
	
		//array_1d<double,3> pressures;

		
		array_1d<double,6>  previous_vel;
		array_1d<double,9>  previous_vel_and_press = ZeroVector(9);
        for(unsigned int iii = 0; iii<3; iii++)
        {
			array_1d<double,3> velocity =  GetGeometry()[iii].FastGetSolutionStepValue(VELOCITY);
			//velocity -= gravity*delta_t; //it was added in the particles, but since we need it here on the rhs instead of inside the velocity, we subtract it in order to add it correctly in the rhs.
			previous_vel(iii*2) = velocity[0];
			previous_vel(iii*2+1) = velocity[1];
			
			previous_vel_and_press(iii*3) = velocity[0];
			previous_vel_and_press(iii*3+1) = velocity[1];
			previous_vel_and_press(iii*3+2) = GetGeometry()[iii].FastGetSolutionStepValue(PRESSURE);
		}
		
		//boost::numeric::ublas::bounded_matrix<double,3,2> msDN_DX;
		//boost::numeric::ublas::bounded_matrix<double,2,2> msD;
		const double one_third = 1.0/3.0;
  		array_1d<double,3> msN; //dimension = number of nodes
		array_1d<double,6> ms_temp; //dimension = number of nodes

		const unsigned int LocalSize = GetGeometry().size()*3; //2 vel + 1 pressure per node
		unsigned int TNumNodes = GetGeometry().size();

		array_1d<double, 3> N;
        // Check sizes and initialize
        if (rRightHandSideVector.size() != LocalSize)
            rRightHandSideVector.resize(LocalSize, false);
        noalias(rRightHandSideVector) = ZeroVector(LocalSize);
        if(rLeftHandSideMatrix.size1() != LocalSize)
			rLeftHandSideMatrix.resize(LocalSize,LocalSize,false);
		noalias(rLeftHandSideMatrix) = ZeroMatrix(LocalSize,LocalSize);	
        // Calculate this element's geometric parameters
        double Area;
        
        boost::numeric::ublas::bounded_matrix<double, 3, 2> DN_DX;

        //boost::numeric::ublas::bounded_matrix<double, TNumNodes*2, 2> N_vel_matrix; //useless to calculate it, it's simply the DN_DX matrix multiplied by 1/3, arranged in a different shape.
        
        GeometryUtils::CalculateGeometryData(this->GetGeometry(), DN_DX, N, Area);
        
        bool has_void_node=false; 
        array_1d<double,3> distances(3);
        for (unsigned int i = 0; i < TNumNodes; i++)
        {
            distances[i] = this->GetGeometry()[i].FastGetSolutionStepValue(DISTANCE);
            if (distances[i]>1.5) 
				has_void_node=true;
        }
        
        
		const double viscosity_air = rCurrentProcessInfo[VISCOSITY_AIR]; //para clindro en rey 100:  0.0005 * 2.0*Area;
		const double viscosity_water = rCurrentProcessInfo[VISCOSITY_WATER];
		
		const double density_air = rCurrentProcessInfo[DENSITY_AIR]; //para clindro en rey 100:  0.0005 * 2.0*Area;
		const double density_water = rCurrentProcessInfo[DENSITY_WATER];
		
		
		

		//KRATOS_WATCH(D_matrix)
		
		if (split_element==false || has_void_node)
		{
			boost::numeric::ublas::bounded_matrix<double, 3, 3> Laplacian_matrix;
			boost::numeric::ublas::bounded_matrix<double, 3, 6 > D_matrix; //(gradient)
			noalias(D_matrix) = ZeroMatrix(3,6);	
			boost::numeric::ublas::bounded_matrix<double, 3, 6 > G_matrix; //(gradient)
			noalias(G_matrix) = ZeroMatrix(3,6);

			boost::numeric::ublas::bounded_matrix<double, 9, 9 > Mass_matrix; //2 vel + 1 pressure per node
			noalias(Mass_matrix) = ZeroMatrix(LocalSize,LocalSize);	
			
			//we start by calculating the non-enriched functions.
			for (unsigned int i = 0; i < 3; i++)
			{
				for (unsigned int j = 0; j < 3; j++)
				{
					D_matrix(i, (j*2) ) = DN_DX(j,0)*Area*one_third;     
					D_matrix(i, (j*2+1) ) =  DN_DX(j,1)*Area*one_third;
					G_matrix(i, (j*2) ) =  - DN_DX(i,0)*Area*one_third;     
					G_matrix(i, (j*2+1) ) = - DN_DX(i,1)*Area*one_third;
				}
			}
			
			
			
			//viscous term:
			double viscosity = 0.0;
			
			if (distances(0)<0.0)
				viscosity = viscosity_water;
			else
				viscosity = viscosity_air;
			
			this->AddViscousTerm(rLeftHandSideMatrix,DN_DX, (viscosity*Area) );
			
			double element_mean_distance=(distances[0]+distances[1]+distances[2])*0.333333333;
			const double water_fraction = 0.5*(1.0-(sin(3.14159*element_mean_distance*0.5)));
			double density = density_water*(water_fraction)+density_air*(1.0-water_fraction);
			if (element_mean_distance>0) density=density_air;
			else density=density_water;
			const double Weight = one_third * Area * density;
			for (unsigned int i=0; i<TNumNodes ; i++)
			{
				Mass_matrix (i*3,i*3) = Weight;
				Mass_matrix (i*3+1,i*3+1) = Weight;
			}	
			
			//PRESSURE TERMS:
			/////////////////// calculating tau
			double AdvVelNorm = sqrt(pow(this->GetGeometry()[0].FastGetSolutionStepValue(VELOCITY_X),2)+pow(this->GetGeometry()[0].FastGetSolutionStepValue(VELOCITY_Y),2));
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
			const double   TauOne = 1.0 / ( ( 1.0 / delta_t + 4.0 * viscosity / (density * mElemSize * mElemSize) + 2.0 * AdvVelNorm / mElemSize) );
			Laplacian_matrix = - prod(DN_DX,trans(DN_DX))*Area*TauOne/density;  // B B^T  (standard laplacian, now we must add the contribution of the condensed new nodes.
			/*
			///DIFFERENT STABILIZATION: - (consistent - lumped mass matrix).
			const double one_sixt=1.0/6.0;
			const double one_twelfth=1.0/12.0;
			Laplacian_matrix(0,0)=-one_sixt*Area;
			Laplacian_matrix(0,1)=one_twelfth*Area;
			Laplacian_matrix(0,2)=one_twelfth*Area;
			Laplacian_matrix(1,0)=one_twelfth*Area;
			Laplacian_matrix(1,1)=-one_sixt*Area;
			Laplacian_matrix(1,2)=one_twelfth*Area;
			Laplacian_matrix(2,0)=one_twelfth*Area;
			Laplacian_matrix(2,1)=one_twelfth*Area;
			Laplacian_matrix(2,2)=-one_sixt*Area;
			Laplacian_matrix *= 1.0 * delta_t/(density * mElemSize * mElemSize);  // 0.5= 0.7^2
			*/

			for (unsigned int i = 0; i < 3; i++)
			{
				for (unsigned int j = 0; j < 3; j++)
				{	
					rLeftHandSideMatrix(i*3+2, j*3+0 ) -= D_matrix(i,j*2);//  + DN_DX(i,0)*TauOne*one_third*Area/delta_t;     
					rLeftHandSideMatrix(i*3+2, j*3+1 ) -= D_matrix(i,j*2+1);// + DN_DX(i,1)*TauOne*one_third*Area/delta_t;
					
					rLeftHandSideMatrix(j*3+0, i*3+2 ) -=  D_matrix(i,j*2);     
					rLeftHandSideMatrix(j*3+1, i*3+2 ) -= D_matrix(i,j*2+1);
					
					rLeftHandSideMatrix(i*3+2, j*3+2 ) += Laplacian_matrix(i,j);
				}
			}
			
			array_1d<double, 2 > rhs_stab = ZeroVector(2);
			for (unsigned int i = 0; i < 3; i++)
			{
				array_1d<double,3>& node_press_proj = GetGeometry()[i].FastGetSolutionStepValue(PRESS_PROJ);
				rhs_stab[0] += node_press_proj(0)*one_third;
				rhs_stab[1] += node_press_proj(1)*one_third;
			}
				
			const bool use_press_proj=true;
			for (unsigned int i = 0; i < 3; i++)
			{
				rRightHandSideVector(i*3+0) += one_third*Area*gravity(0)*density;
				rRightHandSideVector(i*3+1) += one_third*Area*gravity(1)*density;

				if (use_press_proj)
						rRightHandSideVector(i*3+2) -= TauOne*Area*(DN_DX(i,0)*rhs_stab(0)+DN_DX(i,1)*rhs_stab(1));
				else
					rRightHandSideVector(i*3+2) -= TauOne*(DN_DX(i,0)*(gravity(0))+DN_DX(i,1)*(gravity(1)))*Area;
			}

			noalias(rRightHandSideVector) += prod((Mass_matrix/delta_t),previous_vel_and_press);
			noalias(rLeftHandSideMatrix) += Mass_matrix/delta_t;   	
		} 
		else //split element:
		{
			rLeftHandSideMatrix=ZeroMatrix(9,9);
			
			//TOOLS NEEDED TO FIND THE ENRICHMENT SHAPE FUNCTIONS
			//get position of the cut surface
			array_1d<double,6>  local_previous_vel;
			array_1d<double,3>  densities(1);
			boost::numeric::ublas::bounded_matrix<double,3, 4> Nenriched;
			array_1d<double,3>  volumes(1);
			boost::numeric::ublas::bounded_matrix<double,3, 2 > coords;
			boost::numeric::ublas::bounded_matrix<double,3, 3 > Ngauss;
			array_1d<double,3>  signs(1);
			std::vector< Matrix > gauss_gradients(3);
			//std::vector< Matrix > gauss_gradients_in_local_axis(3);
			//boost::numeric::ublas::bounded_matrix<double,2, 2 > RotationMatrix;
			//boost::numeric::ublas::bounded_matrix<double, 3, 2> DN_DX_in_local_axis=ZeroMatrix(3,2);
			
			const int  enrich_velocity_dofs=2; //number of enrichments for the velocity field
			const int  enrich_pressure_dofs=2; //number of enrichments for the pressure field
			const int  enrich_pressure_offset=0; //the enrichment utility returns several (4) shape functions. we will avoid the first 2 and use the third and fourth
			
			boost::numeric::ublas::bounded_matrix<double, 3, 3> Laplacian_matrix =ZeroMatrix(3,3); //standard pressures. we will keep these dof so the size remains			
			boost::numeric::ublas::bounded_matrix<double, 3, 6+enrich_velocity_dofs > D_matrix ; //(divergence) //this matrix will increase its size since we need more  
			noalias(D_matrix) = ZeroMatrix(3,6+enrich_velocity_dofs);	
			boost::numeric::ublas::bounded_matrix<double, 3, 6+enrich_velocity_dofs > G_matrix; //(gradient)
			noalias(G_matrix) = ZeroMatrix(3,6+enrich_velocity_dofs);

			boost::numeric::ublas::bounded_matrix<double, 9+enrich_velocity_dofs+enrich_pressure_dofs, 9+enrich_velocity_dofs+enrich_pressure_dofs > Momentum_matrix; //2 vel + 1 pressure per node   plus 2 enriched pressures and 2 enriched velocities
			noalias(Momentum_matrix) = ZeroMatrix(9+enrich_velocity_dofs+enrich_pressure_dofs,9+enrich_velocity_dofs+enrich_pressure_dofs);	

			boost::numeric::ublas::bounded_matrix<double,(2),2 > rRotationMatrix=ZeroMatrix(2,2);
			boost::numeric::ublas::bounded_matrix<double,(2+1), 2 > rRotatedPoints=ZeroMatrix(3,2);
			boost::numeric::ublas::bounded_matrix<double, (2+1), 2 > DN_DX_in_local_axis=ZeroMatrix(3,2);
			
			
			//fill coordinates
		   
			//unsigned int single_triangle_node = 1;
			for (unsigned int i = 0; i < TNumNodes; i++)
			{
				const array_1d<double, 3 > & xyz = this->GetGeometry()[i].Coordinates();
				volumes[i] = 0.0;
				//KRATOS_WATCH(distances(i));
				for (unsigned int j = 0; j < 2; j++)
					coords(i, j) = xyz[j];
			}
			
			//getting needed information to calculate everything in local coordinates. We will rotate everything. 
			//As this is being written, isotropic material properties are being used, so this is a waste of time. 
			//however it might be useful if different properties along the normal and tangential axis have to be used
			CalculateRotationParameters(
										coords, 
										distances,
										rRotationMatrix, 
										rRotatedPoints,
										DN_DX_in_local_axis);
										
			//DN_DX_in_local_axis=DN_DX;
			//rRotationMatrix=identity_matrix<double>(2);							
			//rRotatedPoints=coords;
			rRotationMatrix=trans(rRotationMatrix);
			for (unsigned int i = 0; i < 3; i++) //3 nodes
			{
				
				boost::numeric::ublas::bounded_matrix<double, (2), 1 >	global_node_vel;
				for (unsigned int k = 0; k < 2; k++) //xy	
					global_node_vel(k,0) = previous_vel(i*2+k);
				boost::numeric::ublas::bounded_matrix<double, (2), 1 >	local_node_vel = prod(rRotationMatrix,global_node_vel);
				for (unsigned int k = 0; k < 2; k++) //xy					
					local_previous_vel(i*2+k)=local_node_vel(k,0);
			}
			boost::numeric::ublas::bounded_matrix<double, (2), 1 > gravity_2d;
			for (unsigned int k = 0; k < 2; k++) //xy
				gravity_2d(k,0)=gravity(k);
			boost::numeric::ublas::bounded_matrix<double, (2), 1 >	local_gravity=prod(rRotationMatrix,gravity_2d);

			
			for (unsigned int i = 0; i < 3; i++)
				gauss_gradients[i].resize(4, 2, false);  //3 values of the 3 shape functions, and derivates in (xy) direction).

			unsigned int ndivisions = EnrichmentUtilities::CalculateEnrichedShapeFuncionsExtended(rRotatedPoints, DN_DX_in_local_axis, distances, volumes, Ngauss, signs, gauss_gradients, Nenriched);

			//interface info:
			double interface_area=0.0;
			array_1d<double,2>  normal=ZeroVector(2);
			array_1d<double,3>  Ninterface=ZeroVector(3);
			this->CalculateInterfaceNormal(coords, distances, normal,interface_area, Ninterface);
			
			//we start by calculating the non-enriched functions.
			for (unsigned int i = 0; i < 3; i++)
			{
				for (unsigned int j = 0; j < 3; j++)
				{
					D_matrix(i, (j*2) ) = DN_DX_in_local_axis(j,0)*Area*one_third;     
					D_matrix(i, (j*2+1) ) =  DN_DX_in_local_axis(j,1)*Area*one_third;
					G_matrix(i, (j*2) ) = - DN_DX_in_local_axis(i,0)*Area*one_third;     
					G_matrix(i, (j*2+1) ) = - DN_DX_in_local_axis(i,1)*Area*one_third;
				}
			}
			//the rest must be done inside the partitions loop
						
			double Weight=0.0;
			double mass=0.0;
			array_1d<double,3>  viscosities;
			for (unsigned int i=0;i<ndivisions;i++) //we go over the three partitions;
			{
				if (signs(i)>0.0)
				{
					Weight += volumes(i)*viscosity_air;
					mass += volumes(i)*density_air;
					viscosities(i) = viscosity_air; 
					densities(i)=density_air;
				}
				else
				{
					Weight += volumes(i)*viscosity_water;
					mass += volumes(i)*density_water;
					viscosities(i) = viscosity_water;
					densities(i)=density_water;
				}
			}
			
			//this->AddViscousTerm(rLeftHandSideMatrix,DN_DX, (Weight) ); //for an  averaged viscosity ( no enrichments)
			this->AddViscousTerm(Momentum_matrix,  //extended matrix, with the velocity enrichments
								DN_DX_in_local_axis,
								distances,
								gauss_gradients, 
								viscosities,
								signs,
								volumes,
								ndivisions);
			
			
			const double element_viscosity=Weight/Area;
			const double element_density= mass/Area;

			array_1d<double,3>  lumped_mass = ZeroVector(3);
			double condensed_dof_mass1=0.0;
			double condensed_dof_mass2=0.0;
			for (unsigned int i = 0; i < 3; i++)  //partitions
			{
				//KRATOS_WATCH(enrich_lhs(i));
				for (unsigned int j = 0; j < 3; j++) //local node (or shape function)
				{
					lumped_mass(j) +=  volumes(i)*Ngauss(i,j)*densities(i);	
					Momentum_matrix(3*j,3*j) += volumes(i)*Ngauss(i,j)*densities(i) / delta_t;	 
					Momentum_matrix(3*j+1,3*j+1) += volumes(i)*Ngauss(i,j)*densities(i) / delta_t;	 
				}
				

				if (enrich_velocity_dofs!=0)
				{
					condensed_dof_mass1 += volumes(i)*Nenriched(i,3)*densities(i);
					condensed_dof_mass2 += volumes(i)*Nenriched(i,3)*densities(i);
					   
					Momentum_matrix(9+0,9+0) += volumes(i)*Nenriched(i,3)*densities(i) / delta_t;	 //degrees of freedom 10 and 11 are the enrichment velocities. degrees of freedom 12 and 13 are the enrichment pressures
					Momentum_matrix(9+1,9+1) += volumes(i)*Nenriched(i,3)*densities(i) / delta_t;	 
				}
			}
			
			
			//PRESSURE TERMS:
			/////////////////// calculating tau
			double AdvVelNorm = sqrt(pow(this->GetGeometry()[0].FastGetSolutionStepValue(VELOCITY_X),2)+pow(this->GetGeometry()[0].FastGetSolutionStepValue(VELOCITY_Y),2));
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
			const double   TauOne = 1.00 / (( 1.0 / delta_t + 4.0 * element_viscosity / (mElemSize * mElemSize * element_density ) + 2.0 * AdvVelNorm / mElemSize) );
			
			////////////////////
			//TauOne=1.0;
			array_1d<double,3> partition_densities;
			array_1d<double,3> inv_densities;
			double Density=0.0; //resetting to zero
			for (unsigned int i = 0; i!=ndivisions; i++) //the 3 partitions of the triangle
			{
				if (signs(i)>0.0)
					partition_densities(i)=rCurrentProcessInfo[DENSITY_AIR];
				else
					partition_densities(i)=rCurrentProcessInfo[DENSITY_WATER];			
				inv_densities(i)=1.0/partition_densities(i);	
				Density+=partition_densities(i)*volumes(i);
			}
			Density*=1.0/Area;			
			
			Laplacian_matrix =  - prod(DN_DX_in_local_axis,trans(DN_DX_in_local_axis)) *(volumes(0)*inv_densities(0)+volumes(1)*inv_densities(1)+volumes(2)*inv_densities(2));

			//and now the rest of the things needed
			boost::numeric::ublas::bounded_matrix<double, enrich_pressure_dofs, enrich_pressure_dofs > Laplacian_enrich;
			noalias(Laplacian_enrich) = ZeroMatrix(enrich_pressure_dofs,enrich_pressure_dofs);	
			boost::numeric::ublas::bounded_matrix<double, enrich_pressure_dofs, 3 > mixed_Laplacian;
			noalias(mixed_Laplacian) = ZeroMatrix(enrich_pressure_dofs,3);	
			//To calculate D* =int (N*DN_DX_enrich).  we must use the N at the gauss point of the partition (returned by the enrichment util) and multiply by the area. 
			//notice that this first row is only the upper part of the mixed D. we also need the lower one for the jump 
			boost::numeric::ublas::bounded_matrix<double, enrich_pressure_dofs, 6 + enrich_velocity_dofs> D_matrix_mixed;
			noalias(D_matrix_mixed) = ZeroMatrix(enrich_pressure_dofs,6 + enrich_velocity_dofs);	
			boost::numeric::ublas::bounded_matrix<double, enrich_pressure_dofs, 6 + enrich_velocity_dofs > G_matrix_mixed;
			noalias(G_matrix_mixed) = ZeroMatrix(enrich_pressure_dofs,6 + enrich_velocity_dofs);	

			//the matrices we need at the end will be:
			//boost::numeric::ublas::bounded_matrix<double, 2, 2 > D_mod;
						
			//double rhs_enrich = 0.0; //the rhs term of the enriched dof 
			array_1d<double,( enrich_pressure_dofs + enrich_velocity_dofs )> rhs_enrich;
			noalias(rhs_enrich) = ZeroVector(enrich_pressure_dofs + enrich_velocity_dofs);	
			
				//			KRATOS_WATCH(Ngauss)
				//KRATOS_WATCH(Nenriched)
			
			for (unsigned int i = 0; i < ndivisions; i++)  //we go through the 3 partitions of the domain
			{
				
				for (unsigned int j = 0; j < enrich_pressure_dofs; j++) //we go through the 4 enrichments
					for (unsigned int k = 0; k < enrich_pressure_dofs; k++) //we go through the 4 enrichments
						Laplacian_enrich(j,k) -= ((gauss_gradients[i](j+enrich_pressure_offset,0)*gauss_gradients[i](k+enrich_pressure_offset,0))+(gauss_gradients[i](j+enrich_pressure_offset,1)*gauss_gradients[i](k+enrich_pressure_offset,1)))*volumes(i)*inv_densities(i);

				//and the 'mixed laplacians' (standard shape functions * enrichments)
				
				
				array_1d<double,3> previous_velocity_at_gauss_point=ZeroVector(3); //needed for the rhs
				
				//first we take care of the standard shape functions of the velocity
				for (unsigned int j = 0; j < 3; j++) //we go through the 3 standard shape functions
				{
					for (unsigned int k = 0; k < enrich_pressure_dofs; k++) //we go through the 3 enrichments
					{
						mixed_Laplacian(k,j) -= (DN_DX_in_local_axis(j,0)*(gauss_gradients[i](k+enrich_pressure_offset,0))+ DN_DX_in_local_axis(j,1)*(gauss_gradients[i](k+enrich_pressure_offset,1)))*volumes(i)*inv_densities(i);
						D_matrix_mixed(k,j*2) += DN_DX_in_local_axis(j,0) *volumes(i)*Nenriched(i,k+enrich_pressure_offset);
						D_matrix_mixed(k,j*2+1) += DN_DX_in_local_axis(j,1) *volumes(i)*Nenriched(i,k+enrich_pressure_offset);
						
						G_matrix_mixed(k,j*2) -= gauss_gradients[i](k+enrich_pressure_offset,0) *volumes(i)*Ngauss(i,j);
						G_matrix_mixed(k,j*2+1) -= gauss_gradients[i](k+enrich_pressure_offset,1) *volumes(i)*Ngauss(i,j);
					}
					previous_velocity_at_gauss_point(0)+=Ngauss(i,j)*local_previous_vel(j*2+0);
					previous_velocity_at_gauss_point(1)+=Ngauss(i,j)*local_previous_vel(j*2+1);
				}
				
				
				//we go through the remaining, new velocity dofs. (one x and one y component
				if (enrich_velocity_dofs!=0)
				{
					for (unsigned int k = 0; k < enrich_pressure_dofs; k++) //we go through the 3 enrichments
					{
						D_matrix_mixed(k,6+0) += gauss_gradients[i](3,0) *volumes(i)*Nenriched(i,k+enrich_pressure_offset);
						D_matrix_mixed(k,6+1) += gauss_gradients[i](3,1) *volumes(i)*Nenriched(i,k+enrich_pressure_offset);
						G_matrix_mixed(k,6+0) -= gauss_gradients[i](k+enrich_pressure_offset,0) *volumes(i)*Nenriched(i,3);
						G_matrix_mixed(k,6+1) -= gauss_gradients[i](k+enrich_pressure_offset,1) *volumes(i)*Nenriched(i,3);
					}
					
					
					//also the standard D and G matrices had to be expanded.

						for (unsigned int k = 0; k < 3; k++) //we go through the 3 standard pressures
						{
							D_matrix(k,6+0) += gauss_gradients[i](3,0) *volumes(i)*Ngauss(i,k);
							D_matrix(k,6+1) += gauss_gradients[i](3,1) *volumes(i)*Ngauss(i,k);
							G_matrix(k,6+0) -= DN_DX_in_local_axis(k,0) *volumes(i)*Nenriched(i,3);
							G_matrix(k,6+1) -= DN_DX_in_local_axis(k,1) *volumes(i)*Nenriched(i,3);
						}
				}
				
				for (unsigned int k = 0; k < enrich_pressure_dofs; k++) //we go through the 4 enrichments
					rhs_enrich(k+enrich_velocity_dofs) -= (gauss_gradients[i](k+enrich_pressure_offset,0)*(local_gravity(0,0))+gauss_gradients[i](k+enrich_pressure_offset,1)*(local_gravity(1,0)))*volumes(i);
				
				
			}
			//KRATOS_WATCH(Laplacian_enrich(0,0));
			//TauOne *=0.01;
			Laplacian_enrich*=TauOne;
			Laplacian_matrix*=TauOne;
			mixed_Laplacian*=TauOne;
			rhs_enrich *=TauOne;

			boost::numeric::ublas::bounded_matrix<double, 1, 9 > penalty_vector=ZeroMatrix(1,9);
		
			array_1d<double,3> standard_velocity_at_interface=ZeroVector(3);
			array_1d<double,3> negative_side_only_velocity_at_interface=ZeroVector(3);
			int negative_nodes=0;
			double sum_n_negative=0.0;
			array_1d<double,3> positive_side_only_velocity_at_interface=ZeroVector(3);
			int positive_nodes=0;
			double sum_n_positive=0.0;
			for (unsigned int i=0; i!=3;i++) //nodes
			{
				if (distances(i)<0.0)
				{
					negative_nodes++;
					sum_n_negative+=Ninterface(i);
					for (unsigned int k = 0; k < 2; k++) //xy
						negative_side_only_velocity_at_interface(k) += local_previous_vel(2*i+k)*Ninterface(i);
				}
				else
				{
					positive_nodes++;
					sum_n_positive+=Ninterface(i);
					for (unsigned int k = 0; k < 2; k++) //xy
						positive_side_only_velocity_at_interface(k) += local_previous_vel(2*i+k)*Ninterface(i);
				}
				for (unsigned int k = 0; k < 2; k++) //xy
					standard_velocity_at_interface(k) += local_previous_vel(2*i+k)*Ninterface(i);
			}	
			
			negative_side_only_velocity_at_interface /= sum_n_negative;
			positive_side_only_velocity_at_interface /= sum_n_positive;
			//we will force the velocity of the previous time step at the interface to match the one of the negative side, as it the other had negligible viscosity
			//the V* shape funcion has value = 1.0 at the interface, so that means the N* will be direcly the difference between what we want (the negative side vel) and the complete one:
			array_1d<double,3> N_vel_star = ZeroVector(3); 
			for (unsigned int i=0; i!=2;i++) //TDim 
				N_vel_star(i) = negative_side_only_velocity_at_interface(i)-positive_side_only_velocity_at_interface(i);
			N_vel_star *= 0.5; //should be 0.5, but might become unstable	
			
			
			
			if (enrich_velocity_dofs!=0)
			{
				rhs_enrich(0) = local_gravity(0,0)*condensed_dof_mass1;
				rhs_enrich(1) = local_gravity(1,0)*condensed_dof_mass2;
				rhs_enrich(0) += N_vel_star(0) * condensed_dof_mass1/delta_t; //we are using the difference in the Y component
				rhs_enrich(1) += N_vel_star(1) * condensed_dof_mass2/delta_t;
				
			}
						
			
			boost::numeric::ublas::bounded_matrix<double, enrich_velocity_dofs+enrich_pressure_dofs, 9 > condensed_rows; //Vx1,Vy1,p1,Vx2,...
			boost::numeric::ublas::bounded_matrix<double, enrich_velocity_dofs+enrich_pressure_dofs, 9 > condensed_columns; //Vx1,Vy1,p1,Vx2,...
			boost::numeric::ublas::bounded_matrix<double, enrich_velocity_dofs+enrich_pressure_dofs, enrich_velocity_dofs+enrich_pressure_dofs > condensed_block; //Vx1,Vy1,p1,Vx2,...
			

			
			for (unsigned int i = 0; i <3; i++)  //we go through the 3 nodes (standard velocity dof + standard pressure dof)
			{
				
				//enriched pressure dof	
				for (unsigned int k = 0; k < enrich_pressure_dofs; k++)
				{
					//enriched
					condensed_rows(enrich_velocity_dofs+k,i*3+0)= - D_matrix_mixed(k,i*2+0);//+mass_stabilization_terms(i*2+0);
					condensed_rows(enrich_velocity_dofs+k,i*3+1)= - D_matrix_mixed(k,i*2+1);//+mass_stabilization_terms(i*2+1);
					condensed_rows(enrich_velocity_dofs+k,i*3+2)= mixed_Laplacian(k,i);
					
					condensed_columns(enrich_velocity_dofs+k,i*3+0)= - D_matrix_mixed(k,i*2+0);
					condensed_columns(enrich_velocity_dofs+k,i*3+1)= - D_matrix_mixed(k,i*2+1);
					condensed_columns(enrich_velocity_dofs+k,i*3+2)= mixed_Laplacian(k,i);	
				}
				
				
				
				for (unsigned int k = 0; k < enrich_velocity_dofs; k++) //we go through the 3 enrichments
				{
					condensed_columns(k,i*3+0)=Momentum_matrix(9+k,i*3+0); //add the viscosity matrix
					condensed_columns(k,i*3+1)=Momentum_matrix(9+k,i*3+1); //add the viscosity matrix
					
					condensed_rows(k,i*3+0)=Momentum_matrix(9+k,i*3+0); //add the viscosity matrix
					condensed_rows(k,i*3+1)=Momentum_matrix(9+k,i*3+1); //add the viscosity matrix
					
					
					///WARNING, WHEN THE MATRIX IS REARRANGED, the condensed rows have the gradient of the pressure and the columns have the divergence. that is the reason for the mixed indexes.
					///if the gradient was not integrated by parts, then G matrix should be used in the columns instead of the rows, unlike the case of standard velocity DOFs
					condensed_rows(k,i*3+2)= -G_matrix(i, 6+k);
					condensed_columns(k,i*3+2)= -G_matrix(i, 6+k);
				}
				
			}

			
			//now the condensed block terms:
			//the condensed block has 4 submatrices:    1 [ K+M*] [ G*+]  2 
			//											3 [ D *+] [ L* ]  4
												
			//first block
			for (unsigned int i = 0; i < enrich_velocity_dofs; i++) //we go through the 3 enrichments
				for (unsigned int k = 0; k < enrich_velocity_dofs; k++) //we go through the 3 enrichments
					condensed_block(i,k)=Momentum_matrix(9+i,9+k);
			//second block
			for (unsigned int i = 0; i < enrich_velocity_dofs; i++) //we go through the 3 enrichments
				for (unsigned int k = 0; k < enrich_pressure_dofs; k++) //we go through the 3 enrichments
					condensed_block(i,k+enrich_velocity_dofs)= -G_matrix_mixed(k,6+i);		// in  this case, we are in the gradient side and we should use the gradient if we do not integrate it by parts.		
			//third block
			for (unsigned int i = 0; i < enrich_pressure_dofs; i++) //we go through the 3 enrichments
				for (unsigned int k = 0; k < enrich_velocity_dofs; k++) //we go through the 3 enrichments
					condensed_block(i+enrich_velocity_dofs,k)= -G_matrix_mixed(i,6+k);		// in  this case, we are in the divergence side and we should use the gradient if we want to integrate the divergence it by parts.
			
			//fourth block
			for (unsigned int i = 0; i < enrich_pressure_dofs; i++) //we go through the 3 enrichments
				for (unsigned int k = 0; k < enrich_pressure_dofs; k++) //we go through the 3 enrichments
					condensed_block(i+enrich_velocity_dofs,k+enrich_velocity_dofs)=Laplacian_enrich(i,k);		//
					
					
					
			//inverting
			boost::numeric::ublas::bounded_matrix<double, enrich_velocity_dofs+enrich_pressure_dofs , enrich_velocity_dofs+enrich_pressure_dofs  > inverse_enrichments;
			this->InvertMatrix(condensed_block,inverse_enrichments);
			
			//condensing
			boost::numeric::ublas::bounded_matrix<double, 9 , enrich_pressure_dofs+enrich_velocity_dofs  > temp_matrix;
			temp_matrix = prod(trans(condensed_columns),inverse_enrichments);
			rLeftHandSideMatrix -=  prod(temp_matrix,condensed_rows);
			noalias(rRightHandSideVector) -= prod(temp_matrix,rhs_enrich);

			//standard part
			for (unsigned int i = 0; i < 3; i++)
			{
				for (unsigned int j = 0; j < 3; j++)
				{
					rLeftHandSideMatrix(i*3+2, j*3+0 ) -= D_matrix(i,j*2);//  + DN_DX(i,0)*TauOne*one_third*Area/delta_t;     
					rLeftHandSideMatrix(i*3+2, j*3+1 ) -= D_matrix(i,j*2+1);// + DN_DX(i,1)*TauOne*one_third*Area/delta_t;
					
					rLeftHandSideMatrix(j*3+0, i*3+2 ) -=  D_matrix(i,j*2);     
					rLeftHandSideMatrix(j*3+1, i*3+2 ) -= D_matrix(i,j*2+1);
					
					rLeftHandSideMatrix(i*3+2, j*3+2 ) += Laplacian_matrix(i,j);
				}
			}

			
			//and finally the RHS of the velocity plus the RHS on the pressure due to the stabilzation terms:
			for (unsigned int i = 0; i < 3; i++)
			{
				rRightHandSideVector(i*3+0) +=local_gravity(0,0)*lumped_mass(i);
				rRightHandSideVector(i*3+1) +=local_gravity(1,0)*lumped_mass(i);
				
				rRightHandSideVector(i*3+0) += lumped_mass(i)* local_previous_vel(i*2) / delta_t;
				rRightHandSideVector(i*3+1) += lumped_mass(i)* local_previous_vel(i*2+1) / delta_t;
				
				rRightHandSideVector(i*3+2) -= TauOne*(DN_DX_in_local_axis(i,0)*local_gravity(0,0)+DN_DX_in_local_axis(i,1)*local_gravity(1,0))*Area;
			}
			
			for (unsigned int i = 0; i < 9; i++) 
				for (unsigned int j = 0; j < 9; j++)
					rLeftHandSideMatrix(i,j) += Momentum_matrix(i,j);
			

			
			//having created the matrix in local coordiantes, we rotate it into the new system
			//we rotate it by blocks:
			boost::numeric::ublas::bounded_matrix<double, 9 , 9 > LocalLeftHandSideMatrix=rLeftHandSideMatrix;
			array_1d<double,9>  LocalRightHandSideVector = rRightHandSideVector;
			
			boost::numeric::ublas::bounded_matrix<double, 3 , 3 > L_rotation_matrix=ZeroMatrix(3,3);
			for (unsigned int k = 0; k < 2; k++) 
				for (unsigned int l = 0; l < 2; l++)
					L_rotation_matrix(k,l)=rRotationMatrix(k,l);
			L_rotation_matrix(2,2)=1.0;
					
									
			for (unsigned int i = 0; i < 3; i++) 
			{
				//the momentum(mass)block has to be rotated with L^T*K*L
				
				for (unsigned int j = 0; j < 3; j++)
				{
					//rotationg as a whole (simpler)
					boost::numeric::ublas::bounded_matrix<double, 3 , 3 > LocalBlock;
					//saving unrotated local block
					for (unsigned int k = 0; k < 3; k++) 
						for (unsigned int l = 0; l < 3; l++)
							LocalBlock(k,l)=LocalLeftHandSideMatrix(i*3+k,j*3+l);	
					//rotationg	
					boost::numeric::ublas::bounded_matrix<double, 3 , 3 > tempBlock=prod(LocalBlock,L_rotation_matrix);
					boost::numeric::ublas::bounded_matrix<double, 3 , 3 > GlobalBlock=prod(trans(L_rotation_matrix),tempBlock);
					//saving info into the global matrix
					for (unsigned int k = 0; k < 3; k++) 
						for (unsigned int l = 0; l < 3; l++)
							rLeftHandSideMatrix(i*3+k,j*3+l)=GlobalBlock(k,l);

				}
				//now the RHS.
				boost::numeric::ublas::bounded_matrix<double, 3 , 1 >  LocalRHS;
				for (unsigned int k = 0; k < 3; k++) 
					LocalRHS(k,0) = LocalRightHandSideVector(i*3+k);	
				boost::numeric::ublas::bounded_matrix<double, 3 , 1 > GlobalRHS = prod(trans(L_rotation_matrix),LocalRHS);
				for (unsigned int k = 0; k < 3; k++) 
					rRightHandSideVector(i*3+k) = GlobalRHS(k,0);
				
			}
		}
		
		//******************************
		//finally, for both split and non split elements, since we use a residualbased scheme
		noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix,previous_vel_and_press);
		
		
		
		
		///////////FOR EMPTY ELEMENTS!!!!!!!!!!!!!!
		const int number_of_particles_in_elem = this->GetValue(NUMBER_OF_FLUID_PARTICLES);
		if( (number_of_particles_in_elem==0) )
		{
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
			
			rLeftHandSideMatrix = ZeroMatrix(9,9);
			for (unsigned int i = 0; i < 3; i++)
			{
				if (GetGeometry()[i].FastGetSolutionStepValue(YP)<0.001)
				{
					rLeftHandSideMatrix(i*3+0,i*3+0) = 1.0 * Area/delta_t ;
					rLeftHandSideMatrix(i*3+1,i*3+1) = 1.0 * Area/delta_t ;
					rLeftHandSideMatrix(i*3+2,i*3+2) = -1.0 * Area /(mElemSize*mElemSize)*delta_t;
				}
			}
			rLeftHandSideMatrix*=0.000001;
			noalias(rRightHandSideVector) = -prod(rLeftHandSideMatrix,previous_vel_and_press);
		 }
		KRATOS_CATCH("");
	}
	
	
	
	
	//*****************************************************************************
	//***************************************************************************
	

		void MonolithicPFEM22D::CalculateViscousRHS(ProcessInfo& CurrentProcessInfo)
	{
		KRATOS_TRY
		
		
		if ((this->GetValue(IS_INACTIVE))==false) //elements can be inactive to add temporary walls. fractional velocity is integrated by parts so walls are seen as having zero velocity without doing anything 
				{
					//ProcessInfo& CurrentProcessInfo = mr_model_part.GetProcessInfo();
					//double viscosity = CurrentProcessInfo[VISCOSITY];
					double viscosity_air = CurrentProcessInfo[VISCOSITY_AIR]; // * (1.0-theta) ;
					double viscosity_water = CurrentProcessInfo[VISCOSITY_WATER]; // * (1.0-theta) ;
					double mDENSITY_AIR = CurrentProcessInfo[DENSITY_AIR]; // * (1.0-theta) ;
					double mDENSITY_WATER = CurrentProcessInfo[DENSITY_WATER]; // * (1.0-theta) ;
					double delta_t = CurrentProcessInfo[DELTA_TIME];
					//array_1d<double,3> & gravity= CurrentProcessInfo[GRAVITY];
					//double x_force  = gravity(0);
					//double y_force  = gravity(1);	
					const array_1d<double,3> zero3 = ZeroVector(3);
					const double nodal_weight = 1.0/ (1.0 + double (2) );
					
					
					
					boost::numeric::ublas::bounded_matrix<double, (2-1)*6, (2-1)*6 > Viscosity_matrix = ZeroMatrix((2-1)*6, (2-1)*6); 
					boost::numeric::ublas::bounded_matrix<double, (2+1), 3 > velocities = ZeroMatrix((2+1), 3);
										
					boost::numeric::ublas::bounded_matrix<double, (2+1), 2 > DN_DX;
					array_1d<double, (2+1) > N;
					array_1d<double, 2 > vel_gauss;
					double Area;
					Geometry<Node<3> >& geom = this->GetGeometry();
					GeometryUtils::CalculateGeometryData(geom, DN_DX, N, Area);

					
					if ((this->GetValue(SPLIT_ELEMENT))==true)
					{
						//TOOLS NEEDED TO FIND THE ENRICHMENT SHAPE FUNCTIONS
						//get position of the cut surface
						//const unsigned int LocalSize = geom.size()*2;
						unsigned int TNumNodes = geom.size();
						array_1d<double,(2+1)> distances;
						boost::numeric::ublas::bounded_matrix<double,3*(2-1), 2> Nenriched;
						array_1d<double,(3*(2-1))> volumes;
						array_1d<double,(3*(2-1))> densities;
						array_1d<double,(3*(2-1))> inv_densities;
						boost::numeric::ublas::bounded_matrix<double,(2+1), 2 > coords;
						boost::numeric::ublas::bounded_matrix<double, 3*(2-1), (2+1) > Ngauss;
						array_1d<double,(3*(2-1))> signs;
						std::vector< Matrix > gauss_gradients(3*(2-1));

						//fill coordinates
					   
						//unsigned int single_triangle_node;
						array_1d<unsigned int, 2+1 > fixed_nodes; //unordered : i position in the array might not correspond to i node of the element
						array_1d<bool, 2+1 > is_node_fixed; //i position belongs to i node of the element
						unsigned int number_of_fixed_nodes=0;
						bool boundary_element=false;
						
						for (unsigned int i = 0; i < TNumNodes; i++)
						{
							const array_1d<double, 3 > & xyz = geom[i].Coordinates();
							const array_1d<double, 3 > & velocity = geom[i].FastGetSolutionStepValue(VELOCITY);
							
							volumes[i] = 0.0;
							distances[i] = geom[i].FastGetSolutionStepValue(DISTANCE);
							for (unsigned int j = 0; j < (2); j++)
							{
								coords(i, j) = xyz[j];
								velocities(i,j) = velocity[j];
							}
							
							//to find out stress when we have cut elements:
							if (2==2)
							{	
								if (geom[i].IsFixed(FRACT_VEL_X) && geom[i].IsFixed(FRACT_VEL_Y))
								{
									fixed_nodes[number_of_fixed_nodes]=i;
									is_node_fixed[i]=true;
									number_of_fixed_nodes++;
								}
								else
									is_node_fixed[i]=false;
							}
							else // (2==3)
							{	
								if (geom[i].IsFixed(FRACT_VEL_X) && geom[i].IsFixed(FRACT_VEL_Y) && geom[i].IsFixed(FRACT_VEL_Z) )
								{
									fixed_nodes[number_of_fixed_nodes]=i;
									number_of_fixed_nodes++;
									is_node_fixed[i]=true;
								}
								else
									is_node_fixed[i]=false;
							}
							
						}

						for (unsigned int i = 0; i < 3*(2-1) ; i++)
							gauss_gradients[i].resize(2, (2), false);  //2 values of the 2 shape functions, and derivates in xy(z) direction).
						//calling the enrichment function
						unsigned int ndivisions;
							ndivisions = EnrichmentUtilities::CalculateEnrichedShapeFuncions(coords, DN_DX, distances, volumes, Ngauss, signs, gauss_gradients, Nenriched); //, face_gauss_N, face_gauss_Nenriched, face_Area, face_n  ,single_triangle_node);

						for (unsigned int i = 0; i < ndivisions; i++)
						{
							//in case we've changed the distance function inside the enrichmentutility:
							///geom[i].FastGetSolutionStepValue(DISTANCE)=distances[i];
							
							//we are in a loop so we also change the sign of each partition
							if (signs[i]>0.0)
							{
								densities(i)= mDENSITY_AIR;
							}
							else
							{
								densities(i)= mDENSITY_WATER;
							}
						}


						//now the viscous forces:

						double Weight=0.0;
						double mass=0.0;
						for (unsigned int i=0;i<ndivisions;i++) //we go over the three partitions;
						{
							mass += volumes(i)*densities(i);
							if (signs(i)>0.0)
							{
								Weight += volumes(i)*viscosity_air;
								
							}
							else
							{
								Weight += volumes(i)*viscosity_water;
							}
						}
						
						AddViscousTerm(Viscosity_matrix, DN_DX,
										Weight);
						
						const double element_viscosity=Weight/Area; 
						const double element_density=mass/Area;
						double fourier= element_viscosity/element_density*delta_t/pow((0.6*this->GetValue(MEAN_SIZE)),2);
						//KRATOS_WATCH(fourier)
						double theta = 0.2*fourier - 0.1;
						if (theta<0.0) theta=0.0;
						else if (theta>1.0) theta=1.0;
						///******************************
						Viscosity_matrix*=(1.0-theta);
						///******************************


						
						double plane_point_distance=1.0;
						double fixed_face_area_or_lenght=0.0;
						array_1d<double, 3 > boundary_stress;
						if (number_of_fixed_nodes==2) //it means we have cutted elements!
						{
							boundary_element=true;
							array_1d<double, 3 > normal;
							unsigned int free_node=0;
							if (2==2)
							{
								fixed_face_area_or_lenght = fabs(sqrt(pow((geom[fixed_nodes[1]].Y()-geom[fixed_nodes[0]].Y()),2 ) + pow( (geom[fixed_nodes[1]].X()-geom[fixed_nodes[0]].X() ),2 ) ) );
								normal[0] = geom[fixed_nodes[1]].Y()-geom[fixed_nodes[0]].Y();
								normal[1] = - ( geom[fixed_nodes[1]].X()-geom[fixed_nodes[0]].X() );
								normal[2] = 0.0;
								normal /= sqrt(normal[0]*normal[0]+normal[1]*normal[1]);
								
								if (fixed_nodes[0]==0)
								{
									if (fixed_nodes[1]==1)
										free_node=2;
									else
										free_node=1;
								}
								else
									free_node=0;
								
								//the plane is composed by the unit normal and any of the points of fixed nodes. we will use fixed_nodes[0];	
								plane_point_distance = inner_prod( (geom[free_node].Coordinates()-geom[fixed_nodes[0]].Coordinates()) , normal);
								//boundary_stress = geom[free_node].FastGetSolutionStepValue(VELOCITY)*viscosity/plane_point_distance;
							}
							else //(2==3)
							{
								//the area is obtained from the crossproduct of the 2 vertices:
								normal = MathUtils<double>::CrossProduct((geom[fixed_nodes[1]].Coordinates()-geom[fixed_nodes[0]].Coordinates()),(geom[fixed_nodes[2]].Coordinates()-geom[fixed_nodes[0]].Coordinates()) );
								fixed_face_area_or_lenght = 0.5 * sqrt( pow(normal[0],2) + pow(normal[1],2) + pow(normal[2],2) );
								normal /= 2.0 * fixed_face_area_or_lenght;  //this way it is a unit vector. now we must find the distance from the plane generated by the triangles to the free node:
								
								//fixed_face_area_or_lenght = fabs(fixed_face_area_or_lenght);
								
								for (unsigned int j=0; j!=(2+1); j++)
								{
									if (is_node_fixed[j]==false)
									{
										free_node=j;
										break;
									}
								}
								
								//the plane is composed by the unit normal and any of the points of fixed nodes. we will use fixed_nodes[0];
								plane_point_distance = inner_prod( (geom[free_node].Coordinates()-geom[fixed_nodes[0]].Coordinates()) , normal);
								//boundary_stress = geom[free_node].FastGetSolutionStepValue(VELOCITY)*viscosity/plane_point_distance;
							}
							
							boundary_stress = - geom[free_node].FastGetSolutionStepValue(VELOCITY)/(pow(plane_point_distance,2)); //have in mind we could not add the viscosity yet since it changes at each node!
							//KRATOS_WATCH(plane_point_distance)
							//KRATOS_WATCH(boundary_stress)
							//KRATOS_WATCH(fixed_face_area_or_lenght)
							
						}


						double current_nodal_mass=0.0;
						
						for (unsigned int j=0; j!=(2+1); j++)
						{
							//geom[j].SetLock();
							
							if (boundary_element==true && is_node_fixed[j]==true)
							{
								double density;
								//double viscosity;
								if (distances[j]<0.0)
								{
									density = mDENSITY_WATER;
									//viscosity =viscosity_water;
								}
								else
								{
									density = mDENSITY_AIR;
									//viscosity=viscosity_air;
								}
								geom[j].FastGetSolutionStepValue(NODAL_AREA) += fixed_face_area_or_lenght*0.5;
								geom[j].FastGetSolutionStepValue(NODAL_MASS) += fixed_face_area_or_lenght*0.5*density;
								//array_1d<double, 3 > & current_rhs = geom[j].FastGetSolutionStepValue(RHS);
								//current_rhs += boundary_stress*fixed_face_area_or_lenght*0.5*1.0*viscosity;
								//NO RHS in these nodes!
							}
							
							else if (is_node_fixed[j]==false)//normal viscosity procedure:
							{
								geom[j].FastGetSolutionStepValue(NODAL_AREA) += Area * nodal_weight;
							
								array_1d<double, 3 > & current_rhs = geom[j].FastGetSolutionStepValue(RHS);
								for (unsigned int i=0;i!=(2+1);i++) //neighbour node
									for (unsigned int k=0;k!=(2);k++) //k component of local stress
										for (unsigned int l=0;l!=(2);l++) //affected by l component in i neighbour node
											current_rhs[k] += Viscosity_matrix(j*2+k,i*2+l)*velocities(i,l);
								
								current_nodal_mass=0.0;
								for (unsigned int i=0;i!=ndivisions;i++)
									current_nodal_mass += volumes(i)*Ngauss(i,j)*densities(i);
								geom[j].FastGetSolutionStepValue(NODAL_MASS) += current_nodal_mass;
								geom[j].FastGetSolutionStepValue(TAU)=theta;
							}
							//geom[j].UnSetLock();
						}
					}
					else //uncut(normal) elem
					{
						double density;
						double viscosity;
						if (geom[0].FastGetSolutionStepValue(DISTANCE)<0.0)
						{
							viscosity=viscosity_water;
							density = mDENSITY_WATER;
						}
						else
						{
							viscosity=viscosity_air;
							density = mDENSITY_AIR;
						}
						
						AddViscousTerm(Viscosity_matrix, DN_DX,
										Area*viscosity); //, Ngauss,
						
						
						double fourier= viscosity/density * delta_t/pow((0.6*this->GetValue(MEAN_SIZE)),2);
						double theta = 0.2*fourier - 0.1;
						if (theta<0.0) theta=0.0;
						else if (theta>1.0) theta=1.0;
						///******************************
						Viscosity_matrix*=(1.0-theta);
						///******************************
						
						
						array_1d<unsigned int, 2+1 > fixed_nodes; //unordered : i position in the array might not correspond to i node of the element
						array_1d<bool, 2+1 > is_node_fixed; //i position belongs to i node of the element
						unsigned int number_of_fixed_nodes=0;
						bool boundary_element=false;
						
						
						for (unsigned int i = 0; i < (2+1); i++)
						{
							const array_1d<double, 3 > & velocity = geom[i].FastGetSolutionStepValue(VELOCITY);
							for (unsigned int j = 0; j < (2); j++)
								velocities(i,j) = velocity[j];
							

							if (geom[i].IsFixed(VELOCITY) && geom[i].IsFixed(VELOCITY))
							{
								fixed_nodes[number_of_fixed_nodes]=i;
								is_node_fixed[i]=true;
								number_of_fixed_nodes++;
							}
							else
								is_node_fixed[i]=false;
							
						}
						
						double plane_point_distance=1.0;
						double fixed_face_area_or_lenght=0.0;
						array_1d<double, 3 > boundary_stress;
						if (number_of_fixed_nodes==2) //it means we have cutted elements!
						{
							boundary_element=true;
							array_1d<double, 3 > normal;
							unsigned int free_node=0;
	
							fixed_face_area_or_lenght = fabs(sqrt(pow((geom[fixed_nodes[1]].Y()-geom[fixed_nodes[0]].Y()),2 ) + pow( (geom[fixed_nodes[1]].X()-geom[fixed_nodes[0]].X() ),2 ) ) );
							normal[0] = geom[fixed_nodes[1]].Y()-geom[fixed_nodes[0]].Y();
							normal[1] = - ( geom[fixed_nodes[1]].X()-geom[fixed_nodes[0]].X() );
							normal[2] = 0.0;
							normal /= sqrt(normal[0]*normal[0]+normal[1]*normal[1]);
							
							if (fixed_nodes[0]==0)
							{
								if (fixed_nodes[1]==1)
									free_node=2;
								else
									free_node=1;
							}
							else
								free_node=0;
							
							//the plane is composed by the unit normal and any of the points of fixed nodes. we will use fixed_nodes[0];	
							plane_point_distance = inner_prod( (geom[free_node].Coordinates()-geom[fixed_nodes[0]].Coordinates()) , normal);
							//boundary_stress = geom[free_node].FastGetSolutionStepValue(VELOCITY)*viscosity/plane_point_distance;
							
							
							boundary_stress = - geom[free_node].FastGetSolutionStepValue(VELOCITY)*viscosity/(pow(plane_point_distance,2));
							//KRATOS_WATCH(plane_point_distance)
							//KRATOS_WATCH(boundary_stress)
							//KRATOS_WATCH(fixed_face_area_or_lenght)
							
						}
						
						for (unsigned int j=0; j!=(2+1); j++)
						{
							//geom[j].SetLock();
							if (boundary_element==true && is_node_fixed[j]==true)
							{
								geom[j].FastGetSolutionStepValue(NODAL_AREA) += fixed_face_area_or_lenght*0.5;
								geom[j].FastGetSolutionStepValue(NODAL_MASS) += fixed_face_area_or_lenght*0.5*density;
								//array_1d<double, 3 > & current_rhs = geom[j].FastGetSolutionStepValue(RHS);
								//current_rhs += boundary_stress*fixed_face_area_or_lenght*0.5*1.0;

								
								
							}
							else if (is_node_fixed[j]==false)//normal viscosity procedure:
							{
								geom[j].FastGetSolutionStepValue(NODAL_AREA) += Area*nodal_weight;
														
								array_1d<double, 3 > & current_rhs = geom[j].FastGetSolutionStepValue(RHS);
								for (unsigned int i=0;i!=(2+1);i++) //neighbour node
									for (unsigned int k=0;k!=(2);k++) //k component of local stress
										for (unsigned int l=0;l!=(2);l++) //affected by l component in i neighbour node
											current_rhs[k] += Viscosity_matrix(j*2+k,i*2+l)*velocities(i,l);
								
								geom[j].FastGetSolutionStepValue(NODAL_MASS) += Area*density*nodal_weight;
								geom[j].FastGetSolutionStepValue(TAU) = theta;
							}
							//geom[j].UnSetLock();
						}
					}
				}
		
		KRATOS_CATCH("");
	}
	
	//*****************************************************************************
	//***************************************************************************
	

	void MonolithicPFEM22D::CalculatePressureProjection(ProcessInfo& CurrentProcessInfo)
	{
		KRATOS_TRY
		
		
		const double mDENSITY_AIR = CurrentProcessInfo[DENSITY_AIR]; // * (1.0-theta) ;
		const double mDENSITY_WATER = CurrentProcessInfo[DENSITY_WATER]; // * (1.0-theta) ;
		
				
		double Area;
		Geometry<Node<3> >& geom = this->GetGeometry();
		boost::numeric::ublas::bounded_matrix<double, (2+1), 2 > DN_DX;
		array_1d<double, (2+1) > N;
		GeometryUtils::CalculateGeometryData(geom, DN_DX, N, Area);
		const double mass_factor = 1.0/ (1.0 + double (2) );
		
		const int number_of_particles_in_elem = this->GetValue(NUMBER_OF_FLUID_PARTICLES);
		
		if( (number_of_particles_in_elem>0))
		{
					array_1d<double,3>  pressures = ZeroVector(3); //to calculate the deformation Gradient F. Dimension = velocity dofs
					boost::numeric::ublas::bounded_matrix<double,3, 3 > coords; //coordinates of the nodes
					bool has_negative_node=false;
					bool has_positive_node=false;
					double element_mean_distance=0.0;
					
					for(unsigned int iii = 0; iii<3; iii++)
					{
						//saving everything
						pressures(iii) = GetGeometry()[iii].FastGetSolutionStepValue(PRESSURE);

						//				
						const array_1d<double, 3 > & xyz = this->GetGeometry()[iii].Coordinates();
						for (unsigned int j = 0; j < 3; j++)
							coords(iii, j) = xyz[j];
						//
						if (this->GetGeometry()[iii].FastGetSolutionStepValue(DISTANCE)<0.0)
							has_negative_node=true;
						else
							has_positive_node=true;		
							
						element_mean_distance+=  this->GetGeometry()[iii].FastGetSolutionStepValue(DISTANCE)*0.3333333333;
					}
		
					bool split_element=false;
					if (has_positive_node && has_negative_node)
						split_element=true;
			
					if (split_element==false) //we only calculate if we have a pure solid/fluid element to make things easier
					{

						//const double water_fraction = -( element_mean_distance - 1.0) *0.5; 
						const double water_fraction = 0.5*(1.0-(sin(3.14159*element_mean_distance*0.5)));
						double density = mDENSITY_WATER*(water_fraction)+mDENSITY_AIR*(1.0-water_fraction);
						//if (has_positive_node) density = mDENSITY_AIR;
						//else density = mDENSITY_WATER ;
						
						
						
						boost::numeric::ublas::bounded_matrix<double, (2+1), (2-1)*6 > G_matrix; //(gradient)
						noalias(G_matrix) = ZeroMatrix((2+1), (2-1)*6);	
						
						

						for (unsigned int i = 0; i < (2+1); i++) //loop in shape functions (row)
						{
							for (unsigned int j = 0; j < (2+1) ; j++) //loop in shape functions (column)
							{
								for (unsigned int k = 0; k < (2) ; k++) //x,y,(z)
								{
									G_matrix(i, (j*2)+k ) = DN_DX(i,k)*Area*mass_factor; //mass_factor=(1/3 in 2d, 1/4 in 3d)
								}
							}
						}
						G_matrix /=density;

						//now we save the data:
						for (unsigned int i=0; i!=(2+1); i++) //loop around the nodes of the element to add contribution to node i
						{
							array_1d<double, 3 > & current_press_proj = geom[i].FastGetSolutionStepValue(PRESS_PROJ);
								
							for (unsigned int j = 0; j < (2+1) ; j++) //loop around the nodes of the element to add contribution of node j to node i
							{
								
								for (unsigned int k = 0; k < (2) ; k++) //x,y,(z)
								{
									current_press_proj[k] += G_matrix(j,i*(2)+k)*(pressures(j));///-old_pressures(j)); //gamma=0!
								}
							}
							
							geom[i].FastGetSolutionStepValue(NODAL_AREA) += Area*mass_factor;
							
						}

					} //closing the useful elements
					else// we add some small addition to the area so that we do not have a division by zero.
					{
						for (unsigned int i=0; i!=(2+1); i++) //loop around the nodes of the element to add contribution to node i
						{
							geom[i].FastGetSolutionStepValue(NODAL_AREA) += Area*mass_factor*0.000000000001;							
						}
					}
					
		} //closing the if(is_inactive==false)
		else// we add some small addition to the area so that we do not have a division by zero.
		{
			for (unsigned int i=0; i!=(2+1); i++) //loop around the nodes of the element to add contribution to node i
			{
				geom[i].FastGetSolutionStepValue(NODAL_AREA) += Area*mass_factor*0.000000000001;							
			}
		}
		
		KRATOS_CATCH("");
	}

	//*****************************************************************************
	//***************************************************************************

	void MonolithicPFEM22D::InitializeSolutionStep(ProcessInfo& CurrentProcessInfo)
	{
		KRATOS_TRY

		KRATOS_CATCH("");
	}
	
	//************************************************************************************
	//************************************************************************************
	void MonolithicPFEM22D::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo)
	{
		const SizeType NumNodes = 3;
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
		}
	}

	//************************************************************************************
	//************************************************************************************
	void MonolithicPFEM22D::GetDofList(DofsVectorType& ElementalDofList,ProcessInfo& CurrentProcessInfo)
	{
		const SizeType NumNodes = 3;
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
		}
	}
	
	//************************************************************************************
	
	//************************************************************************************
	
	//non partitioned elements using MatrixType
	void MonolithicPFEM22D::AddViscousTerm(MatrixType& OutputMatrix,
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
		C_matrix*=Weight;
		
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
		
	}
	

	//using bounded_matrix, for the explict step:
	void MonolithicPFEM22D::AddViscousTerm(boost::numeric::ublas::bounded_matrix<double, (2-1)*6, (2-1)*6 >& rDampMatrix,
                         boost::numeric::ublas::bounded_matrix<double, (2+1), 2 >& rShapeDeriv,
                         const double Weight)
	{

		
		boost::numeric::ublas::bounded_matrix<double, (2-1)*6,(2-1)*3 > B_matrix = ZeroMatrix((2-1)*6,(2-1)*3);
		for (unsigned int i=0; i!=(2+1); i++) //i node
		{
			for (unsigned int j=0; j!=(2); j++) //x,y,z
				B_matrix(i*(2)+j,j)=rShapeDeriv(i,j);
			
			//useful for both 2d and 3d:	
			//relating 12 and 21 stresses
			B_matrix(i*(2)+0,2)=rShapeDeriv(i,1);
			B_matrix(i*(2)+1,2)=rShapeDeriv(i,0);
			

		}
		
		int counter=0;
		boost::numeric::ublas::bounded_matrix<double, (2-1)*3, (2-1)*3 > C_matrix = ZeroMatrix((2-1)*3,(2-1)*3);
		
		for (unsigned int i=0; i!=(2); i++)
		{
			C_matrix(counter,counter)=2.0;
			counter++;
		}
		for (unsigned int i=0; i!=((2-2)*2+1); i++)
		{
			C_matrix(counter,counter)=1.0;
			counter++;
		}
		
		C_matrix*= -Weight;
		
		boost::numeric::ublas::bounded_matrix<double, (2-1)*3 , (2-1)*6  > temp_matrix = prod(C_matrix,trans(B_matrix));
		rDampMatrix = prod(B_matrix, temp_matrix );
		
	}
	

	
		//************************************************************************************
	//with rotation and returnging also the row that will be enriched
	void MonolithicPFEM22D::AddViscousTerm(boost::numeric::ublas::bounded_matrix<double, 13, 13 > & output,
						  boost::numeric::ublas::bounded_matrix<double, (2+1), 2 >& rShapeDeriv,
						  array_1d<double,3>&  distances,
                          std::vector< Matrix >& gauss_gradients, 
						  array_1d<double,3>&  viscosities,
						  array_1d<double,3>&  signs,
						  array_1d<double,3>&  volumes ,
						  const unsigned int ndivisions)
	{
		
		boost::numeric::ublas::bounded_matrix<double, 8, 8 > ExtendedDampMatrix=ZeroMatrix(8,8);
		//boost::numeric::ublas::bounded_matrix<double, 8, 8 > rExtendedDampMatrix= ZeroMatrix(8,8);

		boost::numeric::ublas::bounded_matrix<double, 8,3 > B_matrix = ZeroMatrix(8,(2-1)*3);

		
		int counter=0;
		boost::numeric::ublas::bounded_matrix<double, (2-1)*3, (2-1)*3 > C_matrix = ZeroMatrix((2-1)*3,(2-1)*3);
			
		for (unsigned int i=0; i!=(2); i++)
		{
			C_matrix(counter,counter)=2.0;
			counter++;
		}
		for (unsigned int i=0; i!=((2-2)*2+1); i++)
		{
			C_matrix(counter,counter)=1.0;
			counter++;
		}
		
		//now the enriched part:
		//we have to construct (ndivisions) rTempDampMatrix and asseble add them to rDampMatrix
		for (unsigned int division=0; division!=ndivisions; division++)
		{
			//standard shape functions:
			for (unsigned int i=0; i!=(2+1); i++) //i node
			{

				for (unsigned int j=0; j!=(2); j++) //x,y,z
					B_matrix(i*(2)+j,j)=rShapeDeriv(i,j);
				//useful for both 2d and 3d:	
				//relating 12 and 21 stresses
				B_matrix(i*(2)+0,2)=rShapeDeriv(i,1);
				B_matrix(i*(2)+1,2)=rShapeDeriv(i,0);
			}
			
			for (unsigned int j=0; j!=(2); j++) //x,y,z
				B_matrix(3*(2)+j,j)= gauss_gradients[division](3,j);
			//useful for both 2d and 3d:	
			//relating 12 and 21 stresses
			B_matrix(3*(2)+0,2)=gauss_gradients[division](3,1);
			B_matrix(3*(2)+1,2)=gauss_gradients[division](3,0);

			
			boost::numeric::ublas::bounded_matrix<double, (2-1)*3 , 8  > temp_matrix = prod(C_matrix,trans(B_matrix));
			ExtendedDampMatrix += viscosities(division)*volumes(division)*prod(B_matrix, temp_matrix );
		}
		
		//now we put it all toghether in the big matrix:
		for (unsigned int i=0; i!=(4); i++) //3 nodes + 1dof in the new virtual node
			for (unsigned int j=0; j!=(4); j++) //3 nodes + 1dof in the new virtual node
				for (unsigned int k=0; k!=(2); k++) //x,y,(z)
					for (unsigned int l=0; l!=(2); l++) //x,y,(z)
						 output(i*3+k,j*3+l) += ExtendedDampMatrix(i*2+k,j*2+l);
						 
		//KRATOS_WATCH(ExtendedDampMatrix)
	}
	
	
	void MonolithicPFEM22D::CalculateInterfaceNormal(boost::numeric::ublas::bounded_matrix<double, 3, 2 >& rPoints, array_1d<double,3>&  rDistances, array_1d<double,2>&  normal, double & interface_area, array_1d<double,3>&  Ninterface)
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
	
	
	inline void MonolithicPFEM22D::CalculatePosition(const bounded_matrix<double, 3, 3 > & coordinates,
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
        /*
        if (N[0] >= 0.0 && N[1] >= 0.0 && N[2] >= 0.0 && N[0] <= 1.0 && N[1] <= 1.0 && N[2] <= 1.0) //if the xc yc is inside the triangle return true
            return true;

        return false;
        */
    }

    inline double MonolithicPFEM22D::CalculateVol(const double x0, const double y0,
                                      const double x1, const double y1,
                                      const double x2, const double y2
                                     )
    {
        return 0.5 * ((x1 - x0)*(y2 - y0)- (y1 - y0)*(x2 - x0));
    }
    
    inline double MonolithicPFEM22D::CalculateVolume2D(
			const bounded_matrix<double, 3, 3 > & coordinates)
		{
			double x10 = coordinates(1,0) - coordinates(0,0);
			double y10 = coordinates(1,1) - coordinates(0,1);

			double x20 = coordinates(2,0) - coordinates(0,0);
			double y20 = coordinates(2,1) - coordinates (0,1);
			double detJ = x10 * y20-y10 * x20;
			return 0.5*detJ;
		}
    
template<class T>
bool MonolithicPFEM22D::InvertMatrix(const T& input, T& inverse)
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

    void MonolithicPFEM22D::CalculateRotationParameters(
			boost::numeric::ublas::bounded_matrix<double,(2+1), 2 >& rOriginalPoints, 
			array_1d<double,(2+1)>& rDistances,
            boost::numeric::ublas::bounded_matrix<double,(2),2 >& rRotationMatrix, 
            boost::numeric::ublas::bounded_matrix<double,(2+1), 2 >& rRotatedPoints,
            boost::numeric::ublas::bounded_matrix<double, (2+1), 2 > & DN_DX_in_local_axis)

        {
        KRATOS_TRY
	
		const double one_third=1.0/3.0;
		boost::numeric::ublas::bounded_matrix<double,3,2> aux_points; //for auxiliary nodes 4(between 1 and 2) ,5(between 2 and 3) ,6 (between 3 and 1)		
		
		
		double most_common_sign=0; //the side of the cut in which two nodes are found (same sign) will be the ones that remains unchanged when builing the discontinuity
		double Area;//area of the complete element
		Area = CalculateVolume2D( rOriginalPoints );
		array_1d<bool,3> cut_edges;
		array_1d<double,3> aux_nodes_relative_locations;
		boost::numeric::ublas::bounded_matrix<int,3,2> aux_nodes_father_nodes;

        //to begin with we must check whether our element is cut or not by the interfase.
        if( (rDistances(0)*rDistances(1))>0.0 && (rDistances(0)*rDistances(2))>0.0 ) //it means that this element IS NOT cut by the interfase. we must return data of a normal, non-enriched element
		{
			return;
		}
		
		//else //we must create the enrichement, it can be in 2 or 3 parts. we'll start with 3 always.

			
		array_1d<double, 3> exact_distance = rDistances;
		array_1d<double, 3> abs_distance = ZeroVector(3);
		
		
		//KRATOS_WATCH("one element IS in the intefase")
		if ((rDistances(0)*rDistances(1))<0.0) //edge 12 is cut
			cut_edges[0]=true;
		else
			cut_edges[0]=false;
		if ((rDistances(1)*rDistances(2))<0.0) //edge 23 is cut. 
			cut_edges[1]=true;
		else
			cut_edges[1]=false;
		if ((rDistances(2)*rDistances(0))<0.0) //edge 23 is cut. 
			cut_edges[2]=true;
		else
			cut_edges[2]=false;
		
		
		
		//'TRICK' TO AVOID HAVING THE INTERFASE TOO CLOSE TO THE NODES:
		//since we cannot collapse node because we have to contemplate the possibility of discontinuities, we will move a little the intefase so that it is not that close.
		const double unsigned_distance0=fabs(rDistances(0));
		const double unsigned_distance1=fabs(rDistances(1));
		const double unsigned_distance2=fabs(rDistances(2));
		//we begin by finding the largest distance:
		double longest_distance=fabs(unsigned_distance0);
		if (unsigned_distance1>longest_distance)
			longest_distance=unsigned_distance1;
		if (unsigned_distance2>longest_distance)
			longest_distance=unsigned_distance2;
		//Now we set a maximum relative distance
		const double tolerable_distance =longest_distance*0.001;	// (1/1,000,000 seems to have good results)
		//and now we check if a distance is too small:
		if (unsigned_distance0<tolerable_distance)
			rDistances[0]=tolerable_distance*(rDistances[0]/fabs(rDistances[0]));
		if (unsigned_distance1<tolerable_distance)
			rDistances[1]=tolerable_distance*(rDistances[1]/fabs(rDistances[1]));
		if (unsigned_distance2<tolerable_distance)
			rDistances[2]=tolerable_distance*(rDistances[2]/fabs(rDistances[2]));
		//END OF TRICK. REMEMBER TO OVERWRITE THE DISTANCE VARIABLE IN THE ELEMENT IN CASE THESE LINES HAVE MODIFIED THEM (distances)
		 
		 
		//for (int jj = 0; jj < 3; jj++)
		//	KRATOS_WATCH(rDistances(jj));
		for (unsigned int i=0; i<3; i++) //we go over the 3 edges:
		{
			int edge_begin_node=i;
			int edge_end_node=i+1;
			if (edge_end_node==3) edge_end_node=0; //it's a triangle, so node 3 is actually node 0
			
			if(cut_edges(i)==true)
			{
				aux_nodes_relative_locations(i)=fabs(rDistances(edge_end_node)/(rDistances(edge_end_node)-rDistances(edge_begin_node) ) ) ; //position in 'natural' coordinates of edge 12, 1 when it passes over node 1. (it is over the edge 01)
				aux_nodes_father_nodes(i,0)=edge_begin_node;
				aux_nodes_father_nodes(i,1)=edge_end_node;
			}
			else
			{
				if(fabs(rDistances(edge_end_node))>fabs(rDistances(edge_begin_node))) //if edge is not cut, we collapse the aux node into the node which has the highest absolute value to have "nicer" (less "slivery") subelements
				{
					aux_nodes_relative_locations(i)=0.0;
					aux_nodes_father_nodes(i,0)=edge_end_node;
					aux_nodes_father_nodes(i,1)=edge_end_node;
				}
				else
				{
					aux_nodes_relative_locations(i)=1.0;
					aux_nodes_father_nodes(i,0)=edge_begin_node;
					aux_nodes_father_nodes(i,1)=edge_begin_node;
				}
			}
			
			//and we save the coordinate of the new aux nodes:
			for (unsigned int j=0;j<2;j++)	//x,y coordinates
				aux_points(i,j)= rOriginalPoints(edge_begin_node,j) * aux_nodes_relative_locations(i) + rOriginalPoints(edge_end_node,j) * (1.0- aux_nodes_relative_locations(i));
		}
		
		//having gathered all the points location in local coordinates, we will proceed to move everything to a local system. the reference point will be the first point of aux_points.
		double x_reference=aux_points(0,0);
		double y_reference=aux_points(0,1);
		double cosinus=0.0;
		double sinus=0.0;
		if (cut_edges[0]==false) //then the segment is defined by aux_points 1 and 2. so the values used to initialize were wrong. changing them:
		{
			x_reference=aux_points(1,0);
			y_reference=aux_points(1,1);
			const double one_over_interfase_lenght = 1.0/sqrt( pow((aux_points(2,0) - x_reference),2) + pow((aux_points(2,1) - y_reference),2));
			cosinus = (aux_points(2,0) - x_reference)*one_over_interfase_lenght;
			sinus = - (aux_points(2,1) - y_reference)*one_over_interfase_lenght; //WARNING; THERE IS A MINUS IN FRONT TO GET THE INVERSE ROTATION (FROM REFERENCE TO LOCAL)
		}
		else if(cut_edges[1]==true) //edges 1 and 0 are cut
		{
			const double one_over_interfase_lenght = 1.0/sqrt( pow((aux_points(1,0) - x_reference),2) + pow((aux_points(1,1) - y_reference),2));
			cosinus = (aux_points(1,0) - x_reference)*one_over_interfase_lenght;
			sinus = - (aux_points(1,1) - y_reference)*one_over_interfase_lenght; //WARNING; THERE IS A MINUS IN FRONT TO GET THE INVERSE ROTATION (FROM REFERENCE TO LOCAL)
		}
		else //edges 2 and 0 are cut
		{
			const double one_over_interfase_lenght = 1.0/sqrt( pow((aux_points(2,0) - x_reference),2) + pow((aux_points(2,1) - y_reference),2));
			cosinus = (aux_points(2,0) - x_reference)*one_over_interfase_lenght;
			sinus = - (aux_points(2,1) - y_reference)*one_over_interfase_lenght; //WARNING; THERE IS A MINUS IN FRONT TO GET THE INVERSE ROTATION (FROM REFERENCE TO LOCAL)
		}
			
		for (unsigned int i=0; i<3; i++) //we go over the 3 nodes and 3 aux nodes to move them to the new coordinates:
		{
			rRotatedPoints(i,0)= cosinus * (rOriginalPoints(i,0)-x_reference) - sinus * (rOriginalPoints(i,1)-y_reference);
			rRotatedPoints(i,1)= cosinus * (rOriginalPoints(i,1)-y_reference) + sinus * (rOriginalPoints(i,0)-x_reference);
			
			//and we directly change the aux points, anyway they are used only locally so it's fine:
			double aux_x_coord=aux_points(i,0);
			aux_points(i,0)= cosinus * (aux_x_coord-x_reference) - sinus * (aux_points(i,1)-y_reference);
			aux_points(i,1)= cosinus * (aux_points(i,1)-y_reference) + sinus * (aux_x_coord-x_reference);
		}
		
		rRotationMatrix(0,0)=cosinus;
		rRotationMatrix(0,1)= sinus;
		rRotationMatrix(1,0)= -sinus;
		rRotationMatrix(1,1)=cosinus;
		
		double temp_area;
		//to calculate the new rigidity matrix in local coordinates, the element will need the derivated in the rotated axis and the rotation matrix:
		CalculateGeometryData(rRotatedPoints, DN_DX_in_local_axis, temp_area);
		KRATOS_CATCH("")
	}
	
	inline void MonolithicPFEM22D::CalculateGeometryData(
			const bounded_matrix<double, 3, 3 > & coordinates,
			boost::numeric::ublas::bounded_matrix<double,3,2>& DN_DX,
			double& Area)
		{
			double x10 = coordinates(1,0) - coordinates(0,0);
			double y10 = coordinates(1,1) - coordinates(0,1);

			double x20 = coordinates(2,0) - coordinates(0,0);
			double y20 = coordinates(2,1) - coordinates (0,1);

			//Jacobian is calculated:
			//  |dx/dxi  dx/deta|	|x1-x0   x2-x0|
			//J=|				|=	|			  |
			//  |dy/dxi  dy/deta|	|y1-y0   y2-y0|


			double detJ = x10 * y20-y10 * x20;

			DN_DX(0,0) = -y20 + y10;
			DN_DX(0,1) = x20 - x10;
			DN_DX(1,0) =  y20	   ;
			DN_DX(1,1) = -x20     ;
			DN_DX(2,0) = -y10	   ;
			DN_DX(2,1) = x10	   ;

			DN_DX /= detJ;

			Area = 0.5*detJ;
		}

	
	//************************************************************************************

} // Namespace Kratos
