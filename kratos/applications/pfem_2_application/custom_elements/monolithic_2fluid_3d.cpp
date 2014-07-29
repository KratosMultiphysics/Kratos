//   
//   Project Name:        Kratos       
//   Last modified by:    $Author: it's me! $
//   Date:                $Date: 2008-08-08 23:58:38 $
//   Revision:            $Revision: 1.0 $
//
//

// Project includes 
#include "includes/define.h"
#include "custom_elements/monolithic_2fluid_3d.h"
#include "pfem_2_application.h"
#include "utilities/math_utils.h"
#include "utilities/geometry_utilities.h" 
#include "utilities/enrichment_utilities.h"

#include "processes/process.h"
#include "includes/node.h"
//#include "includes/element.h"
#include "includes/model_part.h"



namespace Kratos
{
	


	//************************************************************************************
	//************************************************************************************
	MonolithicPFEM23D::MonolithicPFEM23D(IndexType NewId, GeometryType::Pointer pGeometry)
		: Element(NewId, pGeometry)
	{		
		//DO NOT ADD DOFS HERE!!!
	}

	//************************************************************************************
	//************************************************************************************
	MonolithicPFEM23D::MonolithicPFEM23D(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
		: Element(NewId, pGeometry, pProperties)
	{

	}

	Element::Pointer MonolithicPFEM23D::Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const
	{
		return Element::Pointer(new MonolithicPFEM23D(NewId, GetGeometry().Create(ThisNodes), pProperties));
	}

	MonolithicPFEM23D::~MonolithicPFEM23D()
	{
	}

	//************************************************************************************
	//************************************************************************************
	
	
	void MonolithicPFEM23D::AddExplicitContribution(ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY;

		switch ( rCurrentProcessInfo[FRACTIONAL_STEP] )
		{
			default:
			{
				KRATOS_ERROR(std::logic_error,"Unexpected value for FRACTIONAL_STEP index: ",rCurrentProcessInfo[FRACTIONAL_STEP]);
			}
		}

		KRATOS_CATCH("");
	}
	
	

	

	
	
	
	
	
	//************************************************************************************
	//************************************************************************************

	void MonolithicPFEM23D::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY

	    double delta_t = rCurrentProcessInfo[DELTA_TIME];	
	    array_1d<double, 3 > gravity = rCurrentProcessInfo[GRAVITY];	
		//KRATOS_WATCH(gradient_discontinuity);
		const bool & split_element = (this->GetValue(SPLIT_ELEMENT));
		

		double theta=1.0;
	
		
		array_1d<double,12>  previous_vel;
		array_1d<double,16>  previous_vel_and_press;
		array_1d<double, 4 > temperatures;
		array_1d<double,4> distances;
			
        for(unsigned int iii = 0; iii<4; iii++)
        {
			temperatures(iii) = GetGeometry()[iii].FastGetSolutionStepValue(TEMPERATURE);
			distances[iii] = this->GetGeometry()[iii].FastGetSolutionStepValue(DISTANCE);
			array_1d<double,3>& vel = GetGeometry()[iii].FastGetSolutionStepValue(VELOCITY);
			for(unsigned int j = 0; j<3; j++)
			{
				previous_vel(iii*3+j) = vel[j];
				previous_vel_and_press(iii*4+j) = vel[j];
			}
			previous_vel_and_press(iii*4+3) = this->GetGeometry()[iii].FastGetSolutionStepValue(PRESSURE);

		}
		
		const double one_quarter = 0.25;
  		array_1d<double,4> msN; //dimension = number of nodes

		const unsigned int LocalSize = GetGeometry().size()*4;
		unsigned int TNumNodes = GetGeometry().size();

		array_1d<double, 4> N;
        // Check sizes and initialize
        if (rRightHandSideVector.size() != LocalSize)
            rRightHandSideVector.resize(LocalSize, false);
        noalias(rRightHandSideVector) = ZeroVector(LocalSize);
        if(rLeftHandSideMatrix.size1() != LocalSize)
			rLeftHandSideMatrix.resize(LocalSize,LocalSize,false);
		noalias(rLeftHandSideMatrix) = ZeroMatrix(LocalSize,LocalSize);	
		//KRATOS_WATCH(LocalSize)
        // Calculate this element's geometric parameters
        double Area;
        
        boost::numeric::ublas::bounded_matrix<double, 4, 3> DN_DX;

        //boost::numeric::ublas::bounded_matrix<double, TNumNodes*2, 2> N_vel_matrix; //useless to calculate it, it's simply the DN_DX matrix multiplied by 1/3, arranged in a different shape.
		GeometryUtils::CalculateGeometryData(this->GetGeometry(), DN_DX, N, Area);
  
		const double viscosity_air = rCurrentProcessInfo[VISCOSITY_AIR]; //para clindro en rey 100:  0.0005 * 2.0*Area;
		const double viscosity_water = rCurrentProcessInfo[VISCOSITY_WATER];
		const double density_air = rCurrentProcessInfo[DENSITY_AIR]; //para clindro en rey 100:  0.0005 * 2.0*Area;
		const double density_water = rCurrentProcessInfo[DENSITY_WATER];
		


		
		
		
		if (split_element==false)
		//if(true==true)
		{
			//viscous term:
			//double element_temperature=0.0;
			//double density=0.0;
			double element_mean_distance = 0.25*(distances(0)+distances(1)+distances(2)+distances(3));
			const double water_fraction = 0.5*(1.0-(sin(3.14159*element_mean_distance*0.5)));
			double density = density_water*(water_fraction)+density_air*(1.0-water_fraction);
			
			double viscosity=0.0;
			/*
			for (unsigned int j = 0; j < 4; j++) //we go through the 3 standard shape functions
			{
				element_temperature += temperatures[j]*one_quarter;
			}
			*/	
			boost::numeric::ublas::bounded_matrix<double, 4, 12 > D_matrix; //(gradient)
			noalias(D_matrix) = ZeroMatrix(4,12);	
			boost::numeric::ublas::bounded_matrix<double, 4, 4 > Laplacian_matrix;
			noalias(Laplacian_matrix) = ZeroMatrix(4,4);	
			boost::numeric::ublas::bounded_matrix<double, 16, 16 > Mass_matrix;
			noalias(Mass_matrix) = ZeroMatrix(LocalSize,LocalSize);	
			for (unsigned int i = 0; i < 4; i++)
			{
				for (unsigned int j = 0; j < 4; j++)
				{
					for (unsigned int k = 0; k < 3; k++)
						D_matrix(i, (j*3+k) ) =  DN_DX(j,k)*Area*one_quarter; //(i,k) -> standard
						//D_matrix(i, (j*3+k) ) =  -DN_DX(i,k)*Area*one_quarter; //(i,k) -> integrated by parts ( or G)
				}
			}
			
			
			double viscosity_for_tau=0.0;;
			if (distances(0)<0.0)
			{
				density=density_water;
				viscosity=viscosity_water;
				/*
				if (element_temperature>100.0)
						viscosity=viscosity_water*(1.0 - (element_temperature-100.0)*0.1 );
				if (viscosity<1.0)
						viscosity=1.0;
				*/
				//this->AddViscousTerm(rLeftHandSideMatrix,DN_DX, (viscosity*Area) );
				viscosity_for_tau=viscosity_water*0.1;
				
			}
			else
			{
				density = density_air;//this->CalculateAirDensity(element_temperature);
				viscosity = viscosity_air;//this->CalculateAirViscosity(element_temperature);
				//this->AddViscousTerm(rLeftHandSideMatrix,DN_DX, (viscosity*Area) );
				viscosity_for_tau=viscosity_air;
			}
			this->AddViscousTerm(rLeftHandSideMatrix,DN_DX, (viscosity*Area) );
			
			
			
			//CALCULATING TAU (only if we are not near the interfase and we are in step 1 or more)
			double TauOne = 0.0; 
			double sqElemSize=0.0;
			{
				//double Viscosity = rCurrentProcessInfo[VISCOSITY];;
				double AdvVelNorm = sqrt(pow(this->GetGeometry()[0].FastGetSolutionStepValue(VELOCITY_X),0)+pow(this->GetGeometry()[0].FastGetSolutionStepValue(VELOCITY_Y),0)+pow(this->GetGeometry()[0].FastGetSolutionStepValue(VELOCITY_Z),0));
				//const double TimeFactor = rCurrentProcessInfo.GetValue(DYNAMIC_TAU);
				double mElemSize;
				array_1d<double,3> Edge(3,0.0);
				Edge = this->GetGeometry()[1].Coordinates() - this->GetGeometry()[0].Coordinates();
				sqElemSize = Edge[0]*Edge[0];
				for (SizeType d = 1; d < 3; d++)
					sqElemSize += Edge[d]*Edge[d];
				//
				for (SizeType i = 2; i < 4; i++)
					for(SizeType j = 0; j < i; j++)
					{
						Edge = this->GetGeometry()[i].Coordinates() - this->GetGeometry()[j].Coordinates();
						double Length = Edge[0]*Edge[0];
						for (SizeType d = 1; d < 3; d++)
							Length += Edge[d]*Edge[d];
						if (Length < sqElemSize) sqElemSize = Length;
					}
				mElemSize=sqrt(sqElemSize);
				//actually we're currently using the kinematic viscosity, so no need to divide the second term by the viscosity, it was already done. Great!
				TauOne =  1.0 / ( ( 1.0/ delta_t + 4.0 * viscosity_for_tau / (mElemSize * mElemSize * density) + 2.0 * AdvVelNorm / mElemSize) );
			}
		
			/*	
			double fourier=0.0;
			if (distances(0)<0.0) fourier = viscosity_water/density_water *delta_t/pow((0.6*this->GetValue(MEAN_SIZE)),2);
			else fourier = viscosity_air/density_air *delta_t/pow((0.6*this->GetValue(MEAN_SIZE)),2);
			theta = 0.2*fourier - 0.1;
			if (theta<0.0) theta=0.0;
			else if (theta>1.0) theta=1.0;
			*/
				
			const double Weight = one_quarter * Area * density;
			
			for (unsigned int i=0; i<TNumNodes ; i++)
				for (unsigned int j=0; j<3 ; j++)
					Mass_matrix (i*4+j,i*4+j) = Weight;
			
			Laplacian_matrix =  -prod(DN_DX,trans(DN_DX))*Area/density;
			Laplacian_matrix*=TauOne;
			
			for (unsigned int i = 0; i < 4; i++)
			{
				for (unsigned int j = 0; j < 4; j++)
				{	
					for (unsigned int k = 0; k < 3; k++)
					{
						rLeftHandSideMatrix(i*4+3, j*4+k ) -= D_matrix(i,j*3+k);//  + DN_DX(i,0)*TauOne*one_third*Area/delta_t;     
						rLeftHandSideMatrix(j*4+k, i*4+3 ) -=  D_matrix(i,j*3+k);
					}
					rLeftHandSideMatrix(i*4+3, j*4+3 ) += Laplacian_matrix(i,j);
				}
			}
			
			
			//and finally the RHS of the velocity plus the RHS on the pressure due to the stabilzation terms:
			array_1d<double,3> mean_velocity = ZeroVector(3);
			double divergence_n = 0.0;
			for (unsigned int i = 0; i < 4; i++)
			{
				mean_velocity += GetGeometry()[i].FastGetSolutionStepValue(VELOCITY)*one_quarter/delta_t;
				divergence_n += one_quarter*Area*(DN_DX(i,0)*GetGeometry()[i].FastGetSolutionStepValue(VELOCITY_X)+DN_DX(i,1)*GetGeometry()[i].FastGetSolutionStepValue(VELOCITY_Y));
			}	
			for (unsigned int i = 0; i < 4; i++)
			{
				rRightHandSideVector(i*4+0) += one_quarter*Area*gravity(0)*density;
				rRightHandSideVector(i*4+1) += one_quarter*Area*gravity(1)*density;
				rRightHandSideVector(i*4+2) += one_quarter*Area*gravity(2)*density;
			
				rRightHandSideVector(i*4+3) -= TauOne*(DN_DX(i,0)*(gravity(0))+DN_DX(i,1)*(gravity(1))+DN_DX(i,2)*(gravity(2)))*Area;

			}
			
			
			noalias(rRightHandSideVector) += prod((Mass_matrix/delta_t),previous_vel_and_press);
		
			noalias(rLeftHandSideMatrix) += Mass_matrix/delta_t;   
			
			
		} 
		else //split element:
		{
			//TOOLS NEEDED TO FIND THE ENRICHMENT SHAPE FUNCTIONS
			//get position of the cut surface
			//KRATOS_ERROR(std::logic_error, "IMPLICIT STEP FIRST STEP NOT YET IMPLEMENTED IN 3D.. USE LOCALFINALVELOCITY", "");
			//TOOLS NEEDED TO FIND THE ENRICHMENT SHAPE FUNCTIONS
			//get position of the cut surface
			array_1d<double,6>  densities(6);
			array_1d<double,6>  viscosities(6);

			boost::numeric::ublas::bounded_matrix<double,6, 2> Nenriched;
			array_1d<double,6>  volumes(6);
			boost::numeric::ublas::bounded_matrix<double,4, 3 > coords;
			boost::numeric::ublas::bounded_matrix<double,6, 4 > Ngauss;
			array_1d<double,6>  signs(6);
			std::vector< Matrix > gauss_gradients(6);
			//fill coordinates
		   	const int  enrich_velocity_dofs=3;
			const int  enrich_pressure_dofs=2;
			const int  enrich_pressure_offset=0;
			
			boost::numeric::ublas::bounded_matrix<double, 4, 4> Laplacian_matrix =ZeroMatrix(4,4); //standard pressures. we will keep these dof so the size remains			
			boost::numeric::ublas::bounded_matrix<double, 4, 12+enrich_velocity_dofs > D_matrix ; //(divergence) //this matrix will increase its size since we need more  
			noalias(D_matrix) = ZeroMatrix(4,12+enrich_velocity_dofs);	
			boost::numeric::ublas::bounded_matrix<double, 4, 12+enrich_velocity_dofs > G_matrix; //(gradient)
			noalias(G_matrix) = ZeroMatrix(4,12+enrich_velocity_dofs);

			boost::numeric::ublas::bounded_matrix<double, 16+enrich_velocity_dofs+enrich_pressure_dofs, 16+enrich_velocity_dofs+enrich_pressure_dofs > Momentum_matrix; //2 vel + 1 pressure per node   plus 2 enriched pressures and 2 enriched velocities
			noalias(Momentum_matrix) = ZeroMatrix(16+enrich_velocity_dofs+enrich_pressure_dofs,16+enrich_velocity_dofs+enrich_pressure_dofs);	
		   
			//unsigned int single_triangle_node = 1;
			for (unsigned int i = 0; i < TNumNodes; i++)
			{
				const array_1d<double, 3 > & xyz = this->GetGeometry()[i].Coordinates();
				//volumes[i] = 0.0;
				//KRATOS_WATCH(distances(i));
				for (unsigned int j = 0; j < 3; j++)
					coords(i, j) = xyz[j];
			}
			for (unsigned int i = 0; i < 6; i++)
				gauss_gradients[i].resize(2, 3, false);  //2 values of the 2 shape functions, and derivates in (xyz) direction).
				
			unsigned int ndivisions = EnrichmentUtilities::CalculateEnrichedShapeFuncions(coords, DN_DX, distances, volumes, Ngauss, signs, gauss_gradients, Nenriched);
			
			
						
			double Weight=0.0;
			double Weight_for_tau=0.0;
			double mass=0.0;
			for (unsigned int i=0;i<ndivisions;i++) //we go over the three partitions;
			{
				//double partition_temperature=0.0;
				//for (unsigned int j = 0; j < 4; j++)
				//	partition_temperature += temperatures[j]*Ngauss(i,j);
				if (signs(i)>0.0)
				{
					double viscosity = viscosity_air;//CalculateAirViscosity(partition_temperature);
					densities(i) = density_air;//CalculateAirDensity(partition_temperature);
					viscosities(i) = viscosity_air;//CalculateAirDensity(partition_temperature);
					Weight += volumes(i)*viscosity;
					Weight_for_tau += volumes(i)*viscosity;
					mass += volumes(i)*densities(i);
				}
				else
				{
					double viscosity=viscosity_water;
					//if (partition_temperature>100.0)
					//	viscosity=viscosity_water;//*(1.0 - (partition_temperature-100.0)*0.1 );
					//if (viscosity<1.0) //(crude oil viscosity)
					//	viscosity=1.0;
					//viscosities(i)=viscosity;
					Weight += volumes(i)*viscosity;
					Weight_for_tau += volumes(i)*viscosity*0.1;
					densities(i)=density_water;
					viscosities(i) = viscosity_water;//CalculateAirDensity(partition_temperature);
					mass += volumes(i)*density_water;
				}
			}
			
			for (unsigned int i = 0; i < 4; i++)
			{
				for (unsigned int j = 0; j < 4; j++)
				{
					for (unsigned int k = 0; k < 3; k++)
						D_matrix(i, (j*3+k) ) =  DN_DX(j,k)*Area*one_quarter; //(i,k) -> standard
						//D_matrix(i, (j*3+k) ) =  -DN_DX(i,k)*Area*one_quarter; //(i,k) -> integrated by parts ( or G)
				}
			}
			
			this->AddViscousTerm(Momentum_matrix,
								DN_DX,
								distances,
								gauss_gradients, 
								viscosities,
								signs,
								volumes,
								ndivisions);
			
			//this->AddViscousTerm(rLeftHandSideMatrix,DN_DX, (Weight) );

			
			const double element_viscosity=Weight/Area;
			const double element_viscosity_for_tau=Weight_for_tau/Area;
			const double element_density= mass/Area;
			/*
			double fourier= element_viscosity/element_density *delta_t/pow((0.6*this->GetValue(MEAN_SIZE)),2);
			theta = 0.2*fourier - 0.1;
			if (theta<0.0) theta=0.0;
			else if (theta>1.0) theta=1.0;
			*/
			//CALCULATING TAU (only if we are not near the interfase and we are in step 1 or more)
			double TauOne = 0.0; 
			double sqElemSize=0.0;
			
			{
				//double Viscosity = rCurrentProcessInfo[VISCOSITY];;
				double AdvVelNorm = sqrt(pow(this->GetGeometry()[0].FastGetSolutionStepValue(VELOCITY_X),0)+pow(this->GetGeometry()[0].FastGetSolutionStepValue(VELOCITY_Y),0)+pow(this->GetGeometry()[0].FastGetSolutionStepValue(VELOCITY_Z),0));
				//const double TimeFactor = rCurrentProcessInfo.GetValue(DYNAMIC_TAU);
				double mElemSize;
				array_1d<double,3> Edge(3,0.0);
				Edge = this->GetGeometry()[1].Coordinates() - this->GetGeometry()[0].Coordinates();
				sqElemSize = Edge[0]*Edge[0];
				for (SizeType d = 1; d < 3; d++)
					sqElemSize += Edge[d]*Edge[d];
				//
				for (SizeType i = 2; i < 4; i++)
					for(SizeType j = 0; j < i; j++)
					{
						Edge = this->GetGeometry()[i].Coordinates() - this->GetGeometry()[j].Coordinates();
						double Length = Edge[0]*Edge[0];
						for (SizeType d = 1; d < 3; d++)
							Length += Edge[d]*Edge[d];
						if (Length < sqElemSize) sqElemSize = Length;
					}
				mElemSize=sqrt(sqElemSize);
				//actually we're currently using the kinematic viscosity, so no need to divide the second term by the viscosity, it was already done. Great!
				TauOne =  1.0 / ( ( 1.0/ delta_t + 4.0 * element_viscosity_for_tau / (mElemSize * mElemSize * element_density) + 2.0 * AdvVelNorm / mElemSize) );
			}
		
			
			
			array_1d<double,4>  lumped_mass = ZeroVector(4);
			double condensed_dof_mass1=0.0;
			double condensed_dof_mass2=0.0;
			double condensed_dof_mass3=0.0;
			double Laplacian_matrix_factor=0.0;
			for (unsigned int i = 0; i < ndivisions; i++) //partitions 
			{
				//KRATOS_WATCH(enrich_lhs(i));
				for (unsigned int j = 0; j < 4 ; j++) // local node (or shape function)
				{
					
					lumped_mass(j) +=  volumes(i)*Ngauss(i,j)*densities(i);	 
					Momentum_matrix(4*j,4*j) += volumes(i)*Ngauss(i,j)*densities(i) / delta_t;	 
					Momentum_matrix(4*j+1,4*j+1) += volumes(i)*Ngauss(i,j)*densities(i) / delta_t;	 
					Momentum_matrix(4*j+2,4*j+2) += volumes(i)*Ngauss(i,j)*densities(i) / delta_t;	 
				}

				
				if (enrich_velocity_dofs!=0)
				{
					condensed_dof_mass1 += volumes(i)*Nenriched(i,0)*densities(i);
					condensed_dof_mass2 += volumes(i)*Nenriched(i,0)*densities(i);
					condensed_dof_mass3 += volumes(i)*Nenriched(i,0)*densities(i);
					   
					Momentum_matrix(16+0,16+0) += volumes(i)*Nenriched(i,0)*densities(i) / delta_t;	 //degrees of freedom 10 and 11 are the enrichment velocities. degrees of freedom 12 and 13 are the enrichment pressures
					Momentum_matrix(16+1,16+1) += volumes(i)*Nenriched(i,0)*densities(i) / delta_t;	 
					Momentum_matrix(16+2,16+2) += volumes(i)*Nenriched(i,0)*densities(i) / delta_t;	 
				}
				
				Laplacian_matrix_factor +=volumes(i)/densities(i);
				
			}
			
			
			
			Laplacian_matrix =  -prod(DN_DX,trans(DN_DX))*Laplacian_matrix_factor;
			Laplacian_matrix*=TauOne;
			
			
			
			boost::numeric::ublas::bounded_matrix<double, enrich_pressure_dofs, enrich_pressure_dofs > Laplacian_enrich;
			noalias(Laplacian_enrich) = ZeroMatrix(enrich_pressure_dofs,enrich_pressure_dofs);	
			boost::numeric::ublas::bounded_matrix<double, enrich_pressure_dofs, 4 > mixed_Laplacian;
			noalias(mixed_Laplacian) = ZeroMatrix(enrich_pressure_dofs,4);	
			//To calculate D* =int (N*DN_DX_enrich).  we must use the N at the gauss point of the partition (returned by the enrichment util) and multiply by the area. 
			//notice that this first row is only the upper part of the mixed D. we also need the lower one for the jump 
			boost::numeric::ublas::bounded_matrix<double, enrich_pressure_dofs, 12 + enrich_velocity_dofs> D_matrix_mixed;
			noalias(D_matrix_mixed) = ZeroMatrix(enrich_pressure_dofs,12 + enrich_velocity_dofs);	
			boost::numeric::ublas::bounded_matrix<double, enrich_pressure_dofs, 12 + enrich_velocity_dofs > G_matrix_mixed;
			noalias(G_matrix_mixed) = ZeroMatrix(enrich_pressure_dofs,12 + enrich_velocity_dofs);	

			//the matrices we need at the end will be:
			//boost::numeric::ublas::bounded_matrix<double, 2, 2 > D_mod;
			
			//array_1d<double,6> mass_stabilization_terms=ZeroVector(6);
			
			//double rhs_enrich = 0.0; //the rhs term of the enriched dof 
			array_1d<double,( enrich_pressure_dofs + enrich_velocity_dofs )> rhs_enrich;
			noalias(rhs_enrich) = ZeroVector(enrich_pressure_dofs + enrich_velocity_dofs);	
			
				//			KRATOS_WATCH(Ngauss)
				//KRATOS_WATCH(Nenriched)
			
			for (unsigned int i = 0; i < ndivisions; i++)  //we go through the 3 partitions of the domain
			{
				
				//KRATOS_WATCH(gauss_gradients[i])

				//double PartitionTauOne = 1.0 / ( ( 1.0/ delta_t + 4.0 * viscosities(i) / (mElemSize * mElemSize*densities(i)) + 2.0 * AdvVelNorm / mElemSize) );
				
				for (unsigned int j = 0; j < enrich_pressure_dofs; j++) //we go through the 4 enrichments
					for (unsigned int k = 0; k < enrich_pressure_dofs; k++) //we go through the 4 enrichments
						Laplacian_enrich(j,k) -= ((gauss_gradients[i](j+enrich_pressure_offset,0)*gauss_gradients[i](k+enrich_pressure_offset,0))+(gauss_gradients[i](j+enrich_pressure_offset,1)*gauss_gradients[i](k+enrich_pressure_offset,1))+(gauss_gradients[i](j+enrich_pressure_offset,2)*gauss_gradients[i](k+enrich_pressure_offset,2)))*volumes(i)/densities(i);

				//and the 'mixed laplacians' (standard shape functions * enrichments)
				
				
				//array_1d<double,3> previous_velocity_at_gauss_point=ZeroVector(3); //needed for the rhs
				
				//first we take care of the standard shape functions of the velocity
				for (unsigned int j = 0; j < 4; j++) //we go through the 4 standard shape functions
				{
					for (unsigned int k = 0; k < enrich_pressure_dofs; k++) //we go through the 3 enrichments
					{
						mixed_Laplacian(k,j) -= (DN_DX(j,0)*(gauss_gradients[i](k+enrich_pressure_offset,0)) + DN_DX(j,1)*(gauss_gradients[i](k+enrich_pressure_offset,1))+ DN_DX(j,2)*(gauss_gradients[i](k+enrich_pressure_offset,2)))*volumes(i)/densities(i);
						D_matrix_mixed(k,j*3) += DN_DX(j,0) *volumes(i)*Nenriched(i,k+enrich_pressure_offset);
						D_matrix_mixed(k,j*3+1) += DN_DX(j,1) *volumes(i)*Nenriched(i,k+enrich_pressure_offset);
						D_matrix_mixed(k,j*3+2) += DN_DX(j,2) *volumes(i)*Nenriched(i,k+enrich_pressure_offset);
						
						G_matrix_mixed(k,j*3) -= gauss_gradients[i](k+enrich_pressure_offset,0) *volumes(i)*Ngauss(i,j);
						G_matrix_mixed(k,j*3+1) -= gauss_gradients[i](k+enrich_pressure_offset,1) *volumes(i)*Ngauss(i,j);
						G_matrix_mixed(k,j*3+2) -= gauss_gradients[i](k+enrich_pressure_offset,2) *volumes(i)*Ngauss(i,j);
					}
					//previous_velocity_at_gauss_point(0)+=Ngauss(i,j)*local_previous_vel(j*2+0);
					//previous_velocity_at_gauss_point(1)+=Ngauss(i,j)*local_previous_vel(j*2+1);
				}
				
				
				//we go through the remaining, new velocity dofs. (one x and one y component
				if (enrich_velocity_dofs!=0)
				{
					for (unsigned int k = 0; k < enrich_pressure_dofs; k++) //we go through the 3 enrichments
					{
						D_matrix_mixed(k,12+0) += gauss_gradients[i](0,0) *volumes(i)*Nenriched(i,k+enrich_pressure_offset);
						D_matrix_mixed(k,12+1) += gauss_gradients[i](0,1) *volumes(i)*Nenriched(i,k+enrich_pressure_offset);
						D_matrix_mixed(k,12+2) += gauss_gradients[i](0,2) *volumes(i)*Nenriched(i,k+enrich_pressure_offset);
						
						G_matrix_mixed(k,12+0) -= gauss_gradients[i](k+enrich_pressure_offset,0) *volumes(i)*Nenriched(i,0);
						G_matrix_mixed(k,12+1) -= gauss_gradients[i](k+enrich_pressure_offset,1) *volumes(i)*Nenriched(i,0);
						G_matrix_mixed(k,12+2) -= gauss_gradients[i](k+enrich_pressure_offset,2) *volumes(i)*Nenriched(i,0);


					}
					
					
					//also the standard D and G matrices had to be expanded.

						for (unsigned int k = 0; k < 4; k++) //we go through the 4 standard pressures
						{
							D_matrix(k,12+0) += gauss_gradients[i](0,0) *volumes(i)*Ngauss(i,k);
							D_matrix(k,12+1) += gauss_gradients[i](0,1) *volumes(i)*Ngauss(i,k);
							D_matrix(k,12+2) += gauss_gradients[i](0,2) *volumes(i)*Ngauss(i,k);

							G_matrix(k,12+0) -= DN_DX(k,0) *volumes(i)*Nenriched(i,0);
							G_matrix(k,12+1) -= DN_DX(k,1) *volumes(i)*Nenriched(i,0);
							G_matrix(k,12+2) -= DN_DX(k,2) *volumes(i)*Nenriched(i,0);
						}
				}
				
				
				for (unsigned int k = 0; k < enrich_pressure_dofs; k++) //we go through the 4 enrichments
					rhs_enrich(k+enrich_velocity_dofs) -= (gauss_gradients[i](k+enrich_pressure_offset,0)*(gravity(0))+gauss_gradients[i](k+enrich_pressure_offset,1)*(gravity(1))+gauss_gradients[i](k+enrich_pressure_offset,2)*(gravity(2)))*volumes(i);
	
			}
			
			Laplacian_enrich*=TauOne;
			mixed_Laplacian*=TauOne;
			rhs_enrich *=TauOne;
			
			
			if (enrich_velocity_dofs!=0)
			{
				rhs_enrich(0) = gravity(0)*condensed_dof_mass1;
				rhs_enrich(1) = gravity(1)*condensed_dof_mass2;
				rhs_enrich(2) = gravity(2)*condensed_dof_mass3;
			}
			
			boost::numeric::ublas::bounded_matrix<double, enrich_velocity_dofs+enrich_pressure_dofs, 16 > condensed_rows; //Vx1,Vy1,p1,Vx2,...
			boost::numeric::ublas::bounded_matrix<double, enrich_velocity_dofs+enrich_pressure_dofs, 16 > condensed_columns; //Vx1,Vy1,p1,Vx2,...
			boost::numeric::ublas::bounded_matrix<double, enrich_velocity_dofs+enrich_pressure_dofs, enrich_velocity_dofs+enrich_pressure_dofs > condensed_block; //Vx1,Vy1,p1,Vx2,...
			
			for (unsigned int i = 0; i <4; i++)  //we go through the 4 nodes (standard velocity dof + standard pressure dof)
			{
				//enriched pressure dof	
				for (unsigned int k = 0; k < enrich_pressure_dofs; k++)
				{
					
					//enriched
					condensed_rows(enrich_velocity_dofs+k,i*4+0)= - D_matrix_mixed(k,i*3+0);//+mass_stabilization_terms(i*2+0);
					condensed_rows(enrich_velocity_dofs+k,i*4+1)= - D_matrix_mixed(k,i*3+1);//+mass_stabilization_terms(i*2+1);
					condensed_rows(enrich_velocity_dofs+k,i*4+2)= - D_matrix_mixed(k,i*3+2);//+mass_stabilization_terms(i*2+1);
					condensed_rows(enrich_velocity_dofs+k,i*4+3)= mixed_Laplacian(k,i);
					
					condensed_columns(enrich_velocity_dofs+k,i*4+0)= - D_matrix_mixed(k,i*3+0);
					condensed_columns(enrich_velocity_dofs+k,i*4+1)= - D_matrix_mixed(k,i*3+1);
					condensed_columns(enrich_velocity_dofs+k,i*4+2)= - D_matrix_mixed(k,i*3+2);
					condensed_columns(enrich_velocity_dofs+k,i*4+3)= mixed_Laplacian(k,i);
				}
				
				
				
				for (unsigned int k = 0; k < enrich_velocity_dofs; k++) //we go through the 3 enrichments
				{
					condensed_columns(k,i*4+0)=Momentum_matrix(16+k,i*4+0); //add the viscosity matrix
					condensed_columns(k,i*4+1)=Momentum_matrix(16+k,i*4+1); //add the viscosity matrix
					condensed_columns(k,i*4+2)=Momentum_matrix(16+k,i*4+2); //add the viscosity matrix

					condensed_rows(k,i*4+0)=Momentum_matrix(16+k,i*4+0); //add the viscosity matrix
					condensed_rows(k,i*4+1)=Momentum_matrix(16+k,i*4+1); //add the viscosity matrix
					condensed_rows(k,i*4+2)=Momentum_matrix(16+k,i*4+2); //add the viscosity matrix
					
					///WARNING, WHEN THE MATRIX IS REARRANGED, the condensed rows have the gradient of the pressure and the columns have the divergence. that is the reason for the mixed indexes.
					///if the gradient was not integrated by parts, then G matrix should be used in the columns instead of the rows, unlike the case of standard velocity DOFs
					condensed_rows(k,i*4+3)= -G_matrix(i, 12+k);
					condensed_columns(k,i*4+3)= -G_matrix(i, 12+k);
				}
				
			}
			//KRATOS_WATCH(Momentum_matrix)
			
			
			//now the condensed block terms:
			//the condensed block has 4 submatrices:    1 [ K+M*] [ G*+]  2 
			//											3 [ D *+] [ L* ]  4
												
			//first block
			for (unsigned int i = 0; i < enrich_velocity_dofs; i++) //we go through the 3 enrichments
				for (unsigned int k = 0; k < enrich_velocity_dofs; k++) //we go through the 3 enrichments
					condensed_block(i,k)=Momentum_matrix(16+i,16+k);
			//second block
			for (unsigned int i = 0; i < enrich_velocity_dofs; i++) //we go through the 3 enrichments
				for (unsigned int k = 0; k < enrich_pressure_dofs; k++) //we go through the 3 enrichments
					condensed_block(i,k+enrich_velocity_dofs)= -G_matrix_mixed(k,12+i);		// in  this case, we are in the gradient side and we should use the gradient if we do not integrate it by parts.		
			//third block
			for (unsigned int i = 0; i < enrich_pressure_dofs; i++) //we go through the 3 enrichments
				for (unsigned int k = 0; k < enrich_velocity_dofs; k++) //we go through the 3 enrichments
					condensed_block(i+enrich_velocity_dofs,k)= -G_matrix_mixed(i,12+k);		// in  this case, we are in the divergence side and we should use the gradient if we want to integrate the divergence it by parts.
			
			//fourth block
			for (unsigned int i = 0; i < enrich_pressure_dofs; i++) //we go through the 3 enrichments
				for (unsigned int k = 0; k < enrich_pressure_dofs; k++) //we go through the 3 enrichments
					condensed_block(i+enrich_velocity_dofs,k+enrich_velocity_dofs)=Laplacian_enrich(i,k);		//
					
			boost::numeric::ublas::bounded_matrix<double, enrich_velocity_dofs+enrich_pressure_dofs , enrich_velocity_dofs+enrich_pressure_dofs  > inverse_enrichments;
			this->InvertMatrix(condensed_block,inverse_enrichments);
			//condensing
			boost::numeric::ublas::bounded_matrix<double, 16 , enrich_pressure_dofs+enrich_velocity_dofs  > temp_matrix;
			temp_matrix = prod(trans(condensed_columns),inverse_enrichments);
			rLeftHandSideMatrix -=  prod(temp_matrix,condensed_rows);
			noalias(rRightHandSideVector) -= prod(temp_matrix,rhs_enrich);
			
			
			
			for (unsigned int i = 0; i < 4; i++)
			{
				for (unsigned int j = 0; j < 4; j++)
				{	
					for (unsigned int k = 0; k < 3; k++)
					{
						rLeftHandSideMatrix(i*4+3, j*4+k ) -= D_matrix(i,j*3+k);//  + DN_DX(i,0)*TauOne*one_third*Area/delta_t;     
						rLeftHandSideMatrix(j*4+k, i*4+3 ) -=  D_matrix(i,j*3+k);
					}
					rLeftHandSideMatrix(i*4+3, j*4+3 ) += Laplacian_matrix(i,j);
				}
			}
			
			
			//and finally the RHS of the velocity plus the RHS on the pressure due to the stabilzation terms:
			array_1d<double,3> mean_velocity = ZeroVector(3);
			double divergence_n = 0.0;
			for (unsigned int i = 0; i < 4; i++)
			{
				mean_velocity += GetGeometry()[i].FastGetSolutionStepValue(VELOCITY)*one_quarter/delta_t;
				divergence_n += one_quarter*Area*(DN_DX(i,0)*GetGeometry()[i].FastGetSolutionStepValue(VELOCITY_X)+DN_DX(i,1)*GetGeometry()[i].FastGetSolutionStepValue(VELOCITY_Y));
			}	
			for (unsigned int i = 0; i < 4; i++)
			{
				rRightHandSideVector(i*4+0) += lumped_mass(i)*gravity(0);
				rRightHandSideVector(i*4+1) += lumped_mass(i)*gravity(1);
				rRightHandSideVector(i*4+2) += lumped_mass(i)*gravity(2);
				
				rRightHandSideVector(i*4+0) += lumped_mass(i)*previous_vel_and_press(i*4+0)/delta_t;
				rRightHandSideVector(i*4+1) += lumped_mass(i)*previous_vel_and_press(i*4+1)/delta_t;
				rRightHandSideVector(i*4+2) += lumped_mass(i)*previous_vel_and_press(i*4+2)/delta_t;
			
				rRightHandSideVector(i*4+3) -= TauOne*(DN_DX(i,0)*(gravity(0))+DN_DX(i,1)*(gravity(1))+DN_DX(i,2)*(gravity(2)))*Area;

			}
			
			for (unsigned int i = 0; i < 16; i++) 
				for (unsigned int j = 0; j < 16; j++)
					rLeftHandSideMatrix(i,j) += Momentum_matrix(i,j);
			
		}
		
		
       
		
		//we must actually use  a fraction of the rigidity matrix. since the Lefthandside only has it so far, we multiply it by 0.5:
		//rLeftHandSideMatrix *= theta;
	
        //finally we add to the righthandside BOTH the mass matrix and the rigidity matrix:

		
		noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix,previous_vel_and_press);
		
		//KRATOS_WATCH("ART")
		//KRATOS_WATCH(rRightHandSideVector)
		//KRATOS_WATCH(rLeftHandSideMatrix)
		
		KRATOS_CATCH("");
	}
	
	
	
	//*****************************************************************************
	//***************************************************************************

	
	void MonolithicPFEM23D::CalculateLocalThermalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY

		//rCurrentProcessInfo[DELTA_TIME]=1.0;
		//ProcessInfo& CurrentProcessInfo = mr_model_part.GetProcessInfo();
		double Density = 0.0 ; //used for stabilization;
		double delta_t = rCurrentProcessInfo[DELTA_TIME];
		Vector densities(3);
		Vector partition_densities(3);
		//KRATOS_WATCH(delta_t);
		//bool neighbour_of_splitted_element=false;

		double enrich_rhs=0.0;
		double extra_enrich_rhs=0.0;

		double input_jump_value=0.0;
		//this->SetValue(PRESS_DISCONTINUITY,input_jump_value);
		//double jump_value=-input_jump_value*delta_t;
		
		const bool & split_element = (this->GetValue(SPLIT_ELEMENT));
		const double & available_element_oxygen = this->GetValue(AVAILABLE_AIR_VOLUME);

		
		//boost::numeric::ublas::bounded_matrix<double,3,2> msDN_DX;
		//boost::numeric::ublas::bounded_matrix<double,2,2> msD;
		const double one_quarter = 0.25;
  		array_1d<double,4> msN; //dimension = number of nodes
		array_1d<double,4> ms_temp; //dimension = number of nodes

		const unsigned int LocalSize = GetGeometry().size();
		unsigned int TNumNodes = GetGeometry().size();
		//boost::numeric::ublas::bounded_matrix<double, 3,1 > enrich_rhs;

		//getting data for the given geometry
//        const unsigned int LocalSize = (2 + 1) * TNumNodes;

		array_1d<double, 4> N;
        // Check sizes and initialize
        if (rRightHandSideVector.size() != LocalSize)
            rRightHandSideVector.resize(LocalSize, false);
        noalias(rRightHandSideVector) = ZeroVector(LocalSize);
        if(rLeftHandSideMatrix.size1() != LocalSize)
			rLeftHandSideMatrix.resize(LocalSize,LocalSize,false);
		noalias(rLeftHandSideMatrix) = ZeroMatrix(LocalSize,LocalSize);	
        // Calculate this element's geometric parameters
        double Area;
        
        boost::numeric::ublas::bounded_matrix<double, 4, 3> DN_DX;
        boost::numeric::ublas::bounded_matrix<double, 4, 4 > Laplacian_matrix;
        boost::numeric::ublas::bounded_matrix<double, 4, 4 > Mass_matrix=ZeroMatrix(4,4);
        //boost::numeric::ublas::bounded_matrix<double, TNumNodes*2, 2> N_vel_matrix; //useless to calculate it, it's simply the DN_DX matrix multiplied by 1/3, arranged in a different shape.
       
        Geometry<Node<3> >& geom = this->GetGeometry();
         
        GeometryUtils::CalculateGeometryData(geom, DN_DX, N, Area);
        //KRATOS_WATCH(Area);

        array_1d<double, 4 > old_temperatures;
        
        for(unsigned int iii = 0; iii<4; iii++)
        {
			old_temperatures(iii) = GetGeometry()[iii].FastGetSolutionStepValue(TEMPERATURE_OLD_IT);
			//if ((GetGeometry()[iii].GetValue(SPLIT_ELEMENT))==true)
			//	neighbour_of_splitted_element=true;
		}

        if (split_element)
        {

			
			//TOOLS NEEDED TO FIND THE ENRICHMENT SHAPE FUNCTIONS
			//get position of the cut surface
			array_1d<double,4> distances;
			boost::numeric::ublas::bounded_matrix<double,6, 2> Nenriched;
			array_1d<double,6>  volumes(6);
			boost::numeric::ublas::bounded_matrix<double,4, 3 > coords;
			boost::numeric::ublas::bounded_matrix<double,6, 4 > Ngauss;
			array_1d<double,6>  signs(6);
			std::vector< Matrix > gauss_gradients(6);
			//array_1d<double,3> Ngauss_in_interfase;
			//array_1d<double,3> Nenriched_in_interfase;
			//double interfase_area;
			//fill coordinates
		   

			
			//unsigned int single_triangle_node = 1;
			for (unsigned int i = 0; i < TNumNodes; i++)
			{
				const array_1d<double, 3 > & xyz = this->GetGeometry()[i].Coordinates();
				distances[i] = this->GetGeometry()[i].FastGetSolutionStepValue(DISTANCE);
				//KRATOS_WATCH(distances(i));
				for (unsigned int j = 0; j < 3; j++)
					coords(i, j) = xyz[j];
			}

			for (unsigned int i = 0; i < 6; i++)
				gauss_gradients[i].resize(2, 3, false);  //2 values of the 2 shape functions, and derivates in (xyz) direction).
			unsigned int ndivisions = EnrichmentUtilities::CalculateEnrichedShapeFuncions(coords, DN_DX, distances, volumes, Ngauss, signs, gauss_gradients,Nenriched); //, Ngauss_in_interfase,Nenriched_in_interfase,interfase_area);
			
			//in case we've changed the distance function: CANNOT: (PARALLEL)
			//for (unsigned int i = 0; i < TNumNodes; i++)
			//	this->GetGeometry()[i].FastGetSolutionStepValue(DISTANCE)=distances[i];
			
			
			
			array_1d<double,6> conductivities;
			array_1d<double,6> volumetric_heat_capacities;
			array_1d<double,6> densities;
			Density=0.0; //resetting to zero
			double air_partition_density=0.0;
			/*
			for (unsigned int i = 0; i!=ndivisions; i++) //the 3 partitions of the triangle
			{
				//double viscosity=0.0;
				//double volumetric_heat_capacity=0.0;
				double partition_temperature=0.0;
				for (unsigned int j = 0; j < 4; j++) //we go through the 3 standard shape functions
				{
						partition_temperature += old_temperatures[j]*Ngauss(i,j);
				}
				if (signs(i)>0)
				{
					densities[i] = this->CalculateAirDensity(partition_temperature);
					air_partition_density=densities[i];
					conductivities[i] =  this->CalculateAirConductivity(partition_temperature);
					volumetric_heat_capacities[i] = rCurrentProcessInfo[SPECIFIC_HEAT_CAPACITY_AIR]*densities[i];
				}
				else
				{
					densities[i] = rCurrentProcessInfo[DENSITY_WATER];
					conductivities[i] =  rCurrentProcessInfo[PERMEABILITY_WATER];
					volumetric_heat_capacities[i] = rCurrentProcessInfo[SPECIFIC_HEAT_CAPACITY_WATER]*densities[i] ; //around 1670 for plastics
				}

			}
			*/ 
			Density*=1.0/Area;

			double enriched_dof_heat_capacity = 0.0;
			for (unsigned int i = 0; i!=ndivisions; i++) //the 3 partitions of the triangle
			{
				enriched_dof_heat_capacity += volumes(i)*Nenriched(i,0)*volumetric_heat_capacities[i];
				for (unsigned int j = 0; j < 4; j++) //we go through the 3 standard shape functions
				{		
					Mass_matrix(j,j) += volumes(i)*Ngauss(i,j)*volumetric_heat_capacities[i];	 
				}
			}
			

			Laplacian_matrix =  prod(DN_DX,trans(DN_DX));
			double coefficient=0.0;
			for (unsigned int i=0; i<ndivisions ; i++)
					coefficient += volumes[i]*conductivities[i];
			Laplacian_matrix *= coefficient; 

			//and now the rest of the things needed
			
			//boost::numeric::ublas::bounded_matrix<double, 2, 2> DN_DX_enrich;
			//boost::numeric::ublas::bounded_matrix<double, 2, 2> DN_DX_enrich_negative;
			boost::numeric::ublas::bounded_matrix<double, 2, 2 > Laplacian_enrich;
			noalias(Laplacian_enrich) = ZeroMatrix(2,2);	
			//double inv_Laplacian_enrich_weighted;
			//double inv_Laplacian_enrich;
			boost::numeric::ublas::bounded_matrix<double, 1, 4 > mixed_Laplacian;
			noalias(mixed_Laplacian) = ZeroMatrix(1,4);	
			
			//boost::numeric::ublas::bounded_matrix<double, 1, 3 > mixed_Laplacian_jump;
			//noalias(mixed_Laplacian_jump) = ZeroMatrix(1,3);

			
				
			
			//the matrices we need at the end will be:
			//boost::numeric::ublas::bounded_matrix<double, 2, 2 > D_mod;

			for (unsigned int i = 0; i < ndivisions; i++)  //we go through the 3 partitions of the domain
			{
				Laplacian_enrich(0,0) += (pow(gauss_gradients[i](0,0),2)+pow(gauss_gradients[i](0,1),2)+pow(gauss_gradients[i](0,2),2))*volumes(i)*conductivities(i);
				//Laplacian_enrich(0,1) -= ((gauss_gradients[i](0,0)*gauss_gradients[i](1,0))+(gauss_gradients[i](0,1)*gauss_gradients[i](1,1)))*volumes(i)*inv_densities(i);
				//Laplacian_enrich(1,1) -= (pow(gauss_gradients[i](1,0),2)+pow(gauss_gradients[i](1,1),2))*volumes(i)*inv_densities(i);
				//Laplacian_enrich(1,0) -= ((gauss_gradients[i](0,0)*gauss_gradients[i](1,0))+(gauss_gradients[i](0,1)*gauss_gradients[i](1,1)))*volumes(i)*inv_densities(i);
				//and the 'mixed laplacians' (standard shape functions * enrichments)
				for (unsigned int j = 0; j < 4; j++) //we go through the 3 standard shape functions
				{
					mixed_Laplacian(0,j) += (DN_DX(j,0)*(gauss_gradients[i](0,0))+ DN_DX(j,1)*(gauss_gradients[i](0,1))+ DN_DX(j,2)*(gauss_gradients[i](0,2)) )*volumes(i)*conductivities(i);
					//mixed_Laplacian_jump(0,j) -= DN_DX(j,0)*(gauss_gradients[i](1,0)*volumes(i)*inv_densities(i))+ DN_DX(j,1)*(gauss_gradients[i](1,1)*volumes(i)*inv_densities(i));
				}
				
				
			}
			
			//the 'jump' degree of freedom will not be condensed since it's imposed
			//therefore we'll only have to calculate the inverse of the reduced matrix of the gradient jump enrich function. which is a 1x1 matrix (inverse of the 1,1 position of the enriched laplacian)
			const double inv_Laplacian_enrich_weighted=1.0/(Laplacian_enrich(0,0)+enriched_dof_heat_capacity/delta_t);
			//this->SetValue(INV_LAPLACIAN_ENRICH,(inv_Laplacian_enrich_weighted/(delta_t+TauOne)));
			
			double input_heat = 1e3*1800.0*available_element_oxygen*air_partition_density/delta_t;
			//since we're imposing this DOF, we'll do:
			
			for (unsigned int i = 0; i < 4; i++)
			{
				for (unsigned int j = 0; j < ndivisions; j++)
				{
					//standard shape functions addition
					//rRightHandSideVector(i) += Ngauss_in_interfase(i)*input_heat;
					rRightHandSideVector(i) += Ngauss(j,i)*input_heat/double(ndivisions);
					//we replace the single gauss point at the iterface for 1 at each partition:
					//enriched shape function addition
					//rRightHandSideVector(i) -= mixed_Laplacian(0,i)*inv_Laplacian_enrich_weighted*Nenriched_in_interfase(0)*input_heat;
					rRightHandSideVector(i) -= mixed_Laplacian(0,i)*inv_Laplacian_enrich_weighted*Nenriched(j,0)*input_heat/double(ndivisions);
					//we replace the single gauss point at the iterface for 1 at each partition:
					
				}
			}
			
			/*
			//finally we calculate the T* in the previous time step to guarantee that the temperature is the same as in the node
			//Tcomplete=NT+N*T* ---> T* = T* = (Tcomplete - NT) /N* 
			double interfase_temperature_from_negative_side=0.0; //the one that we'd like to have
			double interfase_temperature_complete=0.0;
			double sum_N_only_negative=0.0;
			for (unsigned int i = 0; i < 3; i++)
			{
				interfase_temperature_complete+=old_temperatures[i]*Ngauss_in_interfase(i);
				
				if (distances(i)<0.0)
				{
					interfase_temperature_from_negative_side+=old_temperatures[i]*Ngauss_in_interfase(i);
					sum_N_only_negative+=Ngauss_in_interfase(i);
				
				}
			}
			interfase_temperature_from_negative_side /= sum_N_only_negative;
			
			const double Tstar = ( interfase_temperature_from_negative_side - interfase_temperature_complete ) / Nenriched_in_interfase(0);
			for (unsigned int i = 0; i < 3; i++)
			{
				rRightHandSideVector(i) -= mixed_Laplacian(0,i)*inv_Laplacian_enrich_weighted*Nenriched_in_interfase(0)*Tstar*enriched_dof_heat_capacity/delta_t;
			}
			*/	
			//when there's no jump there's nothing in the RHS to condense. but when we add it, it is equal to -Laplacian_enrich*jump_value due to the elimination of this DOF.
			//and just like the LHS, we must condense it and substrac it (actually it's (-)*(-)=+ )
			//enrich_rhs = -Laplacian_enrich(0,1)*jump_value;
				
			
			
			//for (unsigned int i = 0; i < 3; i++)
			//{
			//	(this->GetValue(ENRICH_LHS_ROW))(i)=mixed_Laplacian(0,i)*(delta_t+TauOne);
			//	//KRATOS_WATCH(mixed_Laplacian(0,i));
			//}	
			
			//now we must condense these matrixes into the original system (D_matrix, Laplacian_matrix)
			Laplacian_matrix -= inv_Laplacian_enrich_weighted * prod(trans(mixed_Laplacian),mixed_Laplacian);
				
		}
		/* 
		else //uncut(normal, not interfase) element
        {
			double element_conductivity=0.0;
			double element_temperature=0.0;
			double element_density=0.0;
			double element_volumetric_heat_capacity=0.0;
			double input_heat=0.0;
			for (unsigned int j = 0; j < 4; j++) //we go through the 4 standard shape functions
			{
				element_temperature += old_temperatures[j]*one_quarter;
			}
			if (geom[0].FastGetSolutionStepValue(DISTANCE)>0.0)
			{
				element_density = this->CalculateAirDensity(element_temperature);
				element_conductivity =  this->CalculateAirConductivity(element_temperature);
				element_volumetric_heat_capacity = rCurrentProcessInfo[SPECIFIC_HEAT_CAPACITY_AIR]*element_density;
				if (element_conductivity<0.0)
					KRATOS_WATCH(element_conductivity);
				
				input_heat = 1e3*1800.0*available_element_oxygen*element_density/delta_t;		
			}
			else
			{
				element_density = rCurrentProcessInfo[DENSITY_WATER];
				element_conductivity =  rCurrentProcessInfo[PERMEABILITY_WATER];
				element_volumetric_heat_capacity = rCurrentProcessInfo[SPECIFIC_HEAT_CAPACITY_WATER]*element_density ; //around 1670 for plastics
			}
			
			Laplacian_matrix =  prod(DN_DX,trans(DN_DX))*Area*element_conductivity;  // B B^T  (standard laplacian, now we must add the contribution of the condensed new nodes.
			
			const double Weight = one_quarter * Area * element_volumetric_heat_capacity;
			
			for (unsigned int i=0; i<LocalSize ; i++)
				Mass_matrix (i,i) = Weight;
				
					//now the loads 	
			for (unsigned int i=0; i<LocalSize ; i++)
				rRightHandSideVector[i] += input_heat*one_quarter;	
				
			
		}		
        */

        //done, now we save it into the LHS
        noalias(rLeftHandSideMatrix) = Laplacian_matrix+ Mass_matrix/(delta_t);///gamma=1
        for (unsigned int i=0; i<LocalSize ; i++)
				rRightHandSideVector[i] += Mass_matrix(i,i)*old_temperatures[i]/(delta_t);
				

        
        //enrich_rhs+=extra_enrich_rhs;
        //to recover the gradient jump later, we must save the RHS of the enrichment, that is:
		//this->SetValue(ENRICH_RHS,enrich_rhs);

        
        ///if it is an inactive element we must turn it off. But to avoid problems in case we've turned off all the elements that contribute to a node, we simply do:
		if (this->GetValue(IS_INACTIVE)==true) //it means we must not add it!
		{
			rLeftHandSideMatrix*=1.0e-6;
			rRightHandSideVector*=1.0e-6;
		}
        
        
		//subtracting the dirichlet term
		// RHS -= LHS*DUMMY_UNKNOWNs		
		for(unsigned int iii = 0; iii<LocalSize; iii++)
			ms_temp[iii] = GetGeometry()[iii].FastGetSolutionStepValue(TEMPERATURE);
		noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix,ms_temp);
		
		KRATOS_CATCH("");
	}
	

	

	//************************************************************************************
	//************************************************************************************
	// this subroutine calculates the nodal contributions for the explicit steps of the 
	// fractional step procedure
	void MonolithicPFEM23D::InitializeSolutionStep(ProcessInfo& CurrentProcessInfo)
	{
		KRATOS_TRY

		KRATOS_CATCH("");
	}
	
	//************************************************************************************
	//************************************************************************************


	void MonolithicPFEM23D::ThermalEquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo)
	{
		unsigned int number_of_nodes = GetGeometry().PointsNumber();
		if(rResult.size() != number_of_nodes)
			rResult.resize(number_of_nodes,false);	

		for (unsigned int i=0;i<number_of_nodes;i++)
				rResult[i] = GetGeometry()[i].GetDof(TEMPERATURE).EquationId();
	}
	
	//************************************************************************************
	void MonolithicPFEM23D::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo)
	{
		const SizeType NumNodes = 4;
		const SizeType LocalSize = 16;
		GeometryType& rGeom = this->GetGeometry();

		SizeType LocalIndex = 0;

		if (rResult.size() != LocalSize)
			rResult.resize(LocalSize, false);

		for (SizeType i = 0; i < NumNodes; ++i)
		{
			rResult[LocalIndex++] = rGeom[i].GetDof(VELOCITY_X).EquationId();
			rResult[LocalIndex++] = rGeom[i].GetDof(VELOCITY_Y).EquationId();
			rResult[LocalIndex++] = rGeom[i].GetDof(VELOCITY_Z).EquationId();
			rResult[LocalIndex++] = rGeom[i].GetDof(PRESSURE).EquationId();
		}
	}

	//************************************************************************************
	//************************************************************************************
	
	void MonolithicPFEM23D::GetThermalDofList(DofsVectorType& ElementalDofList,ProcessInfo& CurrentProcessInfo)
	{
		unsigned int number_of_nodes = GetGeometry().PointsNumber();
		if(ElementalDofList.size() != number_of_nodes)
			ElementalDofList.resize(number_of_nodes);	

		for (unsigned int i=0;i<number_of_nodes;i++)
			ElementalDofList[i] = GetGeometry()[i].pGetDof(TEMPERATURE);

	}
	
	//************************************************************************************
	//************************************************************************************
	
	void MonolithicPFEM23D::GetDofList(DofsVectorType& ElementalDofList,ProcessInfo& CurrentProcessInfo)
	{
		const SizeType NumNodes = 4;
		const SizeType LocalSize = 16;
		GeometryType& rGeom = this->GetGeometry();

		if (ElementalDofList.size() != LocalSize)
			ElementalDofList.resize(LocalSize);

		SizeType LocalIndex = 0;

		for (SizeType i = 0; i < NumNodes; ++i)
		{
			ElementalDofList[LocalIndex++] = rGeom[i].pGetDof(VELOCITY_X);
			ElementalDofList[LocalIndex++] = rGeom[i].pGetDof(VELOCITY_Y);
			ElementalDofList[LocalIndex++] = rGeom[i].pGetDof(VELOCITY_Z);
			ElementalDofList[LocalIndex++] = rGeom[i].pGetDof(PRESSURE);
		}
		
	}
	

	//************************************************************************************
	//************************************************************************************
	
	void MonolithicPFEM23D::CalculatePressureProjection(ProcessInfo& CurrentProcessInfo)
	{
		KRATOS_TRY
		
		double Area;
		Geometry<Node<3> >& geom = this->GetGeometry();
		boost::numeric::ublas::bounded_matrix<double, (3+1), 3 > DN_DX;
		array_1d<double, (3+1) > N;
		GeometryUtils::CalculateGeometryData(geom, DN_DX, N, Area);
		const double mass_factor = 1.0/ (1.0 + double (3) );
		
		const int number_of_particles_in_elem = this->GetValue(NUMBER_OF_PARTICLES);
		
		if( (number_of_particles_in_elem>0))
		{
			
					const double density_air = CurrentProcessInfo[DENSITY_AIR]; //para clindro en rey 100:  0.0005 * 2.0*Area;
					const double density_water = CurrentProcessInfo[DENSITY_WATER];
			
					array_1d<double,4>  pressures = ZeroVector(4); //to calculate the deformation Gradient F. Dimension = velocity dofs
					boost::numeric::ublas::bounded_matrix<double,4, 3 > coords; //coordinates of the nodes
					bool has_negative_node=false;
					bool has_positive_node=false;
					double element_mean_distance=0.0;
					for(unsigned int iii = 0; iii<4; iii++)
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
							
						element_mean_distance+=0.25*this->GetGeometry()[iii].FastGetSolutionStepValue(DISTANCE);
					}
		
					bool split_element=false;
					if (has_positive_node && has_negative_node)
						split_element=true;
			
					if (split_element==false) //we only calculate if we have a pure fluid element
					{

						
						boost::numeric::ublas::bounded_matrix<double, (3+1), (3-1)*6 > G_matrix; //(gradient)
						noalias(G_matrix) = ZeroMatrix((3+1), (3-1)*6);	
						
						//const double water_fraction = -( element_mean_distance - 1.0) *0.5; 
						const double water_fraction = 0.5*(1.0-(sin(3.14159*element_mean_distance*0.5)));
						double density = density_water*(water_fraction)+density_air*(1.0-water_fraction);

						for (unsigned int i = 0; i < (3+1); i++) //loop in shape functions (row)
						{
							for (unsigned int j = 0; j < (3+1) ; j++) //loop in shape functions (column)
							{
								for (unsigned int k = 0; k < (3) ; k++) //x,y,(z)
								{
									G_matrix(i, (j*3)+k ) = DN_DX(i,k)*Area*mass_factor; //mass_factor=(1/3 in 2d, 1/4 in 3d)
								}
							}
						}
						G_matrix /=density;

						//now we save the data:
						for (unsigned int i=0; i!=(3+1); i++) //loop around the nodes of the element to add contribution to node i
						{
							array_1d<double, 3 > & current_press_proj = geom[i].FastGetSolutionStepValue(PRESS_PROJ);
								
							for (unsigned int j = 0; j < (3+1) ; j++) //loop around the nodes of the element to add contribution of node j to node i
							{
								
								for (unsigned int k = 0; k < (3) ; k++) //x,y,(z)
								{
									current_press_proj[k] += G_matrix(j,i*(3)+k)*(pressures(j));///-old_pressures(j)); //gamma=0!
								}
							}
							
							geom[i].FastGetSolutionStepValue(NODAL_AREA) += Area*mass_factor;
							
						}

					} //closing the useful elements
					else// we add some small addition to the area so that we do not have a division by zero.
					{
						for (unsigned int i=0; i!=(3+1); i++) //loop around the nodes of the element to add contribution to node i
						{
							geom[i].FastGetSolutionStepValue(NODAL_AREA) += Area*mass_factor*0.00000001;							
						}
					}
					
		} //closing the if(is_inactive==false)
		else// we add some small addition to the area so that we do not have a division by zero.
		{
			for (unsigned int i=0; i!=(3+1); i++) //loop around the nodes of the element to add contribution to node i
			{
				geom[i].FastGetSolutionStepValue(NODAL_AREA) += Area*mass_factor*0.00000001;							
			}
		}
		
		KRATOS_CATCH("");
	}
	
	//with bounded matrix, constant coefficient
	void MonolithicPFEM23D::AddViscousTerm(boost::numeric::ublas::bounded_matrix<double, 12, 12 >& rDampMatrix,
                         boost::numeric::ublas::bounded_matrix<double, 4, 3 >& rShapeDeriv,
                         const double Weight)
	{

		
		boost::numeric::ublas::bounded_matrix<double, 12,6 > B_matrix = ZeroMatrix(12,6);
		for (unsigned int i=0; i!=(4); i++) //i node
		{
			for (unsigned int j=0; j!=(3); j++) //x,y,z
				B_matrix(i*(3)+j,j)=rShapeDeriv(i,j);
			
			//useful for both 2d and 3d:	
			//relating 12 and 21 stresses
			B_matrix(i*(3)+0,3)=rShapeDeriv(i,1);
			B_matrix(i*(3)+1,3)=rShapeDeriv(i,0);
			

			B_matrix(i*(3)+1,4)=rShapeDeriv(i,2);
			B_matrix(i*(3)+2,4)=rShapeDeriv(i,1);
				
			B_matrix(i*(3)+0,5)=rShapeDeriv(i,2);
			B_matrix(i*(3)+2,5)=rShapeDeriv(i,0);

		}
		
		int counter=0;
		boost::numeric::ublas::bounded_matrix<double, 6, 6 > C_matrix = ZeroMatrix(6,6);
		
		for (unsigned int i=0; i!=(3); i++)
		{
			C_matrix(counter,counter)=2.0;
			counter++;
		}
		for (unsigned int i=0; i!=(3); i++)
		{
			C_matrix(counter,counter)=1.0;
			counter++;
		}
		
		C_matrix*= Weight;
		
		boost::numeric::ublas::bounded_matrix<double, 6 , 12  > temp_matrix = prod(C_matrix,trans(B_matrix));
		rDampMatrix += prod(B_matrix, temp_matrix );
		
	}
	
	void MonolithicPFEM23D::AddViscousTerm(boost::numeric::ublas::bounded_matrix<double, 21, 21 > & output,
										  boost::numeric::ublas::bounded_matrix<double, (4), 3 >& rShapeDeriv,
										  array_1d<double,4>&  distances,
										  std::vector< Matrix >& gauss_gradients, 
										  array_1d<double,6>&  viscosities,
										  array_1d<double,6>&  signs,
										  array_1d<double,6>&  volumes ,
										  const unsigned int ndivisions)
	{		
		boost::numeric::ublas::bounded_matrix<double, 15, 15 > ExtendedDampMatrix=ZeroMatrix(15,15);
		//boost::numeric::ublas::bounded_matrix<double, 8, 8 > rExtendedDampMatrix= ZeroMatrix(8,8);

		boost::numeric::ublas::bounded_matrix<double, 15,6 > B_matrix = ZeroMatrix(15,6);

		
		int counter=0;
		boost::numeric::ublas::bounded_matrix<double, 6, 6 > C_matrix = ZeroMatrix(6,6);
			
		for (unsigned int i=0; i!=(3); i++)
		{
			C_matrix(counter,counter)=2.0;
			counter++;
		}
		for (unsigned int i=0; i!=(3); i++)
		{
			C_matrix(counter,counter)=1.0;
			counter++;
		}
		
		//now the enriched part:
		//we have to construct (ndivisions) rTempDampMatrix and asseble add them to rDampMatrix
		for (unsigned int division=0; division!=ndivisions; division++)
		{
			//standard shape functions:
			for (unsigned int i=0; i!=(4); i++) //i node
			{
				for (unsigned int j=0; j!=(3); j++) //x,y,z
					B_matrix(i*(3)+j,j)=rShapeDeriv(i,j);
				
				//useful for both 2d and 3d:	
				//relating 12 and 21 stresses
				B_matrix(i*(3)+0,3)=rShapeDeriv(i,1);
				B_matrix(i*(3)+1,3)=rShapeDeriv(i,0);
				

				B_matrix(i*(3)+1,4)=rShapeDeriv(i,2);
				B_matrix(i*(3)+2,4)=rShapeDeriv(i,1);
					
				B_matrix(i*(3)+0,5)=rShapeDeriv(i,2);
				B_matrix(i*(3)+2,5)=rShapeDeriv(i,0);
			}
			
			
			for (unsigned int j=0; j!=(3); j++) //x,y,z
				B_matrix(3*(2)+j,j)= gauss_gradients[division](0,j);
			
			//useful for both 2d and 3d:	
			//relating 12 and 21 stresses
			B_matrix(4*(3)+0,3)=gauss_gradients[division](0,1);
			B_matrix(4*(3)+1,3)=gauss_gradients[division](0,0);
			
			B_matrix(4*(3)+1,4)=gauss_gradients[division](0,2);
			B_matrix(4*(3)+2,4)=gauss_gradients[division](0,1);
					
			B_matrix(4*(3)+0,5)=gauss_gradients[division](0,2);
			B_matrix(4*(3)+2,5)=gauss_gradients[division](0,0);

			
			boost::numeric::ublas::bounded_matrix<double, 6 , 15  > temp_matrix = prod(C_matrix,trans(B_matrix));
			ExtendedDampMatrix += viscosities(division)*volumes(division)*prod(B_matrix, temp_matrix );
		}
		
		//now we put it all toghether in the big matrix:
		for (unsigned int i=0; i!=(5); i++) //4 nodes + 1dof in the new virtual node
			for (unsigned int j=0; j!=(5); j++) //4 nodes + 1dof in the new virtual node
				for (unsigned int k=0; k!=(3); k++) //x,y,(z)
					for (unsigned int l=0; l!=(3); l++) //x,y,(z)
						 output(i*4+k,j*4+l) += ExtendedDampMatrix(i*3+k,j*3+l);
						 
		//KRATOS_WATCH(ExtendedDampMatrix)
	}	
	
	
	//with matrixtype, constant coefficient
	void MonolithicPFEM23D::AddViscousTerm(MatrixType& rDampMatrix,
                         const boost::numeric::ublas::bounded_matrix<double, 4, 3 >& rShapeDeriv,
                         const double Weight)
	{

		
		boost::numeric::ublas::bounded_matrix<double, 12,6 > B_matrix = ZeroMatrix(12,6);
		for (unsigned int i=0; i!=(4); i++) //i node
		{
			for (unsigned int j=0; j!=(3); j++) //x,y,z
				B_matrix(i*(3)+j,j)=rShapeDeriv(i,j);
			
			//useful for both 2d and 3d:	
			//relating 12 and 21 stresses
			B_matrix(i*(3)+0,3)=rShapeDeriv(i,1);
			B_matrix(i*(3)+1,3)=rShapeDeriv(i,0);
			

			B_matrix(i*(3)+1,4)=rShapeDeriv(i,2);
			B_matrix(i*(3)+2,4)=rShapeDeriv(i,1);
				
			B_matrix(i*(3)+0,5)=rShapeDeriv(i,2);
			B_matrix(i*(3)+2,5)=rShapeDeriv(i,0);

		}
		
		int counter=0;
		boost::numeric::ublas::bounded_matrix<double, 6, 6 > C_matrix = ZeroMatrix(6,6);
		
		for (unsigned int i=0; i!=(3); i++)
		{
			C_matrix(counter,counter)=2.0;
			counter++;
		}
		for (unsigned int i=0; i!=(3); i++)
		{
			C_matrix(counter,counter)=1.0;
			counter++;
		}
		
		C_matrix*= Weight;
		
		boost::numeric::ublas::bounded_matrix<double, 6 , 12  > temp_matrix = prod(C_matrix,trans(B_matrix));
		boost::numeric::ublas::bounded_matrix<double, 12 , 12  > viscosity_matrix = prod(B_matrix, temp_matrix );
		for (unsigned int i=0; i!=4; i++) //i node
		{
			for (unsigned int j=0; j!=4; j++) //j neighbour
			{
				for (unsigned int k=0; k!=3; k++) //xyz
					for (unsigned int l=0; l!=3; l++) //xyz
						rDampMatrix(i*4+k,j*4+l)+=viscosity_matrix(i*3+k,j*3+l); 
			}
		}		
		//KRATOS_WATCH(C_matrix)
		//KRATOS_WATCH(B_matrix)
		//KRATOS_WATCH(rDampMatrix)
		//		KRATOS_ERROR(std::logic_error, "IMPLICIT STEP FIRST STEP NOT YET IMPLEMENTED IN 3D.. USE LOCALFINALVELOCITY", "");

		
	}
	    
	    
	template<class T>
	bool MonolithicPFEM23D::InvertMatrix(const T& input, T& inverse)
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

	

} // Namespace Kratos
