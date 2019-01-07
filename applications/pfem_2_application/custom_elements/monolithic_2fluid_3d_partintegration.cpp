//   
//   Project Name:        Kratos
//   Last modified by:    $Author: it's me! $
//   Date:                $Date: 2008-08-08 23:58:38 $
//   Revision:            $Revision: 1.0 $
//
//

// Project includes
#include "includes/define.h"
#include "custom_elements/monolithic_2fluid_3d_partintegration.h"
#include "pfem_2_application.h"
#include "utilities/math_utils.h"
#include "utilities/geometry_utilities.h"
#include "utilities/enrichment_utilities.h"

#include "processes/process.h"
#include "includes/node.h"
//#include "includes/element.h"
#include "includes/model_part.h"
#include "includes/cfd_variables.h"


namespace Kratos
{



	//************************************************************************************
	//************************************************************************************
	MonolithicAutoSlipPFEM23D::MonolithicAutoSlipPFEM23D(IndexType NewId, GeometryType::Pointer pGeometry)
		: Element(NewId, pGeometry)
	{
		//DO NOT ADD DOFS HERE!!!
	}

	//************************************************************************************
	//************************************************************************************
	MonolithicAutoSlipPFEM23D::MonolithicAutoSlipPFEM23D(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
		: Element(NewId, pGeometry, pProperties)
	{

	}

	Element::Pointer MonolithicAutoSlipPFEM23D::Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const
	{
		return Element::Pointer(new MonolithicAutoSlipPFEM23D(NewId, GetGeometry().Create(ThisNodes), pProperties));
	}

	MonolithicAutoSlipPFEM23D::~MonolithicAutoSlipPFEM23D()
	{
	}

	//************************************************************************************
	//************************************************************************************


	void MonolithicAutoSlipPFEM23D::AddExplicitContribution(ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY;

		switch ( rCurrentProcessInfo[FRACTIONAL_STEP] )
		{
			case 10: //calculating pressure projection. notthing returned. saving data in PRESS_PROJ, PRESS_PROJ_NO_RO , NODAL_MASS and NODAL_AREA
			{
				this->CalculatePressureProjection(rCurrentProcessInfo);
				break;
			}
			default:
			{
				KRATOS_THROW_ERROR(std::logic_error,"Unexpected value for FRACTIONAL_STEP index: ",rCurrentProcessInfo[FRACTIONAL_STEP]);
			}
		}

		KRATOS_CATCH("");
	}

	//************************************************************************************
	//************************************************************************************

	void MonolithicAutoSlipPFEM23D::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY

		const unsigned int TDim=3;

	    const double delta_t = rCurrentProcessInfo[DELTA_TIME];

	    double volume_correction = rCurrentProcessInfo[VOLUME_CORRECTION];

	    const bool use_press_proj= rCurrentProcessInfo[USE_PRESS_PROJ] != 0;

		array_1d<double,TDim*(TDim+1)>  previous_vel;
		array_1d<double,(TDim+1)*(TDim+1)>  previous_vel_and_press;
		array_1d<double,(TDim+1)*(TDim+1)>  current_vel_and_press;
		array_1d<double,(TDim+1)> distances;

        for(unsigned int iii = 0; iii<(TDim+1); iii++)
        {
			distances[iii] = this->GetGeometry()[iii].FastGetSolutionStepValue(DISTANCE);
			array_1d<double,3>& vel = GetGeometry()[iii].FastGetSolutionStepValue(VELOCITY,1);
			array_1d<double,3>& current_vel = GetGeometry()[iii].FastGetSolutionStepValue(VELOCITY);
			for(unsigned int j = 0; j<TDim; j++)
			{
				previous_vel(iii*TDim+j) = vel[j];
				previous_vel_and_press(iii*(TDim+1)+j) = vel[j];
				current_vel_and_press(iii*(TDim+1)+j) = current_vel[j];

			}
			previous_vel_and_press(iii*(TDim+1)+TDim) = this->GetGeometry()[iii].FastGetSolutionStepValue(PRESSURE);
			current_vel_and_press(iii*(TDim+1)+TDim) = this->GetGeometry()[iii].FastGetSolutionStepValue(PRESSURE);

		}

		const double factor = 1.0/(1.0+double(TDim));
  		array_1d<double,TDim+1> msN; //dimension = number of nodes

		const unsigned int LocalSize = GetGeometry().size()*(TDim+1);
		unsigned int TNumNodes = GetGeometry().size();

		array_1d<double, TDim+1> N;
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

        BoundedMatrix<double, TDim+1, TDim> DN_DX;

        //BoundedMatrix<double, TNumNodes*2, 2> N_vel_matrix; //useless to calculate it, it's simply the DN_DX matrix multiplied by 1/3, arranged in a different shape.
		GeometryUtils::CalculateGeometryData(this->GetGeometry(), DN_DX, N, Area);

		const double viscosity_air = rCurrentProcessInfo[VISCOSITY_AIR];
		const double viscosity_water = rCurrentProcessInfo[VISCOSITY_WATER];
		const double density_air = rCurrentProcessInfo[DENSITY_AIR];
		const double density_water = rCurrentProcessInfo[DENSITY_WATER];




		bool split_element=false;
		for(unsigned int iii = 0; iii<(TDim); iii++)
		{
			if (distances[0]*distances[iii+1]<0.0)
			     split_element=true;
		}


		if (split_element==false)
		//if(true==true)
		{
			BoundedMatrix<double, TDim+1, (TDim+1)*TDim > D_matrix; //(gradient)
			noalias(D_matrix) = ZeroMatrix(TDim+1, (TDim+1)*TDim);
			BoundedMatrix<double, TDim+1, TDim+1 > Laplacian_matrix;
			noalias(Laplacian_matrix) = ZeroMatrix(TDim+1,TDim+1);
			BoundedMatrix<double, (TDim+1)*(TDim+1), (TDim+1)*(TDim+1) > Mass_matrix;
			noalias(Mass_matrix) = ZeroMatrix((TDim+1)*(TDim+1), (TDim+1)*(TDim+1));
			for (unsigned int i = 0; i < TDim+1; i++)
			{
				for (unsigned int j = 0; j < TDim+1; j++)
				{
					for (unsigned int k = 0; k < TDim; k++)
						//D_matrix(i, (j*TDim+k) ) =  DN_DX(j,k)*Area*factor; //(i,k) -> standard
						D_matrix(i, (j*3+k) ) =  -DN_DX(i,k)*Area*factor; //(i,k) -> integrated by parts ( or G)
				}
			}


			//double element_mean_distance = 0.25*(distances(0)+distances(1)+distances(2)+distances(3));
			//const double water_fraction = 0.5*(1.0-(sin(3.14159*element_mean_distance*0.5)));
			//double density = density_water*(water_fraction)+density_air*(1.0-water_fraction);
			double density=density_air;
			//double viscosity=viscosity_air;

			double viscosity_for_tau=0.0;;
			if (distances(0)<0.0)
			{
				density=density_water;
				//viscosity=viscosity_water;
				viscosity_for_tau=viscosity_water;

				if (volume_correction>0.1)
				  volume_correction=0.1;
				if (volume_correction<-0.1)
				  volume_correction=-0.1;

			}
			else
			{
				density = density_air;//this->CalculateAirDensity(element_temperature);
				//viscosity = viscosity_air;//this->CalculateAirViscosity(element_temperature);
				viscosity_for_tau=viscosity_air;
				volume_correction=0.0;
			}
			this->AddViscousTerm(rLeftHandSideMatrix,DN_DX, viscosity_for_tau, Area );



			//CALCULATING TAU (only if we are not near the interfase and we are in step 1 or more)
			double TauOne = 0.0;
			double sqElemSize=0.0;
			{
				//double Viscosity = rCurrentProcessInfo[VISCOSITY];;
				//double AdvVelNorm = sqrt(pow(this->GetGeometry()[0].FastGetSolutionStepValue(VELOCITY_X),2)+pow(this->GetGeometry()[0].FastGetSolutionStepValue(VELOCITY_Y),2)+pow(this->GetGeometry()[0].FastGetSolutionStepValue(VELOCITY_Z),2));
				//const double TimeFactor = rCurrentProcessInfo.GetValue(DYNAMIC_TAU);
				double mElemSize;
				array_1d<double,3> Edge(3,0.0);
				Edge = this->GetGeometry()[1].Coordinates() - this->GetGeometry()[0].Coordinates();
				sqElemSize = Edge[0]*Edge[0];
				for (SizeType d = 1; d < TDim; d++)
					sqElemSize += Edge[d]*Edge[d];
				//
				for (SizeType i = 2; i < TDim+1 ; i++)
					for(SizeType j = 0; j < i; j++)
					{
						Edge = this->GetGeometry()[i].Coordinates() - this->GetGeometry()[j].Coordinates();
						double Length = Edge[0]*Edge[0];
						for (SizeType d = 1; d < TDim ; d++)
							Length += Edge[d]*Edge[d];
						if (Length < sqElemSize) sqElemSize = Length;
					}
				mElemSize=sqrt(sqElemSize);
				//actually we're currently using the kinematic viscosity, so no need to divide the second term by the viscosity, it was already done. Great!
				TauOne =  1.0 / ( 1.0/delta_t + (4.0 * viscosity_for_tau / (mElemSize * mElemSize * density)));// + 2.0 * AdvVelNorm / mElemSize );
			}

			/*
			double fourier=0.0;
			if (distances(0)<0.0) fourier = viscosity_water/density_water *delta_t/pow((0.6*this->GetValue(MEAN_SIZE)),2);
			else fourier = viscosity_air/density_air *delta_t/pow((0.6*this->GetValue(MEAN_SIZE)),2);
			theta = 0.2*fourier - 0.1;
			if (theta<0.0) theta=0.0;
			else if (theta>1.0) theta=1.0;
			*/

			const double Weight = factor * Area * density;

			for (unsigned int i=0; i<TNumNodes ; i++)
				for (unsigned int j=0; j<TDim ; j++)
					Mass_matrix (i*(TDim+1)+j,i*(TDim+1)+j) = Weight;

			Laplacian_matrix =  -prod(DN_DX,trans(DN_DX))*Area/density;
			Laplacian_matrix*=TauOne;

			for (unsigned int i = 0; i < (TDim+1); i++)
			{
				for (unsigned int j = 0; j < (TDim+1); j++)
				{
					for (unsigned int k = 0; k < TDim; k++)
					{
						rLeftHandSideMatrix(i*(TDim+1)+TDim, j*(TDim+1)+k ) -= D_matrix(i,j*TDim+k);//  + DN_DX(i,0)*TauOne*one_third*Area/delta_t;
						rLeftHandSideMatrix(j*(TDim+1)+k, i*(TDim+1)+TDim ) -=  D_matrix(i,j*TDim+k);
					}
					rLeftHandSideMatrix(i*(TDim+1)+TDim, j*(TDim+1)+TDim ) += Laplacian_matrix(i,j);
				}
			}


			//array_1d<double,3> mean_velocity = ZeroVector(3);
			//for (unsigned int i = 0; i < (TDim+1); i++)
			//	mean_velocity += GetGeometry()[i].FastGetSolutionStepValue(VELOCITY)*one_quarter/delta_t;

			//and finally the RHS of the velocity plus the RHS on the pressure due to the stabilzation terms:
			array_1d<double, 3 > rhs_stab = ZeroVector(3);
			for (unsigned int i = 0; i < (TDim+1); i++)
			{
				const array_1d<double,3>& node_press_proj = GetGeometry()[i].FastGetSolutionStepValue(PRESS_PROJ);
				for (unsigned int j = 0; j < TDim; j++)
					rhs_stab[j] += node_press_proj(j)*factor;
			}

			//const bool use_press_proj=true;

			for (unsigned int i = 0; i < (TDim+1); i++)
			{
				const array_1d<double,3>& body_force = GetGeometry()[i].FastGetSolutionStepValue(BODY_FORCE);
				for (unsigned int j = 0; j < TDim; j++)
					rRightHandSideVector(i*(TDim+1)+j) += factor*Area*body_force(j)*density;


				if (use_press_proj)
					for (unsigned int j = 0; j < TDim; j++)
						rRightHandSideVector(i*(TDim+1)+TDim) -= TauOne*Area*(DN_DX(i,j)*rhs_stab(j) + fabs(DN_DX(i,j)) * volume_correction *0.1) ;
				else
					for (unsigned int j = 0; j < TDim; j++)
						rRightHandSideVector(i*(TDim+1)+TDim) -= TauOne*Area*(DN_DX(i,j)*body_force(j) + fabs(DN_DX(i,j)) * volume_correction *0.1) ;

			}


			noalias(rRightHandSideVector) += prod((Mass_matrix/delta_t),previous_vel_and_press);

			noalias(rLeftHandSideMatrix) += Mass_matrix/delta_t;


		}
		else //split element:
		{
			const array_1d<double,3>& body_force_0 = GetGeometry()[0].FastGetSolutionStepValue(BODY_FORCE);


			//TOOLS NEEDED TO FIND THE ENRICHMENT SHAPE FUNCTIONS
			//get position of the cut surface

			array_1d<double,3*(TDim-1)>  viscosities(0);
			array_1d<double,3*(TDim-1)>  densities(0);

			BoundedMatrix<double,3*(TDim-1), 2> Nenriched;
			array_1d<double,3*(TDim-1)>  volumes(0);
			BoundedMatrix<double,TDim+1, TDim > coords;
			BoundedMatrix<double,3*(TDim-1), TDim+1 > Ngauss;
			array_1d<double,3*(TDim-1)>  signs(0);
			std::vector< Matrix > gauss_gradients(3*(TDim-1));

		   	const unsigned int  enrich_velocity_dofs=0;//TDim;
			const unsigned int  enrich_pressure_dofs=1;
			const unsigned int  enrich_pressure_offset=0;

			BoundedMatrix<double, TDim+1, TDim+1> Laplacian_matrix =ZeroMatrix(TDim+1, TDim+1); //standard pressures. we will keep these dof so the size remains
			BoundedMatrix<double, TDim+1, (TDim+1)*TDim+enrich_velocity_dofs > D_matrix ; //(divergence) //this matrix will increase its size since we need more
			noalias(D_matrix) = ZeroMatrix(TDim+1,(TDim+1)*TDim+enrich_velocity_dofs);
			BoundedMatrix<double,TDim+1,(TDim+1)*TDim+enrich_velocity_dofs > G_matrix; //(gradient)
			noalias(G_matrix) = ZeroMatrix(TDim+1,(TDim+1)*TDim+enrich_velocity_dofs);

			BoundedMatrix<double, (TDim+1)*(TDim+1)+enrich_velocity_dofs+enrich_pressure_dofs, (TDim+1)*(TDim+1)+enrich_velocity_dofs+enrich_pressure_dofs > Momentum_matrix; //2 vel + 1 pressure per node   plus 2 enriched pressures and 2 enriched velocities
			noalias(Momentum_matrix) = ZeroMatrix((TDim+1)*(TDim+1)+enrich_velocity_dofs+enrich_pressure_dofs,(TDim+1)*(TDim+1)+enrich_velocity_dofs+enrich_pressure_dofs);

			//fill coordinates
			for (unsigned int i = 0; i < TNumNodes; i++)
			{
				const array_1d<double, 3 > & xyz = this->GetGeometry()[i].Coordinates();
				//volumes[i] = 0.0;
				//KRATOS_WATCH(distances(i));
				for (unsigned int j = 0; j < TDim; j++)
					coords(i, j) = xyz[j];
			}
			for (unsigned int i = 0; i < 3*(TDim-1); i++)
				gauss_gradients[i].resize(2, TDim, false);  //2 values of the 2 shape functions, and derivates in (xyz) direction).

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
					Weight_for_tau += volumes(i)*viscosity;
					densities(i)=density_water;
					viscosities(i) = viscosity_water;//CalculateAirDensity(partition_temperature);
					mass += volumes(i)*density_water;
				}
			}

			double element_viscosity_for_tau=Weight_for_tau/Area;
			const double element_density= mass/Area;

			for (unsigned int i = 0; i < TDim+1; i++)
			{
				for (unsigned int j = 0; j < TDim+1; j++)
				{
					for (unsigned int k = 0; k < TDim; k++)
						//D_matrix(i, (j*TDim+k) ) =  DN_DX(j,k)*Area*factor; //(i,k) -> standard
						D_matrix(i, (j*TDim+k) ) =  -DN_DX(i,k)*Area*factor; //(i,k) -> integrated by parts ( or G)
				}
			}

			/*
			this->AddViscousTerm(Momentum_matrix,
								DN_DX,
								distances,
								gauss_gradients,
								viscosities,
								signs,
								volumes,
								ndivisions);
			*/
			this->AddViscousTerm(rLeftHandSideMatrix,DN_DX, element_viscosity_for_tau,Area );


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
				//double AdvVelNorm = sqrt(pow(this->GetGeometry()[0].FastGetSolutionStepValue(VELOCITY_X),2)+pow(this->GetGeometry()[0].FastGetSolutionStepValue(VELOCITY_Y),2)+pow(this->GetGeometry()[0].FastGetSolutionStepValue(VELOCITY_Z),2));
				//const double TimeFactor = rCurrentProcessInfo.GetValue(DYNAMIC_TAU);
				double mElemSize;
				array_1d<double,3> Edge(3,0.0);
				Edge = this->GetGeometry()[1].Coordinates() - this->GetGeometry()[0].Coordinates();
				sqElemSize = Edge[0]*Edge[0];
				for (SizeType d = 1; d < TDim; d++)
					sqElemSize += Edge[d]*Edge[d];
				//
				for (SizeType i = 2; i < (TDim+1); i++)
					for(SizeType j = 0; j < i; j++)
					{
						Edge = this->GetGeometry()[i].Coordinates() - this->GetGeometry()[j].Coordinates();
						double Length = Edge[0]*Edge[0];
						for (SizeType d = 1; d < TDim; d++)
							Length += Edge[d]*Edge[d];
						if (Length < sqElemSize) sqElemSize = Length;
					}
				mElemSize=sqrt(sqElemSize);
				//actually we're currently using the kinematic viscosity, so no need to divide the second term by the viscosity, it was already done. Great!
				TauOne =  0.01 / ( 1.0 / delta_t +  (4.0 * element_viscosity_for_tau / (mElemSize * mElemSize * element_density)));//+ 2.0 * AdvVelNorm / mElemSize );
			}



			array_1d<double,TDim+1>  lumped_mass = ZeroVector(TDim+1);
			double condensed_dof_mass1=0.0;
			double condensed_dof_mass2=0.0;
			double condensed_dof_mass3=0.0;
			double Laplacian_matrix_factor=0.0;
			for (unsigned int i = 0; i < ndivisions; i++) //partitions
			{
				//KRATOS_WATCH(enrich_lhs(i));
				for (unsigned int j = 0; j < TDim+1 ; j++) // local node (or shape function)
				{

					lumped_mass(j) +=  volumes(i)*Ngauss(i,j)*densities(i);
					Momentum_matrix((TDim+1)*j,(TDim+1)*j) += volumes(i)*Ngauss(i,j)*densities(i) / delta_t;
					Momentum_matrix((TDim+1)*j+1,(TDim+1)*j+1) += volumes(i)*Ngauss(i,j)*densities(i) / delta_t;
					Momentum_matrix((TDim+1)*j+2,(TDim+1)*j+2) += volumes(i)*Ngauss(i,j)*densities(i) / delta_t;
				}


				if (enrich_velocity_dofs!=0)
				{
					condensed_dof_mass1 += volumes(i)*Nenriched(i,0)*densities(i);
					condensed_dof_mass2 += volumes(i)*Nenriched(i,0)*densities(i);
					condensed_dof_mass3 += volumes(i)*Nenriched(i,0)*densities(i);

					Momentum_matrix((TDim+1)*(TDim+1)+0,(TDim+1)*(TDim+1)+0) += volumes(i)*Nenriched(i,0)*densities(i) / delta_t;	 //degrees of freedom 10 and 11 are the enrichment velocities. degrees of freedom 12 and 13 are the enrichment pressures
					Momentum_matrix((TDim+1)*(TDim+1)+1,(TDim+1)*(TDim+1)+1) += volumes(i)*Nenriched(i,0)*densities(i) / delta_t;
					Momentum_matrix((TDim+1)*(TDim+1)+2,(TDim+1)*(TDim+1)+2) += volumes(i)*Nenriched(i,0)*densities(i) / delta_t;
				}

				Laplacian_matrix_factor +=volumes(i)/densities(i);

			}



			Laplacian_matrix =  -prod(DN_DX,trans(DN_DX))*Laplacian_matrix_factor;
			Laplacian_matrix*=TauOne;



			BoundedMatrix<double, enrich_pressure_dofs, enrich_pressure_dofs > Laplacian_enrich;
			noalias(Laplacian_enrich) = ZeroMatrix(enrich_pressure_dofs,enrich_pressure_dofs);
			BoundedMatrix<double, enrich_pressure_dofs, (TDim+1) > mixed_Laplacian;
			noalias(mixed_Laplacian) = ZeroMatrix(enrich_pressure_dofs,(TDim+1));
			//To calculate D* =int (N*DN_DX_enrich).  we must use the N at the gauss point of the partition (returned by the enrichment util) and multiply by the area.
			//notice that this first row is only the upper part of the mixed D. we also need the lower one for the jump
			BoundedMatrix<double, enrich_pressure_dofs, (TDim+1)*TDim + enrich_velocity_dofs> D_matrix_mixed;
			noalias(D_matrix_mixed) = ZeroMatrix(enrich_pressure_dofs,(TDim+1)*TDim + enrich_velocity_dofs);
			BoundedMatrix<double, enrich_pressure_dofs, (TDim+1)*TDim + enrich_velocity_dofs > G_matrix_mixed;
			noalias(G_matrix_mixed) = ZeroMatrix(enrich_pressure_dofs,(TDim+1)*TDim + enrich_velocity_dofs);

			//the matrices we need at the end will be:
			//BoundedMatrix<double, 2, 2 > D_mod;

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
					    for (unsigned int l = 0; l < TDim; l++)
							Laplacian_enrich(j,k) -= (gauss_gradients[i](j+enrich_pressure_offset,l)*gauss_gradients[i](k+enrich_pressure_offset,l))*volumes(i)/densities(i);

				//and the 'mixed laplacians' (standard shape functions * enrichments)
				//first we take care of the standard shape functions of the velocity
				for (unsigned int j = 0; j < TDim+1; j++) //we go through the 4 standard shape functions
				{
					for (unsigned int k = 0; k < enrich_pressure_dofs; k++) //we go through the 3 enrichments
					{
						for (unsigned int l = 0; l < TDim; l++)
						{
							mixed_Laplacian(k,j) -= (DN_DX(j,l)*(gauss_gradients[i](k+enrich_pressure_offset,l)))*volumes(i)/densities(i);
							D_matrix_mixed(k,j*TDim+l) += DN_DX(j,l) *volumes(i)*Nenriched(i,k+enrich_pressure_offset);
							G_matrix_mixed(k,j*TDim+l) -= gauss_gradients[i](k+enrich_pressure_offset,l) *volumes(i)*Ngauss(i,j);
						}
					}
				}


				//we go through the remaining, new velocity dofs. (one x and one y component
				if (enrich_velocity_dofs!=0)
				{
					for (unsigned int k = 0; k < enrich_pressure_dofs; k++) //we go through the 3 enrichments
					{
						for (unsigned int l = 0; l < TDim; l++)
						{
							D_matrix_mixed(k,(TDim+1)*TDim+l) += gauss_gradients[i](0,l) *volumes(i)*Nenriched(i,k+enrich_pressure_offset);
							G_matrix_mixed(k,(TDim+1)*TDim+l) -= gauss_gradients[i](k+enrich_pressure_offset,l) *volumes(i)*Nenriched(i,0);
						}
					}


					//also the standard D and G matrices had to be expanded.

						for (unsigned int k = 0; k < TDim; k++) //we go through the 4 standard pressures
						{
							for (unsigned int l = 0; l < TDim; l++)
							{
								D_matrix(k,(TDim+1)*TDim+l) += gauss_gradients[i](0,l) *volumes(i)*Ngauss(i,k);
								G_matrix(k,(TDim+1)*TDim+l) -= DN_DX(k,l) *volumes(i)*Nenriched(i,0);
							}
						}
				}


				for (unsigned int k = 0; k < enrich_pressure_dofs; k++) //we go through the k enrichments
					for (unsigned int l = 0; l < TDim; l++)
						rhs_enrich(k+enrich_velocity_dofs) -= (gauss_gradients[i](k+enrich_pressure_offset,l)*body_force_0(l))*volumes(i);

			}

			Laplacian_enrich*=TauOne;
			mixed_Laplacian*=TauOne;
			rhs_enrich *=TauOne;


			//adding the term int(q* V . n)
			//array_1d<double,3>  boundary_integral=ZeroVector(3);
			//ComputeBoundaryIntegral(coords, distances, boundary_integral);
			//array_1d<double,3> mean_vel = 0.25 * (GetGeometry()[0].FastGetSolutionStepValue(VELOCITY)+GetGeometry()[1].FastGetSolutionStepValue(VELOCITY)+GetGeometry()[2].FastGetSolutionStepValue(VELOCITY)+GetGeometry()[3].FastGetSolutionStepValue(VELOCITY));
			//const double boundary_term = boundary_integral[0]*mean_vel[0]+boundary_integral[1]*mean_vel[1]+boundary_integral[2]*mean_vel[2];
			//rhs_enrich[enrich_velocity_dofs] += 1.0*boundary_term;

			if (enrich_velocity_dofs!=0)
			{
				rhs_enrich(0) = body_force_0(0)*condensed_dof_mass1;
				rhs_enrich(1) = body_force_0(1)*condensed_dof_mass2;
				rhs_enrich(2) = body_force_0(2)*condensed_dof_mass3;
			}

			BoundedMatrix<double, enrich_velocity_dofs+enrich_pressure_dofs, (TDim+1)*(TDim+1) > condensed_rows; //Vx1,Vy1,p1,Vx2,...
			BoundedMatrix<double, enrich_velocity_dofs+enrich_pressure_dofs, (TDim+1)*(TDim+1) > condensed_columns; //Vx1,Vy1,p1,Vx2,...
			BoundedMatrix<double, enrich_velocity_dofs+enrich_pressure_dofs, enrich_velocity_dofs+enrich_pressure_dofs > condensed_block; //Vx1,Vy1,p1,Vx2,...

			for (unsigned int i = 0; i <TDim+1; i++)  //we go through the 4 nodes (standard velocity dof + standard pressure dof)
			{
				//enriched pressure dof
				for (unsigned int k = 0; k < enrich_pressure_dofs; k++)
				{

					//enriched
					for (unsigned int l = 0; l < TDim; l++)
						condensed_rows(enrich_velocity_dofs+k,i*(TDim+1)+l)= - D_matrix_mixed(k,i*TDim+l);//+mass_stabilization_terms(i*2+1);
					condensed_rows(enrich_velocity_dofs+k,i*(TDim+1)+TDim)= mixed_Laplacian(k,i);

					for (unsigned int l = 0; l < TDim; l++)
						condensed_columns(enrich_velocity_dofs+k,i*(TDim+1)+l)= - D_matrix_mixed(k,i*TDim+l);
					condensed_columns(enrich_velocity_dofs+k,i*(TDim+1)+TDim)= mixed_Laplacian(k,i);
				}



				for (unsigned int k = 0; k < enrich_velocity_dofs; k++) //we go through the 3 enrichments
				{
					for (unsigned int l = 0; l < TDim; l++)
					{
						condensed_columns(k,i*(TDim+1)+l)=Momentum_matrix((TDim+1)*(TDim+1)+k,i*(TDim+1)+l); //add the viscosity matrix
						condensed_rows(k,i*(TDim+1)+l)=Momentum_matrix((TDim+1)*(TDim+1)+k,i*(TDim+1)+l); //add the viscosity matrix
					}

					///WARNING, WHEN THE MATRIX IS REARRANGED, the condensed rows have the gradient of the pressure and the columns have the divergence. that is the reason for the mixed indexes.
					///if the gradient was not integrated by parts, then G matrix should be used in the columns instead of the rows, unlike the case of standard velocity DOFs
					condensed_rows(k,i*(TDim+1)+TDim)= -G_matrix(i, (TDim+1)*TDim+k);
					condensed_columns(k,i*(TDim+1)+TDim)= -G_matrix(i, (TDim+1)*TDim+k);
				}

			}
			//KRATOS_WATCH(Momentum_matrix)


			//now the condensed block terms:
			//the condensed block has 4 submatrices:    1 [ K+M*] [ G*+]  2
			//											3 [ D *+] [ L* ]  4

			//first block
			for (unsigned int i = 0; i < enrich_velocity_dofs; i++) //we go through the 3 enrichments
				for (unsigned int k = 0; k < enrich_velocity_dofs; k++) //we go through the 3 enrichments
					condensed_block(i,k)=Momentum_matrix((TDim+1)*(TDim+1)+i,(TDim+1)*(TDim+1)+k);
			//second block
			for (unsigned int i = 0; i < enrich_velocity_dofs; i++) //we go through the 3 enrichments
				for (unsigned int k = 0; k < enrich_pressure_dofs; k++) //we go through the 3 enrichments
					condensed_block(i,k+enrich_velocity_dofs)= -G_matrix_mixed(k,(TDim+1)*TDim+i);		// in  this case, we are in the gradient side and we should use the gradient if we do not integrate it by parts.
			//third block
			for (unsigned int i = 0; i < enrich_pressure_dofs; i++) //we go through the 3 enrichments
				for (unsigned int k = 0; k < enrich_velocity_dofs; k++) //we go through the 3 enrichments
					condensed_block(i+enrich_velocity_dofs,k)= -G_matrix_mixed(i,(TDim+1)*TDim+k);		// in  this case, we are in the divergence side and we should use the gradient if we want to integrate the divergence it by parts.

			//fourth block
			for (unsigned int i = 0; i < enrich_pressure_dofs; i++) //we go through the 3 enrichments
				for (unsigned int k = 0; k < enrich_pressure_dofs; k++) //we go through the 3 enrichments
					condensed_block(i+enrich_velocity_dofs,k+enrich_velocity_dofs)=Laplacian_enrich(i,k);		//

			BoundedMatrix<double, enrich_velocity_dofs+enrich_pressure_dofs , enrich_velocity_dofs+enrich_pressure_dofs  > inverse_enrichments;
			this->InvertMatrix(condensed_block,inverse_enrichments);
			//condensing
			BoundedMatrix<double, (TDim+1)*(TDim+1) , enrich_pressure_dofs+enrich_velocity_dofs  > temp_matrix;
			temp_matrix = prod(trans(condensed_columns),inverse_enrichments);
			rLeftHandSideMatrix -=  prod(temp_matrix,condensed_rows);
			noalias(rRightHandSideVector) -= prod(temp_matrix,rhs_enrich);



			for (unsigned int i = 0; i < TDim+1; i++)
			{
				for (unsigned int j = 0; j < TDim+1; j++)
				{
					for (unsigned int k = 0; k < TDim; k++)
					{
						rLeftHandSideMatrix(i*(TDim+1)+TDim, j*(TDim+1)+k ) -= D_matrix(i,j*TDim+k);//  + DN_DX(i,0)*TauOne*one_third*Area/delta_t;
						rLeftHandSideMatrix(j*(TDim+1)+k, i*(TDim+1)+TDim ) -=  D_matrix(i,j*TDim+k);
					}
					rLeftHandSideMatrix(i*(TDim+1)+TDim, j*(TDim+1)+TDim ) += Laplacian_matrix(i,j);
				}
			}


			//array_1d<double,3> mean_velocity = ZeroVector(3);
			//for (unsigned int i = 0; i < 4; i++)
			//	mean_velocity += GetGeometry()[i].FastGetSolutionStepValue(VELOCITY)*factor;

			//and finally the RHS of the velocity plus the RHS on the pressure due to the stabilzation terms:
			for (unsigned int i = 0; i < TDim+1; i++)
			{
				const array_1d<double,3>& body_force = GetGeometry()[i].FastGetSolutionStepValue(BODY_FORCE);
				for (unsigned int l = 0; l < TDim; l++)
				{
					rRightHandSideVector(i*(TDim+1)+l) += lumped_mass(i)*body_force(l);
					rRightHandSideVector(i*(TDim+1)+l) += lumped_mass(i)*previous_vel_and_press(i*(TDim+1)+l)/delta_t;
					rRightHandSideVector(i*(TDim+1)+TDim) -= TauOne*Area*(DN_DX(i,l)*body_force(l));
				}


			}

			for (unsigned int i = 0; i < (TDim+1)*(TDim+1); i++)
				for (unsigned int j = 0; j < (TDim+1)*(TDim+1); j++)
					rLeftHandSideMatrix(i,j) += Momentum_matrix(i,j);

		}




		//we must actually use  a fraction of the rigidity matrix. since the Lefthandside only has it so far, we multiply it by 0.5:
		//rLeftHandSideMatrix *= theta;

        //finally we add to the righthandside BOTH the mass matrix and the rigidity matrix:


		noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix,current_vel_and_press);

		//KRATOS_WATCH("ART")
		//KRATOS_WATCH(rRightHandSideVector)
		//KRATOS_WATCH(rLeftHandSideMatrix)

		KRATOS_CATCH("");
	}



	//************************************************************************************
	//************************************************************************************
	// this subroutine calculates the nodal contributions for the explicit steps of the
	// fractional step procedure
	void MonolithicAutoSlipPFEM23D::InitializeSolutionStep(ProcessInfo& CurrentProcessInfo)
	{
		KRATOS_TRY

		KRATOS_CATCH("");
	}

	//************************************************************************************
	//************************************************************************************

	void MonolithicAutoSlipPFEM23D::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo)
	{
		const int TDim = 3;

		const SizeType NumNodes = TDim+1;
		const SizeType LocalSize = (TDim+1)*(TDim+1);
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

	void MonolithicAutoSlipPFEM23D::GetDofList(DofsVectorType& ElementalDofList,ProcessInfo& CurrentProcessInfo)
	{
		const int TDim = 3;

		const SizeType NumNodes = TDim+1;
		const SizeType LocalSize = (TDim+1)*(TDim+1);
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

	void MonolithicAutoSlipPFEM23D::CalculatePressureProjection(ProcessInfo& CurrentProcessInfo)
	{
			KRATOS_TRY

		const unsigned int TDim=3;

		double Area;
		Geometry<Node<3> >& geom = this->GetGeometry();
		BoundedMatrix<double, (TDim+1), TDim > DN_DX;
		array_1d<double, (TDim+1) > N;
		GeometryUtils::CalculateGeometryData(geom, DN_DX, N, Area);
		const double factor = 1.0/ (1.0 + double (TDim) );

		//const int number_of_particles_in_elem = this->GetValue(NUMBER_OF_FLUID_PARTICLES);

		//if( (number_of_particles_in_elem>0))
		if (true)
		{

					const double density_air = CurrentProcessInfo[DENSITY_AIR]; //para clindro en rey 100:  0.0005 * 2.0*Area;
					const double density_water = CurrentProcessInfo[DENSITY_WATER];

					array_1d<double,TDim+1>  pressures = ZeroVector(TDim+1); //to calculate the deformation Gradient F. Dimension = velocity dofs
					BoundedMatrix<double,TDim+1, TDim > coords; //coordinates of the nodes
					bool has_negative_node=false;
					bool has_positive_node=false;
					double element_mean_distance=0.0;
					for(unsigned int iii = 0; iii<(TDim+1); iii++)
					{
						//saving everything
						pressures(iii) = GetGeometry()[iii].FastGetSolutionStepValue(PRESSURE);

						//
						const array_1d<double, 3 > & xyz = this->GetGeometry()[iii].Coordinates();
						for (unsigned int j = 0; j < TDim; j++)
							coords(iii, j) = xyz[j];
						//
						if (this->GetGeometry()[iii].FastGetSolutionStepValue(DISTANCE)<0.0)
							has_negative_node=true;
						else
							has_positive_node=true;

						element_mean_distance+=factor*this->GetGeometry()[iii].FastGetSolutionStepValue(DISTANCE);
					}

					bool split_element=false;
					if (has_positive_node && has_negative_node)
						split_element=true;

					if (split_element==false) //we only calculate if we have a pure fluid element
					{


						BoundedMatrix<double, (TDim+1), (TDim+1)*TDim > G_matrix; //(gradient)
						noalias(G_matrix) = ZeroMatrix((TDim+1), (TDim+1)*TDim);

						//const double water_fraction = -( element_mean_distance - 1.0) *0.5;
						//const double water_fraction = 0.5*(1.0-(sin(3.14159*element_mean_distance*0.5)));
						//double density = density_water*(water_fraction)+density_air*(1.0-water_fraction);
						double density = density_air;
						if (has_negative_node==true)
							density= density_water;

						for (unsigned int i = 0; i < (TDim+1); i++) //loop in shape functions (row)
						{
							for (unsigned int j = 0; j < (TDim+1) ; j++) //loop in shape functions (column)
							{
								for (unsigned int k = 0; k < (TDim) ; k++) //x,y,(z)
								{
									G_matrix(i, (j*3)+k ) = DN_DX(i,k)*Area*factor; //mass_factor=(1/3 in 2d, 1/4 in 3d)
								}
							}
						}
						G_matrix /=density;

						//now we save the data:
						for (unsigned int i=0; i!=(TDim+1); i++) //loop around the nodes of the element to add contribution to node i
						{
							geom[i].SetLock();
							array_1d<double, 3 > & current_press_proj = geom[i].FastGetSolutionStepValue(PRESS_PROJ);

							for (unsigned int j = 0; j < (TDim+1) ; j++) //loop around the nodes of the element to add contribution of node j to node i
							{

								for (unsigned int k = 0; k < (TDim) ; k++) //x,y,(z)
								{
									current_press_proj[k] += G_matrix(j,i*(TDim)+k)*(pressures(j));///-old_pressures(j)); //gamma=0!
								}
							}

							geom[i].FastGetSolutionStepValue(NODAL_AREA) += Area*factor;

							geom[i].UnSetLock();

						}

					} //closing the useful elements
					else// we add some small addition to the area so that we do not have a division by zero.
					{
						for (unsigned int i=0; i!=(TDim+1); i++) //loop around the nodes of the element to add contribution to node i
						{
							geom[i].SetLock();
							geom[i].FastGetSolutionStepValue(NODAL_AREA) += Area*factor*0.00000001;
							geom[i].UnSetLock();
						}
					}

		} //closing the if(is_inactive==false)
		else// we add some small addition to the area so that we do not have a division by zero.
		{
			for (unsigned int i=0; i!=(TDim+1); i++) //loop around the nodes of the element to add contribution to node i
			{
				geom[i].SetLock();
				geom[i].FastGetSolutionStepValue(NODAL_AREA) += Area*factor*0.00000001;
				geom[i].UnSetLock();
			}
		}

		KRATOS_CATCH("");
	}


	void MonolithicAutoSlipPFEM23D::AddViscousTerm(BoundedMatrix<double, 21, 21 > & output,
										  BoundedMatrix<double, (4), 3 >& rShapeDeriv,
										  array_1d<double,4>&  distances,
										  std::vector< Matrix >& gauss_gradients,
										  array_1d<double,6>&  viscosities,
										  array_1d<double,6>&  signs,
										  array_1d<double,6>&  volumes ,
										  const unsigned int ndivisions)
	{
		BoundedMatrix<double, 15, 15 > ExtendedDampMatrix=ZeroMatrix(15,15);
		//BoundedMatrix<double, 8, 8 > rExtendedDampMatrix= ZeroMatrix(8,8);

		BoundedMatrix<double, 15,6 > B_matrix = ZeroMatrix(15,6);


		int counter=0;
		BoundedMatrix<double, 6, 6 > C_matrix = ZeroMatrix(6,6);

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


			BoundedMatrix<double, 6 , 15  > temp_matrix = prod(C_matrix,trans(B_matrix));
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
	void MonolithicAutoSlipPFEM23D::AddViscousTerm(MatrixType& rDampMatrix,
                         const BoundedMatrix<double, 4, 3 >& rShapeDeriv,
                         double Viscosity,const double Area)
	{


		BoundedMatrix<double, 12,6 > B_matrix = ZeroMatrix(12,6);
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
		BoundedMatrix<double, 6, 6 > C_matrix = ZeroMatrix(6,6);

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

		C_matrix*= Viscosity*Area;

		BoundedMatrix<double, 6 , 12  > temp_matrix = prod(C_matrix,trans(B_matrix));
		BoundedMatrix<double, 12 , 12  > viscosity_matrix = prod(B_matrix, temp_matrix );
		for (unsigned int i=0; i!=4; i++) //i node
		{
			for (unsigned int j=0; j!=4; j++) //j neighbour
			{
				for (unsigned int k=0; k!=3; k++) //xyz
					for (unsigned int l=0; l!=3; l++) //xyz
						rDampMatrix(i*4+k,j*4+l)+=viscosity_matrix(i*3+k,j*3+l);
			}
		}
	}

	void MonolithicAutoSlipPFEM23D::ComputeBoundaryIntegral(const BoundedMatrix<double, 4, 3 >& coords,
														const array_1d<double,4>&  distances,
														array_1d<double,3>&  output_integral)
	{
		const int node_i[] = {1, 0, 0, 0};
        const int node_j[] = {2, 3, 1, 2};
        const int node_k[] = {3, 2, 3, 1};

        output_integral = ZeroVector(3);

        std::vector< array_1d<double, 3 > > coords_vector(4);
        for (unsigned int i=0; i!=4; i++) //i node
			for (unsigned int j=0; j!=3; j++) //xyz
				coords_vector[i](j) = coords(i,j);

		 for (unsigned int i=0; i!=4; i++) //i face
		 {


			bool is_split=false;
			if (distances[node_i[i]]*distances[node_j[i]]<0.0 || distances[node_i[i]]*distances[node_k[i]]<0.0)
				   is_split=true;

			array_1d<double,3> normal = ZeroVector(3);
			TriangleNormal(coords_vector[node_i[i]],coords_vector[node_j[i]],coords_vector[node_k[i]], normal);

			if (is_split==true)
			{
				output_integral+=0.3*normal; //approximated
			}
		 }


	}


    void MonolithicAutoSlipPFEM23D::TriangleNormal(const array_1d<double,3>&  node_i, const array_1d<double,3>&  node_j, const array_1d<double,3>&  node_k, array_1d<double,3>&  normal)
    {
		array_1d<double,3> v1,v2;
		v1[0] = node_j[0] - node_i[0];
		v1[1] = node_j[1] - node_i[1];
		v1[2] = node_j[2] - node_i[2];

		v2[0] = node_k[0] - node_i[0];;
		v2[1] = node_k[1] - node_i[1];
		v2[2] = node_k[2] - node_i[2];
		MathUtils<double>::CrossProduct(normal,v1,v2);
		normal *= 0.5;
    }


	template<class T>
	bool MonolithicAutoSlipPFEM23D::InvertMatrix(const T& input, T& inverse)
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
