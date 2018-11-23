//   
//   Project Name:        Kratos
//   Last modified by:    $Author: it's me! $
//   Date:                $Date: 2008-08-08 23:58:38 $
//   Revision:            $Revision: 1.0 $
//
//

// Project includes
#include "includes/define.h"
#include "custom_elements/monolithic_3fluid_3d.h"
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
	Monolithic3FluidPFEM23D::Monolithic3FluidPFEM23D(IndexType NewId, GeometryType::Pointer pGeometry)
		: Element(NewId, pGeometry)
	{
		//DO NOT ADD DOFS HERE!!!
	}

	//************************************************************************************
	//************************************************************************************
	Monolithic3FluidPFEM23D::Monolithic3FluidPFEM23D(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
		: Element(NewId, pGeometry, pProperties)
	{

	}

	Element::Pointer Monolithic3FluidPFEM23D::Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const
	{
		return Element::Pointer(new Monolithic3FluidPFEM23D(NewId, GetGeometry().Create(ThisNodes), pProperties));
	}

	Monolithic3FluidPFEM23D::~Monolithic3FluidPFEM23D()
	{
	}

	//************************************************************************************
	//************************************************************************************


	void Monolithic3FluidPFEM23D::AddExplicitContribution(ProcessInfo& rCurrentProcessInfo)
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
				KRATOS_THROW_ERROR(std::logic_error,"Unexpected value for FRACTIONAL_STEP index: ",rCurrentProcessInfo[FRACTIONAL_STEP]);
			}
		}

		KRATOS_CATCH("");
	}










	//************************************************************************************
	//************************************************************************************

	void Monolithic3FluidPFEM23D::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY

		const unsigned int TDim=3;
		const unsigned int TNumNodes = TDim+1;
		const unsigned int LocalSize = ( TNumNodes ) *(TDim+1);


	    double delta_t = rCurrentProcessInfo[DELTA_TIME];
	    array_1d<double, 3 > gravity = rCurrentProcessInfo[GRAVITY];
		//KRATOS_WATCH(gradient_discontinuity);



		double theta=1.0;


		array_1d<double,TDim*TNumNodes >  previous_vel;
		array_1d<double,LocalSize >  previous_vel_and_press;
		//array_1d<double, 4 > temperatures;
		array_1d<double,TNumNodes> distances;
     	array_1d<double,TNumNodes> corrected_distances;

     	bool has_positive_node_distance = false;
     	bool has_negative_node_distance = false;               //interface between  standard "distance" implies that there are volumes of "air" and "water"
     	bool has_positive_node_corrected_distance = false;     //interface between   "corrected distance" implies that there are volumes of "void" material and the rest is either "air or water". we can't have 2 interfaces
     	bool has_negative_node_corrected_distance = false;


        for(unsigned int iii = 0; iii<TNumNodes; iii++)
        {
			//temperatures(iii) = GetGeometry()[iii].FastGetSolutionStepValue(TEMPERATURE);
			distances[iii] = this->GetGeometry()[iii].FastGetSolutionStepValue(DISTANCE);
			corrected_distances[iii] = this->GetGeometry()[iii].FastGetSolutionStepValue(CORRECTED_DISTANCE);
			array_1d<double,3>& vel = GetGeometry()[iii].FastGetSolutionStepValue(VELOCITY);
			for(unsigned int j = 0; j<TDim; j++)
			{
				previous_vel(iii*TDim+j) = vel[j];
				previous_vel_and_press(iii*(TDim+1)+j) = vel[j];
			}
			previous_vel_and_press(iii*(TNumNodes)+TDim) = this->GetGeometry()[iii].FastGetSolutionStepValue(PRESSURE);

			if(distances[iii]<0.0) //negative nodes are not actually really dsf
				has_negative_node_distance = true;
			else if (corrected_distances[iii]<0.0) //if they are positive on corrected_distance, they are actually "void" nodes, not "air" nodes.
				has_positive_node_distance = true;

			if(corrected_distances[iii]<0.0)
				has_negative_node_corrected_distance = true;
			else
				has_positive_node_corrected_distance = true;


		}


		bool split_element_distance = false;
		if (has_negative_node_distance && has_positive_node_distance)
			split_element_distance = true;

		bool split_element_corrected_distance = false;
		if (has_negative_node_corrected_distance && has_positive_node_corrected_distance)
			split_element_corrected_distance = true;


		const double one_quarter = 0.25;
  		array_1d<double,TNumNodes> msN; //dimension = number of nodes

		array_1d<double, TNumNodes > N;
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

        BoundedMatrix<double, TNumNodes, TDim > DN_DX;

		GeometryUtils::CalculateGeometryData(this->GetGeometry(), DN_DX, N, Area);

        //we start all the phases as the light fluid ones
		double viscosity_positive = rCurrentProcessInfo[VISCOSITY_AIR]; //we start all the phases as the light fluid ones
		double viscosity_negative = rCurrentProcessInfo[VISCOSITY_AIR];
		double density_positive = rCurrentProcessInfo[DENSITY_AIR];
		double density_negative = rCurrentProcessInfo[DENSITY_AIR];

		//bool calculate_non_newtonian_viscosity = false;

		if (has_negative_node_distance==true) //that means that the negative side of materials will certainly be the "WATER" material
		{
			viscosity_negative = rCurrentProcessInfo[VISCOSITY_WATER];
			density_negative = rCurrentProcessInfo[DENSITY_WATER];
			//calculate_non_newtonian_viscosity = true;
		}
		if (has_positive_node_corrected_distance==true)
		{
			viscosity_positive = 0.0001;
			density_positive = 1.0 ;
		}

		if (split_element_distance==false && split_element_corrected_distance==false)
		{
			//double element_mean_distance = 0.25*(distances(0)+distances(1)+distances(2)+distances(3));
			//const double negative_fluid_fraction = 0.5*(1.0-(sin(3.14159*element_mean_distance*0.5)));
			//double density = density_negative*(negative_fluid_fraction)+density_positive*(1.0-negative_fluid_fraction);

			/*
			for (unsigned int j = 0; j < 4; j++) //we go through the 3 standard shape functions
			{
				element_temperature += temperatures[j]*one_quarter;
			}
			*/
			BoundedMatrix<double, TNumNodes , TNumNodes*TDim > D_matrix; //(gradient)
			noalias(D_matrix) = ZeroMatrix(TNumNodes , TNumNodes*TDim );
			BoundedMatrix<double, TNumNodes, TNumNodes > Laplacian_matrix;
			noalias(Laplacian_matrix) = ZeroMatrix(TNumNodes,TNumNodes);
			BoundedMatrix<double, LocalSize, LocalSize > Mass_matrix;
			noalias(Mass_matrix) = ZeroMatrix(LocalSize,LocalSize);
			for (unsigned int i = 0; i < TNumNodes; i++)
			{
				for (unsigned int j = 0; j < TNumNodes; j++)
				{
					for (unsigned int k = 0; k < TDim; k++)
						D_matrix(i, (j*TDim+k) ) =  DN_DX(j,k)*Area*one_quarter; //(i,k) -> standard
						//D_matrix(i, (j*3+k) ) =  -DN_DX(i,k)*Area*one_quarter; //(i,k) -> integrated by parts ( or G)
				}
			}

			double density = 1.0;
			double viscosity = 1.0;
			double viscosity_for_tau=0.0;;
			if (distances(0)<0.0)
			{
				density=density_negative;
				viscosity=viscosity_negative;
				/*
				if (element_temperature>100.0)
						viscosity=viscosity_water*(1.0 - (element_temperature-100.0)*0.1 );
				if (viscosity<1.0)
						viscosity=1.0;
				*/
				//this->AddViscousTerm(rLeftHandSideMatrix,DN_DX, (viscosity*Area) );
				viscosity_for_tau=viscosity_negative;

			}
			else
			{
				density = density_positive;//this->CalculateAirDensity(element_temperature);
				viscosity = viscosity_positive;//this->CalculateAirViscosity(element_temperature);
				//this->AddViscousTerm(rLeftHandSideMatrix,DN_DX, (viscosity*Area) );
				viscosity_for_tau=viscosity_positive;
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
				for (SizeType d = 1; d < TDim; d++)
					sqElemSize += Edge[d]*Edge[d];
				//
				for (SizeType i = 2; i < TNumNodes; i++)
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
				for (unsigned int j=0; j<TDim ; j++)
					Mass_matrix (i*(TDim+1)+j,i*(TDim+1)+j) = Weight;

			Laplacian_matrix =  -prod(DN_DX,trans(DN_DX))*Area/density;
			Laplacian_matrix*=TauOne;

			for (unsigned int i = 0; i < TNumNodes; i++)
			{
				for (unsigned int j = 0; j < TNumNodes; j++)
				{
					for (unsigned int k = 0; k < TDim; k++)
					{
						rLeftHandSideMatrix(i*(TDim+1)+TDim, j*(TDim+1)+k ) -= D_matrix(i,j*TDim+k);//  + DN_DX(i,0)*TauOne*one_third*Area/delta_t;
						rLeftHandSideMatrix(j*(TDim+1)+k, i*(TDim+1)+TDim ) -=  D_matrix(i,j*TDim+k);
					}
					rLeftHandSideMatrix(i*(TDim+1)+TDim, j*(TDim+1)+TDim ) += Laplacian_matrix(i,j);
				}
			}


			//and finally the RHS of the velocity plus the RHS on the pressure due to the stabilzation terms:
			array_1d<double,3> mean_velocity = ZeroVector(3);
			double divergence_n = 0.0;
			for (unsigned int i = 0; i < TNumNodes; i++)
			{
				mean_velocity += GetGeometry()[i].FastGetSolutionStepValue(VELOCITY)*one_quarter/delta_t;
				array_1d<double,3> velocity = GetGeometry()[i].FastGetSolutionStepValue(VELOCITY);
				for (unsigned int j = 0; j < TDim; j++)
					divergence_n += one_quarter*Area*DN_DX(i,j)*velocity(i);
			}

			array_1d<double, TDim > rhs_stab = ZeroVector(TDim);
			for (unsigned int i = 0; i < TNumNodes; i++)
			{
				array_1d<double,3>& node_press_proj = GetGeometry()[i].FastGetSolutionStepValue(PRESS_PROJ);
				for (unsigned int j = 0; j < TDim; j++)
					rhs_stab[j] += node_press_proj(j)*one_quarter;
			}

			const bool use_press_proj=false;
			for (unsigned int i = 0; i < TNumNodes; i++)
			{
				for (unsigned int j = 0; j < TDim; j++)
					rRightHandSideVector(i*(TDim+1)+j) += one_quarter*Area*gravity(j)*density;


				if (use_press_proj)
				{
					for (unsigned int j = 0; j < TDim; j++)
						rRightHandSideVector(i*(TDim+1)+TDim) -= TauOne*Area*(DN_DX(i,j)*rhs_stab(j));
				}
				else
				{
					for (unsigned int j = 0; j < TDim; j++)
						rRightHandSideVector(i*(TDim+1)+TDim) -= TauOne*Area*(DN_DX(i,j)*(gravity(j)));
				}

			}


			noalias(rRightHandSideVector) += prod((Mass_matrix/delta_t),previous_vel_and_press);

			noalias(rLeftHandSideMatrix) += Mass_matrix/delta_t;


		}
		else //split element:
		{
			//KRATOS_WATCH("artt")

			// if the "corrected_distance" is the one defining the interface, then it will be the one dominating because it is more accurate:
			if (split_element_corrected_distance==true)
			{
				distances = corrected_distances;
			}

			//TOOLS NEEDED TO FIND THE ENRICHMENT SHAPE FUNCTIONS
			//get position of the cut surface
			//KRATOS_THROW_ERROR(std::logic_error, "IMPLICIT STEP FIRST STEP NOT YET IMPLEMENTED IN 3D.. USE LOCALFINALVELOCITY", "");
			//TOOLS NEEDED TO FIND THE ENRICHMENT SHAPE FUNCTIONS
			//get position of the cut surface
			array_1d<double,(TDim-1)*3>  densities(0);
			array_1d<double,(TDim-1)*3>  viscosities(0);

			BoundedMatrix<double,(TDim-1)*3, 2> Nenriched;
			array_1d<double,(TDim-1)*3>  volumes(0);
			BoundedMatrix<double,TNumNodes, TDim > coords;
			BoundedMatrix<double,(TDim-1)*3, TNumNodes > Ngauss;
			array_1d<double,(TDim-1)*3>  signs(0);
			std::vector< Matrix > gauss_gradients((TDim-1)*3);
			//fill coordinates
		   	const int  enrich_velocity_dofs=0;
			const int  enrich_pressure_dofs=1;
			const int  enrich_pressure_offset=0;

			BoundedMatrix<double, TNumNodes, TNumNodes> Laplacian_matrix =ZeroMatrix(TNumNodes,TNumNodes); //standard pressures. we will keep these dof so the size remains
			BoundedMatrix<double, TNumNodes, TNumNodes*TDim+enrich_velocity_dofs > D_matrix ; //(divergence) //this matrix will increase its size since we need more
			noalias(D_matrix) = ZeroMatrix(TNumNodes, TNumNodes*TDim+enrich_velocity_dofs);
			BoundedMatrix<double, TNumNodes, TNumNodes*TDim+enrich_velocity_dofs > G_matrix; //(gradient)
			noalias(G_matrix) = ZeroMatrix(TNumNodes, TNumNodes*TDim+enrich_velocity_dofs);

			BoundedMatrix<double, LocalSize+enrich_velocity_dofs+enrich_pressure_dofs, LocalSize+enrich_velocity_dofs+enrich_pressure_dofs > Momentum_matrix; //2 vel + 1 pressure per node   plus 2 enriched pressures and 2 enriched velocities
			noalias(Momentum_matrix) = ZeroMatrix(LocalSize+enrich_velocity_dofs+enrich_pressure_dofs,LocalSize+enrich_velocity_dofs+enrich_pressure_dofs);

			//unsigned int single_triangle_node = 1;
			for (unsigned int i = 0; i < TNumNodes; i++)
			{
				const array_1d<double, 3 > & xyz = this->GetGeometry()[i].Coordinates();
				//volumes[i] = 0.0;
				//KRATOS_WATCH(distances(i));
				for (unsigned int j = 0; j < TDim; j++)
					coords(i, j) = xyz[j];
			}
			for (unsigned int i = 0; i < (TDim-1)*3; i++)
				gauss_gradients[i].resize(2, TDim , false);  //2 values of the 2 shape functions, and derivates in (xyz) direction).

			unsigned int ndivisions = EnrichmentUtilities::CalculateEnrichedShapeFuncions(coords, DN_DX, distances, volumes, Ngauss, signs, gauss_gradients, Nenriched);



			double Weight=0.0;
			double Weight_for_tau=0.0;
			double mass=0.0;
			for (unsigned int i=0;i<ndivisions;i++) //we go over the six(three in 2d) partitions;
			{
				//double partition_temperature=0.0;
				//for (unsigned int j = 0; j < 4; j++)
				//	partition_temperature += temperatures[j]*Ngauss(i,j);
				if (signs(i)>0.0)
				{
					double viscosity = viscosity_positive;//CalculateAirViscosity(partition_temperature);
					densities(i) = density_positive;//CalculateAirDensity(partition_temperature);
					viscosities(i) = viscosity_positive;//CalculateAirDensity(partition_temperature);
					Weight += volumes(i)*viscosity;
					Weight_for_tau += volumes(i)*viscosity;
					mass += volumes(i)*densities(i);
				}
				else
				{
					double viscosity=viscosity_negative;
					//if (partition_temperature>100.0)
					//	viscosity=viscosity_water;//*(1.0 - (partition_temperature-100.0)*0.1 );
					//if (viscosity<1.0) //(crude oil viscosity)
					//	viscosity=1.0;
					//viscosities(i)=viscosity;
					Weight += volumes(i)*viscosity;
					Weight_for_tau += volumes(i)*viscosity;
					densities(i)=density_negative;
					viscosities(i) = viscosity_negative;//CalculateAirDensity(partition_temperature);
					mass += volumes(i)*density_negative;
				}
			}

			for (unsigned int i = 0; i < TNumNodes; i++)
			{
				for (unsigned int j = 0; j < TNumNodes; j++)
				{
					for (unsigned int k = 0; k < TDim; k++)
						D_matrix(i, (j*TDim+k) ) =  DN_DX(j,k)*Area*one_quarter; //(i,k) -> standard
						//D_matrix(i, (j*3+k) ) =  -DN_DX(i,k)*Area*one_quarter; //(i,k) -> integrated by parts ( or G)
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
			this->AddViscousTerm(rLeftHandSideMatrix,DN_DX, (Weight) );


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
				for (SizeType d = 1; d < TDim; d++)
					sqElemSize += Edge[d]*Edge[d];
				//
				for (SizeType i = 2; i < TNumNodes; i++)
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
				TauOne =  0.01 / ( ( 1.0/ delta_t + 4.0 * element_viscosity_for_tau / (mElemSize * mElemSize * element_density) + 2.0 * AdvVelNorm / mElemSize) );
			}



			array_1d<double,TNumNodes>  lumped_mass = ZeroVector(TNumNodes);
			double condensed_dof_mass1=0.0;
			double condensed_dof_mass2=0.0;
			//double condensed_dof_mass3=0.0;
			double Laplacian_matrix_factor=0.0;
			for (unsigned int i = 0; i < ndivisions; i++) //partitions
			{
				//KRATOS_WATCH(enrich_lhs(i));
				for (unsigned int j = 0; j < TNumNodes ; j++) // local node (or shape function)
				{

					lumped_mass(j) +=  volumes(i)*Ngauss(i,j)*densities(i);
					for (unsigned int k = 0; k < TDim ; k++)
						Momentum_matrix((TDim+1)*j+k,(TDim+1)*j+k) += volumes(i)*Ngauss(i,j)*densities(i) / delta_t;
				}


				if (enrich_velocity_dofs!=0)
				{
					condensed_dof_mass1 += volumes(i)*Nenriched(i,0)*densities(i);
					condensed_dof_mass2 += volumes(i)*Nenriched(i,0)*densities(i);
					//condensed_dof_mass3 += volumes(i)*Nenriched(i,0)*densities(i);

					Momentum_matrix(LocalSize+0,LocalSize+0) += volumes(i)*Nenriched(i,0)*densities(i) / delta_t;	 //degrees of freedom 10 and 11 are the enrichment velocities. degrees of freedom 12 and 13 are the enrichment pressures
					Momentum_matrix(LocalSize+1,LocalSize+1) += volumes(i)*Nenriched(i,0)*densities(i) / delta_t;
					//if(enrich_velocity_dofs==3)
					//	Momentum_matrix(LocalSize+2,LocalSize+2) += volumes(i)*Nenriched(i,0)*densities(i) / delta_t;
				}

				Laplacian_matrix_factor +=volumes(i)/densities(i);

			}



			Laplacian_matrix =  -prod(DN_DX,trans(DN_DX))*Laplacian_matrix_factor;
			Laplacian_matrix*=TauOne;



			BoundedMatrix<double, enrich_pressure_dofs, enrich_pressure_dofs > Laplacian_enrich;
			noalias(Laplacian_enrich) = ZeroMatrix(enrich_pressure_dofs,enrich_pressure_dofs);
			BoundedMatrix<double, enrich_pressure_dofs, TNumNodes > mixed_Laplacian;
			noalias(mixed_Laplacian) = ZeroMatrix(enrich_pressure_dofs,TNumNodes);
			//To calculate D* =int (N*DN_DX_enrich).  we must use the N at the gauss point of the partition (returned by the enrichment util) and multiply by the area.
			//notice that this first row is only the upper part of the mixed D. we also need the lower one for the jump
			BoundedMatrix<double, enrich_pressure_dofs, TNumNodes*TDim + enrich_velocity_dofs> D_matrix_mixed;
			noalias(D_matrix_mixed) = ZeroMatrix(enrich_pressure_dofs,TNumNodes*TDim + enrich_velocity_dofs);
			BoundedMatrix<double, enrich_pressure_dofs, TNumNodes*TDim + enrich_velocity_dofs > G_matrix_mixed;
			noalias(G_matrix_mixed) = ZeroMatrix(enrich_pressure_dofs,TNumNodes*TDim + enrich_velocity_dofs);

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
					{
						for (unsigned int l = 0; l < TDim; l++)
							Laplacian_enrich(j,k) -= ((gauss_gradients[i](j+enrich_pressure_offset,l)*gauss_gradients[i](k+enrich_pressure_offset,l)))*volumes(i)/densities(i);
					}
				//and the 'mixed laplacians' (standard shape functions * enrichments)


				//array_1d<double,3> previous_velocity_at_gauss_point=ZeroVector(3); //needed for the rhs

				//first we take care of the standard shape functions of the velocity
				for (unsigned int j = 0; j < TNumNodes; j++) //we go through the 4 standard shape functions
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
					//previous_velocity_at_gauss_point(0)+=Ngauss(i,j)*local_previous_vel(j*2+0);
					//previous_velocity_at_gauss_point(1)+=Ngauss(i,j)*local_previous_vel(j*2+1);
				}


				//we go through the remaining, new velocity dofs. (one x and one y component
				if (enrich_velocity_dofs!=0)
				{
					for (unsigned int k = 0; k < enrich_pressure_dofs; k++) //we go through the 3 enrichments
					{
						for (unsigned int l = 0; l < TDim; l++)
						{
							D_matrix_mixed(k,TDim+TNumNodes+l) += gauss_gradients[i](0,l) *volumes(i)*Nenriched(i,k+enrich_pressure_offset);
							G_matrix_mixed(k,TDim+TNumNodes+l) -= gauss_gradients[i](k+enrich_pressure_offset,l) *volumes(i)*Nenriched(i,0);
						}
					}


					//also the standard D and G matrices had to be expanded.

					for (unsigned int k = 0; k < TNumNodes; k++) //we go through the 4 standard pressures
					{
						for (unsigned int l = 0; l < TDim; l++)
						{
							D_matrix(k,TDim+TNumNodes+l) += gauss_gradients[i](0,l) *volumes(i)*Ngauss(i,k);
							G_matrix(k,TDim+TNumNodes+l) -= DN_DX(k,l) *volumes(i)*Nenriched(i,0);
						}
					}
				}


				for (unsigned int k = 0; k < enrich_pressure_dofs; k++) //we go through the 4 enrichments
				{
					for (unsigned int l = 0; l < TDim; l++)
						rhs_enrich(k+enrich_velocity_dofs) -= (gauss_gradients[i](k+enrich_pressure_offset,l)*(gravity(l)))*volumes(i);
				}
			}

			Laplacian_enrich*=TauOne;
			mixed_Laplacian*=TauOne;
			rhs_enrich *=TauOne;


			if (enrich_velocity_dofs!=0)
			{
				rhs_enrich(0) = gravity(0)*condensed_dof_mass1;
				rhs_enrich(1) = gravity(1)*condensed_dof_mass2;
				//if(enrich_velocity_dofs==3)
				//	rhs_enrich(2) = gravity(2)*condensed_dof_mass3;
			}

			BoundedMatrix<double, enrich_velocity_dofs+enrich_pressure_dofs, LocalSize > condensed_rows; //Vx1,Vy1,p1,Vx2,...
			BoundedMatrix<double, enrich_velocity_dofs+enrich_pressure_dofs, LocalSize > condensed_columns; //Vx1,Vy1,p1,Vx2,...
			BoundedMatrix<double, enrich_velocity_dofs+enrich_pressure_dofs, enrich_velocity_dofs+enrich_pressure_dofs > condensed_block; //Vx1,Vy1,p1,Vx2,...

			for (unsigned int i = 0; i <TNumNodes; i++)  //we go through the 4 nodes (standard velocity dof + standard pressure dof)
			{
				//enriched pressure dof
				for (unsigned int k = 0; k < enrich_pressure_dofs; k++)
				{
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
						condensed_columns(k,i*(TDim+1)+l)=Momentum_matrix(LocalSize+k,i*(TDim+1)+l); //add the viscosity matrix
						condensed_rows(k,i*(TDim+1)+l)=Momentum_matrix(LocalSize+k,i*(TDim+1)+l); //add the viscosity matrix
					}

					///WARNING, WHEN THE MATRIX IS REARRANGED, the condensed rows have the gradient of the pressure and the columns have the divergence. that is the reason for the mixed indexes.
					///if the gradient was not integrated by parts, then G matrix should be used in the columns instead of the rows, unlike the case of standard velocity DOFs
					condensed_rows(k,i*(TDim+1)+TDim)= -G_matrix(i, TNumNodes*TDim+k);
					condensed_columns(k,i*(TDim+1)+TDim)= -G_matrix(i, TNumNodes*TDim+k);
				}

			}
			//KRATOS_WATCH(Momentum_matrix)


			//now the condensed block terms:
			//the condensed block has 4 submatrices:    1 [ K+M*] [ G*+]  2
			//											3 [ D *+] [ L* ]  4

			//first block
			for (unsigned int i = 0; i < enrich_velocity_dofs; i++) //we go through the 3 enrichments
				for (unsigned int k = 0; k < enrich_velocity_dofs; k++) //we go through the 3 enrichments
					condensed_block(i,k)=Momentum_matrix(LocalSize+i,LocalSize+k);
			//second block
			for (unsigned int i = 0; i < enrich_velocity_dofs; i++) //we go through the 3 enrichments
				for (unsigned int k = 0; k < enrich_pressure_dofs; k++) //we go through the 3 enrichments
					condensed_block(i,k+enrich_velocity_dofs)= -G_matrix_mixed(k,TDim*TNumNodes+i);		// in  this case, we are in the gradient side and we should use the gradient if we do not integrate it by parts.
			//third block
			for (unsigned int i = 0; i < enrich_pressure_dofs; i++) //we go through the 3 enrichments
				for (unsigned int k = 0; k < enrich_velocity_dofs; k++) //we go through the 3 enrichments
					condensed_block(i+enrich_velocity_dofs,k)= -G_matrix_mixed(i,TDim*TNumNodes+k);		// in  this case, we are in the divergence side and we should use the gradient if we want to integrate the divergence it by parts.

			//fourth block
			for (unsigned int i = 0; i < enrich_pressure_dofs; i++) //we go through the 3 enrichments
				for (unsigned int k = 0; k < enrich_pressure_dofs; k++) //we go through the 3 enrichments
					condensed_block(i+enrich_velocity_dofs,k+enrich_velocity_dofs)=Laplacian_enrich(i,k);		//

			BoundedMatrix<double, enrich_velocity_dofs+enrich_pressure_dofs , enrich_velocity_dofs+enrich_pressure_dofs  > inverse_enrichments;
			this->InvertMatrix(condensed_block,inverse_enrichments);
			//condensing
			BoundedMatrix<double, LocalSize , enrich_pressure_dofs+enrich_velocity_dofs  > temp_matrix;
			temp_matrix = prod(trans(condensed_columns),inverse_enrichments);
			rLeftHandSideMatrix -=  prod(temp_matrix,condensed_rows);
			noalias(rRightHandSideVector) -= prod(temp_matrix,rhs_enrich);



			for (unsigned int i = 0; i < TNumNodes; i++)
			{
				for (unsigned int j = 0; j < TNumNodes; j++)
				{
					for (unsigned int k = 0; k < TDim; k++)
					{
						rLeftHandSideMatrix(i*(TDim+1)+TDim, j*(TDim+1)+k ) -= D_matrix(i,j*TDim+k);//  + DN_DX(i,0)*TauOne*one_third*Area/delta_t;
						rLeftHandSideMatrix(j*(TDim+1)+k, i*(TDim+1)+TDim ) -=  D_matrix(i,j*TDim+k);
					}
					rLeftHandSideMatrix(i*(TDim+1)+TDim, j*(TDim+1)+TDim ) += Laplacian_matrix(i,j);
				}
			}


			//and finally the RHS of the velocity plus the RHS on the pressure due to the stabilzation terms:
			array_1d<double,3> mean_velocity = ZeroVector(3);
			double divergence_n = 0.0;
			for (unsigned int i = 0; i < TNumNodes; i++)
			{
				mean_velocity += GetGeometry()[i].FastGetSolutionStepValue(VELOCITY)*one_quarter/delta_t;
				array_1d<double,3> velocity = GetGeometry()[i].FastGetSolutionStepValue(VELOCITY);
				for (unsigned int k = 0; k < TDim; k++)
					divergence_n += one_quarter*Area*(DN_DX(i,k)*velocity(k));
			}
			for (unsigned int i = 0; i < (TNumNodes); i++)
			{
				for (unsigned int k = 0; k < TDim; k++)
				{
					rRightHandSideVector(i*(TDim+1)+k) += lumped_mass(i)*gravity(k);
					rRightHandSideVector(i*(TDim+1)+k) += lumped_mass(i)*previous_vel_and_press(i*(TDim+1)+k)/delta_t;
					rRightHandSideVector(i*(TDim+1)+TDim) -= TauOne*Area*(DN_DX(i,k)*(gravity(k)));

				}


			}

			for (unsigned int i = 0; i < LocalSize; i++)
				for (unsigned int j = 0; j < LocalSize; j++)
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




	//************************************************************************************
	//************************************************************************************
	// this subroutine calculates the nodal contributions for the explicit steps of the
	// fractional step procedure
	void Monolithic3FluidPFEM23D::InitializeSolutionStep(ProcessInfo& CurrentProcessInfo)
	{
		KRATOS_TRY

		KRATOS_CATCH("");
	}

	//************************************************************************************
	//************************************************************************************

	void Monolithic3FluidPFEM23D::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo)
	{
		const unsigned int TDim=3;
		const SizeType TNumNodes = TDim+1;
		const SizeType LocalSize = TNumNodes*(TDim+1);
		GeometryType& rGeom = this->GetGeometry();

		SizeType LocalIndex = 0;

		if (rResult.size() != LocalSize)
			rResult.resize(LocalSize, false);

		for (SizeType i = 0; i < TNumNodes; ++i)
		{
			rResult[LocalIndex++] = rGeom[i].GetDof(VELOCITY_X).EquationId();
			rResult[LocalIndex++] = rGeom[i].GetDof(VELOCITY_Y).EquationId();
			rResult[LocalIndex++] = rGeom[i].GetDof(VELOCITY_Z).EquationId();
			rResult[LocalIndex++] = rGeom[i].GetDof(PRESSURE).EquationId();
		}
	}


	//************************************************************************************
	//************************************************************************************

	void Monolithic3FluidPFEM23D::GetDofList(DofsVectorType& ElementalDofList,ProcessInfo& CurrentProcessInfo)
	{
		const unsigned int TDim=3;
		const SizeType TNumNodes = TDim+1;
		const SizeType LocalSize = TNumNodes*(TDim+1);
		GeometryType& rGeom = this->GetGeometry();

		if (ElementalDofList.size() != LocalSize)
			ElementalDofList.resize(LocalSize);

		SizeType LocalIndex = 0;

		for (SizeType i = 0; i < TNumNodes; ++i)
		{
			ElementalDofList[LocalIndex++] = rGeom[i].pGetDof(VELOCITY_X);
			ElementalDofList[LocalIndex++] = rGeom[i].pGetDof(VELOCITY_Y);
			ElementalDofList[LocalIndex++] = rGeom[i].pGetDof(VELOCITY_Z);
			ElementalDofList[LocalIndex++] = rGeom[i].pGetDof(PRESSURE);
		}

	}


	//************************************************************************************
	//************************************************************************************

	void Monolithic3FluidPFEM23D::CalculatePressureProjection(ProcessInfo& CurrentProcessInfo)
	{
		KRATOS_TRY

		double Area;
		const unsigned int TDim=3;
		const unsigned int TNumNodes = TDim+1;
		const unsigned int LocalSize = TNumNodes*(TDim+1);

		Geometry<Node<3> >& geom = this->GetGeometry();
		BoundedMatrix<double, TNumNodes, TDim > DN_DX;
		array_1d<double, TNumNodes > N;
		GeometryUtils::CalculateGeometryData(geom, DN_DX, N, Area);
		const double mass_factor = 1.0/ (double (TNumNodes) );

		const int number_of_particles_in_elem = this->GetValue(NUMBER_OF_FLUID_PARTICLES);

		if( (number_of_particles_in_elem>0))
		{

					const double density_positive = CurrentProcessInfo[DENSITY_AIR]; //para clindro en rey 100:  0.0005 * 2.0*Area;
					const double density_negative = CurrentProcessInfo[DENSITY_WATER];

					array_1d<double,TNumNodes>  pressures = ZeroVector(TNumNodes); //to calculate the deformation Gradient F. Dimension = velocity dofs
					array_1d<double,TNumNodes> distances;
					array_1d<double,TNumNodes> corrected_distances;
					double element_mean_distance=0.0;


					bool has_positive_node_distance = false;
					bool has_negative_node_distance = false;               //interface between  standard "distance" implies that there are volumes of "air" and "water"
					bool has_positive_node_corrected_distance = false;     //interface between   "corrected distance" implies that there are volumes of "void" material and the rest is either "air or water". we can't have 2 interfaces
					bool has_negative_node_corrected_distance = false;


					for(unsigned int iii = 0; iii<TNumNodes; iii++)
					{
						//temperatures(iii) = GetGeometry()[iii].FastGetSolutionStepValue(TEMPERATURE);
						distances[iii] = this->GetGeometry()[iii].FastGetSolutionStepValue(DISTANCE);
						corrected_distances[iii] = this->GetGeometry()[iii].FastGetSolutionStepValue(CORRECTED_DISTANCE);
						pressures(iii) = GetGeometry()[iii].FastGetSolutionStepValue(PRESSURE);

						if(distances[iii]<0.0) //negative nodes are not actually really dsf
							has_negative_node_distance = true;
						else if (corrected_distances[iii]>0.0) //if they are positive on corrected_distance, they are actually "void" nodes, not "air" nodes.
							has_positive_node_distance = true;

						if(corrected_distances[iii]<0.0)
							has_negative_node_corrected_distance = true;
						else
							has_positive_node_corrected_distance = true;

						element_mean_distance+=mass_factor*this->GetGeometry()[iii].FastGetSolutionStepValue(DISTANCE);
					}


					bool split_element_distance = false;
					if (has_negative_node_distance && has_positive_node_distance)
						split_element_distance = true;

					bool split_element_corrected_distance = false;
					if (has_negative_node_corrected_distance && has_positive_node_corrected_distance)
						split_element_corrected_distance = true;


					if (split_element_distance==false && split_element_corrected_distance == false) //we only calculate if we have a pure fluid element
					{


						BoundedMatrix<double, TNumNodes, TNumNodes*TDim > G_matrix; //(gradient)
						noalias(G_matrix) = ZeroMatrix(TNumNodes, TNumNodes*TDim);

						//const double water_fraction = -( element_mean_distance - 1.0) *0.5;
						const double negative_phase_fraction = 0.5*(1.0-(sin(3.14159*element_mean_distance*0.5)));
						double density = density_negative*(negative_phase_fraction)+density_positive*(1.0-negative_phase_fraction);

						for (unsigned int i = 0; i < (TNumNodes); i++) //loop in shape functions (row)
						{
							for (unsigned int j = 0; j < (TNumNodes) ; j++) //loop in shape functions (column)
							{
								for (unsigned int k = 0; k < (TDim) ; k++) //x,y,(z)
								{
									G_matrix(i, (j*TDim)+k ) = DN_DX(i,k)*Area*mass_factor; //mass_factor=(1/3 in 2d, 1/4 in 3d)
								}
							}
						}
						G_matrix /=density;

						//now we save the data:
						for (unsigned int i=0; i!=(TNumNodes); i++) //loop around the nodes of the element to add contribution to node i
						{
							array_1d<double, 3 > & current_press_proj = geom[i].FastGetSolutionStepValue(PRESS_PROJ);

							for (unsigned int j = 0; j < (TNumNodes) ; j++) //loop around the nodes of the element to add contribution of node j to node i
							{

								for (unsigned int k = 0; k < (TDim) ; k++) //x,y,(z)
								{
									current_press_proj[k] += G_matrix(j,i*(TDim)+k)*(pressures(j));///-old_pressures(j)); //gamma=0!
								}
							}

							geom[i].FastGetSolutionStepValue(NODAL_AREA) += Area*mass_factor;

						}

					} //closing the useful elements
					else// we add some small addition to the area so that we do not have a division by zero.
					{
						for (unsigned int i=0; i!=(TNumNodes); i++) //loop around the nodes of the element to add contribution to node i
						{
							geom[i].FastGetSolutionStepValue(NODAL_AREA) += Area*mass_factor*0.00000001;
						}
					}

		} //closing the if(is_inactive==false)
		else// we add some small addition to the area so that we do not have a division by zero.
		{
			for (unsigned int i=0; i!=(TNumNodes); i++) //loop around the nodes of the element to add contribution to node i
			{
				geom[i].FastGetSolutionStepValue(NODAL_AREA) += Area*mass_factor*0.00000001;
			}
		}

		KRATOS_CATCH("");
	}


	void Monolithic3FluidPFEM23D::AddViscousTerm(BoundedMatrix<double, 21, 21 > & output,
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
	void Monolithic3FluidPFEM23D::AddViscousTerm(MatrixType& rDampMatrix,
                         const BoundedMatrix<double, 4, 3 >& rShapeDeriv,
                         const double Weight)
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

		C_matrix*= Weight;

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
		//KRATOS_WATCH(C_matrix)
		//KRATOS_WATCH(B_matrix)
		//KRATOS_WATCH(rDampMatrix)
		//		KRATOS_THROW_ERROR(std::logic_error, "IMPLICIT STEP FIRST STEP NOT YET IMPLEMENTED IN 3D.. USE LOCALFINALVELOCITY", "");


	}


	template<class T>
	bool Monolithic3FluidPFEM23D::InvertMatrix(const T& input, T& inverse)
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
