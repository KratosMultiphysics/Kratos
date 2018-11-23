//   
//   Project Name:        Kratos
//   Last modified by:    $Author: it's me! $
//   Date:                $Date: 2008-08-08 23:58:38 $
//   Revision:            $Revision: 1.0 $
//
//

// Project includes
#include "includes/define.h"
#include "custom_elements/vel_enriched_2fluid_2d.h"
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
	VelocityEnrichedPFEM22D::VelocityEnrichedPFEM22D(IndexType NewId, GeometryType::Pointer pGeometry)
		: Element(NewId, pGeometry)
	{
		//DO NOT ADD DOFS HERE!!!
	}

	//************************************************************************************
	//************************************************************************************
	VelocityEnrichedPFEM22D::VelocityEnrichedPFEM22D(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
		: Element(NewId, pGeometry, pProperties)
	{

	}

	Element::Pointer VelocityEnrichedPFEM22D::Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const
	{
		return Element::Pointer(new VelocityEnrichedPFEM22D(NewId, GetGeometry().Create(ThisNodes), pProperties));
	}

	VelocityEnrichedPFEM22D::~VelocityEnrichedPFEM22D()
	{
	}

	//************************************************************************************
	//************************************************************************************


	void VelocityEnrichedPFEM22D::AddExplicitContribution(ProcessInfo& rCurrentProcessInfo)
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

	void VelocityEnrichedPFEM22D::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY

		if(rCurrentProcessInfo[FRACTIONAL_STEP]<100)
		{

		const unsigned int TDim=2;

	    double delta_t = rCurrentProcessInfo[DELTA_TIME];

	    double volume_correction = rCurrentProcessInfo[VOLUME_CORRECTION];

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

		//in case one node is too close to zero
		for(unsigned int iii = 0; iii<(TDim+1); iii++)
		{
			if(fabs(distances(iii)<0.01))
			{
				if(distances(iii)<0.0)
				     distances(iii)=-0.01;
				else
				     distances(iii)=0.01;
			}
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
						D_matrix(i, (j*TDim+k) ) =  DN_DX(j,k)*Area*factor; //(i,k) -> standard
						//D_matrix(i, (j*3+k) ) =  -DN_DX(i,k)*Area*factor; //(i,k) -> integrated by parts ( or G)
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

				volume_correction=0.0; //we must only increase the mass of the water phase, so we turn this coeff off for the air phase
			}
			this->AddViscousTerm(rLeftHandSideMatrix,DN_DX, (viscosity_for_tau), Area );



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
				TauOne =  1.0 / ( ( 1.0/ delta_t + 4.0 * viscosity_for_tau / (mElemSize * mElemSize * density)));// + 2.0 * AdvVelNorm / mElemSize) );
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

			const bool use_press_proj=true;

			for (unsigned int i = 0; i < (TDim+1); i++)
			{
				const array_1d<double,3>& body_force = GetGeometry()[i].FastGetSolutionStepValue(BODY_FORCE);
				for (unsigned int j = 0; j < TDim; j++)
					rRightHandSideVector(i*(TDim+1)+j) += factor*Area*body_force(j)*density;


				if (use_press_proj)
					for (unsigned int j = 0; j < TDim; j++)
						rRightHandSideVector(i*(TDim+1)+TDim) -= TauOne*Area*(DN_DX(i,j)*rhs_stab(j) + fabs(DN_DX(i,j)) * volume_correction *0.25) ;
				else
					for (unsigned int j = 0; j < TDim; j++)
						rRightHandSideVector(i*(TDim+1)+TDim) -= TauOne*Area*(DN_DX(i,j)*body_force(j) + fabs(DN_DX(i,j)) * volume_correction *0.25) ;

			}


			noalias(rRightHandSideVector) += prod((Mass_matrix/delta_t),previous_vel_and_press);

			noalias(rLeftHandSideMatrix) += Mass_matrix/delta_t;


		}
		else //split element:
		{
			//const array_1d<double,3>& body_force_0 = GetGeometry()[0].FastGetSolutionStepValue(BODY_FORCE);


			//TOOLS NEEDED TO FIND THE ENRICHMENT SHAPE FUNCTIONS
			//get position of the cut surface

			array_1d<double,3*(TDim-1)>  viscosities(0);
			array_1d<double,3*(TDim-1)>  densities(0);
			array_1d<double,3*(TDim-1)>  volumes(0);
			BoundedMatrix<double,TDim+1, TDim > coords;
			BoundedMatrix<double,3*(TDim-1), TDim+1 > Ngauss;
			array_1d<double,3*(TDim-1)>  signs(0);
			std::vector< Matrix > gauss_gradients(3*(TDim-1));



			BoundedMatrix<double,3, 5> Nenriched=ZeroMatrix(3,5);
			BoundedMatrix<double,10, 2> rGradientpositive=ZeroMatrix(10,2);
			BoundedMatrix<double,10, 2> rGradientnegative=ZeroMatrix(10,2);
			BoundedMatrix<int,2, 2> nodes=ZeroMatrix(2,2);
			BoundedMatrix<int,3,3> father_nodes;


		   	const unsigned int  enrich_velocity_dofs=8;//TDim;
			const unsigned int  enrich_pressure_dofs=0; //no pressure for the moment!
			//const unsigned int  enrich_pressure_offset=0;

			BoundedMatrix<double, TDim+1, TDim+1> Laplacian_matrix =ZeroMatrix(TDim+1, TDim+1); //standard pressures. we will keep these dof so the size remains
			BoundedMatrix<double, TDim+1, 14 > D_matrix ; //(divergence) //this matrix will increase its size since we need more
			noalias(D_matrix) = ZeroMatrix(TDim+1,14);
			BoundedMatrix<double,TDim+1,14 > G_matrix; //(gradient)
			noalias(G_matrix) = ZeroMatrix(TDim+1,14);

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

			for (unsigned int i = 0; i < 3; i++)
			{
				gauss_gradients[i].resize(5, 2, false);  //2 values of the 2 shape functions, and derivates in (xyz) direction).
			}

			unsigned int ndivisions = EnrichmentUtilities::CalculateEnrichedShapeFuncionsExtendedmodified(coords, DN_DX, distances, volumes, Ngauss, signs, gauss_gradients, Nenriched, rGradientpositive , rGradientnegative, father_nodes);

			//the original gauss_gradients returned 5 shape functions (just as adding two nodes at the interface)
			//now we uncouple the middle nodes to have two different shape functions.
			//they will be later joned using the scalar pseudo-conductivity J.
			std::vector< Matrix > gauss_gradients_discon(3);
			BoundedMatrix<double,3, 7> Nenriched_discon=ZeroMatrix(3,7);
			for (unsigned int i = 0; i < 3; i++)
			{
				gauss_gradients_discon[i].resize(7, 2, false);
				gauss_gradients_discon[i] = ZeroMatrix(7,2);

				for (unsigned int j=0; j<3 ; j++) //the original 3 shape functions remain the same
				{
					gauss_gradients_discon[i](j,0)=gauss_gradients[i](j,0);
					gauss_gradients_discon[i](j,1)=gauss_gradients[i](j,1);
					Nenriched_discon(i,j)=Nenriched(i,j);
				}
				if(signs(i)>0.0)
				{
					gauss_gradients_discon[i](3,0)=gauss_gradients[i](3,0);
					gauss_gradients_discon[i](3,1)=gauss_gradients[i](3,1);
					gauss_gradients_discon[i](4,0)=gauss_gradients[i](4,0);
					gauss_gradients_discon[i](4,1)=gauss_gradients[i](4,1);
					Nenriched_discon(i,3)=Nenriched(i,3);
					Nenriched_discon(i,4)=Nenriched(i,4);
				}
				else
				{
					gauss_gradients_discon[i](5,0)=gauss_gradients[i](3,0);
					gauss_gradients_discon[i](5,1)=gauss_gradients[i](3,1);
					gauss_gradients_discon[i](6,0)=gauss_gradients[i](4,0);
					gauss_gradients_discon[i](6,1)=gauss_gradients[i](4,1);
					Nenriched_discon(i,5)=Nenriched(i,3);
					Nenriched_discon(i,6)=Nenriched(i,4);
				}
			 }
			////////////////////////////////////////
			int index=0;
			for(int aux=0;aux<3; aux++ )
			{
				if(father_nodes(aux,0)!=father_nodes(aux,1) )
				{
					nodes(index,0)=father_nodes(aux,0);
					nodes(index,1)=father_nodes(aux,1);
					index++;
				}
			}
			double sign_correction=1.0;
			//////////
			array_1d<double,2>  interface_segment=ZeroVector(2);
			array_1d<double,2>  normaledge1=ZeroVector(2);
			array_1d<double,2>  normaledge2=ZeroVector(2);
			int indexf= nodes(0,0);
			int indexs=nodes(0,1);
			double norm=0.0;
			///////////
			interface_segment[0] = (GetGeometry()[indexf].X()-GetGeometry()[indexs].X());
			interface_segment[1] = (GetGeometry()[indexf].Y()-GetGeometry()[indexs].Y());
			norm = sqrt(  pow((interface_segment[0]),2) + pow((interface_segment[1]),2));
			double area1=norm;
			normaledge1(0)= -interface_segment[1]*sign_correction/norm;
			normaledge1(1)= interface_segment[0]*sign_correction/norm;
			/////////////////////////////////////////////////////////////////////////////
			BoundedMatrix<double, 2, 2 > InterfacePoints;
			array_1d<double,2>  auxpoint=ZeroVector(2);
			array_1d<bool,3>  cut_edges;
			double relative_position = fabs(this->GetGeometry()[indexs].FastGetSolutionStepValue(DISTANCE)/(this->GetGeometry()[indexs].FastGetSolutionStepValue(DISTANCE) - this->GetGeometry()[indexf].FastGetSolutionStepValue(DISTANCE)));
			InterfacePoints(0,0) = relative_position*GetGeometry()[indexf].X() +  (1.0-relative_position)*GetGeometry()[indexs].X();
			InterfacePoints(0,1) = relative_position* GetGeometry()[indexf].Y() +  (1.0-relative_position)* GetGeometry()[indexs].Y();

			////////////////
			double areaminus1=0.0;
			double areaplus1=0.0;
			///////////
			if (this->GetGeometry()[indexf].FastGetSolutionStepValue(DISTANCE)>0.0)
			{
				interface_segment[0] = (GetGeometry()[indexf].X()-InterfacePoints(0,0));
				interface_segment[1] = (GetGeometry()[indexf].Y()-InterfacePoints(0,1));
				norm = sqrt(  pow((interface_segment[0]),2) + pow((interface_segment[1]),2));
				areaplus1=norm;
				areaminus1=area1-areaplus1;
			}
			else
			{
				interface_segment[0] = (GetGeometry()[indexf].X()-InterfacePoints(0,0));
				interface_segment[1] = (GetGeometry()[indexf].Y()-InterfacePoints(0,1));
				norm = sqrt(  pow((interface_segment[0]),2) + pow((interface_segment[1]),2));
				areaminus1=norm;
				areaplus1=area1-areaminus1;
			}
			/////////////
			indexf= nodes(1,0);
			indexs=nodes(1,1);
			interface_segment[0] = (GetGeometry()[indexf].X()-GetGeometry()[indexs].X());
			interface_segment[1] = (GetGeometry()[indexf].Y()-GetGeometry()[indexs].Y());
			norm = sqrt(  pow((interface_segment[0]),2) + pow((interface_segment[1]),2));
			//////////
			double area2=norm;
			normaledge2(0)= -interface_segment[1]*sign_correction/norm;
			normaledge2(1)= interface_segment[0]*sign_correction/norm;
			///////////
			relative_position = fabs(this->GetGeometry()[indexs].FastGetSolutionStepValue(DISTANCE)/(this->GetGeometry()[indexs].FastGetSolutionStepValue(DISTANCE) - this->GetGeometry()[indexf].FastGetSolutionStepValue(DISTANCE)));
			InterfacePoints(0,0) = relative_position*GetGeometry()[indexf].X() +  (1.0-relative_position)*GetGeometry()[indexs].X();
			InterfacePoints(0,1) = relative_position* GetGeometry()[indexf].Y() +  (1.0-relative_position)* GetGeometry()[indexs].Y();


			/////////////////////////////////
			double areaminus2=0.0;
			double areaplus2=0.0;
			/////////////////////////////////
			if (this->GetGeometry()[indexf].FastGetSolutionStepValue(DISTANCE)>0.0)
			{
				interface_segment[0] = (GetGeometry()[indexf].X()-InterfacePoints(0,0));
				interface_segment[1] = (GetGeometry()[indexf].Y()-InterfacePoints(0,1));
				norm = sqrt(  pow((interface_segment[0]),2) + pow((interface_segment[1]),2));
				areaplus2=norm;
				areaminus2=area2-areaplus2;
			}
			else
			{
				interface_segment[0] = (GetGeometry()[indexf].X()-InterfacePoints(0,0));
				interface_segment[1] = (GetGeometry()[indexf].Y()-InterfacePoints(0,1));
				norm = sqrt(  pow((interface_segment[0]),2) + pow((interface_segment[1]),2));
				areaminus2=norm;
				areaplus2=area2-areaminus2;
			}

			///////////////////////////////////////////////////////////////////////////////////
			double interface_area=0.0;
			array_1d<double,2>  normal=ZeroVector(2);
			array_1d<double,3>  Ninterface=ZeroVector(3);
			BoundedMatrix<double, 2, 2 > rInterfacePoints;
			rInterfacePoints=ZeroMatrix(2,2);
			this->CalculateInterfaceNormal(coords, distances, normal,interface_area, Ninterface,rInterfacePoints);
			///////////////////////////////////////////////////////////////////////////////////





			double Weight=0.0;
			double Weight_for_tau=0.0;
			double mass=0.0;
			double negative_area=0.0;
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
					negative_area += volumes(i);

				}
			}

			const double element_viscosity_for_tau=Weight_for_tau/Area;
			const double element_density= mass/Area;

             /*
			for (unsigned int i = 0; i < TDim+1; i++)
			{
				for (unsigned int j = 0; j < TDim+1; j++)
				{
					for (unsigned int k = 0; k < TDim; k++)
						D_matrix(i, (j*TDim+k) ) =  DN_DX(j,k)*Area*factor; //(i,k) -> standard
						//D_matrix(i, (j*TDim+k) ) =  -DN_DX(i,k)*Area*factor; //(i,k) -> integrated by parts ( or G)
				}
			}
			*/

			this->AddViscousTerm(Momentum_matrix,
								DN_DX,
								distances,
								gauss_gradients_discon,
								viscosities,
								signs,
								volumes,
								ndivisions);


			KRATOS_WATCH(Momentum_matrix)

			const double  one_w=1.0;
			const double  two_w=1.0;

			//Afer adding the viscosity matrix we add the boundary terms:

			//first the block relating the enrichment dofs with themselves
			Matrix Flux_matrix_enrich_enrich;
			Flux_matrix_enrich_enrich.resize(8,8);
			Flux_matrix_enrich_enrich = ZeroMatrix(8,8);

			//positive side
			Flux_matrix_enrich_enrich(0,0) = one_w *( rGradientpositive(3,0) * normaledge1(0) ) * viscosity_air * areaplus1 * 0.5 ;//x
			Flux_matrix_enrich_enrich(1,1) = one_w *( rGradientpositive(3,1) * normaledge1(1) ) * viscosity_air * areaplus1 * 0.5 ;//y
			Flux_matrix_enrich_enrich(0,2) = one_w *( rGradientpositive(4,0) * normaledge1(0) ) * viscosity_air * areaplus1 * 0.5 ; //x
			Flux_matrix_enrich_enrich(1,3) = one_w *( rGradientpositive(4,1) * normaledge1(1) ) * viscosity_air * areaplus1 * 0.5 ; //y
			Flux_matrix_enrich_enrich(2,0) = two_w *( rGradientpositive(8,0) * normaledge2(0) ) * viscosity_air * areaplus2 * 0.5 ; //y
			Flux_matrix_enrich_enrich(3,1) = two_w *( rGradientpositive(8,1) * normaledge2(1) ) * viscosity_air * areaplus2 * 0.5 ; //y
			Flux_matrix_enrich_enrich(2,2) = two_w *( rGradientpositive(9,0) * normaledge2(0) ) * viscosity_air * areaplus2 * 0.5 ; //x
			Flux_matrix_enrich_enrich(3,3) = two_w *( rGradientpositive(9,1) * normaledge2(1) ) * viscosity_air * areaplus2 * 0.5 ; //y
			//negative side
			Flux_matrix_enrich_enrich(4,4) = one_w *( rGradientnegative(3,0) * normaledge1(0) ) * viscosity_water * areaminus1 * 0.5 ;//x
			Flux_matrix_enrich_enrich(5,5) = one_w *( rGradientnegative(3,1) * normaledge1(1) ) * viscosity_water * areaminus1 * 0.5 ;//y
			Flux_matrix_enrich_enrich(4,6) = one_w *( rGradientnegative(4,0) * normaledge1(0) ) * viscosity_water * areaminus1 * 0.5 ;//x
			Flux_matrix_enrich_enrich(5,7) = one_w *( rGradientnegative(4,1) * normaledge1(1) ) * viscosity_water * areaminus1 * 0.5 ;//y
			Flux_matrix_enrich_enrich(6,4) = two_w *( rGradientnegative(8,0) * normaledge2(0) ) * viscosity_water * areaminus2 * 0.5 ;//x
			Flux_matrix_enrich_enrich(7,5) = two_w *( rGradientnegative(8,1) * normaledge2(1) ) * viscosity_water * areaminus2 * 0.5 ;//y
			Flux_matrix_enrich_enrich(6,6) = two_w *( rGradientnegative(9,0) * normaledge2(0) ) * viscosity_water * areaminus2 * 0.5 ;//x
			Flux_matrix_enrich_enrich(7,7) = two_w *( rGradientnegative(9,1) * normaledge2(1) ) * viscosity_water * areaminus2 * 0.5 ;//y


			BoundedMatrix<double, 8, 6 > Flux_matrix_enrich_standard;
			noalias(Flux_matrix_enrich_standard) = ZeroMatrix(8,6);

			bool added_positive_side1=false;
			bool added_negative_side1=false;

			for (unsigned int i = 0; i < ndivisions; i++)
			{
				if(signs(i)>0.0)
				{
					if (added_positive_side1==false)
					{
						for (unsigned int j = 0; j < 3; j++)
						{
							Flux_matrix_enrich_standard(0,j*2) += one_w * (rGradientpositive(j,0) * 1.0 * normaledge1(0)) * 0.5 * areaplus1 * viscosity_air;
							Flux_matrix_enrich_standard(1,j*2+1) += one_w * (rGradientpositive(j,1) * 1.0 * normaledge1(1)) * 0.5 * areaplus1 * viscosity_air;
							Flux_matrix_enrich_standard(2,j*2) += two_w * (rGradientpositive(j+5,0) * 1.0 * normaledge2(0)) * 0.5 * areaplus2 * viscosity_air ;
							Flux_matrix_enrich_standard(3,j*2+1) += two_w * ( rGradientpositive(j+5,1) * 1.0 * normaledge2(1)) * 0.5 * areaplus2 * viscosity_air ;
						}
						added_positive_side1=true;
					}
				}
				else
				{
					if (added_negative_side1==false)
					{
						for (unsigned int j = 0; j < 3; j++)
						{
							Flux_matrix_enrich_standard(4,j*2) += one_w *(rGradientnegative(j,0) * normaledge1(0)) * 0.5 * areaminus1 * viscosity_water;
							Flux_matrix_enrich_standard(5,j*2+1) += one_w *( rGradientnegative(j,1) * normaledge1(1)) * 0.5 * areaminus1 * viscosity_water;
							Flux_matrix_enrich_standard(6,j*2) += two_w * (rGradientnegative(j+5,0) * normaledge2(0) ) * 0.5 * areaminus2 * viscosity_water;
							Flux_matrix_enrich_standard(7,j*2+1) += two_w * ( rGradientnegative(j+5,1) * normaledge2(1)) * 0.5 * areaminus2 * viscosity_water;
						}
						added_negative_side1=true;
					 }
				 }
			 }



	    //KRATOS_WATCH(normal);
	    //KRATOS_WATCH(interface_area);

	    Matrix Jump;
	    Jump.resize(8,8);
	    Jump = ZeroMatrix(8,8);
		double J = this->GetValue(CONDUCTIVITY);
		KRATOS_WATCH(J);

		for(unsigned int iii = 0; iii<TDim; iii++)
		{
			Jump(iii,iii)= J * interface_area*0.3333333; //0.5;// + J * 0.5 * 10000.0;//0.3333333;
			Jump(iii,2+iii)= J * interface_area*0.1666666;//0;// - J * 0.5* 10000.0;//0.1666666;
			Jump(iii,4+iii)=-J * interface_area*0.3333333;//0.5;// - J * 0.5* 10000.0;//0.3333333;
			Jump(iii,6+iii)=-J * interface_area*0.1666666;//0;// + J * 0.5* 10000.0;//0.1666666;

			Jump(2+iii,iii)= J * interface_area*0.1666666;//0;// + J * 0.5* 10000.0;//0.1666666;
			Jump(2+iii,2+iii)= J * interface_area*0.3333333;//0.5;// - J * 0.5* 10000.0;//0.3333333;
			Jump(2+iii,4+iii)=-J * interface_area*0.1666666;//0;// - J * 0.5* 10000.0;//0.1666666;
			Jump(2+iii,6+iii)=-J * interface_area*0.3333333;//0.5;// + J * 0.5* 10000.0;//0.3333333;

			Jump(4+iii,iii)=-J * interface_area*0.3333333;//0.5;// - J * 0.5* 10000.0;//0.3333333;
			Jump(4+iii,2+iii)=-J * interface_area*0.1666666;//0;// + J * 0.5* 10000.0;//0.1666666;
			Jump(4+iii,4+iii)= J * interface_area*0.3333333;//0.5;// + J * 0.5* 10000.0;//0.3333333;
			Jump(4+iii,6+iii)= J * interface_area*0.1666666;//0;// - J * 0.5* 10000.0;//0.1666666;

			Jump(6+iii,iii)=-J * interface_area*0.1666666;//0;// - J * 0.5* 10000.0;//0.1666666;
			Jump(6+iii,2+iii)=-J * interface_area*0.3333333;//0.5;// + J * 0.5* 10000.0;//0.3333333;
			Jump(6+iii,4+iii)= J * interface_area*0.1666666;//0;// + J * 0.5* 10000.0;//0.1666666;
			Jump(6+iii,6+iii)= J * interface_area*0.3333333;//0.5;//- J * 0.5* 10000.0;//0.3333333;
		}
	    //KRATOS_WATCH(Jump);



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
				TauOne =  1.0 / ( ( 1.0/ delta_t + 4.0 * element_viscosity_for_tau / (mElemSize * mElemSize * element_density)));// + 2.0 * AdvVelNorm / mElemSize) );
			}



			//array_1d<double,TDim+1>  lumped_mass_std_nodes = ZeroVector(TDim+1);
			//array_1d<double,4>  lumped_mass_enrich_nodes = ZeroVector(4);
			array_1d<double,7>  lumped_mass = ZeroVector(7);

			double Laplacian_matrix_factor=0.0;
			for (unsigned int i = 0; i < ndivisions; i++) //partitions
			{
				/*
				for (unsigned int j = 0; j < TDim+1 ; j++) // local node (or replacement shape function)
				{

					lumped_mass(j) +=  volumes(i)*Nenriched(i,j)*densities(i);
					Momentum_matrix((TDim+1)*j,(TDim+1)*j) += volumes(i)*Nenriched(i,j)*densities(i) / delta_t;
					Momentum_matrix((TDim+1)*j+1,(TDim+1)*j+1) += volumes(i)*Nenriched(i,j)*densities(i) / delta_t;
				}


				for (unsigned int j = 0; j < 4 ; j++) // new nodes for the enrichment (4 in total
				{
					lumped_mass_enrich_nodes[j] += volumes(i)*Nenriched(i,j+3)*densities(i); // we have to add 3 because we are moving

					Momentum_matrix((TDim+1)*(j+3),(TDim+1)*(j+3)) += volumes(i)*Nenriched(i,j+3)*densities(i) / delta_t;
					Momentum_matrix((TDim+1)*(j+3)+1,(TDim+1)*(j+3)+1) += volumes(i)*Nenriched(i,j+3)*densities(i) / delta_t;
				}
				*/
				for (unsigned int j = 0; j < 7 ; j++) // local node (or replacement shape function)
				{

					lumped_mass(j) +=  volumes(i)*Nenriched(i,j)*densities(i);
					///unsigned int j_assembly=j*3; //for asembly purpouses, we have to do this since the enrichment nodes do not have pressure!, so they do not have to leave gaps for the pressure
					///if(j>2)
					///	j_assembly=9+(j-3)*2;
					//removing mass
					//Momentum_matrix(j_assembly,j_assembly) += volumes(i)*Nenriched(i,j)*densities(i) / delta_t;
					//Momentum_matrix(j_assembly+1,j_assembly+1) += volumes(i)*Nenriched(i,j)*densities(i) / delta_t;
				}

				Laplacian_matrix_factor +=volumes(i)/densities(i);

			}


			//we are using the standard shape functions and nodes for the pressure, so this remain as the original matrix for standard 3 noded triangles.
			Laplacian_matrix =  -prod(DN_DX,trans(DN_DX))*Laplacian_matrix_factor;
			Laplacian_matrix*=TauOne;

			//since we do not have enrichment functions for the pressure, for the moment the laplacian enrich and all these things dissapear: (commenting lines below)
			//BoundedMatrix<double, enrich_pressure_dofs, enrich_pressure_dofs > Laplacian_enrich;
			//noalias(Laplacian_enrich) = ZeroMatrix(enrich_pressure_dofs,enrich_pressure_dofs);
			//BoundedMatrix<double, enrich_pressure_dofs, (TDim+1) > mixed_Laplacian;
			//noalias(mixed_Laplacian) = ZeroMatrix(enrich_pressure_dofs,(TDim+1));
			////To calculate D* =int (N*DN_DX_enrich).  we must use the N at the gauss point of the partition (returned by the enrichment util) and multiply by the area.
			////notice that this first row is only the upper part of the mixed D. we also need the lower one for the jump
			//BoundedMatrix<double, enrich_pressure_dofs, (TDim+1)*TDim + enrich_velocity_dofs> D_matrix_mixed;
			//noalias(D_matrix_mixed) = ZeroMatrix(enrich_pressure_dofs,(TDim+1)*TDim + enrich_velocity_dofs);
			//BoundedMatrix<double, enrich_pressure_dofs, (TDim+1)*TDim + enrich_velocity_dofs > G_matrix_mixed;
			//noalias(G_matrix_mixed) = ZeroMatrix(enrich_pressure_dofs,(TDim+1)*TDim + enrich_velocity_dofs);

			//double rhs_enrich = 0.0; //the rhs term of the enriched dof
			//array_1d<double,( enrich_pressure_dofs + enrich_velocity_dofs )> rhs_complete;
			//noalias(rhs_enrich) = ZeroVector(enrich_pressure_dofs + enrich_velocity_dofs);

				//			KRATOS_WATCH(Ngauss)
				//KRATOS_WATCH(Nenriched)

			for (unsigned int i = 0; i < ndivisions; i++)  //we go through the 3 partitions of the domain
			{

				////no laplacian enrich since we have no pressure enrichment
				//for (unsigned int j = 0; j < enrich_pressure_dofs; j++) //we go through the 4 enrichments
				//	for (unsigned int k = 0; k < enrich_pressure_dofs; k++) //we go through the 4 enrichments
				//	    for (unsigned int l = 0; l < TDim; l++)
				//			Laplacian_enrich(j,k) -= (gauss_gradients[i](j+enrich_pressure_offset,l)*gauss_gradients[i](k+enrich_pressure_offset,l))*volumes(i)/densities(i);

				////and also no mixed laplacian
				//and the 'mixed laplacians' (standard shape functions * enrichments)
				//first we take care of the standard shape functions of the velocity
				/*
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
				*/


				//we go through all the velocity dofs. (one x and one y component
				{
					////no mixed laplacian
					/*
					for (unsigned int k = 0; k < enrich_pressure_dofs; k++) //we go through the 3 enrichments
					{
						for (unsigned int l = 0; l < TDim; l++)
						{
							D_matrix_mixed(k,(TDim+1)*TDim+l) += gauss_gradients[i](0,l) *volumes(i)*Nenriched(i,k+enrich_pressure_offset);
							G_matrix_mixed(k,(TDim+1)*TDim+l) -= gauss_gradients[i](k+enrich_pressure_offset,l) *volumes(i)*Nenriched(i,0);
						}
					}
					*/


					////also the standard D and G matrices had to be expanded.
					////we have to  loop in the 7 nodes!!!
					///buscar error en esta parte en el código original, k debería ser TDim+1 para ya q son las 3 funciones de forma de la por!!!!!!
					for (unsigned int j = 0; j < 7; j++) //we go through the 7 shape functions
					{
						for (unsigned int k = 0; k < TDim+1; k++) //we go through the 4 standard pressures
						{
							for (unsigned int l = 0; l < TDim; l++) //x,y
							{
								D_matrix(k,j*TDim+l) += gauss_gradients_discon[i](j,l) *volumes(i)*Ngauss(i,k);
								G_matrix(k,j*TDim+l) -= DN_DX(k,l) *volumes(i)*Nenriched(i,j);
							}
						}
					}
				}


				//no pressure enrichment
				//for (unsigned int k = 0; k < enrich_pressure_dofs; k++) //we go through the k enrichments
				//	for (unsigned int l = 0; l < TDim; l++)
				//		rhs_enrich(k+enrich_velocity_dofs) -= (gauss_gradients[i](k+enrich_pressure_offset,l)*body_force_0(l))*volumes(i);

			}

			//Laplacian_enrich*=TauOne;
			//mixed_Laplacian*=TauOne;
			//rhs_enrich *=TauOne;

			//warning, i removed all body forces!!!
			/*
			if (enrich_velocity_dofs!=0)
			{
				rhs_enrich(0) = body_force_0(0)*condensed_dof_mass1;
				rhs_enrich(1) = body_force_0(1)*condensed_dof_mass2;
				//rhs_enrich(2) = body_force_0(2)*condensed_dof_mass3;
			}
			*/
			BoundedMatrix<double, 8+enrich_pressure_dofs, (TDim+1)*(TDim+1) > condensed_rows; //Vx1,Vy1,p1,Vx2,...
			BoundedMatrix<double, 8+enrich_pressure_dofs, (TDim+1)*(TDim+1) > condensed_columns; //Vx1,Vy1,p1,Vx2,...
			BoundedMatrix<double, 8+enrich_pressure_dofs, enrich_velocity_dofs+enrich_pressure_dofs > condensed_block; //Vx1,Vy1,p1,Vx2,...

			for (unsigned int i = 0; i <TDim+1; i++)  //we go through the 4 nodes (standard velocity dof + standard pressure dof)
			{
				//enriched pressure dof
				/*
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
				*/


				for (unsigned int k = 0; k < enrich_velocity_dofs; k++) //we go through the 3 enrichments
				{
					for (unsigned int l = 0; l < TDim; l++)
					{
						condensed_columns(k,i*(TDim+1)+l)=Momentum_matrix((TDim+1)*(TDim+1)+k,i*(TDim+1)+l); //add the viscosity matrix
						condensed_rows(k,i*(TDim+1)+l)=Momentum_matrix((TDim+1)*(TDim+1)+k,i*(TDim+1)+l); //add the viscosity matrix
					}

					///WARNING, WHEN THE MATRIX IS REARRANGED, the condensed rows have the gradient of the pressure and the columns have the divergence. that is the reason for the mixed indexes.
					///if the gradient was not integrated by parts, then G matrix should be used in the columns instead of the rows, unlike the case of standard velocity DOFs
					condensed_rows(k,i*(TDim+1)+TDim)= -D_matrix(i, (TDim+1)*TDim+k);
					condensed_columns(k,i*(TDim+1)+TDim)= -D_matrix(i, (TDim+1)*TDim+k);
				}

			}
			//KRATOS_WATCH(Momentum_matrix)

			//condensed_rows += - Flux_matrix_enrich_standard;

			//now the condensed block terms:
			//the condensed block has 4 submatrices:    1 [ K+M*] [ G*+]  2
			//											3 [ D *+] [ L* ]  4

					//first block
			for (unsigned int i = 0; i < enrich_velocity_dofs; i++) //we go through the 3 enrichments
				for (unsigned int k = 0; k < enrich_velocity_dofs; k++) //we go through the 3 enrichments
					condensed_block(i,k)=Momentum_matrix((TDim+1)*(TDim+1)+i,(TDim+1)*(TDim+1)+k);

			condensed_block += Jump; //- Flux_matrix_enrich_enrich;

			//no pressure enrichments!
			/*
			//second block
			for (unsigned int i = 0; i < enrich_velocity_dofs; i++) //we go through the 3 enrichments
				for (unsigned int k = 0; k < enrich_pressure_dofs; k++) //we go through the 3 enrichments
					condensed_block(i,k+enrich_velocity_dofs)= -D_matrix_mixed(k,(TDim+1)*TDim+i);		// in  this case, we are in the gradient side and we should use the gradient if we do not integrate it by parts.
			//third block
			for (unsigned int i = 0; i < enrich_pressure_dofs; i++) //we go through the 3 enrichments
				for (unsigned int k = 0; k < enrich_velocity_dofs; k++) //we go through the 3 enrichments
					condensed_block(i+enrich_velocity_dofs,k)= -D_matrix_mixed(i,(TDim+1)*TDim+k);		// in  this case, we are in the divergence side and we should use the gradient if we want to integrate the divergence it by parts.

			//fourth block
			for (unsigned int i = 0; i < enrich_pressure_dofs; i++) //we go through the 3 enrichments
				for (unsigned int k = 0; k < enrich_pressure_dofs; k++) //we go through the 3 enrichments
					condensed_block(i+enrich_velocity_dofs,k+enrich_velocity_dofs)=Laplacian_enrich(i,k);		//
			*/

			BoundedMatrix<double, enrich_velocity_dofs+enrich_pressure_dofs , enrich_velocity_dofs+enrich_pressure_dofs  > inverse_enrichments;
			this->InvertMatrix(condensed_block,inverse_enrichments);
			//condensing
			BoundedMatrix<double, (TDim+1)*(TDim+1) , enrich_pressure_dofs+enrich_velocity_dofs  > temp_matrix;
			temp_matrix = prod(trans(condensed_columns),inverse_enrichments);
			rLeftHandSideMatrix -=  prod(temp_matrix,condensed_rows);
			//noalias(rRightHandSideVector) -= prod(temp_matrix,rhs_enrich);



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
			/*
			for (unsigned int i = 0; i < TDim+1; i++)
			{
				const array_1d<double,3>& body_force = GetGeometry()[i].FastGetSolutionStepValue(BODY_FORCE);
				for (unsigned int l = 0; l < TDim; l++)
				{
					//rRightHandSideVector(i*(TDim+1)+l) += lumped_mass(i)*body_force(l);
					//rRightHandSideVector(i*(TDim+1)+l) += lumped_mass(i)*previous_vel_and_press(i*(TDim+1)+l)/delta_t;
					//rRightHandSideVector(i*(TDim+1)+TDim) -= TauOne*Area*(DN_DX(i,l)*body_force(l));
				}


			}
			*/
			for (unsigned int i = 0; i < (TDim+1)*(TDim+1); i++)
				for (unsigned int j = 0; j < (TDim+1)*(TDim+1); j++)
					rLeftHandSideMatrix(i,j) += Momentum_matrix(i,j);

		}




		//we must actually use  a fraction of the rigidity matrix. since the Lefthandside only has it so far, we multiply it by 0.5:
		//rLeftHandSideMatrix *= theta;

        //finally we add to the righthandside BOTH the mass matrix and the rigidity matrix:


		noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix,current_vel_and_press);
		}
		else
		{
					unsigned int TNumNodes = GetGeometry().size();
                    const int TDim=2;
					array_1d<double, TDim+1> N;
					// Check sizes and initialize
					if (rRightHandSideVector.size() != TNumNodes)
						rRightHandSideVector.resize(TNumNodes, false);
					noalias(rRightHandSideVector) = ZeroVector(TNumNodes);
					if(rLeftHandSideMatrix.size1() != TNumNodes)
						rLeftHandSideMatrix.resize(TNumNodes,TNumNodes,false);
					noalias(rLeftHandSideMatrix) = ZeroMatrix(TNumNodes,TNumNodes);

					array_1d<double,(TDim+1)> rhs_x,rhs_y,rhs_z,rhs_d;
					rhs_x = ZeroVector((TDim+1));         //resetting vectors
					rhs_y = ZeroVector((TDim+1));         //resetting vectors
					rhs_z = ZeroVector((TDim+1));         //resetting vectors
					rhs_d = ZeroVector((TDim+1));         //resetting vectors

					const int offset = rCurrentProcessInfo[WATER_PARTICLE_POINTERS_OFFSET];

					int & number_of_particles_in_elem= this->GetValue(NUMBER_OF_FLUID_PARTICLES);
					ParticlePointerVector&  element_particle_pointers =  (this->GetValue(FLUID_PARTICLE_POINTERS));
					Geometry<Node<3> >& geom = this->GetGeometry();

					for (int iii=0; iii<number_of_particles_in_elem ; iii++ )
					{

						PFEM_Particle_Fluid & pparticle = element_particle_pointers[offset+iii];

						if (pparticle.GetEraseFlag()==false)
						{

							array_1d<double,3> & position = pparticle.Coordinates();

							const array_1d<float,3>& velocity = pparticle.GetVelocity();

							const float& particle_distance = pparticle.GetDistance();  // -1 if water, +1 if air

							array_1d<double,TDim+1> N;

							double x0 = geom[0].X();
							double y0 = geom[0].Y();
							double x1 = geom[1].X();
							double y1 = geom[1].Y();
							double x2 = geom[2].X();
							double y2 = geom[2].Y();

							double area = 0.5 * ((x1 - x0)*(y2 - y0)- (y1 - y0)*(x2 - x0));
							double inv_area = 1.0 / area;


							N[0] = 0.5 * ((x1 - position[0])*(y2 - position[1])- (y1 - position[1])*(x2 - position[0])) * inv_area;
							N[1] = 0.5 * ((position[0] - x0)*(y2 - y0)- (position[1] - y0)*(x2 - x0))* inv_area;
							N[2] = 0.5 * ((x1 - x0)*(position[1] - y0)- (y1 - y0)*(position[0] - x0))* inv_area;

							for (int j=0 ; j!=(TDim+1); j++) //going through the 3/4 nodes of the element
							{
								double weight=N(j);
								for (int k=0 ; k!=(TDim+1); k++) //building the mass matrix
									rLeftHandSideMatrix(j,k) += weight*N(k);

								rhs_x[j] += weight * double(velocity[0]);
								rhs_y[j] += weight * double(velocity[1]);
								rhs_z[j] += weight * double(velocity[2]);
								rhs_d[j] += weight * double(particle_distance);

								if(rCurrentProcessInfo[FRACTIONAL_STEP]==100)
									rRightHandSideVector[j] += weight * double(particle_distance);
								if(rCurrentProcessInfo[FRACTIONAL_STEP]==101)
								    rRightHandSideVector[j] += weight * double(velocity[0]);
								if(rCurrentProcessInfo[FRACTIONAL_STEP]==102)
								    rRightHandSideVector[j] += weight * double(velocity[1]);

								//adding also a part with the lumped mass matrix to reduce overshoots and undershoots
								/*
								if(true)
								{
									double this_particle_weight = weight*elem_volume/(double(number_of_particles_in_elem))*0.1; //can be increased or reduced to change the lumped mass contrubtion
									nodes_addedweights[j]+= this_particle_weight;
									nodes_added_distance[j] += this_particle_weight*particle_distance;
									for (int k=0 ; k!=(TDim); k++) //x,y,(z)
									{
										nodes_addedvel[j*3+k] += this_particle_weight * double(velocity[k]);
									}
								}
								*/
							}
						}
					}

					array_1d<double,(TDim+1)> previous_solution = ZeroVector(TDim+1);
					for (int j=0 ; j!=(TDim+1); j++)
					{
						if(rCurrentProcessInfo[FRACTIONAL_STEP]==100)
							previous_solution[j] = this->GetGeometry()[j].FastGetSolutionStepValue(DISTANCE);
						if(rCurrentProcessInfo[FRACTIONAL_STEP]==101)
						    previous_solution[j] = this->GetGeometry()[j].FastGetSolutionStepValue(PROJECTED_VELOCITY_X);
						if(rCurrentProcessInfo[FRACTIONAL_STEP]==102)
						    previous_solution[j] = this->GetGeometry()[j].FastGetSolutionStepValue(PROJECTED_VELOCITY_Y);
					}
					noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix,previous_solution);


		}



		KRATOS_CATCH("");
	}



	//************************************************************************************
	//************************************************************************************
	// this subroutine calculates the nodal contributions for the explicit steps of the
	// fractional step procedure
	void VelocityEnrichedPFEM22D::InitializeSolutionStep(ProcessInfo& CurrentProcessInfo)
	{
		KRATOS_TRY

		KRATOS_CATCH("");
	}

	//************************************************************************************
	//************************************************************************************

	void VelocityEnrichedPFEM22D::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo)
	{
		if(CurrentProcessInfo[FRACTIONAL_STEP]<100 )
		{
		const int TDim = 2;

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
			//rResult[LocalIndex++] = rGeom[i].GetDof(VELOCITY_Z).EquationId();
			rResult[LocalIndex++] = rGeom[i].GetDof(PRESSURE).EquationId();
		}
		}
		else
		{
			const int TDim = 2;

			const SizeType NumNodes = TDim+1;
			const SizeType LocalSize = (TDim+1);
			GeometryType& rGeom = this->GetGeometry();

			SizeType LocalIndex = 0;

			if (rResult.size() != LocalSize)
				rResult.resize(LocalSize, false);

			for (SizeType i = 0; i < NumNodes; ++i)
			{
				if(CurrentProcessInfo[FRACTIONAL_STEP]==100)
					rResult[LocalIndex++] = rGeom[i].GetDof(DISTANCE).EquationId();
				if(CurrentProcessInfo[FRACTIONAL_STEP]==101)
					rResult[LocalIndex++] = rGeom[i].GetDof(PROJECTED_VELOCITY_X).EquationId();
				if(CurrentProcessInfo[FRACTIONAL_STEP]==102)
					rResult[LocalIndex++] = rGeom[i].GetDof(PROJECTED_VELOCITY_Y).EquationId();
			}
		}
	}

	//************************************************************************************
	//************************************************************************************

	void VelocityEnrichedPFEM22D::GetDofList(DofsVectorType& ElementalDofList,ProcessInfo& CurrentProcessInfo)
	{
		if(CurrentProcessInfo[FRACTIONAL_STEP]<100 )
		{
		const int TDim = 2;

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
			//ElementalDofList[LocalIndex++] = rGeom[i].pGetDof(VELOCITY_Z);
			ElementalDofList[LocalIndex++] = rGeom[i].pGetDof(PRESSURE);
		}
		}
		else
		{
			const int TDim = 2;

			const SizeType NumNodes = TDim+1;
			const SizeType LocalSize = (TDim+1);
			GeometryType& rGeom = this->GetGeometry();

			if (ElementalDofList.size() != LocalSize)
				ElementalDofList.resize(LocalSize);

			SizeType LocalIndex = 0;

			for (SizeType i = 0; i < NumNodes; ++i)
			{
				if(CurrentProcessInfo[FRACTIONAL_STEP]==100)
					ElementalDofList[LocalIndex++] = rGeom[i].pGetDof(DISTANCE);
				if(CurrentProcessInfo[FRACTIONAL_STEP]==101)
					ElementalDofList[LocalIndex++] = rGeom[i].pGetDof(PROJECTED_VELOCITY_X);
				if(CurrentProcessInfo[FRACTIONAL_STEP]==102)
					ElementalDofList[LocalIndex++] = rGeom[i].pGetDof(PROJECTED_VELOCITY_Y);
			}
		}


	}


	//************************************************************************************
	//************************************************************************************

	void VelocityEnrichedPFEM22D::CalculatePressureProjection(ProcessInfo& CurrentProcessInfo)
	{
		KRATOS_TRY

		const unsigned int TDim=2;

		double Area;
		Geometry<Node<3> >& geom = this->GetGeometry();
		BoundedMatrix<double, (TDim+1), TDim > DN_DX;
		array_1d<double, (TDim+1) > N;
		GeometryUtils::CalculateGeometryData(geom, DN_DX, N, Area);
		const double factor = 1.0/ (1.0 + double (TDim) );

		//const int number_of_particles_in_elem = this->GetValue(NUMBER_OF_FLUID_PARTICLES);

		//if( (number_of_particles_in_elem>0))
		if(true)
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
									G_matrix(i, (j*2)+k ) = DN_DX(i,k)*Area*factor; //mass_factor=(1/3 in 2d, 1/4 in 3d)
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


	void VelocityEnrichedPFEM22D::AddViscousTerm(BoundedMatrix<double, 17, 17 > & output,
										  BoundedMatrix<double, (3), 2 >& rShapeDeriv,
										  array_1d<double,3>&  distances,
										  std::vector< Matrix >& gauss_gradients,
										  array_1d<double,3>&  viscosities,
										  array_1d<double,3>&  signs,
										  array_1d<double,3>&  volumes ,
										  const unsigned int ndivisions)
	{
		BoundedMatrix<double, 14, 14 > ExtendedDampMatrix=ZeroMatrix(14,14);
		//BoundedMatrix<double, 8, 8 > rExtendedDampMatrix= ZeroMatrix(8,8);

		BoundedMatrix<double, 14,3 > B_matrix = ZeroMatrix(14,3);


		int counter=0;
		BoundedMatrix<double, 3, 3 > C_matrix = ZeroMatrix(3,3);

		for (unsigned int i=0; i!=(2); i++)
		{
			C_matrix(counter,counter)=2.0;
			counter++;
		}
		for (unsigned int i=0; i!=(1); i++)
		{
			C_matrix(counter,counter)=1.0;
			counter++;
		}

		//now the enriched part:
		//we have to construct (ndivisions) rTempDampMatrix and asseble add them to rDampMatrix
		for (unsigned int division=0; division!=ndivisions; division++)
		{
			for (unsigned int i=0; i!=(7); i++) //7 nodes!!! (3 original (with modified shape functions)  + 4 new nodes!!!!!
			{
				for (unsigned int j=0; j!=(2); j++) //x,y,z
					B_matrix(3*(i)+j,j)= gauss_gradients[division](i,j);

				//useful for both 2d and 3d:
				//relating 12 and 21 stresses
				B_matrix(i*(2)+0,2)=gauss_gradients[division](i,1);
				B_matrix(i*(2)+1,2)=gauss_gradients[division](i,0);
			}
			BoundedMatrix<double, 3 , 14  > temp_matrix = prod(C_matrix,trans(B_matrix));
			ExtendedDampMatrix += viscosities(division)*volumes(division)*prod(B_matrix, temp_matrix );
		}

		//now we put it all toghether in the big matrix:
		for (unsigned int i=0; i!=(7); i++) //7 nodes
			for (unsigned int j=0; j!=(7); j++) //7 nodes
				for (unsigned int k=0; k!=(2); k++) //x,y,(z)
					for (unsigned int l=0; l!=(2); l++) //x,y,(z)
					{
						 unsigned int i_assembly=i*3; //for asembly purpouses, we have to do this since the enrichment nodes do not have pressure!, so they do not have to leave gaps for the pressure
						 unsigned int j_assembly=j*3; //for asembly purpouses, we have to do this since the enrichment nodes do not have pressure!
						 if(i>2)
							i_assembly=9+(i-3)*2;
						 if(j>2)
							j_assembly=9+(j-3)*2;
						 output(i_assembly+k,j_assembly+l) += ExtendedDampMatrix(i*2+k,j*2+l);
					}
		//KRATOS_WATCH(ExtendedDampMatrix)
	}


		void VelocityEnrichedPFEM22D::AddViscousTerm(BoundedMatrix<double, 12, 12 > & output,
										  BoundedMatrix<double, (3), 2 >& rShapeDeriv,
										  array_1d<double,3>&  distances,
										  std::vector< Matrix >& gauss_gradients,
										  array_1d<double,3>&  viscosities,
										  array_1d<double,3>&  signs,
										  array_1d<double,3>&  volumes ,
										  const unsigned int ndivisions)
	{
		BoundedMatrix<double, 8, 8 > ExtendedDampMatrix=ZeroMatrix(8,8);
		//BoundedMatrix<double, 8, 8 > rExtendedDampMatrix= ZeroMatrix(8,8);

		BoundedMatrix<double, 8,3 > B_matrix = ZeroMatrix(8,3);


		int counter=0;
		BoundedMatrix<double, 3, 3 > C_matrix = ZeroMatrix(3,3);

		for (unsigned int i=0; i!=(2); i++)
		{
			C_matrix(counter,counter)=2.0;
			counter++;
		}
		for (unsigned int i=0; i!=(1); i++)
		{
			C_matrix(counter,counter)=1.0;
			counter++;
		}

		//now the enriched part:
		//we have to construct (ndivisions) rTempDampMatrix and asseble add them to rDampMatrix
		for (unsigned int division=0; division!=ndivisions; division++)
		{
			//standard shape functions:
			for (unsigned int i=0; i!=(3); i++) //i node
			{
				for (unsigned int j=0; j!=(2); j++) //x,y,z
					B_matrix(i*(2)+j,j)=rShapeDeriv(i,j);

				//useful for both 2d and 3d:
				//relating 12 and 21 stresses
				B_matrix(i*(2)+0,2)=rShapeDeriv(i,1);
				B_matrix(i*(2)+1,2)=rShapeDeriv(i,0);

			}


			for (unsigned int j=0; j!=(2); j++) //x,y,z
				B_matrix(3*(2)+j,j)= gauss_gradients[division](0,j);

			//useful for both 2d and 3d:
			//relating 12 and 21 stresses
			B_matrix(3*(2)+0,2)=gauss_gradients[division](0,1);
			B_matrix(3*(2)+1,2)=gauss_gradients[division](0,0);


			BoundedMatrix<double, 3 , 8  > temp_matrix = prod(C_matrix,trans(B_matrix));
			ExtendedDampMatrix += viscosities(division)*volumes(division)*prod(B_matrix, temp_matrix );
		}

		//now we put it all toghether in the big matrix:
		for (unsigned int i=0; i!=(4); i++) //4 nodes + 1dof in the new virtual node
			for (unsigned int j=0; j!=(4); j++) //4 nodes + 1dof in the new virtual node
				for (unsigned int k=0; k!=(2); k++) //x,y,(z)
					for (unsigned int l=0; l!=(2); l++) //x,y,(z)
						 output(i*3+k,j*3+l) += ExtendedDampMatrix(i*2+k,j*2+l);
		//KRATOS_WATCH(ExtendedDampMatrix)
	}




	//with matrixtype, constant coefficient
	void VelocityEnrichedPFEM22D::AddViscousTerm(MatrixType& rDampMatrix,
                         const BoundedMatrix<double, 3, 2 >& rShapeDeriv,
                         double Viscosity,const double Area)
	{


		BoundedMatrix<double, 6,3 > B_matrix = ZeroMatrix(6,3);
		for (unsigned int i=0; i!=(3); i++) //i node
		{
			for (unsigned int j=0; j!=(2); j++) //x,y,z
				B_matrix(i*(2)+j,j)=rShapeDeriv(i,j);

			//useful for both 2d and 3d:
			//relating 12 and 21 stresses
			B_matrix(i*(2)+0,2)=rShapeDeriv(i,1);
			B_matrix(i*(2)+1,2)=rShapeDeriv(i,0);
		}

		int counter=0;
		BoundedMatrix<double, 3, 3 > C_matrix = ZeroMatrix(3,3);

		for (unsigned int i=0; i!=(2); i++)
		{
			C_matrix(counter,counter)=2.0;
			counter++;
		}
		for (unsigned int i=0; i!=(1); i++)
		{
			C_matrix(counter,counter)=1.0;
			counter++;
		}

		C_matrix*= Viscosity*Area;

		BoundedMatrix<double, 3 , 6  > temp_matrix = prod(C_matrix,trans(B_matrix));
		BoundedMatrix<double, 6 , 6  > viscosity_matrix = prod(B_matrix, temp_matrix );
		for (unsigned int i=0; i!=3; i++) //i node
		{
			for (unsigned int j=0; j!=3; j++) //j neighbour
			{
				for (unsigned int k=0; k!=2; k++) //xyz
					for (unsigned int l=0; l!=2; l++) //xyz
						rDampMatrix(i*3+k,j*3+l)+=viscosity_matrix(i*2+k,j*2+l);
			}
		}
	}


	template<class T>
	bool VelocityEnrichedPFEM22D::InvertMatrix(const T& input, T& inverse)
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







	  void VelocityEnrichedPFEM22D::CalculateInterfaceNormal(BoundedMatrix<double, 3, 2 >& rPoints, array_1d<double,3>&  rDistances, array_1d<double,2>&  normal, double & interface_area, array_1d<double,3>&  Ninterface, BoundedMatrix<double, 2, 2 >& rInterfacePoints)

        {
                double sign_correction=1.0;
                BoundedMatrix<double, 2, 2 > InterfacePoints;
                array_1d<bool,3>  cut_edges;
                array_1d<double,2>  interface_segment=ZeroVector(2);
                if ((rDistances(0)*rDistances(1))<0.0) cut_edges[0]=true;//edge 12 is cut
                else         cut_edges[0]=false;

                if ((rDistances(1)*rDistances(2))<0.0) cut_edges[1]=true;//edge 23 is cut.
                else         cut_edges[1]=false;

                if ((rDistances(2)*rDistances(0))<0.0) cut_edges[2]=true;//edge 13 is cut.
                else         cut_edges[2]=false;

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
            interface_area=norm;
	        rInterfacePoints(0,0)=InterfacePoints(0,0);
	        rInterfacePoints(0,1)=InterfacePoints(0,1);
	        rInterfacePoints(1,0)=InterfacePoints(1,0);
	        rInterfacePoints(1,1)=InterfacePoints(1,1);

            const double x_interface = 0.5*(InterfacePoints(0,0)+InterfacePoints(1,0));
            const double y_interface = 0.5*(InterfacePoints(0,1)+InterfacePoints(1,1));
            CalculatePosition(rPoints, x_interface, y_interface, 0.0, Ninterface );
        }



//************************************************************************************
  //************************************************************************************
   inline void VelocityEnrichedPFEM22D::CalculatePosition(const bounded_matrix<double, 3, 3 > & coordinates, const double xc, const double yc, const double zc, array_1d<double, 3 > & N)

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
            KRATOS_THROW_ERROR(std::logic_error, "element with zero area found", "");
        }
        else
        {
            inv_area = 1.0 / area;
        }
        N[0] = CalculateVol(x1, y1, x2, y2, xc, yc) * inv_area;
        N[1] = CalculateVol(x2, y2, x0, y0, xc, yc) * inv_area;
        N[2] = CalculateVol(x0, y0, x1, y1, xc, yc) * inv_area;

    }

    //************************************************************************************
    //************************************************************************************
    inline double VelocityEnrichedPFEM22D::CalculateVol(const double x0, const double y0,
                                      const double x1, const double y1,
                                      const double x2, const double y2
                                     )

    {
        return 0.5 * ((x1 - x0)*(y2 - y0)- (y1 - y0)*(x2 - x0));
    }







} // Namespace Kratos
