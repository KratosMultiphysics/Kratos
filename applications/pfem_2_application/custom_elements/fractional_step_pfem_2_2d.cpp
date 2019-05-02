//   
//   Project Name:        Kratos
//   Last modified by:    $Author: it's me! $
//   Date:                $Date: 2008-08-08 23:58:38 $
//   Revision:            $Revision: 1.0 $
//
//

// Project includes
#include "includes/define.h"
#include "custom_elements/fractional_step_pfem_2_2d.h"
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
	FractionalStepPFEM22D::FractionalStepPFEM22D(IndexType NewId, GeometryType::Pointer pGeometry)
		: Element(NewId, pGeometry)
	{
		//DO NOT ADD DOFS HERE!!!
	}

	//************************************************************************************
	//************************************************************************************
	FractionalStepPFEM22D::FractionalStepPFEM22D(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
		: Element(NewId, pGeometry, pProperties)
	{

	}

	Element::Pointer FractionalStepPFEM22D::Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const
	{
		return Element::Pointer(new FractionalStepPFEM22D(NewId, GetGeometry().Create(ThisNodes), pProperties));
	}

	FractionalStepPFEM22D::~FractionalStepPFEM22D()
	{
	}

	//************************************************************************************
	//************************************************************************************


	void FractionalStepPFEM22D::AddExplicitContribution(ProcessInfo& rCurrentProcessInfo)
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
			KRATOS_THROW_ERROR(std::logic_error,"Unexpected value for FRACTIONAL_STEP index: ",rCurrentProcessInfo[FRACTIONAL_STEP]);
		}
		}

		KRATOS_CATCH("");
	}

	void FractionalStepPFEM22D::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
													VectorType& rRightHandSideVector,
													ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY;

		switch ( rCurrentProcessInfo[FRACTIONAL_STEP] )
		{
		case 2:
		{
			//KRATOS_WATCH("calculating an element")
			this->CalculateLocalPressureSystem(rLeftHandSideMatrix,rRightHandSideVector,rCurrentProcessInfo);
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

	void FractionalStepPFEM22D::EquationIdVector(EquationIdVectorType& rResult,
												ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY;

		switch ( rCurrentProcessInfo[FRACTIONAL_STEP] )
		{
		case 2:
		{
			this->PressureEquationIdVector(rResult,rCurrentProcessInfo);
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


	void FractionalStepPFEM22D::GetDofList(DofsVectorType& rElementalDofList,
										  ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY;

		switch ( rCurrentProcessInfo[FRACTIONAL_STEP] )
		{
		case 2:
		{
			this->GetPressureDofList(rElementalDofList,rCurrentProcessInfo);
			break;
		}

		default:
		{
			KRATOS_THROW_ERROR(std::logic_error,"Unexpected value for FRACTIONAL_STEP index: ",rCurrentProcessInfo[FRACTIONAL_STEP]);
		}
		}

		KRATOS_CATCH("");
	}




	//*************************************************************************************
	//*************************************************************************************




	void FractionalStepPFEM22D::CalculateLocalPressureSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY

		double delta_t = rCurrentProcessInfo[DELTA_TIME];
		int iteration_number =  rCurrentProcessInfo[NL_ITERATION_NUMBER];
		const double density_air = rCurrentProcessInfo[DENSITY_AIR];
		const double density_water = rCurrentProcessInfo[DENSITY_WATER];

		const unsigned int TDim = 2;

		//const double viscosity_air = rCurrentProcessInfo[VISCOSITY_AIR];
		//const double viscosity_water = rCurrentProcessInfo[VISCOSITY_WATER];

		double enrich_rhs=0.0;
		double extra_enrich_rhs=0.0;

		//double input_jump_value=0.0;
		//this->SetValue(PRESS_DISCONTINUITY,input_jump_value);
		//double jump_value= -input_jump_value*delta_t;
		//if (iteration_number>1) jump_value=0.0;

		//const bool & split_element = false;//(this->GetValue(SPLIT_ELEMENT));


		//const double mass_factor = 1.0/double(1+TDim);

		const unsigned int LocalSize = GetGeometry().size();
		unsigned int TNumNodes = GetGeometry().size();
		double TNumNodes_double = double(TDim+1);

		array_1d<double,TDim+1> msN; //dimension = number of nodes
		array_1d<double,TDim+1> ms_temp; //dimension = number of nodes

		//BoundedMatrix<double, 3,1 > enrich_rhs;

        // Check sizes and initialize
        if (rRightHandSideVector.size() != LocalSize)
            rRightHandSideVector.resize(LocalSize, false);
        noalias(rRightHandSideVector) = ZeroVector(LocalSize);
        if(rLeftHandSideMatrix.size1() != LocalSize)
			rLeftHandSideMatrix.resize(LocalSize,LocalSize,false);
		noalias(rLeftHandSideMatrix) = ZeroMatrix(LocalSize,LocalSize);

        // Calculate this element's geometric parameters
        double Area;
        array_1d<double, TDim+1> N;
        BoundedMatrix<double, TDim+1, TDim> DN_DX;
        BoundedMatrix<double, TDim+1, (TDim+1)*TDim > D_matrix; //(gradient)
        BoundedMatrix<double, TDim+1, TDim+1 > Laplacian_matrix;
        //BoundedMatrix<double, TNumNodes*2, 2> N_vel_matrix; //useless to calculate it, it's simply the DN_DX matrix multiplied by 1/3, arranged in a different shape.

        Geometry<Node<3> >& geom = this->GetGeometry();

        GeometryUtils::CalculateGeometryData(geom, DN_DX, N, Area);
        //KRATOS_WATCH(Area);

        BoundedMatrix<double, (TDim+1)*TDim , 1 > fract_vel;
        array_1d<double, TDim+1 > old_pressures;


        for(unsigned int iii = 0; iii<TDim+1; iii++)
        {
			array_1d<double,3> velocity = GetGeometry()[iii].FastGetSolutionStepValue(VELOCITY);
			fract_vel(iii*(TDim),0) = velocity[0];
			fract_vel(iii*(TDim)+1,0) = velocity[1];
			//fract_vel(iii*(TDim)+2,0) = velocity[2]; //3D
			old_pressures(iii) = GetGeometry()[iii].FastGetSolutionStepValue(PREVIOUS_ITERATION_PRESSURE);
		}

		//if we want to use stabilization
		//const array_1d<double, 3 > & proj0 = GetGeometry()[0].FastGetSolutionStepValue(PRESS_PROJ);
		//const array_1d<double, 3 > & proj1 = GetGeometry()[1].FastGetSolutionStepValue(PRESS_PROJ);
		//const array_1d<double, 3 > & proj2 = GetGeometry()[2].FastGetSolutionStepValue(PRESS_PROJ);


        //we start by calculating the non-enriched functions.
        for (unsigned int i = 0; i < (TDim+1); i++)
		{
			for (unsigned int j = 0; j < (TDim+1); j++)
			{
				D_matrix(i, (j*TDim) ) =  -DN_DX(i,0)*Area*(1.0/TNumNodes_double);     //integrated by parts. so this is actually matrix -G
				D_matrix(i, (j*TDim+1) ) =  -DN_DX(i,1)*Area*(1.0/TNumNodes_double);
				//D_matrix(i, (j*TDim+2) ) =  -DN_DX(i,2)*Area*(1.0/TNumNodes_double); //3D
			}
		}

        //if (split_element)
		if((geom[0].FastGetSolutionStepValue(DISTANCE)*geom[1].FastGetSolutionStepValue(DISTANCE))<0.0 || (geom[1].FastGetSolutionStepValue(DISTANCE)*geom[2].FastGetSolutionStepValue(DISTANCE))<0.0)
        {
			//if it is the first iteration, we reset the gradient discontinuity:
			if (iteration_number==1)
			{
				this->GetValue(GRADIENT_DISCONTINUITY)=0.0;
				this->GetValue(ENRICH_RHS)=0.0;
			}

			//TOOLS NEEDED TO FIND THE ENRICHMENT SHAPE FUNCTIONS
			//get position of the cut surface
			array_1d<double,(TDim+1)> distances;
			array_1d<double,3*(TDim-1)>  densities(0);
			BoundedMatrix<double,3*(TDim-1), 2> Nenriched;
			array_1d<double,3*(TDim-1)>  volumes(0);
			BoundedMatrix<double,TDim+1, TDim > coords;
			BoundedMatrix<double,3*(TDim-1), TDim+1 > Ngauss;
			array_1d<double,3*(TDim-1)>  signs(0);
			std::vector< Matrix > gauss_gradients(3*(TDim-1));
			//fill coordinates


			//unsigned int single_triangle_node = 1;
			for (unsigned int i = 0; i < TNumNodes; i++)
			{
				const array_1d<double, 3 > & xyz = this->GetGeometry()[i].Coordinates();
				volumes[i] = 0.0;
				distances[i] = this->GetGeometry()[i].FastGetSolutionStepValue(DISTANCE);
				//KRATOS_WATCH(distances(i));
				for (unsigned int j = 0; j < TDim; j++)
				{
					coords(i, j) = xyz[j];
				}
			}

			for (unsigned int i = 0; i < 3*(TDim-1); i++)
				gauss_gradients[i].resize(2, TDim, false);  //2 values of the 2 shape functions, and derivates in (xy) direction).
			unsigned int ndivisions = EnrichmentUtilities::CalculateEnrichedShapeFuncions(coords, DN_DX, distances, volumes, Ngauss, signs, gauss_gradients, Nenriched);

			array_1d<double,3*(TDim-1)>  laplacian_coeffs(0);
			for (unsigned int i = 0; i!=ndivisions; i++) //the 3 partitions of the triangle
			{
				if (signs[i]>0)
				{
					densities(i) = density_air;
				}
				else
				{
					densities(i) = density_water;
				}
				laplacian_coeffs(i) = delta_t/(densities(i));
			}

			Laplacian_matrix =  - prod(DN_DX,trans(DN_DX)) *(volumes(0)*laplacian_coeffs(0)+volumes(1)*laplacian_coeffs(1)+volumes(2)*laplacian_coeffs(2)); //add more in 3d. better a loop!

			//and now the rest of the things needed
			BoundedMatrix<double, 1, 1 > Laplacian_enrich;
			noalias(Laplacian_enrich) = ZeroMatrix(1,1);
			BoundedMatrix<double, 1, 3 > mixed_Laplacian;
			noalias(mixed_Laplacian) = ZeroMatrix(1,TDim+1);


			//To calculate D* =int (N*DN_DX_enrich).  we must use the N at the gauss point of the partition (returned by the enrichment util) and multiply by the area.
			//notice that this first row is only the upper part of the mixed D. we also need the lower one for the jump
			BoundedMatrix<double, 1, (TDim+1)*TDim > D_matrix_mixed;
			noalias(D_matrix_mixed) = ZeroMatrix(1,(TDim+1)*TDim);



			//the matrices we need at the end will be:
			//BoundedMatrix<double, 2, 2 > D_mod;

			for (unsigned int i = 0; i < ndivisions; i++)  //we go through the 3 partitions of the domain
			{
				Laplacian_enrich(0,0) -= (pow(gauss_gradients[i](0,0),2)+pow(gauss_gradients[i](0,1),2))*volumes(i)*laplacian_coeffs(i); //add z direction in 3d
				//and the 'mixed laplacians' (standard shape functions * enrichments)
				for (unsigned int j = 0; j < TDim+1; j++) //we go through the 3 standard shape functions
				{
					mixed_Laplacian(0,j) -= DN_DX(j,0)*(gauss_gradients[i](0,0)*volumes(i)*laplacian_coeffs(i))+ DN_DX(j,1)*(gauss_gradients[i](0,1)*volumes(i)*laplacian_coeffs(i));
					//and also the D matrixes

					D_matrix_mixed(0,j*TDim) -= gauss_gradients[i](0,0) *volumes(i)*Ngauss(i,j); //integrated by parts!
					D_matrix_mixed(0,j*TDim+1) -= gauss_gradients[i](0,1) *volumes(i)*Ngauss(i,j);
					//D_matrix_mixed(0,j*2+2) -= gauss_gradients[i](0,2) *volumes(i)*Ngauss(i,j); //3D

					//D_matrix_mixed(0,j*2) += DN_DX(j,0) * volumes(i) * Nenriched(i,0); //D matrix
					//D_matrix_mixed(0,j*2+1) += DN_DX(j,1) * volumes(i) * Nenriched(i,0); //D matrix
				}


			}

			//the 'jump' degree of freedom will not be condensed since it's imposed
			//therefore we'll only have to calculate the inverse of the reduced matrix of the gradient jump enrich function. which is a 1x1 matrix (inverse of the 1,1 position of the enriched laplacian)
			const double inv_Laplacian_enrich_weighted=1.0/Laplacian_enrich(0,0);
			this->SetValue(INV_LAPLACIAN_ENRICH,(inv_Laplacian_enrich_weighted));///(delta_t+TauOne)));


			//since we're imposing this DOF, we'll do:
			//for (unsigned int i = 0; i < 3; i++)
			//	rRightHandSideVector(i) -= mixed_Laplacian_jump(0,i)*jump_value;

			//when there's no jump there's nothing in the RHS to condense. but when we add it, it is equal to -Laplacian_enrich*jump_value due to the elimination of this DOF.
			//and just like the LHS, we must condense it and substrac it (actually it's (-)*(-)=+ )
			//enrich_rhs = -Laplacian_enrich(0,1)*jump_value;

			for (unsigned int i = 0; i < (TDim+1)*TDim; i++)
				extra_enrich_rhs +=D_matrix_mixed(0,i)*fract_vel(i,0) ;


			///gamma=1;
			//for (unsigned int i = 0; i < 3; i++)
			//	enrich_rhs += mixed_Laplacian(0,i) *old_pressures(i);

			for (unsigned int i = 0; i < TDim+1; i++)
				rRightHandSideVector(i) -= mixed_Laplacian(0,i) *inv_Laplacian_enrich_weighted * (enrich_rhs + extra_enrich_rhs);



			for (unsigned int i = 0; i < TDim+1; i++)
			{
				(this->GetValue(ENRICH_LHS_ROW))(i)=mixed_Laplacian(0,i);//*(delta_t+TauOne);
				//KRATOS_WATCH(mixed_Laplacian(0,i));
			}

			//now we must condense these matrixes into the original system (D_matrix, Laplacian_matrix)

			///D_matrix-= inv_Laplacian_enrich_weighted * prod(trans(mixed_Laplacian),D_matrix_mixed); //actually inv_Lap.. has 1/delta_t and mix_Lap has delta_t so at the end we do not write it. but do not forget that they are part of the system. to condense we must use the lefthandside, which is actually Laplacian*delta_t
			///no need to do this ^ cos it is in the rhs, done when added : rRightHandSideVector(i) -= mixed_Laplacian(0,i) *inv_Laplacian_enrich_weighted * (enrich_rhs + extra_enrich_rhs);


			///gamma=1!
			//rRightHandSideVector += prod(Laplacian_matrix,old_pressures)*delta_t;

			Laplacian_matrix -= inv_Laplacian_enrich_weighted * prod(trans(mixed_Laplacian),mixed_Laplacian);

		}

		else //uncut(normal, not interfase) element
        {

			double laplacian_coefficient = 1.0;
			if ((geom[0].FastGetSolutionStepValue(DISTANCE)+geom[1].FastGetSolutionStepValue(DISTANCE)+geom[2].FastGetSolutionStepValue(DISTANCE))<0.0) //that is, water.
			{
				laplacian_coefficient = delta_t/density_water;
			}
			else
			{
				laplacian_coefficient = delta_t/density_air;
			}
			Laplacian_matrix =  -prod(DN_DX,trans(DN_DX))*Area*laplacian_coefficient;  // B B^T  (standard laplacian, now we must add the contribution of the condensed new nodes.
		}


        //done, now we save it into the LHS
        noalias(rLeftHandSideMatrix) = Laplacian_matrix;//*(delta_t+TauOne);///gamma=1
        BoundedMatrix<double, TDim+1,1 > temp_rhs;
        noalias(temp_rhs) = prod(D_matrix,fract_vel);
        //KRATOS_WATCH(temp_rhs);
        rRightHandSideVector(0)  += temp_rhs(0,0);
        rRightHandSideVector(1)  += temp_rhs(1,0);
        rRightHandSideVector(2)  += temp_rhs(2,0);
        //rRightHandSideVector(3)  += temp_rhs(2,0);


        //stabilization contribution and previous pressure:
		if (rCurrentProcessInfo[NL_ITERATION_NUMBER]>1)
		{

			//const array_1d<double, 3 > & proj0 = GetGeometry()[0].FastGetSolutionStepValue(PRESS_PROJ);
			//const array_1d<double, 3 > & proj1 = GetGeometry()[1].FastGetSolutionStepValue(PRESS_PROJ);
			//const array_1d<double, 3 > & proj2 = GetGeometry()[2].FastGetSolutionStepValue(PRESS_PROJ);
			//array_1d<double, 2 > vel_gauss;
			//vel_gauss[0] = N[0] * proj0[0] + N[1] * proj1[0] + N[2] * proj2[0];
			//vel_gauss[1] = N[0] * proj0[1] + N[1] * proj1[1] + N[2] * proj2[1];
			//vel_gauss *= TauOne*Area;
			//noalias(rRightHandSideVector) -= prod(DN_DX, vel_gauss);
			//noalias(rLeftHandSideMatrix) += TauOne/delta_t*Laplacian_matrix;

			rRightHandSideVector += prod(Laplacian_matrix,old_pressures);
		}
		//warning! this should be on only for iteration_number>1!!
		//stabilization contribution:


        enrich_rhs+=extra_enrich_rhs;
        //to recover the gradient jump later, we must save the RHS of the enrichment, that is:
		this->SetValue(ENRICH_RHS,enrich_rhs);


        ///if it is an inactive element we must turn it off. But to avoid problems in case we've turned off all the elements that contribute to a node, we simply do:
		//if (this->GetValue(IS_INACTIVE)==true) //it means we must not add it!
		//{
		//	rLeftHandSideMatrix*=1.0e-6;
		//	rRightHandSideVector*=1.0e-6;
		//}


		//subtracting the dirichlet term
		// RHS -= LHS*DUMMY_UNKNOWNs
		for(unsigned int iii = 0; iii<LocalSize; iii++)
			ms_temp[iii] = GetGeometry()[iii].FastGetSolutionStepValue(PRESSURE);
		noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix,ms_temp);

		KRATOS_CATCH("");
	}


	//************************************************************************************
	//************************************************************************************


		void FractionalStepPFEM22D::CalculateViscousRHS(ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY

					const double density_air = rCurrentProcessInfo[DENSITY_AIR];
					const double density_water = rCurrentProcessInfo[DENSITY_WATER];

					const double viscosity_air = rCurrentProcessInfo[VISCOSITY_AIR];
					const double viscosity_water = rCurrentProcessInfo[VISCOSITY_WATER];

					const unsigned int TDim=2;

					double delta_t = rCurrentProcessInfo[DELTA_TIME];
					//array_1d<double,3> & gravity= CurrentProcessInfo[GRAVITY];
					//double x_force  = gravity(0);
					//double y_force  = gravity(1);
					const array_1d<double,3> zero3 = ZeroVector(3);
					const double nodal_weight = 1.0/ (1.0 + double (TDim) );

					BoundedMatrix<double, (TDim+1), 3 > fluid_velocities = ZeroMatrix((TDim+1), 3);
					BoundedMatrix<double, (TDim+1), TDim > DN_DX;
					array_1d<double, (TDim+1) > N;
					double Area;
					Geometry<Node<3> >& geom = this->GetGeometry();
					GeometryUtils::CalculateGeometryData(geom, DN_DX, N, Area);

					double density =1.0;
					double lhs_coeff = 1.0;
					double rhs_fluid_coeff = 1.0;

					if((geom[0].FastGetSolutionStepValue(DISTANCE)*geom[1].FastGetSolutionStepValue(DISTANCE))<0.0 || (geom[1].FastGetSolutionStepValue(DISTANCE)*geom[2].FastGetSolutionStepValue(DISTANCE))<0.0) //check also 3D
					{
						//TOOLS NEEDED TO FIND THE ENRICHMENT SHAPE FUNCTIONS
						//get position of the cut surface
						//const unsigned int LocalSize = geom.size()*2;
						unsigned int TNumNodes = geom.size();
						array_1d<double,(TDim+1)> distances;
						BoundedMatrix<double,3*(TDim-1), 2> Nenriched;
						array_1d<double,(3*(TDim-1))> volumes;
						array_1d<double,(3*(TDim-1))> densities;
						BoundedMatrix<double,(TDim+1), TDim > coords;
						BoundedMatrix<double, 3*(TDim-1), (TDim+1) > Ngauss;
						array_1d<double,(3*(TDim-1))> signs;
						std::vector< Matrix > gauss_gradients(3*(TDim-1));

						//fill coordinates
						for (unsigned int i = 0; i < TNumNodes; i++)
						{
							const array_1d<double, 3 > & xyz = geom[i].Coordinates();
							const array_1d<double, 3 > & velocity = geom[i].FastGetSolutionStepValue(VELOCITY);

							volumes[i] = 0.0;
							distances[i] = geom[i].FastGetSolutionStepValue(DISTANCE);
							for (unsigned int j = 0; j < (TDim); j++)
							{
								coords(i, j) = xyz[j];
								fluid_velocities(i,j) = velocity[j];
							}
						}

						for (unsigned int i = 0; i < 3*(TDim-1) ; i++)
							gauss_gradients[i].resize(2, TDim, false);  //2 values of the 2 shape functions, and derivates in xy(z) direction).
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
								densities(i)= density_air;
							}
							else
							{
								densities(i)= density_water;

							}
						}

						double Weight = Area * viscosity_air;
						double maximum_Weight = Area * 0.2 * density_air/delta_t * Area; //the last Area is in fact h^2 . Supposing h=sqrt(Area) we get this result.
						                                                       //this way we are ensuring MAXIMUM DIFFUSITIVY in the fluid, using Fo=0.2.
						if (Weight>maximum_Weight)
							Weight = maximum_Weight;

					    BoundedMatrix<double, TDim*(TDim+1), TDim*(TDim+1) > Viscosity_matrix = ZeroMatrix(TDim*(TDim+1), TDim*(TDim+1));
                        this->AddViscousTerm(Viscosity_matrix, DN_DX,
										Weight);

						double lhs_node =0.0;
						double rhs_fluid_node =0.0;
						double rhs_bodyforce_node =0.0;


						for (unsigned int j=0; j!=(TDim+1); j++)
						{

							geom[j].SetLock();

							geom[j].FastGetSolutionStepValue(NODAL_AREA) += Area * nodal_weight;

							array_1d<double, 3 > & current_rhs = geom[j].FastGetSolutionStepValue(RHS);

							lhs_node =0.0;
							rhs_fluid_node =0.0;
						    rhs_bodyforce_node =0.0;

							for (unsigned int i=0;i!=ndivisions;i++)
							{
								rhs_fluid_node     += volumes(i)*Ngauss(i,j)*densities(i)/delta_t;
								lhs_node           += volumes(i)*Ngauss(i,j)*densities(i)/delta_t;
								rhs_bodyforce_node += volumes(i)*Ngauss(i,j)*densities(i);


							}

							array_1d<double, 3 > & bodyforce = geom[j].FastGetSolutionStepValue(BODY_FORCE);

							geom[j].FastGetSolutionStepValue(NODAL_MASS) += lhs_node;
							for (unsigned int k=0;k!=(TDim);k++) //k component of local stress
								current_rhs[k] += (rhs_fluid_node*fluid_velocities(j,k)+bodyforce(k)*rhs_bodyforce_node);

							for (unsigned int i=0;i!=(TDim+1);i++) //neighbour node
									for (unsigned int k=0;k!=(TDim);k++) //k component of local stress
										for (unsigned int l=0;l!=(TDim);l++) //affected by l component in i neighbour node
											current_rhs[k] += Viscosity_matrix(j*TDim+k,i*TDim+l)*fluid_velocities(i,l);

							 geom[j].UnSetLock();
						}
					}

					else
					{
						double viscosity=1.0;

						if ((geom[0].FastGetSolutionStepValue(DISTANCE)+geom[1].FastGetSolutionStepValue(DISTANCE)+geom[2].FastGetSolutionStepValue(DISTANCE))<0.0)
						{

								lhs_coeff = density_water/delta_t;
								rhs_fluid_coeff = density_water/delta_t;
								density = density_water;
								viscosity = viscosity_water;
						}
						else
						{
								lhs_coeff = density_air/delta_t;
								rhs_fluid_coeff = density_air/delta_t;
								density = density_air;
								viscosity = viscosity_air;
						}



						for (unsigned int i = 0; i < (TDim+1); i++)
						{
							const array_1d<double, 3 > & velocity = geom[i].FastGetSolutionStepValue(VELOCITY);
							for (unsigned int j = 0; j < (TDim); j++)
							{
								fluid_velocities(i,j) = velocity[j];
							}

						}

						double Weight = Area * viscosity;
						double maximum_Weight = Area * 0.2 * rhs_fluid_coeff * Area; //the last Area is in fact h^2 . Supposing h=sqrt(Area) we get this result.
						                                                       //this way we are ensuring MAXIMUM DIFFUSITIVY in the fluid, using Fo=0.2.
						if (Weight>maximum_Weight)
							Weight = maximum_Weight;


					    BoundedMatrix<double, TDim*(TDim+1), TDim*(TDim+1) > Viscosity_matrix = ZeroMatrix(TDim*(TDim+1), TDim*(TDim+1));
                        this->AddViscousTerm(Viscosity_matrix, DN_DX,
										Weight);

						for (unsigned int j=0; j!=(TDim+1); j++)
						{

							    geom[j].SetLock();

							    array_1d<double, 3 > & bodyforce = geom[j].FastGetSolutionStepValue(BODY_FORCE);

								geom[j].FastGetSolutionStepValue(NODAL_AREA) += Area*nodal_weight;

								array_1d<double, 3 > & current_rhs = geom[j].FastGetSolutionStepValue(RHS);
								for (unsigned int k=0;k!=(TDim);k++) //k component of local stress
										current_rhs[k] += (rhs_fluid_coeff*fluid_velocities(j,k)+bodyforce(k)*density)*nodal_weight*Area;

								for (unsigned int i=0;i!=(TDim+1);i++) //neighbour node
									for (unsigned int k=0;k!=(TDim);k++) //k component of local stress
										for (unsigned int l=0;l!=(TDim);l++) //affected by l component in i neighbour node
											current_rhs[k] += Viscosity_matrix(j*TDim+k,i*TDim+l)*fluid_velocities(i,l);

								geom[j].FastGetSolutionStepValue(NODAL_MASS) += Area*lhs_coeff*nodal_weight;

							    geom[j].UnSetLock();
						}
					}


		KRATOS_CATCH("");
	}

	//*****************************************************************************
	//***************************************************************************


	void FractionalStepPFEM22D::CalculatePressureProjection(ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY
		            const unsigned int TDim=2;

					const double delta_t = rCurrentProcessInfo[DELTA_TIME];
					const double density_air = rCurrentProcessInfo[DENSITY_AIR]; // * (1.0-theta) ;
					const double density_water = rCurrentProcessInfo[DENSITY_WATER]; // * (1.0-theta) ;

					const array_1d<double,3> zero3 = ZeroVector(3);
					const double mass_factor = 1.0/ (1.0 + double (TDim) );

					BoundedMatrix<double, (TDim+1), TDim > DN_DX;
					array_1d<double, (TDim+1) > N;
					double Area;
					Geometry<Node<3> >& geom = this->GetGeometry();
					//const unsigned int LocalSize = geom.size()*2;
					unsigned int TNumNodes = geom.size();
					//BoundedMatrix<double, (2+1),1 > enrich_rhs;
					GeometryUtils::CalculateGeometryData(geom, DN_DX, N, Area);

					array_1d<double,(TDim+1)> nodal_masses = ZeroVector(TDim+1);

					array_1d<double,(TDim+1)> pressure;
					for (unsigned int i=0; i<(TDim+1);i++)
					{
						pressure(i) = geom[i].FastGetSolutionStepValue(PRESSURE);
					}


					BoundedMatrix<double, (TDim+1), TDim*(TDim+1) > G_matrix_no_ro; //(gradient)
					noalias(G_matrix_no_ro) = ZeroMatrix((TDim+1), TDim*(TDim+1));
					//new (upper) part of D, corresponding the first enrichment function (gradient discontinuity
					BoundedMatrix<double, 1, TDim*(TDim+1) > G_matrix_mixed_no_ro;
					noalias(G_matrix_mixed_no_ro) = ZeroMatrix(1,TDim*(TDim+1));
					//lower part of D, the one corresponding to the second enrichmend function (jump)
					BoundedMatrix<double, 1, TDim*(TDim+1)> G_matrix_mixed_jump_no_ro;
					noalias(G_matrix_mixed_jump_no_ro) = ZeroMatrix(1,TDim*(TDim+1));


					//for the enrichment:
					array_1d<double,(TDim+1)> distances;
					BoundedMatrix<double,3*(TDim-1), 2> Nenriched;
					array_1d<double,(3*(TDim-1))> volumes;
					array_1d<double,(3*(TDim-1))> densities;
					BoundedMatrix<double,(TDim+1), TDim > coords;
					BoundedMatrix<double, 3*(TDim-1), (TDim+1) > Ngauss;
					array_1d<double,(3*(TDim-1))> signs;
					std::vector< Matrix > gauss_gradients((TDim-1)*3);


					if((geom[0].FastGetSolutionStepValue(DISTANCE)*geom[1].FastGetSolutionStepValue(DISTANCE))<0.0 || (geom[1].FastGetSolutionStepValue(DISTANCE)*geom[2].FastGetSolutionStepValue(DISTANCE))<0.0) //do not forget the 3d
					{
						//to begin with we calculate the gradient discontinuity:
						const int iteration_number =  rCurrentProcessInfo[NL_ITERATION_NUMBER];
						if (true)
						{

							double gradient_discontinuity = this->GetValue(ENRICH_RHS) ;
							double gradient_discontinuity_complete_rhs = this->GetValue(GRADIENT_DISCONTINUITY)/((this->GetValue(INV_LAPLACIAN_ENRICH)));

							array_1d<double,(TDim+1)> enrich_lhs;
							enrich_lhs = this->GetValue(ENRICH_LHS_ROW);

							array_1d<double,(TDim+1)> pressure_array;
							array_1d<double,(TDim+1)> previous_iteration_pressure_array;
							for (unsigned int i=0; i<(TDim+1);i++)
							{
								pressure_array(i) = geom[i].FastGetSolutionStepValue(PRESSURE);//-geom[i].GetSolutionStepValue(PRESSURE,1);
								previous_iteration_pressure_array(i) = geom[i].FastGetSolutionStepValue(PREVIOUS_ITERATION_PRESSURE);
							}
							for (unsigned int i = 0; i < (TDim+1); i++)  // node i
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



						//double input_jump_value = this->GetValue(PRESS_DISCONTINUITY);
						//double jump_value=-input_jump_value;
						//if (iteration_number>1) jump_value *=double(iteration_number);
						double gradient_discontinuity; // = this->GetValue(ENRICH_RHS);
						//array_1d<double,(2+1)> & enrich_lhs = this->GetValue(ENRICH_LHS_ROW);

						noalias(G_matrix_no_ro) = ZeroMatrix((TDim+1),(TDim+1)*TDim);
						noalias(G_matrix_mixed_no_ro) = ZeroMatrix(1,(TDim+1)*TDim);
						noalias(G_matrix_mixed_jump_no_ro) = ZeroMatrix(1,(TDim+1)*TDim);



						//TOOLS NEEDED TO FIND THE ENRICHMENT SHAPE FUNCTIONS
						//get position of the cut surface

						for (unsigned int i = 0; i < TNumNodes; i++)
						{
							const array_1d<double, 3 > & xyz = geom[i].Coordinates();
							distances[i] = geom[i].FastGetSolutionStepValue(DISTANCE);
							//KRATOS_WATCH(distances(i));
							for (unsigned int j = 0; j < (TDim); j++)
								coords(i, j) = xyz[j];
						}

						for (unsigned int i = 0; i < ((TDim-1)*3) ; i++)
							gauss_gradients[i].resize(2, TDim, false);  //2 values of the 2 shape functions, and derivates in (xy) direction).
						unsigned int ndivisions;
						ndivisions = EnrichmentUtilities::CalculateEnrichedShapeFuncions(coords, DN_DX, distances, volumes, Ngauss, signs, gauss_gradients, Nenriched); //, face_gauss_N, face_gauss_Nenriched, face_Area, face_n  ,single_triangle_node);


					   	gradient_discontinuity = this->GetValue(GRADIENT_DISCONTINUITY);

						for (unsigned int i = 0; i < ndivisions; i++)
						{
							//geom[i].FastGetSolutionStepValue(DISTANCE)=distances[i];

							if (signs[i]>0)
							{
								densities(i) = density_air;
							}
							else
							{
								densities(i) = density_water;
							}
						}

						///*****************************************************************************************
						for (unsigned int i = 0; i < ndivisions; i++)  // partiton i
						{
							for (unsigned int j = 0; j < (TDim+1); j++) //shape function (Nj)
							{


								for (unsigned int k=0;k!=TDim;k++) //x,y,(z)
								{
									G_matrix_mixed_no_ro(0,j*(TDim)+k) += gauss_gradients[i](0,k) *volumes(i)*Ngauss(i,j);

								}

								nodal_masses(j) += volumes(i)*Ngauss(i,j)*densities(i)*(1.0/delta_t);
							}
						}
						//for the massless g matrix (to calculate the press proj for the stabilization) we do not need to loop the partitions, so we create a new loop
						for (unsigned int j = 0; j < (TDim+1); j++) //shape function (Nj)
						{
							for (unsigned int k=0;k!=TDim;k++) //x,y,(z)
							{
								for (unsigned int l = 0; l < (TDim+1); l++) //shape function derivative (dNl/dk)
								{
									G_matrix_no_ro(l, (j*(TDim)+k) ) = DN_DX(l,k)*Area*mass_factor;
								}
							}
						}
						//now we save the data:
						for (unsigned int i=0; i!=(TDim+1); i++) //loop around the nodes of the element to add contribution to node i
						{

							geom[i].SetLock();

							array_1d<double, 3 > & current_press_proj_no_ro = geom[i].FastGetSolutionStepValue(PRESS_PROJ_NO_RO);
							//array_1d<double, 3 > & current_press_proj = geom[i].FastGetSolutionStepValue(PRESS_PROJ);

							for (unsigned int j = 0; j < (TDim+1) ; j++) //loop around the nodes of the element to add contribution of node j to node i
							{

								for (unsigned int k = 0; k < (TDim) ; k++) //x,y,(z)
								{
									//current_press_proj[k] += reduction_factor*G_matrix(j,i*(2)+k)*(pressure(j));///-old_pressures(j)); //gamma=0!
									current_press_proj_no_ro[k] += G_matrix_no_ro(j,i*(TDim)+k)*(pressure(j));///-old_pressures(j)); //gamma=0!
								}
							}

							for (unsigned int k = 0; k < (TDim) ; k++) //x,y,(z)
							{
								current_press_proj_no_ro[k] += G_matrix_mixed_no_ro(0,TDim*i+k)*gradient_discontinuity;// + G_matrix_mixed_jump_no_ro(0,2*i+k)*jump_value;
							}

							//geom[i].FastGetSolutionStepValue(NODAL_AREA) += Area*mass_factor*reduction_factor;
							geom[i].FastGetSolutionStepValue(NODAL_MASS) += nodal_masses(i);

							geom[i].UnSetLock();

						}

					}
					else //normal (uncut) element
					{
						noalias(G_matrix_no_ro) = ZeroMatrix((TDim+1),(TDim+1)*TDim);
						//noalias(G_matrix) = ZeroMatrix((2+1),(2-1)*6);

						//double density =1.0;
						double lhs_coeff = 1.0;
						//double rhs_coeff = 1.0;

						if ((geom[0].FastGetSolutionStepValue(DISTANCE)+geom[1].FastGetSolutionStepValue(DISTANCE)+geom[2].FastGetSolutionStepValue(DISTANCE))<0.0) //that is, water. do not forget the 3d
						{

							lhs_coeff = density_water/delta_t;
						}
						else
						{
							lhs_coeff = density_air/delta_t;
						}



						for (unsigned int i = 0; i < (TDim+1); i++) //loop in shape functions (row)
						{
							for (unsigned int j = 0; j < (TDim+1) ; j++) //loop in shape functions (column)
							{
								for (unsigned int k = 0; k < (TDim) ; k++) //x,y,(z)
								{
									//G_matrix(i, (j*2)+k ) = DN_DX(i,k)*Area*mass_factor;
									G_matrix_no_ro(i, (j*TDim)+k ) = DN_DX(i,k)*Area*mass_factor;
								}
							}
						}

						//G_matrix /=density;

						//now we save the data:
						for (unsigned int i=0; i!=(TDim+1); i++) //loop around the nodes of the element to add contribution to node i
						{
							geom[i].SetLock();

							array_1d<double, 3 > & current_press_proj_no_ro = geom[i].FastGetSolutionStepValue(PRESS_PROJ_NO_RO);
							//array_1d<double, 3 > & current_press_proj = geom[i].FastGetSolutionStepValue(PRESS_PROJ);

							for (unsigned int j = 0; j < (TDim+1) ; j++) //loop around the nodes of the element to add contribution of node j to node i
							{

								for (unsigned int k = 0; k < (TDim) ; k++) //x,y,(z)
								{
									//current_press_proj[k] += G_matrix(j,i*(2)+k)*(pressure(j));///-old_pressures(j)); //gamma=0!
									current_press_proj_no_ro[k] += G_matrix_no_ro(j,i*(TDim)+k)*(pressure(j));///-old_pressures(j)); //gamma=0!
								}
							}

							//geom[i].FastGetSolutionStepValue(NODAL_AREA) += Area*mass_factor;
							geom[i].FastGetSolutionStepValue(NODAL_MASS) += Area*mass_factor*lhs_coeff; //lhs coeff already includes the delta_t!!!!!!!!!!

							geom[i].UnSetLock();

						}

					} //closing the normal (uncut) element


		KRATOS_CATCH("");
	}

	//*****************************************************************************
	//***************************************************************************



	//************************************************************************************
	//************************************************************************************
	// this subroutine calculates the nodal contributions for the explicit steps of the
	// fractional step procedure
	void FractionalStepPFEM22D::InitializeSolutionStep(ProcessInfo& CurrentProcessInfo)
	{
		KRATOS_TRY
		///WARNING!!!! not calculating the pressure projection since i'm supposing it's being done before!
		KRATOS_CATCH("");
	}

	//************************************************************************************
	//************************************************************************************


	//************************************************************************************
	//************************************************************************************
	void FractionalStepPFEM22D::PressureEquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo)
	{
		unsigned int number_of_nodes = GetGeometry().PointsNumber();
		if(rResult.size() != number_of_nodes)
			rResult.resize(number_of_nodes,false);

		for (unsigned int i=0;i<number_of_nodes;i++)
				rResult[i] = GetGeometry()[i].GetDof(PRESSURE).EquationId();
	}


	//************************************************************************************
	//************************************************************************************
	void FractionalStepPFEM22D::GetPressureDofList(DofsVectorType& ElementalDofList,ProcessInfo& CurrentProcessInfo)
	{
		unsigned int number_of_nodes = GetGeometry().PointsNumber();
		if(ElementalDofList.size() != number_of_nodes)
			ElementalDofList.resize(number_of_nodes);

		for (unsigned int i=0;i<number_of_nodes;i++)
			ElementalDofList[i] = GetGeometry()[i].pGetDof(PRESSURE);

	}


	//************************************************************************************
	//************************************************************************************


	void FractionalStepPFEM22D::AddViscousTerm(BoundedMatrix<double, (2-1)*6, (2-1)*6 >& rDampMatrix,
                         BoundedMatrix<double, (2+1), 2 >& rShapeDeriv,
                         const double Weight)
	{

        const unsigned int TDim=2;

		BoundedMatrix<double, (TDim+1)*TDim,(TDim-1)*3 > B_matrix = ZeroMatrix((TDim+1)*TDim,(TDim-1)*3);
		for (unsigned int i=0; i!=(TDim+1); i++) //i node
		{
			for (unsigned int j=0; j!=(TDim); j++) //x,y,z
				B_matrix(i*(TDim)+j,j)=rShapeDeriv(i,j);

			//useful for both 2d and 3d:
			//relating 12 and 21 stresses
			B_matrix(i*(TDim)+0,2)=rShapeDeriv(i,1);
			B_matrix(i*(TDim)+1,2)=rShapeDeriv(i,0);


		}

		int counter=0;
		BoundedMatrix<double, (TDim-1)*3, (TDim-1)*3 > C_matrix = ZeroMatrix((TDim-1)*3,(TDim-1)*3);

		for (unsigned int i=0; i!=(TDim); i++)
		{
			C_matrix(counter,counter)=2.0;
			counter++;
		}
		for (unsigned int i=0; i!=((TDim-2)*2+1); i++)
		{
			C_matrix(counter,counter)=1.0;
			counter++;
		}

		C_matrix*= -Weight;

		BoundedMatrix<double, (2-1)*3 , (2-1)*6  > temp_matrix = prod(C_matrix,trans(B_matrix));
		rDampMatrix = prod(B_matrix, temp_matrix );

	}



template<class T>
bool FractionalStepPFEM22D::InvertMatrix(const T& input, T& inverse)
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


   //TO PRINT ELEMENTAL VARIABLES (ONLY ONE GAUSS POINT PER ELEMENT)
    void FractionalStepPFEM22D::GetValueOnIntegrationPoints(const Variable<double>& rVariable,
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


	  void FractionalStepPFEM22D::GetValueOnIntegrationPoints(const Variable<array_1d<double, 3 > >& rVariable,
            std::vector<array_1d<double, 3 > >& rValues,
            const ProcessInfo& rCurrentProcessInfo)
	{
		if (rVariable == PRESS_PROJ_NO_RO)
		{
			// Set output vector (for a single integration point)
			rValues.resize(1);
			rValues[0]=this->GetValue(PRESS_PROJ_NO_RO);
		}
		else // Default behaviour (returns elemental data)
		{
			rValues.resize(1);
			//const VMS<Dim,NumNodes>* const_this = static_cast< const VMS<Dim,NumNodes>* >(this);
			//rOutput[0] = const_this->GetValue(rVariable);
		}

	}



} // Namespace Kratos
