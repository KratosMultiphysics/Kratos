//   
//   Project Name:        Kratos       
//   Last modified by:    $Author: it's me! $
//   Date:                $Date: 2008-08-08 23:58:38 $
//   Revision:            $Revision: 1.0 $
//
//

// Project includes 
#include "includes/define.h"
#include "custom_elements/fluid_phase_2d.h"
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
	FluidPhasePFEM22D::FluidPhasePFEM22D(IndexType NewId, GeometryType::Pointer pGeometry)
		: Element(NewId, pGeometry)
	{		
		//DO NOT ADD DOFS HERE!!!
	}

	//************************************************************************************
	//************************************************************************************
	FluidPhasePFEM22D::FluidPhasePFEM22D(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
		: Element(NewId, pGeometry, pProperties)
	{

	}

	Element::Pointer FluidPhasePFEM22D::Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const
	{
		return Element::Pointer(new FluidPhasePFEM22D(NewId, GetGeometry().Create(ThisNodes), pProperties));
	}

	FluidPhasePFEM22D::~FluidPhasePFEM22D()
	{
	}

	//************************************************************************************
	//************************************************************************************
	
	
	void FluidPhasePFEM22D::AddExplicitContribution(ProcessInfo& rCurrentProcessInfo)
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
	
	void FluidPhasePFEM22D::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
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
		case 10:
		{
			//KRATOS_WATCH("calculating an element")
			this->CalculateSolidLocalSystem(rLeftHandSideMatrix,rRightHandSideVector,rCurrentProcessInfo);
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
	
	void FluidPhasePFEM22D::EquationIdVector(EquationIdVectorType& rResult,
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
		case 10:
		{
			this->SolidEquationIdVector(rResult,rCurrentProcessInfo);
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
	

	void FluidPhasePFEM22D::GetDofList(DofsVectorType& rElementalDofList,
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
		case 10:
		{
			this->GetSolidDofList(rElementalDofList,rCurrentProcessInfo);
			break;
		}
		default:
		{
			KRATOS_ERROR(std::logic_error,"Unexpected value for FRACTIONAL_STEP index: ",rCurrentProcessInfo[FRACTIONAL_STEP]);
		}
		}

		KRATOS_CATCH("");
	}
	
	
	
	
	//*************************************************************************************
	//*************************************************************************************
	
	
	

	void FluidPhasePFEM22D::CalculateLocalPressureSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY

		//rCurrentProcessInfo[DELTA_TIME]=1.0;
		//ProcessInfo& CurrentProcessInfo = mr_model_part.GetProcessInfo();
		double Density = 0.0 ; //used for stabilization;
		double delta_t = rCurrentProcessInfo[DELTA_TIME];
		int iteration_number =  rCurrentProcessInfo[NL_ITERATION_NUMBER];
		double mDENSITY_AIR = rCurrentProcessInfo[DENSITY_AIR]; // * (1.0-theta) ;
		double mDENSITY_WATER = rCurrentProcessInfo[DENSITY_WATER]; // * (1.0-theta) ;
		
		const double porosity=rCurrentProcessInfo[POROSITY];
		const double diameter=rCurrentProcessInfo[DIAMETER];
		const double water_visc=0.0001;
		const double water_dens=1000.0;
		double darcy_coeff = 150.0*pow((1.0-porosity),2)/(pow(porosity,3))*water_visc/(pow(diameter,2));//rCurrentProcessInfo[LIN_DARCY_COEF];
		double nonlin_darcy_coeff = 1.75*(1.0-porosity)/(pow(porosity,3))*water_dens/(diameter); ;//rCurrentProcessInfo[NONLIN_DARCY_COEF];

		double enrich_rhs=0.0;
		double extra_enrich_rhs=0.0;

		//double input_jump_value=0.0;
		//this->SetValue(PRESS_DISCONTINUITY,input_jump_value);
		//double jump_value= -input_jump_value*delta_t;
		//if (iteration_number>1) jump_value=0.0;
		
		//const bool & split_element = false;//(this->GetValue(SPLIT_ELEMENT));
		
		//boost::numeric::ublas::bounded_matrix<double,3,2> msDN_DX;
		//boost::numeric::ublas::bounded_matrix<double,2,2> msD;
		const double one_third = 1.0/3.0;
  		array_1d<double,3> msN; //dimension = number of nodes
		array_1d<double,3> ms_temp; //dimension = number of nodes

		const unsigned int LocalSize = GetGeometry().size();
		unsigned int TNumNodes = GetGeometry().size();
		//boost::numeric::ublas::bounded_matrix<double, 3,1 > enrich_rhs;

		//getting data for the given geometry
//        const unsigned int LocalSize = (2 + 1) * TNumNodes;

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
        boost::numeric::ublas::bounded_matrix<double, 3, 6 > D_matrix; //(gradient)
        boost::numeric::ublas::bounded_matrix<double, 3, 3 > Laplacian_matrix;
        //boost::numeric::ublas::bounded_matrix<double, TNumNodes*2, 2> N_vel_matrix; //useless to calculate it, it's simply the DN_DX matrix multiplied by 1/3, arranged in a different shape.
       
        Geometry<Node<3> >& geom = this->GetGeometry();
         
        GeometryUtils::CalculateGeometryData(geom, DN_DX, N, Area);
        //KRATOS_WATCH(Area);

        boost::numeric::ublas::bounded_matrix<double, 6, 1 > fract_vel;
        array_1d<double, 3 > old_pressures;
        array_1d<double, 3 > vectorial_mean_vel=ZeroVector(3);
        
        
        
        for(unsigned int iii = 0; iii<3; iii++)
        {
			array_1d<double,3> velocity = GetGeometry()[iii].FastGetSolutionStepValue(WATER_VELOCITY);
			array_1d<double,3> solid_velocity = GetGeometry()[iii].FastGetSolutionStepValue(VELOCITY);
			fract_vel(iii*2,0) = velocity[0];
			fract_vel(iii*2+1,0) = velocity[1];
			//old_pressures(iii) = GetGeometry()[iii].GetSolutionStepValue(PRESSURE,1);
			old_pressures(iii) = GetGeometry()[iii].FastGetSolutionStepValue(PREVIOUS_ITERATION_PRESSURE);
			//if ((GetGeometry()[iii].GetValue(SPLIT_ELEMENT))==true)
			//	neighbour_of_splitted_element=true;
			vectorial_mean_vel += one_third * (velocity-solid_velocity); //we just need the difference to compute the nonlinear darcy coeff
		}
		
		const double mean_vel = sqrt(pow(vectorial_mean_vel[0],2)+pow(vectorial_mean_vel[1],2)); //do not forget to add the z vel in 3d!
		
		//const array_1d<double, 3 > & proj0 = GetGeometry()[0].FastGetSolutionStepValue(PRESS_PROJ);
		//const array_1d<double, 3 > & proj1 = GetGeometry()[1].FastGetSolutionStepValue(PRESS_PROJ);
		//const array_1d<double, 3 > & proj2 = GetGeometry()[2].FastGetSolutionStepValue(PRESS_PROJ);


        //we start by calculating the non-enriched functions.
        for (unsigned int i = 0; i < 3; i++)
		{
			for (unsigned int j = 0; j < 3; j++)
			{
				D_matrix(i, (j*2) ) =  -DN_DX(i,0)*Area*one_third;     
				D_matrix(i, (j*2+1) ) =  -DN_DX(i,1)*Area*one_third;
			}
		}
		
        //if (split_element)
		if((geom[0].FastGetSolutionStepValue(WATER_DISTANCE)*geom[1].FastGetSolutionStepValue(WATER_DISTANCE))<0.0 || (geom[1].FastGetSolutionStepValue(WATER_DISTANCE)*geom[2].FastGetSolutionStepValue(WATER_DISTANCE))<0.0)
        {
			//if it is the first iteration, we reset the gradient discontinuity:
			if (iteration_number==1)
			{
				this->GetValue(GRADIENT_DISCONTINUITY)=0.0;
				this->GetValue(ENRICH_RHS)=0.0;
			}
			
			//TOOLS NEEDED TO FIND THE ENRICHMENT SHAPE FUNCTIONS
			//get position of the cut surface
			array_1d<double,3> distances;
			array_1d<double,3>  densities(0);
			array_1d<double,3>  darcies(0);
			boost::numeric::ublas::bounded_matrix<double,3, 2> Nenriched;
			array_1d<double,3>  volumes(0);
			boost::numeric::ublas::bounded_matrix<double,3, 2 > coords;
			boost::numeric::ublas::bounded_matrix<double,3, 3 > Ngauss;
			array_1d<double,3>  signs(0);
			std::vector< Matrix > gauss_gradients(3);
			//fill coordinates
		   
			
			//unsigned int single_triangle_node = 1;
			for (unsigned int i = 0; i < TNumNodes; i++)
			{
				const array_1d<double, 3 > & xyz = this->GetGeometry()[i].Coordinates();
				volumes[i] = 0.0;
				distances[i] = this->GetGeometry()[i].FastGetSolutionStepValue(WATER_DISTANCE);
				//KRATOS_WATCH(distances(i));
				for (unsigned int j = 0; j < 2; j++)
				{
					coords(i, j) = xyz[j];
				}
			}

			for (unsigned int i = 0; i < 3; i++)
				gauss_gradients[i].resize(2, 2, false);  //2 values of the 2 shape functions, and derivates in (xy) direction).
			unsigned int ndivisions = EnrichmentUtilities::CalculateEnrichedShapeFuncions(coords, DN_DX, distances, volumes, Ngauss, signs, gauss_gradients, Nenriched);
			
			//in case we've changed the distance function: CANNOT: (PARALLEL)
			//for (unsigned int i = 0; i < TNumNodes; i++)
			//	this->GetGeometry()[i].FastGetSolutionStepValue(DISTANCE)=distances[i];

			if ((geom[0].FastGetSolutionStepValue(DISTANCE)+geom[1].FastGetSolutionStepValue(DISTANCE)+geom[2].FastGetSolutionStepValue(DISTANCE))>0.0) //that is, there is  NO solid inside this element
			{
				nonlin_darcy_coeff=0.0;
				darcy_coeff=0.0;
			}
				
			array_1d<double,3>  laplacian_coeffs(0);
			for (unsigned int i = 0; i!=ndivisions; i++) //the 3 partitions of the triangle
			{
				if (signs[i]>0)
				{
					densities(i) = mDENSITY_AIR;
					darcies(i) = darcy_coeff*0.0;
				}
				else
				{
					densities(i) = mDENSITY_WATER;
					darcies(i) = darcy_coeff + mean_vel*nonlin_darcy_coeff;
				}
				laplacian_coeffs(i) = 1.0/(densities(i)*(1.0/delta_t+darcies(i)));
			}
			
			Laplacian_matrix =  - prod(DN_DX,trans(DN_DX)) *(volumes(0)*laplacian_coeffs(0)+volumes(1)*laplacian_coeffs(1)+volumes(2)*laplacian_coeffs(2));

			//and now the rest of the things needed
			
			//boost::numeric::ublas::bounded_matrix<double, 2, 2> DN_DX_enrich;
			//boost::numeric::ublas::bounded_matrix<double, 2, 2> DN_DX_enrich_negative;
			boost::numeric::ublas::bounded_matrix<double, 1, 1 > Laplacian_enrich;
			noalias(Laplacian_enrich) = ZeroMatrix(1,1);	
			//double inv_Laplacian_enrich_weighted;
			//double inv_Laplacian_enrich;
			boost::numeric::ublas::bounded_matrix<double, 1, 3 > mixed_Laplacian;
			noalias(mixed_Laplacian) = ZeroMatrix(1,3);	
			
			
			//To calculate D* =int (N*DN_DX_enrich).  we must use the N at the gauss point of the partition (returned by the enrichment util) and multiply by the area. 
			//notice that this first row is only the upper part of the mixed D. we also need the lower one for the jump 
			boost::numeric::ublas::bounded_matrix<double, 1, 6 > D_matrix_mixed;
			noalias(D_matrix_mixed) = ZeroMatrix(1,6);	
			
				
			
			//the matrices we need at the end will be:
			//boost::numeric::ublas::bounded_matrix<double, 2, 2 > D_mod;

			for (unsigned int i = 0; i < ndivisions; i++)  //we go through the 3 partitions of the domain
			{
				Laplacian_enrich(0,0) -= (pow(gauss_gradients[i](0,0),2)+pow(gauss_gradients[i](0,1),2))*volumes(i)*laplacian_coeffs(i);
				//and the 'mixed laplacians' (standard shape functions * enrichments)
				for (unsigned int j = 0; j < 3; j++) //we go through the 3 standard shape functions
				{
					mixed_Laplacian(0,j) -= DN_DX(j,0)*(gauss_gradients[i](0,0)*volumes(i)*laplacian_coeffs(i))+ DN_DX(j,1)*(gauss_gradients[i](0,1)*volumes(i)*laplacian_coeffs(i));
					//and also the D matrixes

					//D_matrix_mixed(0,j*2) -= gauss_gradients[i](0,0) *volumes(i)*Ngauss(i,j); //integrated by parts!
					//D_matrix_mixed(0,j*2+1) -= gauss_gradients[i](0,1) *volumes(i)*Ngauss(i,j);
					D_matrix_mixed(0,j*2) += DN_DX(j,0) * volumes(i) * Nenriched(i,0); //D matrix
					D_matrix_mixed(0,j*2+1) += DN_DX(j,1) * volumes(i) * Nenriched(i,0); //D matrix
					
					
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
			
			for (unsigned int i = 0; i < 6; i++)
				extra_enrich_rhs +=D_matrix_mixed(0,i)*fract_vel(i,0) ;
			
			
			///gamma=1;
			//for (unsigned int i = 0; i < 3; i++)
			//	enrich_rhs += mixed_Laplacian(0,i) *old_pressures(i);	
				
			for (unsigned int i = 0; i < 3; i++)
				rRightHandSideVector(i) -= mixed_Laplacian(0,i) *inv_Laplacian_enrich_weighted * (enrich_rhs + extra_enrich_rhs);
			
			
			
			for (unsigned int i = 0; i < 3; i++)
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
			if ((geom[0].FastGetSolutionStepValue(WATER_DISTANCE)+geom[1].FastGetSolutionStepValue(WATER_DISTANCE)+geom[2].FastGetSolutionStepValue(WATER_DISTANCE))<0.0) //that is, there is water inside this element
			{
				if ((geom[0].FastGetSolutionStepValue(DISTANCE)+geom[1].FastGetSolutionStepValue(DISTANCE)+geom[2].FastGetSolutionStepValue(DISTANCE))>0.0) //that is, there is  NO solid inside this element
					laplacian_coefficient = delta_t/rCurrentProcessInfo[DENSITY_WATER];
				else //the water is flowing through a porous media.
					laplacian_coefficient = 1.0/(rCurrentProcessInfo[DENSITY_WATER]*(1.0/delta_t+(darcy_coeff+nonlin_darcy_coeff*mean_vel)));
			}
			else
			{
				if ((geom[0].FastGetSolutionStepValue(DISTANCE)+geom[1].FastGetSolutionStepValue(DISTANCE)+geom[2].FastGetSolutionStepValue(DISTANCE))>0.0) //that is, there is  NO solid inside this element
					laplacian_coefficient = delta_t/rCurrentProcessInfo[DENSITY_AIR];
				else
					laplacian_coefficient = 1.0/(rCurrentProcessInfo[DENSITY_AIR]*(1.0/delta_t+darcy_coeff*0.0)); //air can move freely
			}	
			Laplacian_matrix =  -prod(DN_DX,trans(DN_DX))*Area*laplacian_coefficient;  // B B^T  (standard laplacian, now we must add the contribution of the condensed new nodes.
		}		
        

        //done, now we save it into the LHS
        noalias(rLeftHandSideMatrix) = Laplacian_matrix;//*(delta_t+TauOne);///gamma=1
        boost::numeric::ublas::bounded_matrix<double, 3,1 > temp_rhs;
        noalias(temp_rhs) = prod(D_matrix,fract_vel);
        //KRATOS_WATCH(temp_rhs);
        rRightHandSideVector(0)  += temp_rhs(0,0);
        rRightHandSideVector(1)  += temp_rhs(1,0);
        rRightHandSideVector(2)  += temp_rhs(2,0);

        
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
	
		
		//KRATOS_WATCH(temp_rhs(0,0));
		//KRATOS_WATCH(temp_rhs(1,0));
		//KRATOS_WATCH(temp_rhs(2,0));
        
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
			ms_temp[iii] = GetGeometry()[iii].FastGetSolutionStepValue(WATER_PRESSURE);
		noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix,ms_temp);
		
		KRATOS_CATCH("");
	}
	

	//************************************************************************************
	//************************************************************************************
	

		void FluidPhasePFEM22D::CalculateViscousRHS(ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY
		
					double mDENSITY_AIR = rCurrentProcessInfo[DENSITY_AIR]; // * (1.0-theta) ;
					double mDENSITY_WATER = rCurrentProcessInfo[DENSITY_WATER]; // * (1.0-theta) ;
					const double porosity=rCurrentProcessInfo[POROSITY];
					const double diameter=rCurrentProcessInfo[DIAMETER];
					const double water_visc=0.0001;
					const double water_dens=1000.0;
					double darcy_coeff = 150.0*pow((1.0-porosity),2)/(pow(porosity,3))*water_visc/(pow(diameter,2));//rCurrentProcessInfo[LIN_DARCY_COEF];
					double nonlin_darcy_coeff = 1.75*(1.0-porosity)/(pow(porosity,3))*water_dens/(diameter); ;//rCurrentProcessInfo[NONLIN_DARCY_COEF];
					double delta_t = rCurrentProcessInfo[DELTA_TIME];
					//array_1d<double,3> & gravity= CurrentProcessInfo[GRAVITY];
					//double x_force  = gravity(0);
					//double y_force  = gravity(1);	
					const array_1d<double,3> zero3 = ZeroVector(3);
					const array_1d<double,3> gravity =  rCurrentProcessInfo[GRAVITY];
					const double nodal_weight = 1.0/ (1.0 + double (2) );
					
					boost::numeric::ublas::bounded_matrix<double, (2+1), 3 > fluid_velocities = ZeroMatrix((2+1), 3);	
					boost::numeric::ublas::bounded_matrix<double, (2+1), 3 > solid_velocities = ZeroMatrix((2+1), 3);	
					boost::numeric::ublas::bounded_matrix<double, (2+1), 2 > DN_DX;
					array_1d<double, (2+1) > N;
					double Area;
					Geometry<Node<3> >& geom = this->GetGeometry();
					GeometryUtils::CalculateGeometryData(geom, DN_DX, N, Area);
					
					double density =1.0;
					double lhs_coeff = 1.0;
					double rhs_fluid_coeff = 1.0;
					
					array_1d<double, 3 > vectorial_mean_vel=ZeroVector(3);
					for(unsigned int iii = 0; iii<3; iii++)
					{
						array_1d<double,3> velocity = GetGeometry()[iii].FastGetSolutionStepValue(WATER_VELOCITY);
						array_1d<double,3> solid_velocity = GetGeometry()[iii].FastGetSolutionStepValue(VELOCITY);
						vectorial_mean_vel += nodal_weight * (velocity-solid_velocity); //we just need the difference to compute the nonlinear darcy coeff
					}
					const double mean_vel = sqrt(pow(vectorial_mean_vel[0],2)+pow(vectorial_mean_vel[1],2)); //do not forget to add the z vel in 3d!
					
					
					if((geom[0].FastGetSolutionStepValue(WATER_DISTANCE)*geom[1].FastGetSolutionStepValue(WATER_DISTANCE))<0.0 || (geom[1].FastGetSolutionStepValue(WATER_DISTANCE)*geom[2].FastGetSolutionStepValue(WATER_DISTANCE))<0.0)
					{
						//TOOLS NEEDED TO FIND THE ENRICHMENT SHAPE FUNCTIONS
						//get position of the cut surface
						//const unsigned int LocalSize = geom.size()*2;
						unsigned int TNumNodes = geom.size();
						array_1d<double,(2+1)> distances;
						boost::numeric::ublas::bounded_matrix<double,3*(2-1), 2> Nenriched;
						array_1d<double,(3*(2-1))> volumes;
						array_1d<double,(3*(2-1))> densities;
						array_1d<double,(3*(2-1))> darcies;
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
							const array_1d<double, 3 > & velocity = geom[i].FastGetSolutionStepValue(WATER_VELOCITY);
							const array_1d<double, 3 > & solid_velocity = geom[i].FastGetSolutionStepValue(VELOCITY);

							
							volumes[i] = 0.0;
							distances[i] = geom[i].FastGetSolutionStepValue(WATER_DISTANCE);
							for (unsigned int j = 0; j < (2); j++)
							{
								coords(i, j) = xyz[j];
								fluid_velocities(i,j) = velocity[j];
								solid_velocities(i,j) = solid_velocity[j];
							}
						}

						for (unsigned int i = 0; i < 3*(2-1) ; i++)
							gauss_gradients[i].resize(2, (2), false);  //2 values of the 2 shape functions, and derivates in xy(z) direction).
						//calling the enrichment function
						unsigned int ndivisions;
							ndivisions = EnrichmentUtilities::CalculateEnrichedShapeFuncions(coords, DN_DX, distances, volumes, Ngauss, signs, gauss_gradients, Nenriched); //, face_gauss_N, face_gauss_Nenriched, face_Area, face_n  ,single_triangle_node);

						if ((geom[0].FastGetSolutionStepValue(DISTANCE)+geom[1].FastGetSolutionStepValue(DISTANCE)+geom[2].FastGetSolutionStepValue(DISTANCE))>0.0) //that is, there is  NO solid inside this element
						{
							nonlin_darcy_coeff=0.0;
							darcy_coeff=0.0;
						}
						for (unsigned int i = 0; i < ndivisions; i++)
						{
							//in case we've changed the distance function inside the enrichmentutility:
							///geom[i].FastGetSolutionStepValue(DISTANCE)=distances[i];
							
							//we are in a loop so we also change the sign of each partition
							
							if (signs[i]>0.0)
							{
								densities(i)= mDENSITY_AIR;
								darcies[i] = darcy_coeff*0.0;
							}
							else
							{
								densities(i)= mDENSITY_WATER;
								darcies[i] = darcy_coeff+nonlin_darcy_coeff*mean_vel;
							}
						}

						double rhs_darcy_node=0.0;
						double lhs_node =0.0;
						double rhs_fluid_node =0.0;
						
						
						for (unsigned int j=0; j!=(2+1); j++)
						{
							
							geom[j].FastGetSolutionStepValue(NODAL_AREA) += Area * nodal_weight;
							
							array_1d<double, 3 > & current_rhs = geom[j].FastGetSolutionStepValue(RHS);

							rhs_darcy_node=0.0;
							lhs_node =0.0;
							rhs_fluid_node =0.0;
						
							for (unsigned int i=0;i!=ndivisions;i++)
							{
								rhs_fluid_node += volumes(i)*Ngauss(i,j)*densities(i)/delta_t;
								lhs_node += volumes(i)*Ngauss(i,j)*densities(i)*(1.0/delta_t+darcies(i));
								rhs_darcy_node += volumes(i)*Ngauss(i,j)*densities(i)*darcies(i);

							}
							geom[j].FastGetSolutionStepValue(NODAL_MASS) += lhs_node;
							for (unsigned int k=0;k!=(2);k++) //k component of local stress
								current_rhs[k] += (rhs_darcy_node*solid_velocities(j,k)+rhs_fluid_node*fluid_velocities(j,k)+gravity(k)*rhs_fluid_node*delta_t);

						}
					}
					
					else	
					{	
						if ((geom[0].FastGetSolutionStepValue(WATER_DISTANCE)+geom[1].FastGetSolutionStepValue(WATER_DISTANCE)+geom[2].FastGetSolutionStepValue(WATER_DISTANCE))<0.0) 
						{
							if ((geom[0].FastGetSolutionStepValue(DISTANCE)+geom[1].FastGetSolutionStepValue(DISTANCE)+geom[2].FastGetSolutionStepValue(DISTANCE))>0.0) //that is, there is  NO solid inside this element
							{
								lhs_coeff = rCurrentProcessInfo[DENSITY_WATER]/delta_t;
								rhs_fluid_coeff = rCurrentProcessInfo[DENSITY_WATER]/delta_t;
								darcy_coeff=0.0; //if it is pure water, the darcy coeff goes to zero (no porous media to slow it down
								density = rCurrentProcessInfo[DENSITY_WATER];
							}
							else //the water is flowing through a porous media.
							{
								lhs_coeff = (rCurrentProcessInfo[DENSITY_WATER]*(1.0/delta_t+darcy_coeff+nonlin_darcy_coeff*mean_vel));
								rhs_fluid_coeff = (rCurrentProcessInfo[DENSITY_WATER]*(1.0/delta_t));
								density = rCurrentProcessInfo[DENSITY_WATER];
							}
						}
						else
						{
							if ((geom[0].FastGetSolutionStepValue(DISTANCE)+geom[1].FastGetSolutionStepValue(DISTANCE)+geom[2].FastGetSolutionStepValue(DISTANCE))>0.0) //that is, there is  NO solid inside this element
							{
								lhs_coeff = rCurrentProcessInfo[DENSITY_AIR]/delta_t;
								rhs_fluid_coeff = rCurrentProcessInfo[DENSITY_AIR]/delta_t;
								darcy_coeff=0.0; //the porous media is invisible to the air phase
								density = rCurrentProcessInfo[DENSITY_AIR];
							}
							else
							{
								lhs_coeff = rCurrentProcessInfo[DENSITY_AIR]*(1.0/delta_t+darcy_coeff*0.0);
								rhs_fluid_coeff = rCurrentProcessInfo[DENSITY_AIR]/delta_t;
								darcy_coeff=darcy_coeff*0.0; //the porous media is invisible to the air phase
								density = rCurrentProcessInfo[DENSITY_AIR];
							}
						}
					
	
						
						for (unsigned int i = 0; i < (2+1); i++)
						{
							const array_1d<double, 3 > & velocity = geom[i].FastGetSolutionStepValue(WATER_VELOCITY);
							const array_1d<double, 3 > & solid_velocity = geom[i].FastGetSolutionStepValue(VELOCITY);
							for (unsigned int j = 0; j < (2); j++)
							{
								fluid_velocities(i,j) = velocity[j];
								solid_velocities(i,j) = solid_velocity[j];
							}
							
						}

						for (unsigned int j=0; j!=(2+1); j++)
						{
							
								geom[j].FastGetSolutionStepValue(NODAL_AREA) += Area*nodal_weight;
														
								array_1d<double, 3 > & current_rhs = geom[j].FastGetSolutionStepValue(RHS);
								for (unsigned int k=0;k!=(2);k++) //k component of local stress
										current_rhs[k] += (density*(darcy_coeff+nonlin_darcy_coeff*mean_vel)*solid_velocities(j,k)+rhs_fluid_coeff*fluid_velocities(j,k)+gravity(k)*density)*nodal_weight*Area;
								
								geom[j].FastGetSolutionStepValue(NODAL_MASS) += Area*lhs_coeff*nodal_weight;							
							//geom[j].UnSetLock();
						}
					}
				
		
		KRATOS_CATCH("");
	}
	
	//*****************************************************************************
	//***************************************************************************
	

	void FluidPhasePFEM22D::CalculatePressureProjection(ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY
		

					const double delta_t = rCurrentProcessInfo[DELTA_TIME];
					const double mDENSITY_AIR = rCurrentProcessInfo[DENSITY_AIR]; // * (1.0-theta) ;
					const double mDENSITY_WATER = rCurrentProcessInfo[DENSITY_WATER]; // * (1.0-theta) ;
										const double porosity=rCurrentProcessInfo[POROSITY];
					const double diameter=rCurrentProcessInfo[DIAMETER];
					const double water_visc=0.0001;
					const double water_dens=1000.0;
					double darcy_coeff = 150.0*pow((1.0-porosity),2)/(pow(porosity,3))*water_visc/(pow(diameter,2));//rCurrentProcessInfo[LIN_DARCY_COEF];
					double nonlin_darcy_coeff = 1.75*(1.0-porosity)/(pow(porosity,3))*water_dens/(diameter); ;//rCurrentProcessInfo[NONLIN_DARCY_COEF];
					
					//const double delta_t = rCurrentProcessInfo[DELTA_TIME];
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
					array_1d<double,(3)> vectorial_mean_vel=ZeroVector(3);
					for (unsigned int i=0; i<(2+1);i++)
					{
						pressure(i) = geom[i].FastGetSolutionStepValue(WATER_PRESSURE);
						vectorial_mean_vel += mass_factor*(geom[i].FastGetSolutionStepValue(VELOCITY) - geom[i].FastGetSolutionStepValue(WATER_VELOCITY));
					}
					const double mean_vel = sqrt(pow(vectorial_mean_vel[0],2)+pow(vectorial_mean_vel[1],2)); // do not forget to add the z vel in 3d!
						
					/*
					boost::numeric::ublas::bounded_matrix<double, (2+1), (2-1)*6 > G_matrix; //(gradient)
					noalias(G_matrix) = ZeroMatrix((2+1), (2-1)*6);	
					//new (upper) part of D, corresponding the first enrichment function (gradient discontinuity)
					boost::numeric::ublas::bounded_matrix<double, 1, (2-1)*6 > G_matrix_mixed;
					noalias(G_matrix_mixed) = ZeroMatrix(1,(2-1)*6);	
					//lower part of D, the one corresponding to the second enrichmend function (jump)
					boost::numeric::ublas::bounded_matrix<double, 1, (2-1)*6 > G_matrix_mixed_jump;
					noalias(G_matrix_mixed_jump) = ZeroMatrix(1,(2-1)*6);	
					*/
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
					array_1d<double,(3*(2-1))> darcies;
					boost::numeric::ublas::bounded_matrix<double,(2+1), 2 > coords;
					boost::numeric::ublas::bounded_matrix<double, 3*(2-1), (2+1) > Ngauss;
					array_1d<double,(3*(2-1))> signs;
					std::vector< Matrix > gauss_gradients(((2-1)*3));
					//fill coordinates
				   
					//unsigned int single_triangle_node;

					
			
						
						
					
					
					if((geom[0].FastGetSolutionStepValue(WATER_DISTANCE)*geom[1].FastGetSolutionStepValue(WATER_DISTANCE))<0.0 || (geom[1].FastGetSolutionStepValue(WATER_DISTANCE)*geom[2].FastGetSolutionStepValue(WATER_DISTANCE))<0.0)
					{
						//to begin with we calculate the gradient discontinuity:
						const int iteration_number =  rCurrentProcessInfo[NL_ITERATION_NUMBER];
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
								pressure_array(i) = geom[i].FastGetSolutionStepValue(WATER_PRESSURE);//-geom[i].GetSolutionStepValue(PRESSURE,1);
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
						
						
						
						//double input_jump_value = this->GetValue(PRESS_DISCONTINUITY);
						//double jump_value=-input_jump_value;
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
							distances[i] = geom[i].FastGetSolutionStepValue(WATER_DISTANCE);
							//KRATOS_WATCH(distances(i));
							for (unsigned int j = 0; j < (2); j++)
								coords(i, j) = xyz[j];
						}

						for (unsigned int i = 0; i < ((2-1)*3) ; i++)
							gauss_gradients[i].resize(2, 2, false);  //2 values of the 2 shape functions, and derivates in (xy) direction).
						unsigned int ndivisions;
						ndivisions = EnrichmentUtilities::CalculateEnrichedShapeFuncions(coords, DN_DX, distances, volumes, Ngauss, signs, gauss_gradients, Nenriched); //, face_gauss_N, face_gauss_Nenriched, face_Area, face_n  ,single_triangle_node);


						//in case we've changed the distance function:
						if ((geom[0].FastGetSolutionStepValue(DISTANCE)+geom[1].FastGetSolutionStepValue(DISTANCE)+geom[2].FastGetSolutionStepValue(DISTANCE))>0.0) //that is, there is  NO solid inside this element
						{
						   darcy_coeff=0.0;
						   nonlin_darcy_coeff=0.0;
					   }
						for (unsigned int i = 0; i < ndivisions; i++)
						{
							//geom[i].FastGetSolutionStepValue(DISTANCE)=distances[i];
							
							if (signs[i]>0)
							{
								densities(i) = mDENSITY_AIR;
								darcies(i) = darcy_coeff*0.0;
							}
							else
							{
								densities(i) = mDENSITY_WATER;
								darcies(i) = darcy_coeff+nonlin_darcy_coeff*mean_vel;
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
									//G_matrix_mixed_jump_no_ro(0,j*(2)+k) += gauss_gradients[i](1,k) *volumes(i)*Ngauss(i,j);

									
									//for (unsigned int l = 0; l < (2+1); l++) //shape function derivative (dNl/dk)
									//{
									//	G_matrix(l, (j*(2)+k) ) += DN_DX(l,k)*volumes(i)*Ngauss(i,j)*inv_densities(i);
									//}
									
								}
								
								nodal_masses(j) += volumes(i)*Ngauss(i,j)*densities(i)*(1.0/delta_t+darcies(i));
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
						
						//const double reduction_factor=0.00000001; //the splitted elements should not add anything to the stabilization terms to avoid problems. we keep a small value to avoid problems due to zero nodal area.
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
								current_press_proj_no_ro[k] += G_matrix_mixed_no_ro(0,2*i+k)*gradient_discontinuity;// + G_matrix_mixed_jump_no_ro(0,2*i+k)*jump_value;
							}
							
							//geom[i].FastGetSolutionStepValue(NODAL_AREA) += Area*mass_factor*reduction_factor;
							geom[i].FastGetSolutionStepValue(NODAL_MASS) += nodal_masses(i);
							
						}
						
					}
					else //normal (uncut) element
					{
						noalias(G_matrix_no_ro) = ZeroMatrix((2+1),(2-1)*6);
						//noalias(G_matrix) = ZeroMatrix((2+1),(2-1)*6);
						
						//double density =1.0;
						double lhs_coeff = 1.0;
						//double rhs_coeff = 1.0;

						if ((geom[0].FastGetSolutionStepValue(WATER_DISTANCE)+geom[1].FastGetSolutionStepValue(WATER_DISTANCE)+geom[2].FastGetSolutionStepValue(WATER_DISTANCE))<0.0) //that is, there is water inside this element
						{
							if ((geom[0].FastGetSolutionStepValue(DISTANCE)+geom[1].FastGetSolutionStepValue(DISTANCE)+geom[2].FastGetSolutionStepValue(DISTANCE))>0.0) //that is, there is  NO solid inside this element
							{
								lhs_coeff = rCurrentProcessInfo[DENSITY_WATER]/delta_t;
								//rhs_coeff = rCurrentProcessInfo[DENSITY_WATER]/delta_t
								//darcy_coeff=0.0; //if it is pure water, the darcy coeff goes to zero (no porous media to slow it down
								//density = rCurrentProcessInfo[DENSITY_WATER];
							}
							else //the water is flowing through a porous media.
							{
								lhs_coeff = (rCurrentProcessInfo[DENSITY_WATER]*(1.0/delta_t+darcy_coeff+nonlin_darcy_coeff*mean_vel));
								//rhs_coeff = (rCurrentProcessInfo[DENSITY_WATER]*(1.0/delta_t));
								//density = rCurrentProcessInfo[DENSITY_WATER];
							}
						}
						else
						{
							if ((geom[0].FastGetSolutionStepValue(DISTANCE)+geom[1].FastGetSolutionStepValue(DISTANCE)+geom[2].FastGetSolutionStepValue(DISTANCE))>0.0) //that is, there is  NO solid inside this element
								lhs_coeff = rCurrentProcessInfo[DENSITY_AIR]/delta_t;
							else	
								lhs_coeff = rCurrentProcessInfo[DENSITY_AIR]*(1.0/delta_t+darcy_coeff*0.0);
							//rhs_coeff = rCurrentProcessInfo[DENSITY_AIR]/delta_t;
							//darcy_coeff=0.0; //the porous media is invisible to the air phase
							//density = rCurrentProcessInfo[DENSITY_AIR];
						}
							

						
						for (unsigned int i = 0; i < (2+1); i++) //loop in shape functions (row)
						{
							for (unsigned int j = 0; j < (2+1) ; j++) //loop in shape functions (column)
							{
								for (unsigned int k = 0; k < (2) ; k++) //x,y,(z)
								{
								//G_matrix(i, (j*2)+k ) = DN_DX(i,k)*Area*mass_factor;
								G_matrix_no_ro(i, (j*2)+k ) = DN_DX(i,k)*Area*mass_factor;
								}
							}
						}
						
						//G_matrix /=density;

						//now we save the data:
						for (unsigned int i=0; i!=(2+1); i++) //loop around the nodes of the element to add contribution to node i
						{
							
							array_1d<double, 3 > & current_press_proj_no_ro = geom[i].FastGetSolutionStepValue(PRESS_PROJ_NO_RO);
							//array_1d<double, 3 > & current_press_proj = geom[i].FastGetSolutionStepValue(PRESS_PROJ);
								
							for (unsigned int j = 0; j < (2+1) ; j++) //loop around the nodes of the element to add contribution of node j to node i
							{
								
								for (unsigned int k = 0; k < (2) ; k++) //x,y,(z)
								{
									//current_press_proj[k] += G_matrix(j,i*(2)+k)*(pressure(j));///-old_pressures(j)); //gamma=0!
									current_press_proj_no_ro[k] += G_matrix_no_ro(j,i*(2)+k)*(pressure(j));///-old_pressures(j)); //gamma=0!
								}
							}
							
							//geom[i].FastGetSolutionStepValue(NODAL_AREA) += Area*mass_factor;
							geom[i].FastGetSolutionStepValue(NODAL_MASS) += Area*mass_factor*lhs_coeff;
							
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
	void FluidPhasePFEM22D::InitializeSolutionStep(ProcessInfo& CurrentProcessInfo)
	{
		KRATOS_TRY
		///WARNING!!!! not calculating the pressure projection since i'm supposing it's being done before!
		/* 

		int FractionalStepNumber = CurrentProcessInfo[FRACTIONAL_STEP];

		if (FractionalStepNumber == 5) //calculation of stabilization terms
		{
			boost::numeric::ublas::bounded_matrix<double, 3, 2 > DN_DX;
			array_1d<double, 3 > N;
			array_1d<double, 2 > vel_gauss;
			double Area;
			
			GeometryUtils::CalculateGeometryData(GetGeometry(), DN_DX, N, Area);

			//array_1d<double, 3 > & fv0 = GetGeometry()[0].FastGetSolutionStepValue(FRACT_VEL);
			//array_1d<double, 3 > & w0 = GetGeometry()[0].FastGetSolutionStepValue(MESH_VELOCITY);
			array_1d<double, 3 > & press_proj0 = GetGeometry()[0].FastGetSolutionStepValue(PRESS_PROJ);
			//array_1d<double, 3 > & conv_proj0 = GetGeometry()[0].FastGetSolutionStepValue(CONV_PROJ);
			double p0 = GetGeometry()[0].FastGetSolutionStepValue(PRESSURE);
			//const double rho0 = GetGeometry()[0].FastGetSolutionStepValue(DENSITY);

			//array_1d<double, 3 > & fv1 = GetGeometry()[1].FastGetSolutionStepValue(FRACT_VEL);
			//array_1d<double, 3 > & w1 = GetGeometry()[1].FastGetSolutionStepValue(MESH_VELOCITY);
			array_1d<double, 3 > & press_proj1 = GetGeometry()[1].FastGetSolutionStepValue(PRESS_PROJ);
			//array_1d<double, 3 > & conv_proj1 = GetGeometry()[1].FastGetSolutionStepValue(CONV_PROJ);
			double p1 = GetGeometry()[1].FastGetSolutionStepValue(PRESSURE);
			//const double rho1 = GetGeometry()[1].FastGetSolutionStepValue(DENSITY);

			//array_1d<double, 3 > & fv2 = GetGeometry()[2].FastGetSolutionStepValue(FRACT_VEL);
			//array_1d<double, 3 > & w2 = GetGeometry()[2].FastGetSolutionStepValue(MESH_VELOCITY);
			array_1d<double, 3 > & press_proj2 = GetGeometry()[2].FastGetSolutionStepValue(PRESS_PROJ);
			//array_1d<double, 3 > & conv_proj2 = GetGeometry()[2].FastGetSolutionStepValue(CONV_PROJ);
			double p2 = GetGeometry()[2].FastGetSolutionStepValue(PRESSURE);
			//const double rho2 = GetGeometry()[2].FastGetSolutionStepValue(DENSITY);

			//double density = 0.3333333333333333333333 * (rho0 + rho1 + rho2);

			//calculation of the pressure gradient (saved in vel_gauss)
			//note that here we calculate it "strong"
			vel_gauss[0] = DN_DX(0, 0)*(p0) + DN_DX(1, 0)*(p1) + DN_DX(2, 0)*(p2);
			vel_gauss[1] = DN_DX(0, 1)*(p0) + DN_DX(1, 1)*(p1) + DN_DX(2, 1)*(p2);
			vel_gauss *= Area;

			//press_proj += G*p
			GetGeometry()[0].SetLock();
			press_proj0[0] += N[0] * vel_gauss[0];
			press_proj0[1] += N[0] * vel_gauss[1];
			GetGeometry()[0].UnSetLock();

			GetGeometry()[1].SetLock();
			press_proj1[0] += N[1] * vel_gauss[0];
			press_proj1[1] += N[1] * vel_gauss[1];
			GetGeometry()[1].UnSetLock();

			GetGeometry()[2].SetLock();
			press_proj2[0] += N[2] * vel_gauss[0];
			press_proj2[1] += N[2] * vel_gauss[1];
			GetGeometry()[2].UnSetLock();
		}
        */
		KRATOS_CATCH("");
	}
	
	//************************************************************************************
	//************************************************************************************


	//************************************************************************************
	//************************************************************************************
	void FluidPhasePFEM22D::PressureEquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo)
	{
		unsigned int number_of_nodes = GetGeometry().PointsNumber();
		if(rResult.size() != number_of_nodes)
			rResult.resize(number_of_nodes,false);	

		for (unsigned int i=0;i<number_of_nodes;i++)
				rResult[i] = GetGeometry()[i].GetDof(WATER_PRESSURE).EquationId();
	}
	

	//************************************************************************************
	//************************************************************************************
	void FluidPhasePFEM22D::GetPressureDofList(DofsVectorType& ElementalDofList,ProcessInfo& CurrentProcessInfo)
	{
		unsigned int number_of_nodes = GetGeometry().PointsNumber();
		if(ElementalDofList.size() != number_of_nodes)
			ElementalDofList.resize(number_of_nodes);	

		for (unsigned int i=0;i<number_of_nodes;i++)
			ElementalDofList[i] = GetGeometry()[i].pGetDof(WATER_PRESSURE);

	}
	
	
	//************************************************************************************
	//************************************************************************************
	
	void FluidPhasePFEM22D::CalculateSolidLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY
		array_1d<double,3> gravity= rCurrentProcessInfo[GRAVITY];

		//WE MUST CALCULATE Mass and G(or D) matrixes

	    double delta_t = rCurrentProcessInfo[DELTA_TIME];	
		//KRATOS_WATCH(gradient_discontinuity);
		//const bool & split_element = (this->GetValue(SPLIT_ELEMENT));

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
		const double density_water = rCurrentProcessInfo[DENSITY_WATER]*2.0;
		
		
		

		//KRATOS_WATCH(D_matrix)
		//if(false)
		//if((distances[0]*distances[1])<0.0 || (distances[1]*distances[2])<0.0 || has_void_node)
		if(true)
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
			{
				viscosity = viscosity_water;
				if (( this->GetGeometry()[0].FastGetSolutionStepValue(WATER_DISTANCE)+ this->GetGeometry()[1].FastGetSolutionStepValue(WATER_DISTANCE)+ this->GetGeometry()[2].FastGetSolutionStepValue(WATER_DISTANCE))<0.0) //that is, there is  NO water inside this element
					viscosity*=0.1;
			}	
			else
				viscosity = viscosity_air;
			
			
			this->AddViscousTerm(rLeftHandSideMatrix,DN_DX, (viscosity*Area) );

				
			//mass matrix:
			/*
			double density;
			
			if (distances(0)<0.0)
				density=density_water;
			else
				density=density_air;
			*/
			double element_mean_distance=(distances[0]+distances[1]+distances[2])*0.333333333;
			//const double water_fraction = -( element_mean_distance - 1.0) *0.5; 
			const double water_fraction = 0.5*(1.0-(sin(3.14159*element_mean_distance*0.5)));
			double density = density_water*(water_fraction)+density_air*(1.0-water_fraction);
							
			double fourier=0.0;
			if (distances(0)<0.0) fourier = viscosity_water/density_water *delta_t/pow((0.6*this->GetValue(MEAN_SIZE)),2);
			else fourier = viscosity_air/density_air *delta_t/pow((0.6*this->GetValue(MEAN_SIZE)),2);
			theta = 0.2*fourier - 0.1;
			if (theta<0.0) theta=0.0;
			else if (theta>1.0) theta=1.0;
				
			const double Weight = one_third * Area * density;
			
			for (unsigned int i=0; i<TNumNodes ; i++)
			{
				Mass_matrix (i*3,i*3) = Weight;
				Mass_matrix (i*3+1,i*3+1) = Weight;
			}	
			//rLeftHandSideMatrix += Mass_matrix;
			
			//PRESSURE TERMS:

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
			const double   TauOne = 1.0 / ( ( 1.0 / delta_t + 4.0 * viscosity / (density * mElemSize * mElemSize) + 2.0 * AdvVelNorm / mElemSize) );
			//actually we're currently using the kinematic viscosity, so no need to divide the second term by the viscosity, it was already done. Great!
			//const double TauOne =  - 1.0 / (  1.0/ delta_t + 4.0 * (viscosity/density) / (mElemSize * mElemSize) + 2.0 * AdvVelNorm / mElemSize) ;
			////////////////////
			//KRATOS_WATCH(TauOne);
			
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
			 //KRATOS_WATCH(mElemSize * mElemSize);
			//Laplacian_matrix = - prod(DN_DX,trans(DN_DX))*Area*TauOne;  // B B^T  (standard laplacian, now we must add the contribution of the condensed new nodes.

			for (unsigned int i = 0; i < 3; i++)
			{
				for (unsigned int j = 0; j < 3; j++)
				{	
					/*
					rLeftHandSideMatrix(i*3+2, j*3+0 ) -= D_matrix(i,j*2);//  + DN_DX(i,0)*TauOne*one_third*Area/delta_t;     
					rLeftHandSideMatrix(i*3+2, j*3+1 ) -= D_matrix(i,j*2+1);// + DN_DX(i,1)*TauOne*one_third*Area/delta_t;
					
					rLeftHandSideMatrix(j*3+0, i*3+2 ) -=  D_matrix(i,j*2);     
					rLeftHandSideMatrix(j*3+1, i*3+2 ) -= D_matrix(i,j*2+1);
					*/ 
					rLeftHandSideMatrix(i*3+2, j*3+0 ) -= D_matrix(i,j*2);//  + DN_DX(i,0)*TauOne*one_third*Area/delta_t;     
					rLeftHandSideMatrix(i*3+2, j*3+1 ) -= D_matrix(i,j*2+1);// + DN_DX(i,1)*TauOne*one_third*Area/delta_t;
					
					rLeftHandSideMatrix(j*3+0, i*3+2 ) -=  D_matrix(i,j*2);     
					rLeftHandSideMatrix(j*3+1, i*3+2 ) -= D_matrix(i,j*2+1);
					
					rLeftHandSideMatrix(i*3+2, j*3+2 ) += Laplacian_matrix(i,j);
					
					
					
				}
			}
			
			
			//and finally the RHS of the velocity plus the RHS on the pressure due to the stabilzation terms:
			array_1d<double,3> mean_velocity = ZeroVector(3);
			double divergence_n = 0.0;
			
			for (unsigned int i = 0; i < 3; i++)
			{
				mean_velocity += GetGeometry()[i].FastGetSolutionStepValue(VELOCITY)*one_third/delta_t;
				divergence_n += one_third*Area*(DN_DX(i,0)*GetGeometry()[i].FastGetSolutionStepValue(VELOCITY_X)+DN_DX(i,1)*GetGeometry()[i].FastGetSolutionStepValue(VELOCITY_Y));
			}
			
			
			array_1d<double, 2 > rhs_stab = ZeroVector(2);
				//Vector pressure_rhs_stab(3);
			for (unsigned int i = 0; i < 3; i++)
			{
				array_1d<double,3>& node_press_proj = GetGeometry()[i].FastGetSolutionStepValue(PRESS_PROJ);
				rhs_stab[0] += node_press_proj(0)*one_third;
				rhs_stab[1] += node_press_proj(1)*one_third;
			}
				
			const bool use_press_proj=false;
			
				
			for (unsigned int i = 0; i < 3; i++)
			{
				rRightHandSideVector(i*3+0) += one_third*Area*gravity(0)*density;
				rRightHandSideVector(i*3+1) += one_third*Area*gravity(1)*density;
				//rRightHandSideVector(i*3+2) += divergence_n;
				//rRightHandSideVector(i*3+2) -= delta_t*DN_DX(i,1)*gravity(1)*Area*one_third*0.5;
				//rRightHandSideVector(i*3+2) -= TauOne*(DN_DX(i,0)*(gravity(0)+previous_vel_and_press(3*i+0)/delta_t)+DN_DX(i,1)*(gravity(1)+previous_vel_and_press(3*i+1)/delta_t))*Area;
				//rRightHandSideVector(i*3+2) -= TauOne*(DN_DX(i,0)*(gravity(0)-mean_velocity(0))+DN_DX(i,1)*(gravity(1)-mean_velocity(1)))*Area;
				if (use_press_proj)
						rRightHandSideVector(i*3+2) -= TauOne*Area*(DN_DX(i,0)*rhs_stab(0)+DN_DX(i,1)*rhs_stab(1));
				else
					rRightHandSideVector(i*3+2) -= TauOne*(DN_DX(i,0)*(gravity(0))+DN_DX(i,1)*(gravity(1)))*Area;
				//for (unsigned int j = 0; j < 3; j++)
				//	rRightHandSideVector(i*3+2) -= Laplacian_matrix(i,j)*GetGeometry()[j].FastGetSolutionStepValue(PRESSURE);
				//rRightHandSideVector(i*3+2) -= 10.0*Area*one_third*(DN_DX(i,0)*gravity(0)+DN_DX(i,1)*gravity(1));
				//for (unsigned int j = 0; j < 3; j++)
				//	rRightHandSideVector(i*3+2) += Laplacian_matrix(i,j)*(DN_DX(j,0)*gravity(0)+DN_DX(j,1)*gravity(1));
				/*
				for (unsigned int j = 0; j < 3; j++) //j neighbour node
				{
					rRightHandSideVector(i*3+2) += TauOne*(DN_DX(i,0)*(gravity(0)-previous_vel_and_press(3*j+0)/delta_t)+DN_DX(i,1)*(gravity(1)-previous_vel_and_press(3*j+1))/delta_t)*Area*one_third;
				}
				*/ 
			}
			
			//we must actually use a fraction of the rigidity matrix. since the Lefthandside only has it so far, we multiply it by tita:
			//rLeftHandSideMatrix *= theta;

			noalias(rRightHandSideVector) += prod((Mass_matrix/delta_t),previous_vel_and_press);

			noalias(rLeftHandSideMatrix) += Mass_matrix/delta_t;   

			
		} 
		else //split element:
		{
			//TOOLS NEEDED TO FIND THE ENRICHMENT SHAPE FUNCTIONS
			//get position of the cut surface
			
			
			rLeftHandSideMatrix=ZeroMatrix(9,9);
			
			//TOOLS NEEDED TO FIND THE ENRICHMENT SHAPE FUNCTIONS
			//get position of the cut surface
			array_1d<double,6>  local_previous_vel;
			array_1d<double,3>  densities(3);
			boost::numeric::ublas::bounded_matrix<double,3, 4> Nenriched;
			array_1d<double,3>  volumes(3);
			boost::numeric::ublas::bounded_matrix<double,3, 2 > coords;
			boost::numeric::ublas::bounded_matrix<double,3, 3 > Ngauss;
			array_1d<double,3>  signs(3);
			std::vector< Matrix > gauss_gradients(3);
			//std::vector< Matrix > gauss_gradients_in_local_axis(3);
			//boost::numeric::ublas::bounded_matrix<double,2, 2 > RotationMatrix;
			//boost::numeric::ublas::bounded_matrix<double, 3, 2> DN_DX_in_local_axis=ZeroMatrix(3,2);
			
			const int  enrich_velocity_dofs=2;
			const int  enrich_pressure_dofs=2;
			const int  enrich_pressure_offset=2;
			
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
			
			//getting needed information to calculate everything in local coordinates
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
			//KRATOS_WATCH(local_gravity)
			//KRATOS_WATCH(local_gravity);
			
			for (unsigned int i = 0; i < 3; i++)
				gauss_gradients[i].resize(4, 2, false);  //3 values of the 3 shape functions, and derivates in (xy) direction).

			unsigned int ndivisions = EnrichmentUtilities::CalculateEnrichedShapeFuncionsExtended(rRotatedPoints, DN_DX_in_local_axis, distances, volumes, Ngauss, signs, gauss_gradients, Nenriched);

			
			//we start by calculating the non-enriched functions.
			for (unsigned int i = 0; i < 3; i++)
			{
				for (unsigned int j = 0; j < 3; j++)
				{
					//D_matrix(i, (j*2) ) = DN_DX(j,0)*Area*one_third;     
					//D_matrix(i, (j*2+1) ) =  DN_DX(j,1)*Area*one_third;
					//G_matrix(i, (j*2) ) = - DN_DX(i,0)*Area*one_third;     
					//G_matrix(i, (j*2+1) ) = - DN_DX(i,1)*Area*one_third;
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
			//this->AddViscousTerm(rLeftHandSideMatrix,DN_DX, (Weight) );
			
			this->AddViscousTerm(Momentum_matrix,
								DN_DX,
								distances,
								gauss_gradients, 
								viscosities,
								signs,
								volumes,
								ndivisions);
			
			//this->AddViscousTerm(rLeftHandSideMatrix,DN_DX, viscosity_air, viscosity_water, volumes, signs, Ngauss, Area);
			//boost::numeric::ublas::bounded_matrix<double, 8, 8 > LocalAxisExtendedDampMatrix = ZeroMatrix(8,8);
			
			//this->AddViscousTerm(Momentum_matrix, DN_DX, viscosities, signs, volumes, ndivisions);
			
			
			const double element_viscosity=Weight/Area;
			const double element_density= mass/Area;
			double fourier= element_viscosity/element_density *delta_t/pow((0.6*this->GetValue(MEAN_SIZE)),2);
			theta = 0.2*fourier - 0.1;
			if (theta<0.0) theta=0.0;
			else if (theta>1.0) theta=1.0;
			//element_viscosity=0.0;
		
			/*
			KRATOS_WATCH(Nenriched)
			KRATOS_WATCH(gauss_gradients[0])
			KRATOS_WATCH(gauss_gradients[1])
			KRATOS_WATCH(gauss_gradients[2])
			*/
			array_1d<double,3>  lumped_mass = ZeroVector(3);
			double condensed_dof_mass1=0.0;
			double condensed_dof_mass2=0.0;
			double condensed_dof_mass3=0.0;
			double condensed_dof_mass4=0.0;
			double condensed_dof_mass5=0.0;
			double condensed_dof_mass6=0.0;
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
					//condensed_dof_mass1 += volumes(i)*Nenriched(i,2)*densities(i);
					//condensed_dof_mass2 += volumes(i)*Nenriched(i,3)*densities(i);
					//condensed_dof_mass3 += volumes(i)*Nenriched(i,3)*densities(i);
					condensed_dof_mass1 += volumes(i)*Nenriched(i,0)*densities(i);
					condensed_dof_mass2 += volumes(i)*Nenriched(i,0)*densities(i);
					condensed_dof_mass3 += volumes(i)*Nenriched(i,1)*densities(i);
					condensed_dof_mass4 += volumes(i)*Nenriched(i,1)*densities(i);
					if (enrich_velocity_dofs>4)
					{	
						condensed_dof_mass5 += volumes(i)*Nenriched(i,2)*densities(i);
						condensed_dof_mass6 += volumes(i)*Nenriched(i,3)*densities(i);
					}
					   
					Momentum_matrix(9+0,9+0) += volumes(i)*Nenriched(i,3)*densities(i) / delta_t;	 //degrees of freedom 10 and 11 are the enrichment velocities. degrees of freedom 12 and 13 are the enrichment pressures
					Momentum_matrix(9+1,9+1) += volumes(i)*Nenriched(i,3)*densities(i) / delta_t;	 
					//Momentum_matrix(9+2,9+2) += volumes(i)*Nenriched(i,3)*densities(i) / delta_t;	 
					//Momentum_matrix(9+0,9+0) += volumes(i)*Nenriched(i,0)*densities(i) / delta_t;	 
					//Momentum_matrix(9+1,9+1) += volumes(i)*Nenriched(i,0)*densities(i) / delta_t;	 
					//Momentum_matrix(9+2,9+2) += volumes(i)*Nenriched(i,1)*densities(i) / delta_t;	 
					//Momentum_matrix(9+3,9+3) += volumes(i)*Nenriched(i,1)*densities(i) / delta_t;
					if (enrich_velocity_dofs>4)
					{	 
						Momentum_matrix(9+4,9+4) += volumes(i)*Nenriched(i,2)*densities(i) / delta_t;	 
						Momentum_matrix(9+5,9+5) += volumes(i)*Nenriched(i,3)*densities(i) / delta_t;	 
					}

				}
			}
			
			
			//PRESSURE TERMS:
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
			const double   TauOne = 0.01 / (( 1.0 / delta_t + 4.0 * element_viscosity / (mElemSize * mElemSize * element_density ) + 2.0 * AdvVelNorm / mElemSize) );
			//actually we're currently using the kinematic viscosity, so no need to divide the second term by the viscosity, it was already done. Great!
			
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
			
			//inv_densities(0)=1.0;
			//inv_densities(1)=1.0;
			//inv_densities(2)=1.0;
			
			Laplacian_matrix =  - prod(DN_DX_in_local_axis,trans(DN_DX_in_local_axis)) *(volumes(0)*inv_densities(0)+volumes(1)*inv_densities(1)+volumes(2)*inv_densities(2));

			//and now the rest of the things needed
			
			//boost::numeric::ublas::bounded_matrix<double, 2, 2> DN_DX_enrich;
			//boost::numeric::ublas::bounded_matrix<double, 2, 2> DN_DX_enrich_negative;
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
			
			array_1d<double,6> mass_stabilization_terms=ZeroVector(6);
			
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
						//D_matrix_mixed(k,6+0) += gauss_gradients[i](2,0) *volumes(i)*Nenriched(i,k+enrich_pressure_offset);  //we simply use the shape function number 4 (continuous function, discontinuous gradient)
						//D_matrix_mixed(k,6+1) += gauss_gradients[i](3,0) *volumes(i)*Nenriched(i,k+enrich_pressure_offset);
						//D_matrix_mixed(k,6+2) += gauss_gradients[i](3,1) *volumes(i)*Nenriched(i,k+enrich_pressure_offset);
						
						D_matrix_mixed(k,6+0) += gauss_gradients[i](3,0) *volumes(i)*Nenriched(i,k+enrich_pressure_offset);
						D_matrix_mixed(k,6+1) += gauss_gradients[i](3,1) *volumes(i)*Nenriched(i,k+enrich_pressure_offset);
						//D_matrix_mixed(k,6+2) += gauss_gradients[i](1,0) *volumes(i)*Nenriched(i,k+enrich_pressure_offset);
						//D_matrix_mixed(k,6+3) += gauss_gradients[i](1,1) *volumes(i)*Nenriched(i,k+enrich_pressure_offset);
						//if (enrich_velocity_dofs>4)
						//{
						//	D_matrix_mixed(k,6+4) += gauss_gradients[i](2,0) *volumes(i)*Nenriched(i,k+enrich_pressure_offset);
						//	D_matrix_mixed(k,6+5) += gauss_gradients[i](3,0) *volumes(i)*Nenriched(i,k+enrich_pressure_offset);
						//}
						
						//G_matrix_mixed(k,6+0) -= gauss_gradients[i](k+enrich_pressure_offset,0) *volumes(i)*Nenriched(i,2);
						//G_matrix_mixed(k,6+1) -= gauss_gradients[i](k+enrich_pressure_offset,0) *volumes(i)*Nenriched(i,3);
						//G_matrix_mixed(k,6+2) -= gauss_gradients[i](k+enrich_pressure_offset,1) *volumes(i)*Nenriched(i,3);
						
						G_matrix_mixed(k,6+0) -= gauss_gradients[i](k+enrich_pressure_offset,0) *volumes(i)*Nenriched(i,3);
						G_matrix_mixed(k,6+1) -= gauss_gradients[i](k+enrich_pressure_offset,1) *volumes(i)*Nenriched(i,3);
						//G_matrix_mixed(k,6+2) -= gauss_gradients[i](k+enrich_pressure_offset,0) *volumes(i)*Nenriched(i,1);
						//G_matrix_mixed(k,6+3) -= gauss_gradients[i](k+enrich_pressure_offset,1) *volumes(i)*Nenriched(i,1);
						//if (enrich_velocity_dofs>4)
						//{
						//	G_matrix_mixed(k,6+4) -= gauss_gradients[i](k+enrich_pressure_offset,0) *volumes(i)*Nenriched(i,2);
						//	G_matrix_mixed(k,6+5) -= gauss_gradients[i](k+enrich_pressure_offset,0) *volumes(i)*Nenriched(i,3);
						//}
					}
					
					
					//also the standard D and G matrices had to be expanded.

						for (unsigned int k = 0; k < 3; k++) //we go through the 3 standard pressures
						{
							//D_matrix(k,6+0) += gauss_gradients[i](2,0) *volumes(i)*Ngauss(i,k);  //we simply use the shape function number 4 (continuous function, discontinuous gradient)
							//D_matrix(k,6+1) += gauss_gradients[i](3,0) *volumes(i)*Ngauss(i,k);
							//D_matrix(k,6+2) += gauss_gradients[i](3,1) *volumes(i)*Ngauss(i,k);
							
							D_matrix(k,6+0) += gauss_gradients[i](3,0) *volumes(i)*Ngauss(i,k);
							D_matrix(k,6+1) += gauss_gradients[i](3,1) *volumes(i)*Ngauss(i,k);
							//D_matrix(k,6+2) += gauss_gradients[i](1,0) *volumes(i)*Ngauss(i,k);
							//D_matrix(k,6+3) += gauss_gradients[i](1,1) *volumes(i)*Ngauss(i,k);
							//if (enrich_velocity_dofs>4)
							//{
							//	D_matrix(k,6+4) += gauss_gradients[i](2,0) *volumes(i)*Ngauss(i,k);
							//	D_matrix(k,6+5) += gauss_gradients[i](3,0) *volumes(i)*Ngauss(i,k);
							//}
							
							//G_matrix(k,6+0) -= DN_DX_in_local_axis(k,0) *volumes(i)*Nenriched(i,2);
							//G_matrix(k,6+1) -= DN_DX_in_local_axis(k,0) *volumes(i)*Nenriched(i,3);
							//G_matrix(k,6+2) -= DN_DX_in_local_axis(k,1) *volumes(i)*Nenriched(i,3);
							
							G_matrix(k,6+0) -= DN_DX_in_local_axis(k,0) *volumes(i)*Nenriched(i,3);
							G_matrix(k,6+1) -= DN_DX_in_local_axis(k,1) *volumes(i)*Nenriched(i,3);
							//G_matrix(k,6+2) -= DN_DX_in_local_axis(k,0) *volumes(i)*Nenriched(i,1);
							//G_matrix(k,6+3) -= DN_DX_in_local_axis(k,1) *volumes(i)*Nenriched(i,1);
							//if (enrich_velocity_dofs>4)
							//{
							//	G_matrix(k,6+4) -= DN_DX_in_local_axis(k,0) *volumes(i)*Nenriched(i,2);
							//	G_matrix(k,6+5) -= DN_DX_in_local_axis(k,0) *volumes(i)*Nenriched(i,3);
							//}
						}
				}
				
				//KRATOS_WATCH(volumes(i))
				//KRATOS_WATCH(partition_densities(i))
				//and the rhs:
				for (unsigned int j = 0; j < 3; j++) //we go through the 3 standard shape functions
				{
					//array_1d<double,3> velocity = GetGeometry()[j].FastGetSolutionStepValue(VELOCITY);
					//rhs_enrich -= (gauss_gradients[i](0,0)*(gravity(0)+velocity(0))+gauss_gradients[i](0,1)*(gravity(1)+velocity(1)))*volumes(i)*Ngauss(i,j);
				}
				
				for (unsigned int k = 0; k < enrich_pressure_dofs; k++) //we go through the 4 enrichments
					rhs_enrich(k+enrich_velocity_dofs) -= (gauss_gradients[i](k+enrich_pressure_offset,0)*(local_gravity(0,0))+gauss_gradients[i](k+enrich_pressure_offset,1)*(local_gravity(1,0)))*volumes(i);
				
				
			}
			//KRATOS_WATCH(Laplacian_enrich(0,0));
			//TauOne *=0.01;
			const double continuous_factor=1.0;//1;
			Laplacian_enrich*=TauOne;
			Laplacian_matrix*=TauOne*continuous_factor;
			mixed_Laplacian*=TauOne;
			rhs_enrich *=TauOne;
			mass_stabilization_terms *=TauOne/delta_t;
			
			


			
			
			//KRATOS_WATCH(mass_stabilization_terms)
			
			if (enrich_velocity_dofs!=0)
			{
				//rhs_enrich(0) = local_gravity(0,0)*condensed_dof_mass1;
				//rhs_enrich(1) = local_gravity(0,0)*condensed_dof_mass2;
				//rhs_enrich(2) = local_gravity(1,0)*condensed_dof_mass3;
				rhs_enrich(0) = local_gravity(0,0)*condensed_dof_mass1;
				rhs_enrich(1) = local_gravity(1,0)*condensed_dof_mass2;
				//rhs_enrich(2) = local_gravity(0,0)*condensed_dof_mass3;
				//rhs_enrich(3) = local_gravity(1,0)*condensed_dof_mass4;
				//if (enrich_velocity_dofs>4)
				//{
				//	rhs_enrich(4) = local_gravity(0,0)*condensed_dof_mass5;
				//	rhs_enrich(5) = local_gravity(0,0)*condensed_dof_mass6;
				//	

				//	rhs_enrich(4) += N_vel_star(0) * condensed_dof_mass5/delta_t;
				//	rhs_enrich(5) += N_vel_star(0) * condensed_dof_mass6/delta_t;
				//}
				//rhs_enrich(0) += N_vel_star(0) * condensed_dof_mass1/delta_t; //we are using the difference in the Y component
				//rhs_enrich(1) += N_vel_star(1) * condensed_dof_mass2/delta_t;
				
				
				/*
				if (this->Id() == 44752)
				{
					KRATOS_WATCH(negative_side_only_velocity_at_interface)
					KRATOS_WATCH(standard_velocity_at_interface)
					KRATOS_WATCH(N_vel_star);
					KRATOS_WATCH(Ninterface)
					for (unsigned int i=0; i!=3;i++) //nodes
						KRATOS_WATCH(GetGeometry()[i].FastGetSolutionStepValue(VELOCITY));
				}
				* */
			}
			
			//KRATOS_WATCH(rhs_enrich)
			
			
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
			//KRATOS_WATCH(Momentum_matrix)
			
			
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
					
					
					
			//KRATOS_WATCH(condensed_block)
			//KRATOS_WATCH(condensed_rows)
			//KRATOS_WATCH(condensed_columns)
			//KRATOS_WATCH(Laplacian_enrich)
			
			//KRATOS_WATCH(condensed_row)
			//KRATOS_WATCH(rhs_enrich)
			
			//const double inv_Laplacian_enrich_weighted=1.0/Laplacian_enrich(0,0);
			//KRATOS_WATCH(Laplacian_enrich(0,0))
			boost::numeric::ublas::bounded_matrix<double, enrich_velocity_dofs+enrich_pressure_dofs , enrich_velocity_dofs+enrich_pressure_dofs  > inverse_enrichments;
			
			//this->invert44(condensed_block,  inverse_enrichments);
			//this->invert33(condensed_block,  inverse_enrichments);
			/*
			double determinant = (Laplacian_enrich(0,0)*Laplacian_enrich(1,1)) - (Laplacian_enrich(0,1)*Laplacian_enrich(1,0));
			inverse_enrichments(0,0)=Laplacian_enrich(1,1);
			inverse_enrichments(1,1)=Laplacian_enrich(0,0);
			inverse_enrichments(0,1)=-Laplacian_enrich(0,1);
			inverse_enrichments(1,0)=-Laplacian_enrich(1,0);
			inverse_enrichments /= determinant;			
			*/
			//inverse_enrichments(0,0)=1.0/Laplacian_enrich(0,0);
			this->InvertMatrix(condensed_block,inverse_enrichments);
			//KRATOS_WATCH(Momentum_matrix)
			
			//condensing

			
			//KRATOS_WATCH
			boost::numeric::ublas::bounded_matrix<double, 9 , enrich_pressure_dofs+enrich_velocity_dofs  > temp_matrix;
			temp_matrix = prod(trans(condensed_columns),inverse_enrichments);
			rLeftHandSideMatrix -=  prod(temp_matrix,condensed_rows);
			//boost::numeric::ublas::bounded_matrix<double, 9, 9 > art = prod(trans(condensed_column),condensed_row);
			//KRATOS_WATCH(art)
			
			noalias(rRightHandSideVector) -= prod(temp_matrix,rhs_enrich);
			/*
			for (unsigned int i = 0; i < 9; i++)
			{
				rRightHandSideVector(i) -= prod(condensed_column(0,i) * rhs_enrich;
			}
			*/
			///KRATOS_WATCH(rhs_enrich)
			//KRATOS_WATCH(condensed_column);
			//KRATOS_WATCH(D_matrix)
			//standard part:
			
			for (unsigned int i = 0; i < 3; i++)
			{
				for (unsigned int j = 0; j < 3; j++)
				{
					//rLeftHandSideMatrix(j*3+2, i*3+0 ) -=  D_matrix(i,j*2);     
					//rLeftHandSideMatrix(j*3+2, i*3+1 ) -= D_matrix(i,j*2+1);
					rLeftHandSideMatrix(i*3+2, j*3+0 ) -= D_matrix(i,j*2);//  + DN_DX(i,0)*TauOne*one_third*Area/delta_t;     
					rLeftHandSideMatrix(i*3+2, j*3+1 ) -= D_matrix(i,j*2+1);// + DN_DX(i,1)*TauOne*one_third*Area/delta_t;
					
					rLeftHandSideMatrix(j*3+0, i*3+2 ) -=  D_matrix(i,j*2);     
					rLeftHandSideMatrix(j*3+1, i*3+2 ) -= D_matrix(i,j*2+1);
					
					rLeftHandSideMatrix(i*3+2, j*3+2 ) += Laplacian_matrix(i,j);
				}
			}
			//KRATOS_WATCH(rLeftHandSideMatrix);
			/*
			boost::numeric::ublas::bounded_matrix<double, 10, 10 > my_matrix=ZeroMatrix(10,10);
			boost::numeric::ublas::bounded_matrix<double, 1, 10 > my_vector=ZeroMatrix(1,10);
			for (unsigned int i = 0; i < 3; i++)
			{
				for (unsigned int j = 0; j < 3; j++)
				{
					my_matrix(i*2+0,i*2+0)=Mass_matrix(i*3,i*3);
					my_matrix(i*2+1,i*2+1)=Mass_matrix(i*3,i*3);
					
					my_matrix(6+i,j*2+0)=D_matrix(i,j*2);
					my_matrix(6+i,j*2+1)=D_matrix(i,j*2+1);
					
					my_matrix(j*2+0, 6+i)=D_matrix(i,j*2);
					my_matrix(j*2+1, 6+i)=D_matrix(i,j*2+1);
					
					my_matrix(6+i,6+j)=Laplacian_matrix(i,j);
				}
				
				my_matrix(9,i*2+0) = D_matrix_mixed(0,i*2+0);
				my_matrix(9,i*2+1) = D_matrix_mixed(0,i*2+1);
				
				my_matrix(i*2+0,9) = D_matrix_mixed(0,i*2+0);
				my_matrix(i*2+1,9) = D_matrix_mixed(0,i*2+1);
				
				my_matrix(9,6+i) = mixed_Laplacian(0,i);
				
				my_matrix(6+i,9) = mixed_Laplacian(0,i);
				
				my_vector(0,i*2+0) = gravity(0)*Mass_matrix(i*3,i*3);
				my_vector(0,i*2+1) = gravity(1)*Mass_matrix(i*3,i*3);
				
				my_vector(0,6+i) = - (gravity(0)*DN_DX(i,0)+gravity(1)*DN_DX(i,1))*Area;
				
			}
			
			my_matrix(9,9)= Laplacian_enrich(0,0);
			//my_vector(0,9) = rhs_enrich;
			*/
			//KRATOS_WATCH(my_matrix)
			//KRATOS_WATCH(my_vector)
			
			//and finally the RHS of the velocity plus the RHS on the pressure due to the stabilzation terms:
			for (unsigned int i = 0; i < 3; i++)
			{
				//KRATOS_WATCH(Mass_matrix(i*3,i*3));

				//rRightHandSideVector(i*3+0) += one_third*gravity(0)*((0.25+0.125)*1000.0+0.125*1.0);//Mass_matrix(i*3,i*3)*gravity(0);
				//rRightHandSideVector(i*3+1) += 0.16666*gravity(1)*((0.25+0.125)*1.0+0.125*1000.0);//Mass_matrix(i*3+1,i*3+1)*gravity(1);
				//rRightHandSideVector(i*3+0) +=gravity(0)*Mass_matrix(i*3,i*3);
				//rRightHandSideVector(i*3+1) +=gravity(1)*Mass_matrix(i*3,i*3);
				rRightHandSideVector(i*3+0) +=local_gravity(0,0)*lumped_mass(i);
				rRightHandSideVector(i*3+1) +=local_gravity(1,0)*lumped_mass(i);
				
				rRightHandSideVector(i*3+2) -= TauOne*continuous_factor*(DN_DX_in_local_axis(i,0)*local_gravity(0,0)+DN_DX_in_local_axis(i,1)*local_gravity(1,0))*Area;
				
				
				for (unsigned int j = 0; j < 3; j++) //i node-  j partition k neighbour node
				{
					//rRightHandSideVector(i*3+1) += Nenriched(j,0)*volumes(j)*partition_densities(j)*gravity(1);
					for (unsigned int k = 0; k < 3; k++) //i node-  j partition k neighbour node
					{
						
					//rRightHandSideVector(i*3+2) -= TauOne*(DN_DX(i,0)*(gravity(0)+previous_vel_and_press(3*i+0))+DN_DX(i,1)*(gravity(1)+previous_vel_and_press(3*i+1)))*volumes(j)*partition_densities(j);
					//rRightHandSideVector(i*3+2) -= TauOne*(DN_DX(i,0)*(gravity(0)+previous_vel_and_press(3*k+0))+DN_DX(i,1)*(gravity(1)+previous_vel_and_press(3*k+1)))*volumes(j)*Ngauss(j,k);
					}
				}
				 
			}
			
			
			//KRATOS_WATCH(TauOne)
			
			
			//we must actually use a fraction of the rigidity matrix. since the Lefthandside only has it so far, we multiply it by tita:
			//rLeftHandSideMatrix *= theta;

			//noalias(rRightHandSideVector) += prod((Mass_matrix/delta_t),previous_vel_and_press);
			for (unsigned int j = 0; j < 3; j++) //i node-  j partition k neighbour node
			{
				rRightHandSideVector(j*3+0) += lumped_mass(j)* local_previous_vel(j*2) / delta_t;
				rRightHandSideVector(j*3+1) += lumped_mass(j)* local_previous_vel(j*2+1) / delta_t;
			}
			
			
			for (unsigned int i = 0; i < 9; i++) 
				for (unsigned int j = 0; j < 9; j++)
					rLeftHandSideMatrix(i,j) += Momentum_matrix(i,j);
			



			//KRATOS_WATCH(rLeftHandSideMatrix);
			
			
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
					
					//or by subblocks. (more work)
					/*
					boost::numeric::ublas::bounded_matrix<double, 2 , 2 > LocalKBlock;
					//saving unrotated local block
					for (unsigned int k = 0; k < 2; k++) 
						for (unsigned int l = 0; l < 2; l++)
							LocalKBlock(k,l)=LocalLeftHandSideMatrix(i*3+k,j*3+l);	
					//rotationg	
					boost::numeric::ublas::bounded_matrix<double, 2 , 2 > tempBlock=prod(LocalKBlock,rRotationMatrix);
					boost::numeric::ublas::bounded_matrix<double, 2 , 2 > GlobalKBlock=prod(trans(rRotationMatrix),tempBlock);
					//saving info into the global matrix
					for (unsigned int k = 0; k < 2; k++) 
						for (unsigned int l = 0; l < 2; l++)
							rLeftHandSideMatrix(i*3+k,j*3+l)=GlobalKBlock(k,l);

					//now the Dmatrix block, G matrix block and the RHS.
					boost::numeric::ublas::bounded_matrix<double, 1 , 2 >  LocalDBlock;
					boost::numeric::ublas::bounded_matrix<double, 2 , 1 >  LocalGBlock;
					for (unsigned int k = 0; k < 2; k++) 
					{
						LocalDBlock(0,k) = LocalLeftHandSideMatrix(i*3+2,j*3+k);	
						LocalGBlock(k,0) = LocalLeftHandSideMatrix(i*3+k,j*3+2);	
					}
					//rotating	
					boost::numeric::ublas::bounded_matrix<double, 1 , 2 >  GlobalDBlock = prod(LocalDBlock,rRotationMatrix);
					boost::numeric::ublas::bounded_matrix<double, 2 , 1 > GlobalGBlock = prod(trans(rRotationMatrix),LocalGBlock);
					//saving into the global system
					for (unsigned int k = 0; k < 2; k++) 
					{
						rLeftHandSideMatrix(i*3+2,j*3+k) = GlobalDBlock(0,k);
						rLeftHandSideMatrix(i*3+k,j*3+2) = GlobalGBlock(k,0);
					}
					*/
					
							
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
		
		
		//KRATOS_WATCH(rRightHandSideVector)

		
        

		//subtracting the dirichlet term
		// RHS -= LHS*DUMMY_UNKNOWNs
		
		

		//KRATOS_WATCH((rRightHandSideVector/(one_third*Area))*delta_t);


		noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix,previous_vel_and_press);
		
		//KRATOS_WATCH(rRightHandSideVector)
		//		KRATOS_WATCH(rLeftHandSideMatrix)

		
		
		///////////TO MAKE IT TIDIER LATER!!!!!!!!!!!!!
		const int number_of_particles_in_elem = this->GetValue(NUMBER_OF_PARTICLES);
		
		if( (number_of_particles_in_elem==0) )
		//if(false)
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
			//noalias(rRightHandSideVector) = ZeroVector(9);
			noalias(rRightHandSideVector) = -prod(rLeftHandSideMatrix,previous_vel_and_press);
		 }
		KRATOS_CATCH("");
	}
	
	
	
	

	void FluidPhasePFEM22D::SolidCalculatePressureProjection(ProcessInfo& CurrentProcessInfo)
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


	//************************************************************************************
	//************************************************************************************
	void FluidPhasePFEM22D::SolidEquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo)
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
	void FluidPhasePFEM22D::GetSolidDofList(DofsVectorType& ElementalDofList,ProcessInfo& CurrentProcessInfo)
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
	void FluidPhasePFEM22D::AddViscousTerm(MatrixType& OutputMatrix,
                                       const boost::numeric::ublas::bounded_matrix<double, 3, 2>& rShapeDeriv,
                                       const double Weight)
	{
		/*
		const SizeType NumNodes = this->GetGeometry().PointsNumber();

		const double FourThirds = 4.0 / 3.0;
		const double nTwoThirds = -2.0 / 3.0;
		

		SizeType FirstRow(0),FirstCol(0);

		for (SizeType j = 0; j < NumNodes; ++j)
		{
			for (SizeType i = 0; i < NumNodes; ++i)
			{
				// First Row
				rDampMatrix(FirstRow,FirstCol) += Weight * ( FourThirds * rShapeDeriv(i,0) * rShapeDeriv(j,0) + rShapeDeriv(i,1) * rShapeDeriv(j,1) );
				rDampMatrix(FirstRow,FirstCol+1) += Weight * ( nTwoThirds * rShapeDeriv(i,0) * rShapeDeriv(j,1) + rShapeDeriv(i,1) * rShapeDeriv(j,0) );

				// Second Row
				rDampMatrix(FirstRow+1,FirstCol) += Weight * ( nTwoThirds * rShapeDeriv(i,1) * rShapeDeriv(j,0) + rShapeDeriv(i,0) * rShapeDeriv(j,1) );
				rDampMatrix(FirstRow+1,FirstCol+1) += Weight * ( FourThirds * rShapeDeriv(i,1) * rShapeDeriv(j,1) + rShapeDeriv(i,0) * rShapeDeriv(j,0) );

				// Update Counter
				FirstRow += 2;
			}
			FirstRow = 0;
			FirstCol += 2;
		}
		*/
		
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
		C_matrix(0,0)=2.0;
		C_matrix(1,1)=2.0;
		C_matrix(2,2)=1.0;
		/*
		for (unsigned int i=0; i!=3; i++) //i node
		{
			for (unsigned int j=0; j!=3; j++) //j neighbour
			{
					C_matrix(i*3+0,j*3+0)=2.0;
					C_matrix(i*3+1,j*3+1)=2.0;
					C_matrix(i*3+2,j*3+2)=1.0;
			}
		}
		*/ 
		
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
	
	
	//************************************************************************************
	//with rotation and returnging also the row that will be enriched
	void FluidPhasePFEM22D::AddViscousTerm(boost::numeric::ublas::bounded_matrix<double, 13, 13 > & output,
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
				//if (distances(i)*signs(division)>0.0)
				if (true)
				{
					for (unsigned int j=0; j!=(2); j++) //x,y,z
						B_matrix(i*(2)+j,j)=rShapeDeriv(i,j);
					
					//useful for both 2d and 3d:	
					//relating 12 and 21 stresses
					B_matrix(i*(2)+0,2)=rShapeDeriv(i,1);
					B_matrix(i*(2)+1,2)=rShapeDeriv(i,0);
				}
				else
				{
					for (unsigned int j=0; j!=(2); j++) //x,y,z
						B_matrix(i*(2)+j,j)=0.0;
					
					//useful for both 2d and 3d:	
					//relating 12 and 21 stresses
					B_matrix(i*(2)+0,2)=0.0;
					B_matrix(i*(2)+1,2)=0.0;

				}

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
	

	
    inline double FluidPhasePFEM22D::CalculateVol(const double x0, const double y0,
                                      const double x1, const double y1,
                                      const double x2, const double y2
                                     )
    {
        return 0.5 * ((x1 - x0)*(y2 - y0)- (y1 - y0)*(x2 - x0));
    }
    
    inline double FluidPhasePFEM22D::CalculateVolume2D(
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
bool FluidPhasePFEM22D::InvertMatrix(const T& input, T& inverse)
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

    void FluidPhasePFEM22D::CalculateRotationParameters(
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
	
	inline void FluidPhasePFEM22D::CalculateGeometryData(
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

