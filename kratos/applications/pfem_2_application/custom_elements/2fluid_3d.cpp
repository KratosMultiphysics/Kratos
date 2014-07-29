//   
//   Project Name:        Kratos       
//   Last modified by:    $Author: it's me! $
//   Date:                $Date: 2008-08-08 23:58:38 $
//   Revision:            $Revision: 1.0 $
//
//

// Project includes 
#include "includes/define.h"
#include "custom_elements/2fluid_3d.h"
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
	PFEM23D::PFEM23D(IndexType NewId, GeometryType::Pointer pGeometry)
		: Element(NewId, pGeometry)
	{		
		//DO NOT ADD DOFS HERE!!!
	}

	//************************************************************************************
	//************************************************************************************
	PFEM23D::PFEM23D(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
		: Element(NewId, pGeometry, pProperties)
	{

	}

	Element::Pointer PFEM23D::Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const
	{
		return Element::Pointer(new PFEM23D(NewId, GetGeometry().Create(ThisNodes), pProperties));
	}

	PFEM23D::~PFEM23D()
	{
	}

	//************************************************************************************
	//************************************************************************************
	
	
	void PFEM23D::AddExplicitContribution(ProcessInfo& rCurrentProcessInfo)
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
		
		case 7: //calculating pressure projection. notthing returned. saving data in PRESS_PROJ, PRESS_PROJ_NO_RO , NODAL_MASS and NODAL_AREA
		{
			this->CalculateMassMatrix(rCurrentProcessInfo);
			break;
		}
		
		default:
		{
			KRATOS_ERROR(std::logic_error,"Unexpected value for FRACTIONAL_STEP index: ",rCurrentProcessInfo[FRACTIONAL_STEP]);
		}
		}

		KRATOS_CATCH("");
	}
	
	void PFEM23D::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
													VectorType& rRightHandSideVector,
													ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY;

		switch ( rCurrentProcessInfo[FRACTIONAL_STEP] )
		{
		case 1:
		{
			this->CalculateLocalFractionalVelocitySystem(rLeftHandSideMatrix,rRightHandSideVector,rCurrentProcessInfo);
			break;
		}
		case 2:
		{
			//KRATOS_WATCH("calculating an element")
			this->CalculateLocalPressureSystem(rLeftHandSideMatrix,rRightHandSideVector,rCurrentProcessInfo);
			break;
		}
		case 4:
		{
			this->CalculateLocalFinalVelocitySystem(rLeftHandSideMatrix,rRightHandSideVector,rCurrentProcessInfo);
			break;
			///CAMBIAR ESTO DESPUÉS!!!!!!!!!!!!!
			//KRATOS_ERROR(std::logic_error,"Full solution of end of step velocity is not implemented, see Calculate(VELOCITY)","");
			//break;
		}
		
		case 5:
		{
			this->CalculateLocalThermalSystem(rLeftHandSideMatrix,rRightHandSideVector,rCurrentProcessInfo);
			break;
			///CAMBIAR ESTO DESPUÉS!!!!!!!!!!!!!
			//KRATOS_ERROR(std::logic_error,"Full solution of end of step velocity is not implemented, see Calculate(VELOCITY)","");
			//break;
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
	
	void PFEM23D::EquationIdVector(EquationIdVectorType& rResult,
												ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY;

		switch ( rCurrentProcessInfo[FRACTIONAL_STEP] )
		{
		case 1:
		{
			this->FractionalVelocityEquationIdVector(rResult,rCurrentProcessInfo);
			break;
		}
		case 2:
		{
			this->PressureEquationIdVector(rResult,rCurrentProcessInfo);
			break;
		}
		case 4:
		{
			this->VelocityEquationIdVector(rResult,rCurrentProcessInfo);
			break;
		}
		case 5:
		{
			this->ThermalEquationIdVector(rResult,rCurrentProcessInfo);
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
	

	void PFEM23D::GetDofList(DofsVectorType& rElementalDofList,
										  ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY;

		switch ( rCurrentProcessInfo[FRACTIONAL_STEP] )
		{
		case 1:
		{
			this->GetFractionalVelocityDofList(rElementalDofList,rCurrentProcessInfo);
			break;
		}
		case 2:
		{
			this->GetPressureDofList(rElementalDofList,rCurrentProcessInfo);
			break;
		}
		case 4:
		{
			this->GetVelocityDofList(rElementalDofList,rCurrentProcessInfo);
			break;
		}
		case 5:
		{
			this->GetThermalDofList(rElementalDofList,rCurrentProcessInfo);
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
	
	
	
	void PFEM23D::CalculateLocalFractionalVelocitySystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY
	
		KRATOS_ERROR(std::logic_error, "USE LOCALFINALVELOCITY. SAME PROCEDURE, ONLY ADDING THE VISCOUS TERM", "");
		
		KRATOS_CATCH("");
	}
	
	//************************************************************************************
	//************************************************************************************
	
	// LAPLACIAN with 1/density for volume conservation
	
	void PFEM23D::CalculateLocalPressureSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY

		//rCurrentProcessInfo[DELTA_TIME]=1.0;
		//ProcessInfo& CurrentProcessInfo = mr_model_part.GetProcessInfo();
		const double viscosity_water = rCurrentProcessInfo[VISCOSITY_WATER]; // * (1.0-theta) ;
		double Density = 0.0 ; //used for stabilization;
		double delta_t = rCurrentProcessInfo[DELTA_TIME];
		int iteration_number =  rCurrentProcessInfo[NL_ITERATION_NUMBER];
		Vector densities(4);
		Vector partition_densities(6);
		Vector partition_viscosities(6);
		//KRATOS_WATCH(delta_t);
		//bool neighbour_of_splitted_element=false;

		double enrich_rhs=0.0;
		double extra_enrich_rhs=0.0;

		double input_jump_value=0.0;
		this->SetValue(PRESS_DISCONTINUITY,input_jump_value);
		double jump_value=-input_jump_value*0.5;
		
		const bool & split_element = (this->GetValue(SPLIT_ELEMENT));
		
		//boost::numeric::ublas::bounded_matrix<double,3,2> msDN_DX;
		//boost::numeric::ublas::bounded_matrix<double,2,2> msD;
		const double one_quarter = 1.0/4.0;
		array_1d<double,4> msN; //dimension = number of nodes
		array_1d<double,4> ms_temp; //dimension = number of nodes

		const unsigned int LocalSize = GetGeometry().size();
		unsigned int TNumNodes = GetGeometry().size();
		//boost::numeric::ublas::bounded_matrix<double, 3,1 > enrich_rhs;

		//getting data for the given geometry
//        const unsigned int LocalSize = (TDim + 1) * TNumNodes;

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
        boost::numeric::ublas::bounded_matrix<double, 4, 12 > D_matrix; //(gradient)
        boost::numeric::ublas::bounded_matrix<double, 4, 4 > Laplacian_matrix;
        //boost::numeric::ublas::bounded_matrix<double, TNumNodes*2, TDim> N_vel_matrix; //useless to calculate it, it's simply the DN_DX matrix multiplied by 1/3, arranged in a different shape.
       
        Geometry<Node<3> >& geom = this->GetGeometry();
         
        GeometryUtils::CalculateGeometryData(geom, DN_DX, N, Area);
        //KRATOS_WATCH(Area);

        boost::numeric::ublas::bounded_matrix<double, 12, 1 > fract_vel;
        array_1d<double, 4 > old_pressures;
        array_1d<double, 4 > temperatures;
        
        for(unsigned int iii = 0; iii<4; iii++)
        {
			const array_1d<double, 3 > & velocity = geom[iii].FastGetSolutionStepValue(VELOCITY);
			fract_vel(iii*3,0) = velocity[0];///delta_t;
			fract_vel(iii*3+1,0) = velocity[1];///delta_t;
			fract_vel(iii*3+2,0) = velocity[2];///delta_t;
			//old_pressures(iii) = GetGeometry()[iii].GetSolutionStepValue(PRESSURE,1);
			old_pressures(iii) = GetGeometry()[iii].FastGetSolutionStepValue(PREVIOUS_ITERATION_PRESSURE);
			temperatures(iii) = GetGeometry()[iii].FastGetSolutionStepValue(TEMPERATURE);
		}
		const array_1d<double, 3 > & proj0 = GetGeometry()[0].FastGetSolutionStepValue(PRESS_PROJ);
		const array_1d<double, 3 > & proj1 = GetGeometry()[1].FastGetSolutionStepValue(PRESS_PROJ);
		const array_1d<double, 3 > & proj2 = GetGeometry()[2].FastGetSolutionStepValue(PRESS_PROJ);
		const array_1d<double, 3 > & proj3 = GetGeometry()[3].FastGetSolutionStepValue(PRESS_PROJ);


        //we start by calculating the non-enriched functions.
        for (unsigned int i = 0; i < 4; i++)
		{
			for (unsigned int j = 0; j < 4; j++)
			{
				D_matrix(i, (j*3) ) =  -DN_DX(i,0)*Area*one_quarter;     
				D_matrix(i, (j*3+1) ) =  -DN_DX(i,1)*Area*one_quarter;
				D_matrix(i, (j*3+2) ) =  -DN_DX(i,2)*Area*one_quarter;
			}
		}
		
		//CALCULATING TAU (only if we are not near the interfase and we are in step 1 or more)
		double TauOne = 0.0; 
		double sqElemSize=0.0;
		//if( ((rCurrentProcessInfo[NL_ITERATION_NUMBER])>1) && (split_element==false) && (neighbour_of_splitted_element==false) )
		if (true)
		{
			//double Viscosity = rCurrentProcessInfo[VISCOSITY];;
			//double AdvVelNorm = sqrt(pow(this->GetGeometry()[0].FastGetSolutionStepValue(VELOCITY_X),2)+pow(this->GetGeometry()[0].FastGetSolutionStepValue(VELOCITY_Y),2));
			//const double TimeFactor = rCurrentProcessInfo.GetValue(DYNAMIC_TAU);
			//double mElemSize;
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
			//mElemSize = sqrt(mElemSize);
			//const double   TauOne = 1.0 / (Density * ( TimeFactor / DeltaTime + 4.0 * Viscosity / (mElemSize * mElemSize) + 2.0 * AdvVelNorm / mElemSize) );
			//actually we're currently using the kinematic viscosity, so no need to divide the second term by the viscosity, it was already done. Great!
			//TauOne =  1.0 / ( ( 1.0/ delta_t + 4.0 * Viscosity / (mElemSize * mElemSize) + 2.0 * AdvVelNorm / mElemSize) );
			
		}
		 
		
		

        if (split_element)
        {
			//if it is the first iteration, we reset the gradient discontinuity:
			if (iteration_number==1)
			{
				this->GetValue(GRADIENT_DISCONTINUITY)=0.0;
				this->GetValue(ENRICH_RHS)=0.0;
			}
			
			//TOOLS NEEDED TO FIND THE ENRICHMENT SHAPE FUNCTIONS
			//get position of the cut surface
			array_1d<double,(4)> distances;
			boost::numeric::ublas::bounded_matrix<double,6, 2> Nenriched;
			array_1d<double,(6)> volumes;
			array_1d<double,(6)> densities;
			array_1d<double,(6)> inv_densities;
			array_1d<double,(6)> viscosities;
			boost::numeric::ublas::bounded_matrix<double,(4), 3 > coords;
			boost::numeric::ublas::bounded_matrix<double, 6, 4 > Ngauss;
			array_1d<double,(6)> signs;
			std::vector< Matrix > gauss_gradients(6);
		   

			
			//unsigned int single_triangle_node = 1;
			for (unsigned int i = 0; i < TNumNodes; i++)
			{
				const array_1d<double, 3 > & xyz = this->GetGeometry()[i].Coordinates();
				//volumes[i] = 0.0;
				distances[i] = this->GetGeometry()[i].FastGetSolutionStepValue(DISTANCE);
				//KRATOS_WATCH(distances(i));
				for (unsigned int j = 0; j < 3; j++)
					coords(i, j) = xyz[j];
			}

			for (unsigned int i = 0; i < 6; i++)
				gauss_gradients[i].resize(2, 3, false);  ///WARNING. ITS SIZE IS CURRENTLY only ONE, no JUMP!! //2 values of the 2 shape functions, and derivates in (xyz) directions).
			unsigned int ndivisions = EnrichmentUtilities::CalculateEnrichedShapeFuncions(coords, DN_DX, distances, volumes, Ngauss, signs, gauss_gradients, Nenriched ); // ,face_gauss_N, face_gauss_Nenriched, face_Area, face_n  ,single_triangle_node);
			//KRATOS_WATCH(ndivisions)
			//in case we've changed the distance function:CANNOT: (PARALLEL)
			//for (unsigned int i = 0; i < TNumNodes; i++)
			//	this->GetGeometry()[i].FastGetSolutionStepValue(DISTANCE)=distances[i];
			//KRATOS_WATCH(volumes)
			
			Density=0.0; //resetting to zero
			for (unsigned int i = 0; i!=ndivisions; i++) //the 6 partitions of the thetraedron
			{
				double partition_temperature=0.0;
				for (unsigned int j = 0; j < (4); j++) //shape function (Nj)
					partition_temperature += temperatures(i)*Ngauss(i,j);

				if (signs(i)>0)
				{
					partition_densities(i)=CalculateAirDensity(partition_temperature);
					partition_viscosities(i)=CalculateAirViscosity(partition_temperature);
				}
				else
				{
					partition_densities(i)=rCurrentProcessInfo[DENSITY_WATER];		
					double viscosity=viscosity_water;
					if (partition_temperature>100.0)
						viscosity=viscosity_water*(1.0 - (partition_temperature-100.0)*0.1 );
					if (viscosity<1.0) //(crude oil viscosity)
						viscosity=1.0;	
					partition_viscosities(i)=viscosity;
				}	
				inv_densities(i)=1.0/(partition_densities(i)/delta_t + partition_viscosities(i)/0.04 ) ;	
				Density+=partition_densities(i)*volumes(i);
			}
			Density*=1.0/Area;


			jump_value*=signs(0);
			
			
			Laplacian_matrix =  - prod(DN_DX,trans(DN_DX));
			
			if (true) //did this to delete the coefficient later
			{
				double coefficient=0.0;
				for (unsigned int i=0; i<ndivisions ; i++)

					coefficient += volumes[i]*inv_densities[i];
				Laplacian_matrix *= coefficient;
			}

			//and now the rest of the things needed
			
			//boost::numeric::ublas::bounded_matrix<double, 2, 2> DN_DX_enrich;
			//boost::numeric::ublas::bounded_matrix<double, 2, 2> DN_DX_enrich_negative;
			boost::numeric::ublas::bounded_matrix<double, 2, 2 > Laplacian_enrich;
			noalias(Laplacian_enrich) = ZeroMatrix(2,2);	
			//double inv_Laplacian_enrich_weighted;
			//double inv_Laplacian_enrich;
			boost::numeric::ublas::bounded_matrix<double, 1, 4 > mixed_Laplacian;
			noalias(mixed_Laplacian) = ZeroMatrix(1,4);	
			
			boost::numeric::ublas::bounded_matrix<double, 1, 4 > mixed_Laplacian_jump;
			noalias(mixed_Laplacian_jump) = ZeroMatrix(1,4);
			
			
			//To calculate D* =int (N*DN_DX_enrich).  we must use the N at the gauss point of the partition (returned by the enrichment util) and multiply by the area. 
			//notice that this first row is only the upper part of the mixed D. we also need the lower one for the jump 
			boost::numeric::ublas::bounded_matrix<double, 1, 12 > D_matrix_mixed;
			noalias(D_matrix_mixed) = ZeroMatrix(1,12);	
			//lower part of D, the one corresponding to the second enrichmend function (jump)
			boost::numeric::ublas::bounded_matrix<double, 1, 12 > D_matrix_mixed_jump;
			noalias(D_matrix_mixed_jump) = ZeroMatrix(1,12);	
			
				
			
			//the matrices we need at the end will be:
			//boost::numeric::ublas::bounded_matrix<double, 2, 2 > D_mod;

			for (unsigned int i = 0; i < ndivisions; i++)  //we go through the [ndivisions] partitions of the domain
			{
				Laplacian_enrich(0,0) -= (pow(gauss_gradients[i](0,0),2)+pow(gauss_gradients[i](0,1),2)+pow(gauss_gradients[i](0,2),2))*volumes(i)*inv_densities(i);
				Laplacian_enrich(0,1) -= ((gauss_gradients[i](0,0)*gauss_gradients[i](1,0))+(gauss_gradients[i](0,1)*gauss_gradients[i](1,1))+(gauss_gradients[i](0,2)*gauss_gradients[i](1,2)))*volumes(i)*inv_densities(i);
				Laplacian_enrich(1,1) -= (pow(gauss_gradients[i](1,0),2)+pow(gauss_gradients[i](1,1),2)+pow(gauss_gradients[i](1,2),2))*volumes(i)*inv_densities(i);
				Laplacian_enrich(1,0) -= ((gauss_gradients[i](0,0)*gauss_gradients[i](1,0))+(gauss_gradients[i](0,1)*gauss_gradients[i](1,1))+(gauss_gradients[i](0,2)*gauss_gradients[i](1,2)))*volumes(i)*inv_densities(i);
				//and the 'mixed laplacians' (standard shape functions * enrichments)
				for (unsigned int j = 0; j < 4; j++) //we go through the 4 standard shape functions
				{
					mixed_Laplacian(0,j) -= (DN_DX(j,0)*(gauss_gradients[i](0,0) ) + DN_DX(j,1)*(gauss_gradients[i](0,1)) + DN_DX(j,2)*(gauss_gradients[i](0,2)) ) *volumes(i)*inv_densities(i);
					mixed_Laplacian_jump(0,j) -= ( DN_DX(j,0)*(gauss_gradients[i](1,0)) + DN_DX(j,1)*(gauss_gradients[i](1,1)) + DN_DX(j,2)*(gauss_gradients[i](1,2) ) ) *volumes(i)*inv_densities(i);
					//and also the D matrixes

					D_matrix_mixed(0,j*3) -= gauss_gradients[i](0,0) *volumes(i)*Ngauss(i,j);
					D_matrix_mixed(0,j*3+1) -= gauss_gradients[i](0,1) *volumes(i)*Ngauss(i,j);
					D_matrix_mixed(0,j*3+2) -= gauss_gradients[i](0,2) *volumes(i)*Ngauss(i,j);
					D_matrix_mixed_jump(0,j*3) -= gauss_gradients[i](1,0) *volumes(i)*Ngauss(i,j);
					D_matrix_mixed_jump(0,j*3+1) -= gauss_gradients[i](1,1) *volumes(i)*Ngauss(i,j);
					D_matrix_mixed_jump(0,j*3+2) -= gauss_gradients[i](1,2) *volumes(i)*Ngauss(i,j);
					

				}
				
				
			}
			
			//the 'jump' degree of freedom will not be condensed since it's imposed
			//therefore we'll only have to calculate the inverse of the reduced matrix of the gradient jump enrich function. which is a 1x1 matrix (inverse of the 1,1 position of the enriched laplacian)
			const double inv_Laplacian_enrich_weighted=1.0/Laplacian_enrich(0,0);
			this->SetValue(INV_LAPLACIAN_ENRICH,(inv_Laplacian_enrich_weighted)) ;//(delta_t+TauOne)));
			//KRATOS_WATCH(inv_Laplacian_enrich_weighted)
			
			//since we're imposing this DOF, we'll do:
			//NONO! not defined yet the 3d function!
			//for (unsigned int i = 0; i < 4; i++)
			//	rRightHandSideVector(i) -= mixed_Laplacian_jump(0,i)*jump_value;
				
			//when there's no jump there's nothing in the RHS to condense. but when we add it, it is equal to -Laplacian_enrich*jump_value due to the elimination of this DOF.
			//and just like the LHS, we must condense it and substrac it (actually it's (-)*(-)=+ )
			
			//enrich_rhs = -Laplacian_enrich(0,1)*jump_value;
			//NONO! not defined yet the 3d function!
			for (unsigned int i = 0; i < 12; i++)
				extra_enrich_rhs +=D_matrix_mixed(0,i)*fract_vel(i,0) ;
			
			
			///gamma=1;
			//for (unsigned int i = 0; i < 3; i++)
			//	enrich_rhs += mixed_Laplacian(0,i) *old_pressures(i);	
			
			//NOT ENRICHENING!!	
			//for (unsigned int i = 0; i < 4; i++)
			//	rRightHandSideVector(i) -= mixed_Laplacian(0,i) *inv_Laplacian_enrich_weighted * (enrich_rhs + extra_enrich_rhs);
			//Laplacian_matrix -= inv_Laplacian_enrich_weighted * prod(trans(mixed_Laplacian),mixed_Laplacian);
			
			if((this->GetValue(ENRICH_LHS_ROW_3D)).size()!=4)
				KRATOS_ERROR(std::logic_error, "WRONG SIZE in LOCALSYSTEM", "");
				//(this->GetValue(ENRICH_LHS_ROW_3D)).resize(4);
			for (unsigned int i = 0; i < 4; i++)
			{
				(this->GetValue(ENRICH_LHS_ROW_3D))(i)=mixed_Laplacian(0,i); //*(delta_t+TauOne);
				//KRATOS_WATCH(mixed_Laplacian(0,i));
			}
			//KRATOS_WATCH(this->GetValue(ENRICH_LHS_ROW_3D))	
			
			//now we must condense these matrixes into the original system (D_matrix, Laplacian_matrix)
			
			///D_matrix-= inv_Laplacian_enrich_weighted * prod(trans(mixed_Laplacian),D_matrix_mixed); //actually inv_Lap.. has 1/delta_t and mix_Lap has delta_t so at the end we do not write it. but do not forget that they are part of the system. to condense we must use the lefthandside, which is actually Laplacian*delta_t
			///no need to do this ^ cos it is in the rhs, done when added : rRightHandSideVector(i) -= mixed_Laplacian(0,i) *inv_Laplacian_enrich_weighted * (enrich_rhs + extra_enrich_rhs);
			
			
			///gamma=1!
			//rRightHandSideVector += prod(Laplacian_matrix,old_pressures)*delta_t;
			
			
				
		}
		 
		else //uncut(normal, not interfase) element
        {
			
			double element_temperature= 0.25*(temperatures(0)+temperatures(1)+temperatures(2)+temperatures(3));
			//double density_air = CalculateAirDensity(element_temperature);
			double inv_density;
			if (geom[0].FastGetSolutionStepValue(DISTANCE)<0.0)
			{
				
				double viscosity=viscosity_water;
				if (element_temperature>100.0)
						viscosity=viscosity_water*(1.0 - (element_temperature-100.0)*0.1 );
				if (viscosity<1.0)
						viscosity=1.0;
				inv_density = 1.0/(rCurrentProcessInfo[DENSITY_WATER]/delta_t + viscosity/0.04);
			}
			else
			{
				double viscosity= CalculateAirViscosity(element_temperature);
				inv_density = 1.0/(this->CalculateAirDensity(element_temperature)/delta_t + viscosity/0.04);
			}
			Laplacian_matrix =  -prod(DN_DX,trans(DN_DX))*Area*inv_density;  // B B^T  (standard laplacian, now we must add the contribution of the condensed new nodes.
		}		
        

        //done, now we save it into the LHS
        noalias(rLeftHandSideMatrix) = Laplacian_matrix;//*(delta_t+TauOne);///gamma=1
        boost::numeric::ublas::bounded_matrix<double, 4,1 > temp_rhs;
        noalias(temp_rhs) = prod(D_matrix,fract_vel);
        //KRATOS_WATCH(temp_rhs);
        rRightHandSideVector(0)  += temp_rhs(0,0);
        rRightHandSideVector(1)  += temp_rhs(1,0);
        rRightHandSideVector(2)  += temp_rhs(2,0);
        rRightHandSideVector(3)  += temp_rhs(3,0);

        
        //stabilization contribution:
        array_1d<double, 3 > vel_gauss;
        vel_gauss[0] = N[0] * proj0[0] + N[1] * proj1[0] + N[2] * proj2[0] + N[3] * proj3[0];
		vel_gauss[1] = N[0] * proj0[1] + N[1] * proj1[1] + N[2] * proj2[1] + N[3] * proj3[1];
		vel_gauss[2] = N[0] * proj0[2] + N[1] * proj1[2] + N[2] * proj2[2] + N[3] * proj3[2];
		vel_gauss *= TauOne*Area;
		//KRATOS_WATCH(vel_gauss);
		///gammma=1!
		if (rCurrentProcessInfo[NL_ITERATION_NUMBER]>1)
		{
			noalias(rRightHandSideVector) -= prod(DN_DX, vel_gauss);
			rRightHandSideVector += prod(Laplacian_matrix,old_pressures);//*delta_t;
		}
		
		//KRATOS_WATCH(temp_rhs(0,0));
		//KRATOS_WATCH(temp_rhs(1,0));
		//KRATOS_WATCH(temp_rhs(2,0));
        
        enrich_rhs+=extra_enrich_rhs;
        //to recover the gradient jump later, we must save the RHS of the enrichment, that is:
		this->SetValue(ENRICH_RHS,enrich_rhs);
        
        
        ///if it is an inactive element we must turn it off. But to avoid problems in case we've turned off all the elements that contribute to a node, we simply do:
		if (this->GetValue(IS_INACTIVE)==true) //it means we must not add it!
		{
			rLeftHandSideMatrix*=1.0e-6;
			rRightHandSideVector*=1.0e-6;
		}
        
        
		//subtracting the dirichlet term
		// RHS -= LHS*DUMMY_UNKNOWNs		
		for(unsigned int iii = 0; iii<LocalSize; iii++)
			ms_temp[iii] = GetGeometry()[iii].FastGetSolutionStepValue(PRESSURE);
		noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix,ms_temp);
		
		KRATOS_CATCH("");
	}
	
	//************************************************************************************
	//************************************************************************************
	
	//************************************************************************************
	//************************************************************************************
	
	void PFEM23D::CalculateLocalFinalVelocitySystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY

	    double delta_t = rCurrentProcessInfo[DELTA_TIME];	
		//KRATOS_WATCH(gradient_discontinuity);
		const bool & split_element = (this->GetValue(SPLIT_ELEMENT));

		double theta=1.0;
	
		
		array_1d<double,12>  previous_vel;
		array_1d<double, 4 > temperatures;
		array_1d<double,4> distances;
			
        for(unsigned int iii = 0; iii<4; iii++)
        {
			temperatures(iii) = GetGeometry()[iii].FastGetSolutionStepValue(TEMPERATURE);
			distances[iii] = this->GetGeometry()[iii].FastGetSolutionStepValue(DISTANCE);
			array_1d<double,3>& vel = GetGeometry()[iii].FastGetSolutionStepValue(VELOCITY);
			previous_vel(iii*3) = vel[0];
			previous_vel(iii*3+1) = vel[1];
			previous_vel(iii*3+2) = vel[2];

		}
		
		//boost::numeric::ublas::bounded_matrix<double,3,2> msDN_DX;
		//boost::numeric::ublas::bounded_matrix<double,2,2> msD;
		const double one_quarter = 0.25;
  		array_1d<double,4> msN; //dimension = number of nodes

		const unsigned int LocalSize = GetGeometry().size()*3;
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
        //boost::numeric::ublas::bounded_matrix<double, 4, 12 > D_matrix; //(gradient)
        //noalias(D_matrix) = ZeroMatrix(4,12);	

        boost::numeric::ublas::bounded_matrix<double, 12, 12 > Mass_matrix;
        noalias(Mass_matrix) = ZeroMatrix(LocalSize,LocalSize);	
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
			double element_temperature=0.0;
			double density=0.0;
			double viscosity=0.0;
			for (unsigned int j = 0; j < 4; j++) //we go through the 3 standard shape functions
			{
				element_temperature += temperatures[j]*one_quarter;
			}
				
			if (distances(0)<0.0)
			{
				density=density_water;
				viscosity=viscosity_water;
				if (element_temperature>100.0)
						viscosity=viscosity_water*(1.0 - (element_temperature-100.0)*0.1 );
				if (viscosity<1.0)
						viscosity=1.0;
				this->AddViscousTerm(rLeftHandSideMatrix,DN_DX, (viscosity*Area) );
			}
			else
			{
				density = this->CalculateAirDensity(element_temperature);
				viscosity = this->CalculateAirViscosity(element_temperature);
				this->AddViscousTerm(rLeftHandSideMatrix,DN_DX, (viscosity*Area) );
			}
				
				
				
			double fourier=0.0;
			if (distances(0)<0.0) fourier = viscosity_water/density_water *delta_t/pow((0.6*this->GetValue(MEAN_SIZE)),2);
			else fourier = viscosity_air/density_air *delta_t/pow((0.6*this->GetValue(MEAN_SIZE)),2);
			theta = 0.2*fourier - 0.1;
			if (theta<0.0) theta=0.0;
			else if (theta>1.0) theta=1.0;
				
			const double Weight = one_quarter * Area * density;
			
			for (unsigned int i=0; i<LocalSize ; i++)
				Mass_matrix (i,i) = Weight;
				
			//rLeftHandSideMatrix += Mass_matrix;	
			
		} 
		else //split element:
		{
			//TOOLS NEEDED TO FIND THE ENRICHMENT SHAPE FUNCTIONS
			//get position of the cut surface
			//KRATOS_ERROR(std::logic_error, "IMPLICIT STEP FIRST STEP NOT YET IMPLEMENTED IN 3D.. USE LOCALFINALVELOCITY", "");
			//TOOLS NEEDED TO FIND THE ENRICHMENT SHAPE FUNCTIONS
			//get position of the cut surface
			array_1d<double,6>  densities(6);
			boost::numeric::ublas::bounded_matrix<double,6, 2> Nenriched;
			array_1d<double,6>  volumes(6);
			boost::numeric::ublas::bounded_matrix<double,4, 3 > coords;
			boost::numeric::ublas::bounded_matrix<double,6, 4 > Ngauss;
			array_1d<double,6>  signs(6);
			std::vector< Matrix > gauss_gradients(6);
			//fill coordinates
		   
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
			double mass=0.0;
			for (unsigned int i=0;i<ndivisions;i++) //we go over the three partitions;
			{
				double partition_temperature=0.0;
				for (unsigned int j = 0; j < 4; j++)
					partition_temperature += temperatures[j]*Ngauss(i,j);
				if (signs(i)>0.0)
				{
					double viscosity = CalculateAirViscosity(partition_temperature);
					densities(i) = CalculateAirDensity(partition_temperature);
					Weight += volumes(i)*viscosity;
					mass += volumes(i)*densities(i);
				}
				else
				{
					double viscosity=viscosity_water;
					if (partition_temperature>100.0)
						viscosity=viscosity_water*(1.0 - (partition_temperature-100.0)*0.1 );
					if (viscosity<1.0) //(crude oil viscosity)
						viscosity=1.0;
					//viscosities(i)=viscosity;
					Weight += volumes(i)*viscosity;
					densities(i)=density_water;
					mass += volumes(i)*density_water;
				}
			}
			this->AddViscousTerm(rLeftHandSideMatrix,DN_DX, (Weight) );
			//this->AddViscousTerm(rLeftHandSideMatrix,DN_DX, viscosity_air, viscosity_water, volumes, signs, Ngauss, Area);
			
			const double element_viscosity=Weight/Area;
			const double element_density= mass/Area;
			double fourier= element_viscosity/element_density *delta_t/pow((0.6*this->GetValue(MEAN_SIZE)),2);
			theta = 0.2*fourier - 0.1;
			if (theta<0.0) theta=0.0;
			else if (theta>1.0) theta=1.0;


			
			//KRATOS_WATCH(volumes)
			//KRATOS_WATCH(densities)
			//KRATOS_WATCH(signs)
			//KRATOS_WATCH(Ngauss)
			
			for (unsigned int i = 0; i < 4; i++)  // local node (or shape function)
			{
				//KRATOS_WATCH(enrich_lhs(i));
				for (unsigned int j = 0; j < ndivisions; j++) //partitions
				{
					Mass_matrix(3*i,3*i) +=  volumes(j)*Ngauss(j,i)*densities(j);	 
				}
				Mass_matrix(3*i+1,3*i+1) = Mass_matrix(3*i,3*i);
				Mass_matrix(3*i+2,3*i+2) = Mass_matrix(3*i,3*i);
			}
			
		}
		
		
       
		
		//we must actually use  a fraction of the rigidity matrix. since the Lefthandside only has it so far, we multiply it by 0.5:
		//rLeftHandSideMatrix *= theta;
	
        //finally we add to the righthandside BOTH the mass matrix and the rigidity matrix:
        //noalias(rRightHandSideVector) = prod((Mass_matrix/delta_t),previous_vel);
		for (unsigned int i = 0; i < 12; i++)  // local node (or shape function)
		{
			rRightHandSideVector(i)=Mass_matrix(i,i)*previous_vel(i)/delta_t;
			//KRATOS_WATCH(Mass_matrix(i,i));
			//KRATOS_WATCH(previous_vel(i));
		}
		//KRATOS_WATCH(rRightHandSideVector)
		//KRATOS_WATCH(rLeftHandSideMatrix)
		//KRATOS_WATCH(Mass_matrix)
		
		noalias(rLeftHandSideMatrix) += Mass_matrix/delta_t;   
		
		noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix,previous_vel);
		
		//KRATOS_WATCH("ART")
		//KRATOS_WATCH(rRightHandSideVector)
		//KRATOS_WATCH(rLeftHandSideMatrix)
		
		KRATOS_CATCH("");
	}
	
	
	
	//*****************************************************************************
	//***************************************************************************
	
		//*****************************************************************************
	//***************************************************************************
	
	void PFEM23D::CalculateLocalThermalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
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
	

	//*****************************************************************************
	//***************************************************************************
	

		void PFEM23D::CalculateViscousRHS(ProcessInfo& CurrentProcessInfo)
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
					const array_1d<double,3> zero3 = ZeroVector(3);
					const double nodal_weight = 0.25;			
					
					boost::numeric::ublas::bounded_matrix<double, 12, 12 > Viscosity_matrix = ZeroMatrix(12, 12); 
					boost::numeric::ublas::bounded_matrix<double, 4, 3 > velocities = ZeroMatrix(4, 3);
										
					boost::numeric::ublas::bounded_matrix<double, 4, 3 > DN_DX;
					array_1d<double, 4 > N;
					array_1d<double, 3 > vel_gauss;
					double Area;
					Geometry<Node<3> >& geom = this->GetGeometry();
					GeometryUtils::CalculateGeometryData(geom, DN_DX, N, Area);
					
					if ((this->GetValue(SPLIT_ELEMENT))==true)
					{
						//TOOLS NEEDED TO FIND THE ENRICHMENT SHAPE FUNCTIONS
						//get position of the cut surface
						//const unsigned int LocalSize = geom.size()*2;
						unsigned int TNumNodes = geom.size();
						array_1d<double,4> distances;
						boost::numeric::ublas::bounded_matrix<double,6, 2> Nenriched;
						array_1d<double,(6)> volumes;
						array_1d<double,(6)> densities;
						array_1d<double,(6)> inv_densities;
						boost::numeric::ublas::bounded_matrix<double,4, 3 > coords;
						boost::numeric::ublas::bounded_matrix<double, 6, 4 > Ngauss;
						array_1d<double,(6)> signs;
						std::vector< Matrix > gauss_gradients(6);

						//fill coordinates
					   
						//unsigned int single_triangle_node;
						array_1d<unsigned int, 4 > fixed_nodes; //unordered : i position in the array might not correspond to i node of the element
						array_1d<bool, 4 > is_node_fixed; //i position belongs to i node of the element
						unsigned int number_of_fixed_nodes=0;
						bool boundary_element=false;
						
						for (unsigned int i = 0; i < TNumNodes; i++)
						{
							const array_1d<double, 3 > & xyz = geom[i].Coordinates();
							const array_1d<double, 3 > & velocity = geom[i].FastGetSolutionStepValue(VELOCITY);
							
							volumes[i] = 0.0;
							distances[i] = geom[i].FastGetSolutionStepValue(DISTANCE);
							for (unsigned int j = 0; j < (3); j++)
							{
								coords(i, j) = xyz[j];
								velocities(i,j) = velocity[j];
							}
							
							//to find out stress when we have cut elements:
							if (geom[i].IsFixed(VELOCITY_X) && geom[i].IsFixed(VELOCITY_Y) && geom[i].IsFixed(VELOCITY_Z) )
							{
								fixed_nodes[number_of_fixed_nodes]=i;
								number_of_fixed_nodes++;
								is_node_fixed[i]=true;
							}
								else
									is_node_fixed[i]=false;
							
							
						}

						for (unsigned int i = 0; i < 6 ; i++)
							gauss_gradients[i].resize(2, (3), false);  //2 values of the 2 shape functions, and derivates in xy(z) direction).
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
						
						/*
						AddViscousTerm(Viscosity_matrix, DN_DX,
										viscosity_air, viscosity_water,
										volumes, signs, ndivisions); //, Ngauss,ç
						*/
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
						if (number_of_fixed_nodes==3) //it means we have a surface element
						{
							boundary_element=true;
							array_1d<double, 3 > normal;
							unsigned int free_node=0;

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
							
							
							boundary_stress = - geom[free_node].FastGetSolutionStepValue(VELOCITY)/(pow(plane_point_distance,2)); //have in mind we could not add the viscosity yet since it changes at each node!
							//KRATOS_WATCH(plane_point_distance)
							//KRATOS_WATCH(boundary_stress)
							//KRATOS_WATCH(fixed_face_area_or_lenght)
							
						}


						double current_nodal_mass=0.0;
						
						for (unsigned int j=0; j!=(4); j++)
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
								//NO RHS in these nodes!
								//array_1d<double, 3 > & current_rhs = geom[j].FastGetSolutionStepValue(RHS);
								//current_rhs += boundary_stress*fixed_face_area_or_lenght*0.5*1.0*viscosity;
								
							}
							
							else if (is_node_fixed[j]==false)//normal viscosity procedure:
							{
								geom[j].FastGetSolutionStepValue(NODAL_AREA) += Area * nodal_weight;
							
								array_1d<double, 3 > & current_rhs = geom[j].FastGetSolutionStepValue(RHS);
								for (unsigned int i=0;i!=(4);i++) //neighbour node
									for (unsigned int k=0;k!=(3);k++) //k component of local stress
										for (unsigned int l=0;l!=(3);l++) //affected by l component in i neighbour node
											current_rhs[k] += Viscosity_matrix(j*3+k,i*3+l)*velocities(i,l);
								
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
						
						
						array_1d<unsigned int, 4 > fixed_nodes; //unordered : i position in the array might not correspond to i node of the element
						array_1d<bool,4 > is_node_fixed; //i position belongs to i node of the element
						unsigned int number_of_fixed_nodes=0;
						bool boundary_element=false;
						
						
						for (unsigned int i = 0; i < 4; i++)
						{
							const array_1d<double, 3 > & velocity = geom[i].FastGetSolutionStepValue(VELOCITY);
							for (unsigned int j = 0; j < (3); j++)
								velocities(i,j) = velocity[j];
							

							if (geom[i].IsFixed(VELOCITY_X) && geom[i].IsFixed(VELOCITY_Y) && geom[i].IsFixed(VELOCITY_Z) )
							{
								fixed_nodes[number_of_fixed_nodes]=i;
								number_of_fixed_nodes++;
								is_node_fixed[i]=true;
							}
							else
								is_node_fixed[i]=false;

						}
						
						double plane_point_distance=1.0;
						double fixed_face_area_or_lenght=0.0;
						array_1d<double, 3 > boundary_stress;
						if (number_of_fixed_nodes==3) //it means we have cutted elements!
						{
							boundary_element=true;
							array_1d<double, 3 > normal;
							unsigned int free_node=0;

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

							
							boundary_stress = - geom[free_node].FastGetSolutionStepValue(VELOCITY)*viscosity/(pow(plane_point_distance,2));
							//KRATOS_WATCH(plane_point_distance)
							//KRATOS_WATCH(boundary_stress)
							//KRATOS_WATCH(fixed_face_area_or_lenght)
							
						}
						
						for (unsigned int j=0; j!=(4); j++)
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
								for (unsigned int i=0;i!=(4);i++) //neighbour node
									for (unsigned int k=0;k!=(3);k++) //k component of local stress
										for (unsigned int l=0;l!=(3);l++) //affected by l component in i neighbour node
											current_rhs[k] += Viscosity_matrix(j*3+k,i*3+l)*velocities(i,l);
								
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
	

	void PFEM23D::CalculatePressureProjection(ProcessInfo& CurrentProcessInfo)
	{
		KRATOS_TRY
		
		
		if ((this->GetValue(IS_INACTIVE))==false) //elements can be inactive to add temporary walls. fractional velocity is integrated by parts so walls are seen as having zero velocity without doing anything 
				{
					const double mDENSITY_AIR = CurrentProcessInfo[DENSITY_AIR]; // * (1.0-theta) ;
					const double mDENSITY_WATER = CurrentProcessInfo[DENSITY_WATER]; // * (1.0-theta) ;
					const double delta_t = CurrentProcessInfo[DELTA_TIME]; // * (1.0-theta) ;
					const double mINV_DENSITY_AIR = 1.0/mDENSITY_AIR;
					const double mINV_DENSITY_WATER = 1.0/mDENSITY_WATER;
					const double viscosity_water = CurrentProcessInfo[VISCOSITY_WATER]; // * (1.0-theta) ;
					
					//const double delta_t = CurrentProcessInfo[DELTA_TIME];	
					const array_1d<double,3> zero3 = ZeroVector(3);
					const double mass_factor = 0.25;
					
					boost::numeric::ublas::bounded_matrix<double, 4, 3 > DN_DX;
					array_1d<double, 4 > N;
					array_1d<double, 3 > vel_gauss;
					
					double Area;
					Geometry<Node<3> >& geom = this->GetGeometry();
					//const unsigned int LocalSize = geom.size()*3;
					unsigned int TNumNodes = geom.size();
					//boost::numeric::ublas::bounded_matrix<double, (2+1),1 > enrich_rhs;
					GeometryUtils::CalculateGeometryData(geom, DN_DX, N, Area);
					
					array_1d<double,(4)> nodal_masses = ZeroVector(4);
					array_1d<double,(4)> nodal_viscosities = ZeroVector(4);

					array_1d<double,(4)> pressure;
					array_1d<double,(2+1)> temperatures;
					for (unsigned int i=0; i<(4);i++)
					{
						pressure(i) = geom[i].FastGetSolutionStepValue(PRESSURE);
						temperatures(i) = geom[i].FastGetSolutionStepValue(TEMPERATURE);
					}
						
					
					boost::numeric::ublas::bounded_matrix<double, 4, 12 > G_matrix; //(gradient)
					noalias(G_matrix) = ZeroMatrix(4, 12);	
					//new (upper) part of D, corresponding the first enrichment function (gradient discontinuity)
					boost::numeric::ublas::bounded_matrix<double, 1, 12 > G_matrix_mixed;
					noalias(G_matrix_mixed) = ZeroMatrix(1,12);	
					//lower part of D, the one corresponding to the second enrichmend function (jump)
					boost::numeric::ublas::bounded_matrix<double, 1, 12 > G_matrix_mixed_jump;
					noalias(G_matrix_mixed_jump) = ZeroMatrix(1,12);	
					
					boost::numeric::ublas::bounded_matrix<double, 4, 12 > G_matrix_no_ro; //(gradient)
					noalias(G_matrix_no_ro) = ZeroMatrix(4,12);	
					//new (upper) part of D, corresponding the first enrichment function (gradient discontinuity
					boost::numeric::ublas::bounded_matrix<double, 1, 12 > G_matrix_mixed_no_ro;
					noalias(G_matrix_mixed_no_ro) = ZeroMatrix(1,12);	
					//lower part of D, the one corresponding to the second enrichmend function (jump)
					boost::numeric::ublas::bounded_matrix<double, 1, 12> G_matrix_mixed_jump_no_ro;
					noalias(G_matrix_mixed_jump_no_ro) = ZeroMatrix(1,12);		
					
					
					//for the enrichment:
					array_1d<double,(4)> distances;
					boost::numeric::ublas::bounded_matrix<double,6, 2> Nenriched;
					array_1d<double,(6)> volumes;
					array_1d<double,(6)> densities;
					array_1d<double,(6)> inv_densities;
					array_1d<double,(6)> viscosities;
					boost::numeric::ublas::bounded_matrix<double,(4), 3 > coords;
					boost::numeric::ublas::bounded_matrix<double, 6, 4 > Ngauss;
					array_1d<double,(6)> signs;
					std::vector< Matrix > gauss_gradients(6);
					//fill coordinates
				   
					unsigned int single_triangle_node;

					
					double sqElemSize=0.0;
					//if( ((rCurrentProcessInfo[NL_ITERATION_NUMBER])>1) && (split_element==false) && (neighbour_of_splitted_element==false) )
					if (true)
					{
						//double Viscosity = rCurrentProcessInfo[VISCOSITY];;
						//double AdvVelNorm = sqrt(pow(this->GetGeometry()[0].FastGetSolutionStepValue(VELOCITY_X),2)+pow(this->GetGeometry()[0].FastGetSolutionStepValue(VELOCITY_Y),2));
						//const double TimeFactor = rCurrentProcessInfo.GetValue(DYNAMIC_TAU);
						//double mElemSize;
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
						//mElemSize = sqrt(mElemSize);
						//const double   TauOne = 1.0 / (Density * ( TimeFactor / DeltaTime + 4.0 * Viscosity / (mElemSize * mElemSize) + 2.0 * AdvVelNorm / mElemSize) );
						//actually we're currently using the kinematic viscosity, so no need to divide the second term by the viscosity, it was already done. Great!
						//TauOne =  1.0 / ( ( 1.0/ delta_t + 4.0 * Viscosity / (mElemSize * mElemSize) + 2.0 * AdvVelNorm / mElemSize) );
						
					}
									
						
					
					
					if ((this->GetValue(SPLIT_ELEMENT))==true)
					{
						//to begin with we calculate the gradient discontinuity:
						//if (true) 
						//{
							const int iteration_number =  CurrentProcessInfo[NL_ITERATION_NUMBER];
							double gradient_discontinuity = this->GetValue(ENRICH_RHS) ;
							double gradient_discontinuity_complete_rhs = this->GetValue(GRADIENT_DISCONTINUITY)/((this->GetValue(INV_LAPLACIAN_ENRICH)));
							
							if((this->GetValue(ENRICH_LHS_ROW_3D)).size()!=4)
							{
								KRATOS_ERROR(std::logic_error, "WRONG SIZE in PRESSPROJ", "");
								(this->GetValue(ENRICH_LHS_ROW_3D)).resize(4);
								(this->GetValue(ENRICH_LHS_ROW_3D))=ZeroVector(4);
							}
							//KRATOS_WATCH(this->GetValue(ENRICH_LHS_ROW_3D))
							//KRATOS_WATCH(this->GetValue(ENRICH_RHS))
							
							array_1d<double,(4)> enrich_lhs;
							enrich_lhs = this->GetValue(ENRICH_LHS_ROW_3D);


							array_1d<double,(4)> pressure_array;
							array_1d<double,(4)> previous_iteration_pressure_array;
							for (unsigned int i=0; i<(4);i++)
							{
								pressure_array(i) = geom[i].FastGetSolutionStepValue(PRESSURE);//-geom[i].GetSolutionStepValue(PRESSURE,1);
								previous_iteration_pressure_array(i) = geom[i].FastGetSolutionStepValue(PREVIOUS_ITERATION_PRESSURE);
							}
							for (unsigned int i = 0; i < (4); i++)  // node i
							{ 
								gradient_discontinuity -= pressure_array(i)*enrich_lhs(i);
								if (iteration_number>1)
									gradient_discontinuity_complete_rhs += previous_iteration_pressure_array(i)*enrich_lhs(i);
							}
							//gradient_discontinuity_complete_rhs *= delta_t;
							gradient_discontinuity += gradient_discontinuity_complete_rhs;
							gradient_discontinuity *= this->GetValue(INV_LAPLACIAN_ENRICH);
							//KRATOS_WATCH(this->GetValue(INV_LAPLACIAN_ENRICH))
							this->GetValue(GRADIENT_DISCONTINUITY) = gradient_discontinuity;
						
						


						double input_jump_value = this->GetValue(PRESS_DISCONTINUITY);
						double jump_value=-input_jump_value*0.5;
												//double gradient_discontinuity; // = this->GetValue(ENRICH_RHS);
						//array_1d<double,(2+1)> & enrich_lhs = this->GetValue(ENRICH_LHS_ROW);
					
	
						//TOOLS NEEDED TO FIND THE ENRICHMENT SHAPE FUNCTIONS
						//get position of the cut surface

						for (unsigned int i = 0; i < TNumNodes; i++)
						{
							const array_1d<double, 3 > & xyz = geom[i].Coordinates();
							distances[i] = geom[i].FastGetSolutionStepValue(DISTANCE);
							//KRATOS_WATCH(distances(i));
							for (unsigned int j = 0; j < (3); j++)
								coords(i, j) = xyz[j];
						}

						for (unsigned int i = 0; i < (6) ; i++)
							gauss_gradients[i].resize(2, 3, false);  //2 values of the 2 shape functions, and derivates in (xyz) direction).
						unsigned int ndivisions;
						single_triangle_node=1;
						ndivisions = EnrichmentUtilities::CalculateEnrichedShapeFuncions(coords, DN_DX, distances, volumes, Ngauss, signs, gauss_gradients, Nenriched); //, face_gauss_N, face_gauss_Nenriched, face_Area, face_n  ,single_triangle_node);


						//in case we've changed the distance function:
						for (unsigned int i = 0; i < ndivisions; i++)
						{
							//geom[i].FastGetSolutionStepValue(DISTANCE)=distances[i];
							double partition_temperature=0.0;
							for (unsigned int j = 0; j < (4); j++) //shape function (Nj)
								partition_temperature += temperatures(i)*Ngauss(i,j);
							if (signs[i]>0)
							{
								densities(i) = CalculateAirDensity(partition_temperature);
								viscosities(i) = this->CalculateAirViscosity(partition_temperature);
								inv_densities(i) = 1.0/densities(i);
							}
							else
							{
								densities(i) = mDENSITY_WATER;
								inv_densities(i) = mINV_DENSITY_WATER;
								double viscosity=viscosity_water;
								if (partition_temperature>100.0)
									viscosity=viscosity_water*(1.0 - (partition_temperature-100.0)*0.1 );
								if (viscosity<1.0) //(crude oil viscosity)
									viscosity=1.0;
									
								viscosities(i) = viscosity;
							}
						}

						//KRATOS_WATCH(gradient_discontinuity);
						//KRATOS_WATCH(jump_value)
						if (ndivisions>6)
							KRATOS_WATCH("ERROR!!!!!!!!!! FOUND TOO MANY PARTITIONS");
						//jump_value*=signs((single_triangle_node-1));

						//gradient_discontinuity = this->GetValue(GRADIENT_DISCONTINUITY);

						///*****************************************************************************************
						for (unsigned int i = 0; i < ndivisions; i++)  // partiton i
						{ 
							for (unsigned int j = 0; j < (4); j++) //shape function (Nj)
							{
								for (unsigned int k=0;k!=3;k++) //x,y,(z)
								{
									G_matrix_mixed(0,j*(3)+k) += gauss_gradients[i](0,k) *volumes(i)*Ngauss(i,j)*inv_densities(i);
									G_matrix_mixed_no_ro(0,j*(3)+k) += gauss_gradients[i](0,k) *volumes(i)*Ngauss(i,j);
									

									G_matrix_mixed_jump(0,j*(2)+k) += gauss_gradients[i](1,k) *volumes(i)*Ngauss(i,j)*inv_densities(i);
									G_matrix_mixed_jump_no_ro(0,j*(2)+k) += gauss_gradients[i](1,k) *volumes(i)*Ngauss(i,j);

									
									for (unsigned int l = 0; l < (4); l++) //shape function derivative (dNl/dk)
									{
										G_matrix(l, (j*(3)+k) ) += DN_DX(l,k)*volumes(i)*Ngauss(i,j)*inv_densities(i);
									}
									
								}
								
								nodal_masses(j) += volumes(i)*Ngauss(i,j)*densities(i);
								nodal_viscosities(j) += volumes(i)*Ngauss(i,j)*viscosities(i);
							}
						}
						//for the massless g matrix (to calculate the press proj for the stabilization) we do not need to loop the partitions, so we create a new loop
						for (unsigned int j = 0; j < (4); j++) //shape function (Nj)
						{
							for (unsigned int k=0;k!=3;k++) //x,y,(z)
							{
								for (unsigned int l = 0; l < (4); l++) //shape function derivative (dNl/dk)
								{
									G_matrix_no_ro(l, (j*(3)+k) ) = DN_DX(l,k)*Area*mass_factor;
								}
							}
						}
								
						//gradient_discontinuity = this->GetValue(GRADIENT_DISCONTINUITY);
							
						//now we save the data:
						for (unsigned int i=0; i!=(4); i++) //loop around the nodes of the element to add contribution to node i
						{
							
							array_1d<double, 3 > & current_press_proj_no_ro = geom[i].FastGetSolutionStepValue(PRESS_PROJ_NO_RO);
							array_1d<double, 3 > & current_press_proj = geom[i].FastGetSolutionStepValue(PRESS_PROJ);
								
							for (unsigned int j = 0; j < (4) ; j++) //loop around the nodes of the element to add contribution of node j to node i
							{
								
								for (unsigned int k = 0; k < (3) ; k++) //x,y,(z)
								{
									current_press_proj[k] += G_matrix(j,i*(3)+k)*(pressure(j));///-old_pressures(j)); //gamma=0!
									current_press_proj_no_ro[k] += G_matrix_no_ro(j,i*(3)+k)*(pressure(j));///-old_pressures(j)); //gamma=0!
								}
							}
							//NOT ENRICHENING!!
							/*
							for (unsigned int k = 0; k < (3) ; k++) //x,y,(z)
							{
								current_press_proj[k] += G_matrix_mixed(0,3*i+k)*gradient_discontinuity;// + G_matrix_mixed_jump(0,3*i+k)*jump_value;
								current_press_proj_no_ro[k] += G_matrix_mixed_no_ro(0,3*i+k)*gradient_discontinuity;// + G_matrix_mixed_jump_no_ro(0,3*i+k)*jump_value;
							}
							*/
							geom[i].FastGetSolutionStepValue(NODAL_AREA) += Area*mass_factor;
							geom[i].FastGetSolutionStepValue(NODAL_MASS) += nodal_masses(i)/delta_t+nodal_viscosities(i)/0.04;
							
						}
						
					}
					else //normal (uncut) element
					{
						
						double density;
						double viscosity;
						double element_temperature= 0.25*(temperatures(0)+temperatures(1)+temperatures(2)+temperatures(3));
						if (geom[0].FastGetSolutionStepValue(DISTANCE)<0.0)
						{
							density = mDENSITY_WATER;
							viscosity=viscosity_water;
							if (element_temperature>100.0)
								viscosity=viscosity_water*(1.0 - (element_temperature-100.0)*0.1 );
							if (viscosity<1.0) //(crude oil viscosity)
								viscosity=1.0;
						}
						else
						{
							
							density = CalculateAirDensity(element_temperature);
							viscosity = this->CalculateAirViscosity(element_temperature);

						}
							

						
						for (unsigned int i = 0; i < (4); i++) //loop in shape functions (row)
						{
							for (unsigned int j = 0; j < (4) ; j++) //loop in shape functions (column)
							{
								for (unsigned int k = 0; k < (3) ; k++) //x,y,(z)
								{
								G_matrix(i, (j*3)+k ) = DN_DX(i,k)*Area*mass_factor;
								G_matrix_no_ro(i, (j*3)+k ) = DN_DX(i,k)*Area*mass_factor;
								}
							}
						}
						
						G_matrix /=density;

						//now we save the data:
						for (unsigned int i=0; i!=(4); i++) //loop around the nodes of the element to add contribution to node i
						{
							
							array_1d<double, 3 > & current_press_proj_no_ro = geom[i].FastGetSolutionStepValue(PRESS_PROJ_NO_RO);
							array_1d<double, 3 > & current_press_proj = geom[i].FastGetSolutionStepValue(PRESS_PROJ);
								
							for (unsigned int j = 0; j < (4) ; j++) //loop around the nodes of the element to add contribution of node j to node i
							{
								
								for (unsigned int k = 0; k < (3) ; k++) //x,y,(z)
								{
									current_press_proj[k] += G_matrix(j,i*(3)+k)*(pressure(j));///-old_pressures(j)); //gamma=0!
									current_press_proj_no_ro[k] += G_matrix_no_ro(j,i*(3)+k)*(pressure(j));///-old_pressures(j)); //gamma=0!
								}
							}
							
							geom[i].FastGetSolutionStepValue(NODAL_AREA) += Area*mass_factor;
							geom[i].FastGetSolutionStepValue(NODAL_MASS) += Area*mass_factor*(density/delta_t+viscosity/0.04);
							
						}

					} //closing the normal (uncut) element
					
				} //closing the if(is_inactive==false)

		
		KRATOS_CATCH("");
	}

	//*****************************************************************************
	//***************************************************************************
	
	//*****************************************************************************
	//***************************************************************************
	

	void PFEM23D::CalculateMassMatrix(ProcessInfo& CurrentProcessInfo)
	{
		KRATOS_TRY
		
		
		if ((this->GetValue(IS_INACTIVE))==false) //elements can be inactive to add temporary walls. fractional velocity is integrated by parts so walls are seen as having zero velocity without doing anything 
				{
					const double mDENSITY_AIR = CurrentProcessInfo[DENSITY_AIR]; // * (1.0-theta) ;
					const double mDENSITY_WATER = CurrentProcessInfo[DENSITY_WATER]; // * (1.0-theta) ;
					const double delta_t = CurrentProcessInfo[DELTA_TIME]; // * (1.0-theta) ;
					const double mINV_DENSITY_AIR = 1.0/mDENSITY_AIR;
					const double mINV_DENSITY_WATER = 1.0/mDENSITY_WATER;
					const double viscosity_water = CurrentProcessInfo[VISCOSITY_WATER]; // * (1.0-theta) ;
					
					//const double delta_t = CurrentProcessInfo[DELTA_TIME];	
					const array_1d<double,3> zero3 = ZeroVector(3);
					const double mass_factor = 0.25;
					
					boost::numeric::ublas::bounded_matrix<double, 4, 3 > DN_DX;
					array_1d<double, 4 > N;
					array_1d<double, 3 > vel_gauss;
					
					double Area;
					Geometry<Node<3> >& geom = this->GetGeometry();
					//const unsigned int LocalSize = geom.size()*3;
					unsigned int TNumNodes = geom.size();
					//boost::numeric::ublas::bounded_matrix<double, (2+1),1 > enrich_rhs;
					GeometryUtils::CalculateGeometryData(geom, DN_DX, N, Area);
					
					array_1d<double,(4)> nodal_masses = ZeroVector(4);
					array_1d<double,(4)> nodal_viscosities = ZeroVector(4);

					array_1d<double,(4)> pressure;
					array_1d<double,(2+1)> temperatures;
					for (unsigned int i=0; i<(4);i++)
					{
						pressure(i) = geom[i].FastGetSolutionStepValue(PRESSURE);
						temperatures(i) = geom[i].FastGetSolutionStepValue(TEMPERATURE);
					}
						
					//for the enrichment:
					array_1d<double,(4)> distances;
					boost::numeric::ublas::bounded_matrix<double,6, 2> Nenriched;
					array_1d<double,(6)> volumes;
					array_1d<double,(6)> densities;
					array_1d<double,(6)> inv_densities;
					array_1d<double,(6)> viscosities;
					boost::numeric::ublas::bounded_matrix<double,(4), 3 > coords;
					boost::numeric::ublas::bounded_matrix<double, 6, 4 > Ngauss;
					array_1d<double,(6)> signs;
					std::vector< Matrix > gauss_gradients(6);
					//fill coordinates
				   
					unsigned int single_triangle_node;

					double sqElemSize=0.0;
					//if( ((rCurrentProcessInfo[NL_ITERATION_NUMBER])>1) && (split_element==false) && (neighbour_of_splitted_element==false) )
					if (true)
					{
						//double Viscosity = rCurrentProcessInfo[VISCOSITY];;
						//double AdvVelNorm = sqrt(pow(this->GetGeometry()[0].FastGetSolutionStepValue(VELOCITY_X),2)+pow(this->GetGeometry()[0].FastGetSolutionStepValue(VELOCITY_Y),2));
						//const double TimeFactor = rCurrentProcessInfo.GetValue(DYNAMIC_TAU);
						//double mElemSize;
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
						//mElemSize = sqrt(mElemSize);
						//const double   TauOne = 1.0 / (Density * ( TimeFactor / DeltaTime + 4.0 * Viscosity / (mElemSize * mElemSize) + 2.0 * AdvVelNorm / mElemSize) );
						//actually we're currently using the kinematic viscosity, so no need to divide the second term by the viscosity, it was already done. Great!
						//TauOne =  1.0 / ( ( 1.0/ delta_t + 4.0 * Viscosity / (mElemSize * mElemSize) + 2.0 * AdvVelNorm / mElemSize) );
						
					}
									
						
					
					
					if ((this->GetValue(SPLIT_ELEMENT))==true)
					{

						//TOOLS NEEDED TO FIND THE ENRICHMENT SHAPE FUNCTIONS
						//get position of the cut surface

						for (unsigned int i = 0; i < TNumNodes; i++)
						{
							const array_1d<double, 3 > & xyz = geom[i].Coordinates();
							distances[i] = geom[i].FastGetSolutionStepValue(DISTANCE);
							//KRATOS_WATCH(distances(i));
							for (unsigned int j = 0; j < (3); j++)
								coords(i, j) = xyz[j];
						}

						for (unsigned int i = 0; i < (6) ; i++)
							gauss_gradients[i].resize(2, 3, false);  //2 values of the 2 shape functions, and derivates in (xyz) direction).
						unsigned int ndivisions;
						single_triangle_node=1;
						ndivisions = EnrichmentUtilities::CalculateEnrichedShapeFuncions(coords, DN_DX, distances, volumes, Ngauss, signs, gauss_gradients, Nenriched); //, face_gauss_N, face_gauss_Nenriched, face_Area, face_n  ,single_triangle_node);


						//in case we've changed the distance function:
						for (unsigned int i = 0; i < ndivisions; i++)
						{
							//geom[i].FastGetSolutionStepValue(DISTANCE)=distances[i];
							double partition_temperature=0.0;
							for (unsigned int j = 0; j < (4); j++) //shape function (Nj)
								partition_temperature += temperatures(i)*Ngauss(i,j);
							if (signs[i]>0)
							{
								densities(i) = CalculateAirDensity(partition_temperature);
								viscosities(i) = this->CalculateAirViscosity(partition_temperature);
								inv_densities(i) = 1.0/densities(i);
							}
							else
							{
								densities(i) = mDENSITY_WATER;
								inv_densities(i) = mINV_DENSITY_WATER;
								double viscosity=viscosity_water;
								if (partition_temperature>100.0)
									viscosity=viscosity_water*(1.0 - (partition_temperature-100.0)*0.1 );
								if (viscosity<1.0) //(crude oil viscosity)
									viscosity=1.0;
									
								viscosities(i) = viscosity;
							}
						}

						//KRATOS_WATCH(gradient_discontinuity);
						//KRATOS_WATCH(jump_value)
						if (ndivisions>6)
							KRATOS_WATCH("ERROR!!!!!!!!!! FOUND TOO MANY PARTITIONS");
						//jump_value*=signs((single_triangle_node-1));

						//gradient_discontinuity = this->GetValue(GRADIENT_DISCONTINUITY);

						///*****************************************************************************************
						for (unsigned int i = 0; i < ndivisions; i++)  // partiton i
						{ 
							for (unsigned int j = 0; j < (4); j++) //shape function (Nj)
							{
								nodal_masses(j) += volumes(i)*Ngauss(i,j)*densities(i);
								nodal_viscosities(j) += volumes(i)*Ngauss(i,j)*viscosities(i);
							}
						}
								
						//gradient_discontinuity = this->GetValue(GRADIENT_DISCONTINUITY);
							
						//now we save the data:
						for (unsigned int i=0; i!=(4); i++) //loop around the nodes of the element to add contribution to node i
						{
							geom[i].FastGetSolutionStepValue(NODAL_AREA) += nodal_masses(i)/delta_t+nodal_viscosities(i)/0.04;	
							geom[i].FastGetSolutionStepValue(NODAL_MASS) += nodal_masses(i)/delta_t;
						}
						
					}
					else //normal (uncut) element
					{
						
						double density;
						double viscosity;
						double element_temperature= 0.25*(temperatures(0)+temperatures(1)+temperatures(2)+temperatures(3));
						if (geom[0].FastGetSolutionStepValue(DISTANCE)<0.0)
						{
							density = mDENSITY_WATER;
							viscosity=viscosity_water;
							if (element_temperature>100.0)
								viscosity=viscosity_water*(1.0 - (element_temperature-100.0)*0.1 );
							if (viscosity<1.0) //(crude oil viscosity)
								viscosity=1.0;
						}
						else
						{
							
							density = CalculateAirDensity(element_temperature);
							viscosity = this->CalculateAirViscosity(element_temperature);

						}
							


						//now we save the data:
						for (unsigned int i=0; i!=(4); i++) //loop around the nodes of the element to add contribution to node i
						{
							
							geom[i].FastGetSolutionStepValue(NODAL_AREA) += Area*mass_factor*(density/delta_t+viscosity/0.04);
							geom[i].FastGetSolutionStepValue(NODAL_MASS) += Area*mass_factor*(density/delta_t);	
						}

					} //closing the normal (uncut) element
					
				} //closing the if(is_inactive==false)

		
		KRATOS_CATCH("");
	}

	//*****************************************************************************
	//***************************************************************************	
	


	//************************************************************************************
	//************************************************************************************
	// this subroutine calculates the nodal contributions for the explicit steps of the 
	// fractional step procedure
	void PFEM23D::InitializeSolutionStep(ProcessInfo& CurrentProcessInfo)
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
	void PFEM23D::PressureEquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo)
	{
		unsigned int number_of_nodes = GetGeometry().PointsNumber();
		if(rResult.size() != number_of_nodes)
			rResult.resize(number_of_nodes,false);	

		for (unsigned int i=0;i<number_of_nodes;i++)
				rResult[i] = GetGeometry()[i].GetDof(PRESSURE).EquationId();
	}
	
	void PFEM23D::ThermalEquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo)
	{
		unsigned int number_of_nodes = GetGeometry().PointsNumber();
		if(rResult.size() != number_of_nodes)
			rResult.resize(number_of_nodes,false);	

		for (unsigned int i=0;i<number_of_nodes;i++)
				rResult[i] = GetGeometry()[i].GetDof(TEMPERATURE).EquationId();
	}
	
	//************************************************************************************
	void PFEM23D::VelocityEquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo)
	{
		const SizeType NumNodes = 4;
		const SizeType LocalSize = 12;
		GeometryType& rGeom = this->GetGeometry();

		SizeType LocalIndex = 0;

		if (rResult.size() != LocalSize)
			rResult.resize(LocalSize, false);

		for (SizeType i = 0; i < NumNodes; ++i)
		{
			rResult[LocalIndex++] = rGeom[i].GetDof(VELOCITY_X).EquationId();
			rResult[LocalIndex++] = rGeom[i].GetDof(VELOCITY_Y).EquationId();
			rResult[LocalIndex++] = rGeom[i].GetDof(VELOCITY_Z).EquationId();
		}
	}
		//************************************************************************************
	void PFEM23D::FractionalVelocityEquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo)
	{
		const SizeType NumNodes = 4;
		const SizeType LocalSize = 12;
		GeometryType& rGeom = this->GetGeometry();

		SizeType LocalIndex = 0;

		if (rResult.size() != LocalSize)
			rResult.resize(LocalSize, false);

		for (SizeType i = 0; i < NumNodes; ++i)
		{
			rResult[LocalIndex++] = rGeom[i].GetDof(FRACT_VEL_X).EquationId();
			rResult[LocalIndex++] = rGeom[i].GetDof(FRACT_VEL_Y).EquationId();
			rResult[LocalIndex++] = rGeom[i].GetDof(FRACT_VEL_Z).EquationId();
		}
	}

	//************************************************************************************
	//************************************************************************************
	void PFEM23D::GetPressureDofList(DofsVectorType& ElementalDofList,ProcessInfo& CurrentProcessInfo)
	{
		unsigned int number_of_nodes = GetGeometry().PointsNumber();
		if(ElementalDofList.size() != number_of_nodes)
			ElementalDofList.resize(number_of_nodes);	

		for (unsigned int i=0;i<number_of_nodes;i++)
			ElementalDofList[i] = GetGeometry()[i].pGetDof(PRESSURE);

	}
	
	void PFEM23D::GetThermalDofList(DofsVectorType& ElementalDofList,ProcessInfo& CurrentProcessInfo)
	{
		unsigned int number_of_nodes = GetGeometry().PointsNumber();
		if(ElementalDofList.size() != number_of_nodes)
			ElementalDofList.resize(number_of_nodes);	

		for (unsigned int i=0;i<number_of_nodes;i++)
			ElementalDofList[i] = GetGeometry()[i].pGetDof(TEMPERATURE);

	}

	//************************************************************************************
	void PFEM23D::GetVelocityDofList(DofsVectorType& ElementalDofList,ProcessInfo& CurrentProcessInfo)
	{
		const SizeType NumNodes = 4;
		const SizeType LocalSize = 12;
		GeometryType& rGeom = this->GetGeometry();

		if (ElementalDofList.size() != LocalSize)
			ElementalDofList.resize(LocalSize);

		SizeType LocalIndex = 0;

		for (SizeType i = 0; i < NumNodes; ++i)
		{
			ElementalDofList[LocalIndex++] = rGeom[i].pGetDof(VELOCITY_X);
			ElementalDofList[LocalIndex++] = rGeom[i].pGetDof(VELOCITY_Y);
			ElementalDofList[LocalIndex++] = rGeom[i].pGetDof(VELOCITY_Z);
		}
	}
	
		//************************************************************************************
	void PFEM23D::GetFractionalVelocityDofList(DofsVectorType& ElementalDofList,ProcessInfo& CurrentProcessInfo)
	{
		const SizeType NumNodes = 4;
		const SizeType LocalSize = 12;
		GeometryType& rGeom = this->GetGeometry();

		if (ElementalDofList.size() != LocalSize)
			ElementalDofList.resize(LocalSize);

		SizeType LocalIndex = 0;

		for (SizeType i = 0; i < NumNodes; ++i)
		{
			ElementalDofList[LocalIndex++] = rGeom[i].pGetDof(FRACT_VEL_X);
			ElementalDofList[LocalIndex++] = rGeom[i].pGetDof(FRACT_VEL_Y);
			ElementalDofList[LocalIndex++] = rGeom[i].pGetDof(FRACT_VEL_Z);
		}
	}
	
	//************************************************************************************
	
	//with bounded matrix, constant coefficient
	void PFEM23D::AddViscousTerm(boost::numeric::ublas::bounded_matrix<double, 12, 12 >& rDampMatrix,
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
	
	//with matrixtype, constant coefficient
	void PFEM23D::AddViscousTerm(MatrixType& rDampMatrix,
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
		rDampMatrix = prod(B_matrix, temp_matrix );
		
		//KRATOS_WATCH(C_matrix)
		//KRATOS_WATCH(B_matrix)
		//KRATOS_WATCH(rDampMatrix)
		//		KRATOS_ERROR(std::logic_error, "IMPLICIT STEP FIRST STEP NOT YET IMPLEMENTED IN 3D.. USE LOCALFINALVELOCITY", "");

		
	}
	
	
	double PFEM23D::CalculateAirDensity(double & temperature_in_celsius)
	{
		return 1.0;
		//return (352.96443/(273.15+temperature_in_celsius) ) ;
	}
	
	double PFEM23D::CalculateAirConductivity(double & temperature_in_celsius)
	{
		const double temperature_in_kelvin = 273.15+temperature_in_celsius;
		return (1.5207e-11 * pow(temperature_in_kelvin,3) - 4.8574e-8 * pow(temperature_in_kelvin,2) + 1.0184e-4 *temperature_in_kelvin - 0.00039333);
	}
	
	double PFEM23D::CalculateAirViscosity(double & temperature_in_celsius)
	{
		return (0.0032/1.000);
		//const double temperature_in_kelvin = 273.15+temperature_in_celsius;
		//return (0.0000000000097838 * pow(temperature_in_kelvin,2) + 0.000000046646 *temperature_in_kelvin + 0.0000051948);
		//return 10000.0;
	}
	
	//************************************************************************************

} // Namespace Kratos
