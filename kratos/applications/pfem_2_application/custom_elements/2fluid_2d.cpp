//   
//   Project Name:        Kratos       
//   Last modified by:    $Author: it's me! $
//   Date:                $Date: 2008-08-08 23:58:38 $
//   Revision:            $Revision: 1.0 $
//
//

// Project includes 
#include "includes/define.h"
#include "custom_elements/2fluid_2d.h"
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
	PFEM22D::PFEM22D(IndexType NewId, GeometryType::Pointer pGeometry)
		: Element(NewId, pGeometry)
	{		
		//DO NOT ADD DOFS HERE!!!
	}

	//************************************************************************************
	//************************************************************************************
	PFEM22D::PFEM22D(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
		: Element(NewId, pGeometry, pProperties)
	{

	}

	Element::Pointer PFEM22D::Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const
	{
		return Element::Pointer(new PFEM22D(NewId, GetGeometry().Create(ThisNodes), pProperties));
	}

	PFEM22D::~PFEM22D()
	{
	}

	//************************************************************************************
	//************************************************************************************
	
	
	void PFEM22D::AddExplicitContribution(ProcessInfo& rCurrentProcessInfo)
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
	
	void PFEM22D::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix,
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
		default:
		{
			KRATOS_ERROR(std::logic_error,"Unexpected value for FRACTIONAL_STEP index: ",rCurrentProcessInfo[FRACTIONAL_STEP]);
		}
		}

		KRATOS_CATCH("");
	}
	
	//************************************************************************************
	//************************************************************************************
	
	void PFEM22D::EquationIdVector(EquationIdVectorType& rResult,
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
		default:
		{
			KRATOS_ERROR(std::logic_error,"Unexpected value for FRACTIONAL_STEP index: ",rCurrentProcessInfo[FRACTIONAL_STEP]);
		}
		}

		KRATOS_CATCH("");
	}
	
	//************************************************************************************
	//************************************************************************************
	

	void PFEM22D::GetDofList(DofsVectorType& rElementalDofList,
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
		default:
		{
			KRATOS_ERROR(std::logic_error,"Unexpected value for FRACTIONAL_STEP index: ",rCurrentProcessInfo[FRACTIONAL_STEP]);
		}
		}

		KRATOS_CATCH("");
	}
	
	
	
	
	//*************************************************************************************
	//*************************************************************************************
	
	
	
	void PFEM22D::CalculateLocalFractionalVelocitySystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY
		array_1d<double,3> gravity= rCurrentProcessInfo[GRAVITY];
		//WE MUST CALCULATE Mass and G(or D) matrixes

	    double delta_t = rCurrentProcessInfo[DELTA_TIME];	
		//KRATOS_WATCH(gradient_discontinuity);
		const bool & split_element = (this->GetValue(SPLIT_ELEMENT));
		const int nl_iteration_number =  rCurrentProcessInfo[NL_ITERATION_NUMBER];
		double theta=0.0;

	
	
		//array_1d<double,3> pressures;

		
		array_1d<double,6>  previous_vel;
        for(unsigned int iii = 0; iii<3; iii++)
        {
			previous_vel(iii*2) = GetGeometry()[iii].GetSolutionStepValue(MESH_VELOCITY_X,0);//+x_force*delta_t;
			previous_vel(iii*2+1) = GetGeometry()[iii].GetSolutionStepValue(MESH_VELOCITY_Y,0);//+y_force*delta_t;
			
			if (GetGeometry()[iii].IsFixed(FRACT_VEL_X))
				previous_vel(iii*2) = GetGeometry()[iii].GetSolutionStepValue(VELOCITY_X,0);
			if (GetGeometry()[iii].IsFixed(FRACT_VEL_Y))
				previous_vel(iii*2+1) = GetGeometry()[iii].GetSolutionStepValue(VELOCITY_Y,0);
			
				
		}
		
		//boost::numeric::ublas::bounded_matrix<double,3,2> msDN_DX;
		//boost::numeric::ublas::bounded_matrix<double,2,2> msD;
		const double one_third = 1.0/3.0;
  		array_1d<double,3> msN; //dimension = number of nodes
		array_1d<double,6> ms_temp; //dimension = number of nodes

		const unsigned int LocalSize = GetGeometry().size()*2;
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
        boost::numeric::ublas::bounded_matrix<double, 3, 6 > D_matrix; //(gradient)
        noalias(D_matrix) = ZeroMatrix(3,6);	

        boost::numeric::ublas::bounded_matrix<double, 3*2, 3*2 > Mass_matrix;
        noalias(Mass_matrix) = ZeroMatrix(LocalSize,LocalSize);	
        //boost::numeric::ublas::bounded_matrix<double, TNumNodes*2, 2> N_vel_matrix; //useless to calculate it, it's simply the DN_DX matrix multiplied by 1/3, arranged in a different shape.
        
        GeometryUtils::CalculateGeometryData(this->GetGeometry(), DN_DX, N, Area);
         
        array_1d<double,3> distances;
        for (unsigned int i = 0; i < TNumNodes; i++)
        {
            distances[i] = this->GetGeometry()[i].FastGetSolutionStepValue(DISTANCE);
        }
        
        
		const double viscosity_air = rCurrentProcessInfo[VISCOSITY_AIR]; //para clindro en rey 100:  0.0005 * 2.0*Area;
		const double viscosity_water = rCurrentProcessInfo[VISCOSITY_WATER];
		
		const double density_air = rCurrentProcessInfo[DENSITY_AIR]; //para clindro en rey 100:  0.0005 * 2.0*Area;
		const double density_water = rCurrentProcessInfo[DENSITY_WATER];
		
		if (split_element==false)
		{
			//viscous term:
			if (distances(0)<0.0)
				this->AddViscousTerm(rLeftHandSideMatrix,DN_DX, (viscosity_water*Area) );
			else
				this->AddViscousTerm(rLeftHandSideMatrix,DN_DX, (viscosity_air*Area) );
			
			double fourier=0.0;
			if (distances(0)<0.0) fourier = viscosity_water/density_water *delta_t/pow((0.6*this->GetValue(MEAN_SIZE)),2);
			else fourier = viscosity_air/density_air *delta_t/pow((0.6*this->GetValue(MEAN_SIZE)),2);
			theta = 0.2*fourier - 0.1;
			if (theta<0.0) theta=0.0;
			else if (theta>1.0) theta=1.0;
				
			//mass matrix:
			double density;
			if (distances(0)<0.0)
				density=density_water;
			else
				density=density_air;
				
			const double Weight = one_third * Area * density;
			
			for (unsigned int i=0; i<LocalSize ; i++)
				Mass_matrix (i,i) = Weight;
				
			//rLeftHandSideMatrix += Mass_matrix;	
			
		} 
		else //split element:
		{
			//TOOLS NEEDED TO FIND THE ENRICHMENT SHAPE FUNCTIONS
			//get position of the cut surface
			array_1d<double,3>  densities(3);
			boost::numeric::ublas::bounded_matrix<double,3, 2> Nenriched;
			array_1d<double,3>  volumes(3);
			boost::numeric::ublas::bounded_matrix<double,3, 2 > coords;
			boost::numeric::ublas::bounded_matrix<double,3, 3 > Ngauss;
			array_1d<double,3>  signs(3);
			std::vector< Matrix > gauss_gradients(3);
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
			for (unsigned int i = 0; i < 3; i++)
				gauss_gradients[i].resize(2, 2, false);  //2 values of the 2 shape functions, and derivates in (xy) direction).
			unsigned int ndivisions = EnrichmentUtilities::CalculateEnrichedShapeFuncions(coords, DN_DX, distances, volumes, Ngauss, signs, gauss_gradients, Nenriched);
			
			
			
			
			
			double Weight=0.0;
			double mass=0.0;
			for (unsigned int i=0;i<ndivisions;i++) //we go over the three partitions;
			{
				if (signs(i)>0.0)
				{
					Weight += volumes(i)*viscosity_air;
					mass += volumes(i)*density_air;
				}
				else
				{
					Weight += volumes(i)*viscosity_water;
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
				


			
			for (unsigned int i = 0; i < 3; i++)  //partition
			{
				if (signs(i)<0.0)
					densities(i)=density_water;
				else
					densities(i)=density_air;
			}
			
			
			for (unsigned int i = 0; i < 3; i++)  // local node (or shape function)
			{
				//KRATOS_WATCH(enrich_lhs(i));
				for (unsigned int j = 0; j < ndivisions; j++) //partitions
				{
					Mass_matrix(2*i,2*i) +=  volumes(j)*Ngauss(j,i)*densities(j);	 
				}
				Mass_matrix(2*i+1,2*i+1) = Mass_matrix(2*i,2*i);
			}
			/*
			const double Weight = one_third * Area * density_water;
			
			for (unsigned int i=0; i<LocalSize ; i++)
				Mass_matrix (i,i) = Weight;		
			*/
			
		}
       
		
		//we must actually use 1/2 of the rigidity matrix. since the Lefthandside only has it so far, we multiply it by 0.5:
		//but if it is a first order scheme, then we do not do that:
		if(nl_iteration_number!=0)
			rLeftHandSideMatrix *= 0.5;
		//now we save the mass matrix into the lefthandside
		//noalias(rLeftHandSideMatrix) += Mass_matrix/delta_t;  
	
        //finally we add to the righthandside BOTH the mass matrix and the rigidity matrix:
        //noalias(rRightHandSideVector) = prod(rLeftHandSideMatrix,previous_vel);
        noalias(rRightHandSideVector) = prod((Mass_matrix/delta_t),previous_vel);
        //noalias(rRightHandSideVector) += 
        
        for(unsigned int iii = 0; iii<3; iii++)
        {
			previous_vel(iii*2) = GetGeometry()[iii].GetSolutionStepValue(VELOCITY_X,0);
			previous_vel(iii*2+1) = GetGeometry()[iii].GetSolutionStepValue(VELOCITY_Y,0);
			/*
			if (GetGeometry()[iii].IsFixed(FRACT_VEL_X))
				previous_vel(iii*2) = 0.0;
			if (GetGeometry()[iii].IsFixed(FRACT_VEL_Y))
				previous_vel(iii*2+1) = 0.0;
				*/ 
		}
        
		if(nl_iteration_number!=0)
			rRightHandSideVector += prod(rLeftHandSideMatrix,previous_vel); 			//since now we only add the CHANGE in velocity, for the moment we leave it here, then we substract it
		
		rLeftHandSideMatrix *= theta;
		noalias(rLeftHandSideMatrix) += Mass_matrix/delta_t;   
		
		/*
		for (unsigned int i = 0; i < 3 ; i++)
		{
			rRightHandSideVector(2*i+0) += Mass_matrix(2*i+0,2*i+0)*x_force;
			rRightHandSideVector(2*i+1) += Mass_matrix(2*i+1,2*i+1)*y_force;
		}
        */
        

		//subtracting the dirichlet term
		// RHS -= LHS*DUMMY_UNKNOWNs
		
		for(unsigned int iii = 0; iii<3; iii++)
		{
			ms_temp[iii*2] = GetGeometry()[iii].FastGetSolutionStepValue(FRACT_VEL_X);        
			ms_temp[iii*2+1] = GetGeometry()[iii].FastGetSolutionStepValue(FRACT_VEL_Y);
		}
		noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix,ms_temp);
		
		if (theta>0.05)
			KRATOS_WATCH(theta)
		
		KRATOS_CATCH("");
	}
	
	
	
	//************************************************************************************
	//************************************************************************************

	void PFEM22D::CalculateLocalPressureSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY

		//rCurrentProcessInfo[DELTA_TIME]=1.0;
		//ProcessInfo& CurrentProcessInfo = mr_model_part.GetProcessInfo();
		double Density = 0.0 ; //used for stabilization;
		double delta_t = rCurrentProcessInfo[DELTA_TIME];
		int iteration_number =  rCurrentProcessInfo[NL_ITERATION_NUMBER];
		Vector densities(3);
		Vector partition_densities(3);
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
        
        for(unsigned int iii = 0; iii<3; iii++)
        {
			array_1d<double,3> velocity = GetGeometry()[iii].FastGetSolutionStepValue(VELOCITY);
			fract_vel(iii*2,0) = velocity[0];
			fract_vel(iii*2+1,0) = velocity[1];
			//old_pressures(iii) = GetGeometry()[iii].GetSolutionStepValue(PRESSURE,1);
			old_pressures(iii) = GetGeometry()[iii].FastGetSolutionStepValue(PREVIOUS_ITERATION_PRESSURE);
			//if ((GetGeometry()[iii].GetValue(SPLIT_ELEMENT))==true)
			//	neighbour_of_splitted_element=true;
		}
		const array_1d<double, 3 > & proj0 = GetGeometry()[0].FastGetSolutionStepValue(PRESS_PROJ);
		const array_1d<double, 3 > & proj1 = GetGeometry()[1].FastGetSolutionStepValue(PRESS_PROJ);
		const array_1d<double, 3 > & proj2 = GetGeometry()[2].FastGetSolutionStepValue(PRESS_PROJ);


        //we start by calculating the non-enriched functions.
        for (unsigned int i = 0; i < 3; i++)
		{
			for (unsigned int j = 0; j < 3; j++)
			{
				D_matrix(i, (j*2) ) =  -DN_DX(i,0)*Area*one_third;     
				D_matrix(i, (j*2+1) ) =  -DN_DX(i,1)*Area*one_third;
			}
		}
		
		//CALCULATING TAU (only if we are not near the interfase and we are in step 1 or more)
		double TauOne = 0.0; 
		
		if( ((rCurrentProcessInfo[NL_ITERATION_NUMBER])>1) && (split_element==false) )
		{
			double Viscosity = rCurrentProcessInfo[VISCOSITY];;
			double AdvVelNorm = sqrt(pow(this->GetGeometry()[0].FastGetSolutionStepValue(VELOCITY_X),2)+pow(this->GetGeometry()[0].FastGetSolutionStepValue(VELOCITY_Y),2));
			//const double TimeFactor = rCurrentProcessInfo.GetValue(DYNAMIC_TAU);
			double mElemSize;
			array_1d<double,3> Edge(3,0.0);
			Edge = this->GetGeometry()[1].Coordinates() - this->GetGeometry()[0].Coordinates();
			mElemSize = Edge[0]*Edge[0];
			for (SizeType d = 1; d < 2; d++)
				mElemSize += Edge[d]*Edge[d];
			//
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
			//const double   TauOne = 1.0 / (Density * ( TimeFactor / DeltaTime + 4.0 * Viscosity / (mElemSize * mElemSize) + 2.0 * AdvVelNorm / mElemSize) );
			//actually we're currently using the kinematic viscosity, so no need to divide the second term by the viscosity, it was already done. Great!
			TauOne =  1.0 / ( ( 1.0/ delta_t + 4.0 * Viscosity / (mElemSize * mElemSize) + 2.0 * AdvVelNorm / mElemSize) );
		}
		//TauOne=0.0;
		//KRATOS_WATCH(TauOne)
		
		

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
			array_1d<double,3> distances;
			array_1d<double,3>  densities(3);
			boost::numeric::ublas::bounded_matrix<double,3, 2> Nenriched;
			array_1d<double,3>  volumes(3);
			boost::numeric::ublas::bounded_matrix<double,3, 2 > coords;
			boost::numeric::ublas::bounded_matrix<double,3, 3 > Ngauss;
			array_1d<double,3>  signs(3);
			std::vector< Matrix > gauss_gradients(3);
			//fill coordinates
		   

			
			//unsigned int single_triangle_node = 1;
			for (unsigned int i = 0; i < TNumNodes; i++)
			{
				const array_1d<double, 3 > & xyz = this->GetGeometry()[i].Coordinates();
				volumes[i] = 0.0;
				distances[i] = this->GetGeometry()[i].FastGetSolutionStepValue(DISTANCE);
				//KRATOS_WATCH(distances(i));
				for (unsigned int j = 0; j < 2; j++)
					coords(i, j) = xyz[j];
			}

			for (unsigned int i = 0; i < 3; i++)
				gauss_gradients[i].resize(2, 2, false);  //2 values of the 2 shape functions, and derivates in (xy) direction).
			unsigned int ndivisions = EnrichmentUtilities::CalculateEnrichedShapeFuncions(coords, DN_DX, distances, volumes, Ngauss, signs, gauss_gradients, Nenriched);
			
			//in case we've changed the distance function: CANNOT: (PARALLEL)
			//for (unsigned int i = 0; i < TNumNodes; i++)
			//	this->GetGeometry()[i].FastGetSolutionStepValue(DISTANCE)=distances[i];
			
			
			
			array_1d<double,3> inv_densities;
			Density=0.0; //resetting to zero
			for (unsigned int i = 0; i!=ndivisions; i++) //the 3 partitions of the triangle
			{
				if (signs(i)>0)
					partition_densities(i)=rCurrentProcessInfo[DENSITY_AIR];
				else
					partition_densities(i)=rCurrentProcessInfo[DENSITY_WATER];			
				inv_densities(i)=1.0/partition_densities(i);	
				Density+=partition_densities(i)*volumes(i);
			}
			Density*=1.0/Area;


			jump_value*=signs(0);
			
			
			Laplacian_matrix =  - prod(DN_DX,trans(DN_DX)) *(volumes(0)*inv_densities(0)+volumes(1)*inv_densities(1)+volumes(2)*inv_densities(2));

			//and now the rest of the things needed
			
			//boost::numeric::ublas::bounded_matrix<double, 2, 2> DN_DX_enrich;
			//boost::numeric::ublas::bounded_matrix<double, 2, 2> DN_DX_enrich_negative;
			boost::numeric::ublas::bounded_matrix<double, 2, 2 > Laplacian_enrich;
			noalias(Laplacian_enrich) = ZeroMatrix(2,2);	
			//double inv_Laplacian_enrich_weighted;
			//double inv_Laplacian_enrich;
			boost::numeric::ublas::bounded_matrix<double, 1, 3 > mixed_Laplacian;
			noalias(mixed_Laplacian) = ZeroMatrix(1,3);	
			
			boost::numeric::ublas::bounded_matrix<double, 1, 3 > mixed_Laplacian_jump;
			noalias(mixed_Laplacian_jump) = ZeroMatrix(1,3);
			
			
			//To calculate D* =int (N*DN_DX_enrich).  we must use the N at the gauss point of the partition (returned by the enrichment util) and multiply by the area. 
			//notice that this first row is only the upper part of the mixed D. we also need the lower one for the jump 
			boost::numeric::ublas::bounded_matrix<double, 1, 6 > D_matrix_mixed;
			noalias(D_matrix_mixed) = ZeroMatrix(1,6);	
			//lower part of D, the one corresponding to the second enrichmend function (jump)
			boost::numeric::ublas::bounded_matrix<double, 1, 6 > D_matrix_mixed_jump;
			noalias(D_matrix_mixed_jump) = ZeroMatrix(1,6);	
			
				
			
			//the matrices we need at the end will be:
			//boost::numeric::ublas::bounded_matrix<double, 2, 2 > D_mod;

			for (unsigned int i = 0; i < ndivisions; i++)  //we go through the 3 partitions of the domain
			{
				Laplacian_enrich(0,0) -= (pow(gauss_gradients[i](0,0),2)+pow(gauss_gradients[i](0,1),2))*volumes(i)*inv_densities(i);
				Laplacian_enrich(0,1) -= ((gauss_gradients[i](0,0)*gauss_gradients[i](1,0))+(gauss_gradients[i](0,1)*gauss_gradients[i](1,1)))*volumes(i)*inv_densities(i);
				Laplacian_enrich(1,1) -= (pow(gauss_gradients[i](1,0),2)+pow(gauss_gradients[i](1,1),2))*volumes(i)*inv_densities(i);
				Laplacian_enrich(1,0) -= ((gauss_gradients[i](0,0)*gauss_gradients[i](1,0))+(gauss_gradients[i](0,1)*gauss_gradients[i](1,1)))*volumes(i)*inv_densities(i);
				//and the 'mixed laplacians' (standard shape functions * enrichments)
				for (unsigned int j = 0; j < 3; j++) //we go through the 3 standard shape functions
				{
					mixed_Laplacian(0,j) -= DN_DX(j,0)*(gauss_gradients[i](0,0)*volumes(i)*inv_densities(i))+ DN_DX(j,1)*(gauss_gradients[i](0,1)*volumes(i)*inv_densities(i));
					mixed_Laplacian_jump(0,j) -= DN_DX(j,0)*(gauss_gradients[i](1,0)*volumes(i)*inv_densities(i))+ DN_DX(j,1)*(gauss_gradients[i](1,1)*volumes(i)*inv_densities(i));
					//and also the D matrixes

					D_matrix_mixed(0,j*2) -= gauss_gradients[i](0,0) *volumes(i)*Ngauss(i,j);
					D_matrix_mixed(0,j*2+1) -= gauss_gradients[i](0,1) *volumes(i)*Ngauss(i,j);
					D_matrix_mixed_jump(0,j*2) -= gauss_gradients[i](1,0) *volumes(i)*Ngauss(i,j);
					D_matrix_mixed_jump(0,j*2+1) -= gauss_gradients[i](1,1) *volumes(i)*Ngauss(i,j);
					
				}
				
				
			}
			
			//the 'jump' degree of freedom will not be condensed since it's imposed
			//therefore we'll only have to calculate the inverse of the reduced matrix of the gradient jump enrich function. which is a 1x1 matrix (inverse of the 1,1 position of the enriched laplacian)
			const double inv_Laplacian_enrich_weighted=1.0/Laplacian_enrich(0,0);
			this->SetValue(INV_LAPLACIAN_ENRICH,(inv_Laplacian_enrich_weighted/(delta_t+TauOne)));
			
			
			//since we're imposing this DOF, we'll do:
			for (unsigned int i = 0; i < 3; i++)
				rRightHandSideVector(i) -= mixed_Laplacian_jump(0,i)*jump_value;
				
			//when there's no jump there's nothing in the RHS to condense. but when we add it, it is equal to -Laplacian_enrich*jump_value due to the elimination of this DOF.
			//and just like the LHS, we must condense it and substrac it (actually it's (-)*(-)=+ )
			enrich_rhs = -Laplacian_enrich(0,1)*jump_value;
			
			for (unsigned int i = 0; i < 6; i++)
				extra_enrich_rhs +=D_matrix_mixed(0,i)*fract_vel(i,0) ;
			
			
			///gamma=1;
			//for (unsigned int i = 0; i < 3; i++)
			//	enrich_rhs += mixed_Laplacian(0,i) *old_pressures(i);	
				
			for (unsigned int i = 0; i < 3; i++)
				rRightHandSideVector(i) -= mixed_Laplacian(0,i) *inv_Laplacian_enrich_weighted * (enrich_rhs + extra_enrich_rhs);
			
			
			
			for (unsigned int i = 0; i < 3; i++)
			{
				(this->GetValue(ENRICH_LHS_ROW))(i)=mixed_Laplacian(0,i)*(delta_t+TauOne);
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
			double inv_density;
			if (geom[0].FastGetSolutionStepValue(DISTANCE)<0.0)
				inv_density = 1.0/rCurrentProcessInfo[DENSITY_WATER];
			else
				inv_density = 1.0/rCurrentProcessInfo[DENSITY_AIR];
			Laplacian_matrix =  -prod(DN_DX,trans(DN_DX))*Area*inv_density;  // B B^T  (standard laplacian, now we must add the contribution of the condensed new nodes.
		}		
        

        //done, now we save it into the LHS
        noalias(rLeftHandSideMatrix) = Laplacian_matrix*(delta_t+TauOne);///gamma=1
        boost::numeric::ublas::bounded_matrix<double, 3,1 > temp_rhs;
        noalias(temp_rhs) = prod(D_matrix,fract_vel);
        //KRATOS_WATCH(temp_rhs);
        rRightHandSideVector(0)  += temp_rhs(0,0);
        rRightHandSideVector(1)  += temp_rhs(1,0);
        rRightHandSideVector(2)  += temp_rhs(2,0);

        
        //stabilization contribution:
        array_1d<double, 2 > vel_gauss;
        vel_gauss[0] = N[0] * proj0[0] + N[1] * proj1[0] + N[2] * proj2[0];
		vel_gauss[1] = N[0] * proj0[1] + N[1] * proj1[1] + N[2] * proj2[1];
		vel_gauss *= TauOne*Area;
		//KRATOS_WATCH(vel_gauss);
		///gammma=1!
		if (rCurrentProcessInfo[NL_ITERATION_NUMBER]>1)
		{
			noalias(rRightHandSideVector) -= prod(DN_DX, vel_gauss);
			rRightHandSideVector += prod(Laplacian_matrix,old_pressures)*delta_t;
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
	
	void PFEM22D::CalculateLocalFinalVelocitySystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY
		//array_1d<double,3> gravity= rCurrentProcessInfo[GRAVITY];

		//WE MUST CALCULATE Mass and G(or D) matrixes

	    double delta_t = rCurrentProcessInfo[DELTA_TIME];	
		//KRATOS_WATCH(gradient_discontinuity);
		const bool & split_element = (this->GetValue(SPLIT_ELEMENT));

		double theta=1.0;
	
		//array_1d<double,3> pressures;

		
		array_1d<double,6>  previous_vel;
        for(unsigned int iii = 0; iii<3; iii++)
        {
			array_1d<double,3> velocity =  GetGeometry()[iii].FastGetSolutionStepValue(VELOCITY);
			previous_vel(iii*2) = velocity[0];
			previous_vel(iii*2+1) = velocity[1];
		}
		
		//boost::numeric::ublas::bounded_matrix<double,3,2> msDN_DX;
		//boost::numeric::ublas::bounded_matrix<double,2,2> msD;
		const double one_third = 1.0/3.0;
  		array_1d<double,3> msN; //dimension = number of nodes
		array_1d<double,6> ms_temp; //dimension = number of nodes

		const unsigned int LocalSize = GetGeometry().size()*2;
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
        boost::numeric::ublas::bounded_matrix<double, 3, 6 > D_matrix; //(gradient)
        noalias(D_matrix) = ZeroMatrix(3,6);	

        boost::numeric::ublas::bounded_matrix<double, 3*2, 3*2 > Mass_matrix;
        noalias(Mass_matrix) = ZeroMatrix(LocalSize,LocalSize);	
        //boost::numeric::ublas::bounded_matrix<double, TNumNodes*2, 2> N_vel_matrix; //useless to calculate it, it's simply the DN_DX matrix multiplied by 1/3, arranged in a different shape.
        
        GeometryUtils::CalculateGeometryData(this->GetGeometry(), DN_DX, N, Area);
         
        array_1d<double,3> distances(3);
        for (unsigned int i = 0; i < TNumNodes; i++)
        {
            distances[i] = this->GetGeometry()[i].FastGetSolutionStepValue(DISTANCE);
        }
        
        
		const double viscosity_air = rCurrentProcessInfo[VISCOSITY_AIR]; //para clindro en rey 100:  0.0005 * 2.0*Area;
		const double viscosity_water = rCurrentProcessInfo[VISCOSITY_WATER];
		
		const double density_air = rCurrentProcessInfo[DENSITY_AIR]; //para clindro en rey 100:  0.0005 * 2.0*Area;
		const double density_water = rCurrentProcessInfo[DENSITY_WATER];
		
		if (split_element==false)
		{
			//viscous term:
			if (distances(0)<0.0)
				this->AddViscousTerm(rLeftHandSideMatrix,DN_DX, (viscosity_water*Area) );
			else
				this->AddViscousTerm(rLeftHandSideMatrix,DN_DX, (viscosity_air*Area) );
				
				
			//mass matrix:
			double density;
			if (distances(0)<0.0)
				density=density_water;
			else
				density=density_air;
				
			double fourier=0.0;
			if (distances(0)<0.0) fourier = viscosity_water/density_water *delta_t/pow((0.6*this->GetValue(MEAN_SIZE)),2);
			else fourier = viscosity_air/density_air *delta_t/pow((0.6*this->GetValue(MEAN_SIZE)),2);
			theta = 0.2*fourier - 0.1;
			if (theta<0.0) theta=0.0;
			else if (theta>1.0) theta=1.0;
				
			const double Weight = one_third * Area * density;
			
			for (unsigned int i=0; i<LocalSize ; i++)
				Mass_matrix (i,i) = Weight;
				
			//rLeftHandSideMatrix += Mass_matrix;	
			
		} 
		else //split element:
		{
			//TOOLS NEEDED TO FIND THE ENRICHMENT SHAPE FUNCTIONS
			//get position of the cut surface
			
			//TOOLS NEEDED TO FIND THE ENRICHMENT SHAPE FUNCTIONS
			//get position of the cut surface
			array_1d<double,3>  densities(3);
			boost::numeric::ublas::bounded_matrix<double,3, 2> Nenriched;
			array_1d<double,3>  volumes(3);
			boost::numeric::ublas::bounded_matrix<double,3, 2 > coords;
			boost::numeric::ublas::bounded_matrix<double,3, 3 > Ngauss;
			array_1d<double,3>  signs(3);
			std::vector< Matrix > gauss_gradients(3);
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
			for (unsigned int i = 0; i < 3; i++)
				gauss_gradients[i].resize(2, 2, false);  //2 values of the 2 shape functions, and derivates in (xy) direction).
			unsigned int ndivisions = EnrichmentUtilities::CalculateEnrichedShapeFuncions(coords, DN_DX, distances, volumes, Ngauss, signs, gauss_gradients, Nenriched);
			
			
						
			double Weight=0.0;
			double mass=0.0;
			for (unsigned int i=0;i<ndivisions;i++) //we go over the three partitions;
			{
				if (signs(i)>0.0)
				{
					Weight += volumes(i)*viscosity_air;
					mass += volumes(i)*density_air;
				}
				else
				{
					Weight += volumes(i)*viscosity_water;
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

			
			for (unsigned int i = 0; i < 3; i++)  //partition
			{
				if (signs(i)<0.0)
					densities(i)=density_water;
				else
					densities(i)=density_air;
			}
			
			
			for (unsigned int i = 0; i < 3; i++)  // local node (or shape function)
			{
				//KRATOS_WATCH(enrich_lhs(i));
				for (unsigned int j = 0; j < 3; j++) //partitions
				{
					Mass_matrix(2*i,2*i) +=  volumes(j)*Ngauss(j,i)*densities(j);	 
				}
				Mass_matrix(2*i+1,2*i+1) = Mass_matrix(2*i,2*i);
			}
			
		}
		
		//we must actually use a fraction of the rigidity matrix. since the Lefthandside only has it so far, we multiply it by tita:
		rLeftHandSideMatrix *= theta;

        noalias(rRightHandSideVector) = prod((Mass_matrix/delta_t),previous_vel);

		noalias(rLeftHandSideMatrix) += Mass_matrix/delta_t;   
		
        

		//subtracting the dirichlet term
		// RHS -= LHS*DUMMY_UNKNOWNs
		noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix,previous_vel);
		
		KRATOS_CATCH("");
	}
	
	
	
	//*****************************************************************************
	//***************************************************************************
	

		void PFEM22D::CalculateViscousRHS(ProcessInfo& CurrentProcessInfo)
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
	

	void PFEM22D::CalculatePressureProjection(ProcessInfo& CurrentProcessInfo)
	{
		KRATOS_TRY
		
		
		if ((this->GetValue(IS_INACTIVE))==false) //elements can be inactive to add temporary walls. fractional velocity is integrated by parts so walls are seen as having zero velocity without doing anything 
				{
					const double mDENSITY_AIR = CurrentProcessInfo[DENSITY_AIR]; // * (1.0-theta) ;
					const double mDENSITY_WATER = CurrentProcessInfo[DENSITY_WATER]; // * (1.0-theta) ;
					const double mINV_DENSITY_AIR = 1.0/mDENSITY_AIR;
					const double mINV_DENSITY_WATER = 1.0/mDENSITY_WATER;
					
					//const double delta_t = CurrentProcessInfo[DELTA_TIME];
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

					
			
						
						
					
					
					if ((this->GetValue(SPLIT_ELEMENT))==true)
					{
						//to begin with we calculate the gradient discontinuity:
						if (true) 
						{
							const int iteration_number =  CurrentProcessInfo[NL_ITERATION_NUMBER];
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
						double jump_value=-input_jump_value*0.5;
						double gradient_discontinuity; // = this->GetValue(ENRICH_RHS);
						//array_1d<double,(2+1)> & enrich_lhs = this->GetValue(ENRICH_LHS_ROW);
						
						
						noalias(G_matrix) = ZeroMatrix((2+1), (2-1)*6);	
						noalias(G_matrix_mixed) = ZeroMatrix(1,(2-1)*6);	
						noalias(G_matrix_mixed_jump) = ZeroMatrix(1,(2-1)*6);	
						
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
						single_triangle_node=1;
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

						jump_value*=signs((single_triangle_node-1));

						gradient_discontinuity = this->GetValue(GRADIENT_DISCONTINUITY);

						///*****************************************************************************************
						for (unsigned int i = 0; i < ndivisions; i++)  // partiton i
						{ 
							for (unsigned int j = 0; j < (2+1); j++) //shape function (Nj)
							{
								for (unsigned int k=0;k!=2;k++) //x,y,(z)
								{
									G_matrix_mixed(0,j*(2)+k) += gauss_gradients[i](0,k) *volumes(i)*Ngauss(i,j)*inv_densities(i);
									G_matrix_mixed_no_ro(0,j*(2)+k) += gauss_gradients[i](0,k) *volumes(i)*Ngauss(i,j);
									

									G_matrix_mixed_jump(0,j*(2)+k) += gauss_gradients[i](1,k) *volumes(i)*Ngauss(i,j)*inv_densities(i);
									G_matrix_mixed_jump_no_ro(0,j*(2)+k) += gauss_gradients[i](1,k) *volumes(i)*Ngauss(i,j);

									
									for (unsigned int l = 0; l < (2+1); l++) //shape function derivative (dNl/dk)
									{
										G_matrix(l, (j*(2)+k) ) += DN_DX(l,k)*volumes(i)*Ngauss(i,j)*inv_densities(i);
									}
									
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
								
						gradient_discontinuity = this->GetValue(GRADIENT_DISCONTINUITY);
							
						//now we save the data:
						for (unsigned int i=0; i!=(2+1); i++) //loop around the nodes of the element to add contribution to node i
						{
							
							array_1d<double, 3 > & current_press_proj_no_ro = geom[i].FastGetSolutionStepValue(PRESS_PROJ_NO_RO);
							array_1d<double, 3 > & current_press_proj = geom[i].FastGetSolutionStepValue(PRESS_PROJ);
								
							for (unsigned int j = 0; j < (2+1) ; j++) //loop around the nodes of the element to add contribution of node j to node i
							{
								
								for (unsigned int k = 0; k < (2) ; k++) //x,y,(z)
								{
									current_press_proj[k] += G_matrix(j,i*(2)+k)*(pressure(j));///-old_pressures(j)); //gamma=0!
									current_press_proj_no_ro[k] += G_matrix_no_ro(j,i*(2)+k)*(pressure(j));///-old_pressures(j)); //gamma=0!
								}
							}
							
							for (unsigned int k = 0; k < (2) ; k++) //x,y,(z)
							{
								current_press_proj[k] += G_matrix_mixed(0,2*i+k)*gradient_discontinuity + G_matrix_mixed_jump(0,2*i+k)*jump_value;
								current_press_proj_no_ro[k] += G_matrix_mixed_no_ro(0,2*i+k)*gradient_discontinuity + G_matrix_mixed_jump_no_ro(0,2*i+k)*jump_value;
							}
							
							geom[i].FastGetSolutionStepValue(NODAL_AREA) += Area*mass_factor;
							geom[i].FastGetSolutionStepValue(NODAL_MASS) += nodal_masses(i);
							
						}
						
					}
					else //normal (uncut) element
					{
						noalias(G_matrix_no_ro) = ZeroMatrix((2+1),(2-1)*6);
						noalias(G_matrix) = ZeroMatrix((2+1),(2-1)*6);
						
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
								G_matrix_no_ro(i, (j*2)+k ) = DN_DX(i,k)*Area*mass_factor;
								}
							}
						}
						
						G_matrix /=density;

						//now we save the data:
						for (unsigned int i=0; i!=(2+1); i++) //loop around the nodes of the element to add contribution to node i
						{
							
							array_1d<double, 3 > & current_press_proj_no_ro = geom[i].FastGetSolutionStepValue(PRESS_PROJ_NO_RO);
							array_1d<double, 3 > & current_press_proj = geom[i].FastGetSolutionStepValue(PRESS_PROJ);
								
							for (unsigned int j = 0; j < (2+1) ; j++) //loop around the nodes of the element to add contribution of node j to node i
							{
								
								for (unsigned int k = 0; k < (2) ; k++) //x,y,(z)
								{
									current_press_proj[k] += G_matrix(j,i*(2)+k)*(pressure(j));///-old_pressures(j)); //gamma=0!
									current_press_proj_no_ro[k] += G_matrix_no_ro(j,i*(2)+k)*(pressure(j));///-old_pressures(j)); //gamma=0!
								}
							}
							
							geom[i].FastGetSolutionStepValue(NODAL_AREA) += Area*mass_factor;
							geom[i].FastGetSolutionStepValue(NODAL_MASS) += Area*mass_factor*density;
							
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
	void PFEM22D::InitializeSolutionStep(ProcessInfo& CurrentProcessInfo)
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
	void PFEM22D::PressureEquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo)
	{
		unsigned int number_of_nodes = GetGeometry().PointsNumber();
		if(rResult.size() != number_of_nodes)
			rResult.resize(number_of_nodes,false);	

		for (unsigned int i=0;i<number_of_nodes;i++)
				rResult[i] = GetGeometry()[i].GetDof(PRESSURE).EquationId();
	}
	
	//************************************************************************************
	void PFEM22D::VelocityEquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo)
	{
		const SizeType NumNodes = 3;
		const SizeType LocalSize = 6;
		GeometryType& rGeom = this->GetGeometry();

		SizeType LocalIndex = 0;

		if (rResult.size() != LocalSize)
			rResult.resize(LocalSize, false);

		for (SizeType i = 0; i < NumNodes; ++i)
		{
			rResult[LocalIndex++] = rGeom[i].GetDof(VELOCITY_X).EquationId();
			rResult[LocalIndex++] = rGeom[i].GetDof(VELOCITY_Y).EquationId();
		}
	}
		//************************************************************************************
	void PFEM22D::FractionalVelocityEquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo)
	{
		const SizeType NumNodes = 3;
		const SizeType LocalSize = 6;
		GeometryType& rGeom = this->GetGeometry();

		SizeType LocalIndex = 0;

		if (rResult.size() != LocalSize)
			rResult.resize(LocalSize, false);

		for (SizeType i = 0; i < NumNodes; ++i)
		{
			rResult[LocalIndex++] = rGeom[i].GetDof(FRACT_VEL_X).EquationId();
			rResult[LocalIndex++] = rGeom[i].GetDof(FRACT_VEL_Y).EquationId();
		}
	}

	//************************************************************************************
	//************************************************************************************
	void PFEM22D::GetPressureDofList(DofsVectorType& ElementalDofList,ProcessInfo& CurrentProcessInfo)
	{
		unsigned int number_of_nodes = GetGeometry().PointsNumber();
		if(ElementalDofList.size() != number_of_nodes)
			ElementalDofList.resize(number_of_nodes);	

		for (unsigned int i=0;i<number_of_nodes;i++)
			ElementalDofList[i] = GetGeometry()[i].pGetDof(PRESSURE);

	}

	//************************************************************************************
	void PFEM22D::GetVelocityDofList(DofsVectorType& ElementalDofList,ProcessInfo& CurrentProcessInfo)
	{
		const SizeType NumNodes = 3;
		const SizeType LocalSize = 6;
		GeometryType& rGeom = this->GetGeometry();

		if (ElementalDofList.size() != LocalSize)
			ElementalDofList.resize(LocalSize);

		SizeType LocalIndex = 0;

		for (SizeType i = 0; i < NumNodes; ++i)
		{
			ElementalDofList[LocalIndex++] = rGeom[i].pGetDof(VELOCITY_X);
			ElementalDofList[LocalIndex++] = rGeom[i].pGetDof(VELOCITY_Y);
		}
	}
	
		//************************************************************************************
	void PFEM22D::GetFractionalVelocityDofList(DofsVectorType& ElementalDofList,ProcessInfo& CurrentProcessInfo)
	{
		const SizeType NumNodes = 3;
		const SizeType LocalSize = 6;
		GeometryType& rGeom = this->GetGeometry();

		if (ElementalDofList.size() != LocalSize)
			ElementalDofList.resize(LocalSize);

		SizeType LocalIndex = 0;

		for (SizeType i = 0; i < NumNodes; ++i)
		{
			ElementalDofList[LocalIndex++] = rGeom[i].pGetDof(FRACT_VEL_X);
			ElementalDofList[LocalIndex++] = rGeom[i].pGetDof(FRACT_VEL_Y);
		}
	}
	
	//************************************************************************************
	
	//non partitioned elements using MatrixType
	void PFEM22D::AddViscousTerm(MatrixType& rDampingMatrix,
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
				rDampingMatrix(FirstRow,FirstCol) += Weight * ( FourThirds * rShapeDeriv(i,0) * rShapeDeriv(j,0) + rShapeDeriv(i,1) * rShapeDeriv(j,1) );
				rDampingMatrix(FirstRow,FirstCol+1) += Weight * ( nTwoThirds * rShapeDeriv(i,0) * rShapeDeriv(j,1) + rShapeDeriv(i,1) * rShapeDeriv(j,0) );

				// Second Row
				rDampingMatrix(FirstRow+1,FirstCol) += Weight * ( nTwoThirds * rShapeDeriv(i,1) * rShapeDeriv(j,0) + rShapeDeriv(i,0) * rShapeDeriv(j,1) );
				rDampingMatrix(FirstRow+1,FirstCol+1) += Weight * ( FourThirds * rShapeDeriv(i,1) * rShapeDeriv(j,1) + rShapeDeriv(i,0) * rShapeDeriv(j,0) );

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
		rDampingMatrix = prod(B_matrix, temp_matrix );
		
	}
	
	//using MatriType (for the implicit step)
	void PFEM22D::AddViscousTerm(MatrixType& rDampingMatrix,
                                       const boost::numeric::ublas::bounded_matrix<double, 3, 2>& rShapeDeriv,
                                       const double viscosity_air,
                                       const double viscosity_water,
                                       array_1d<double,3>  volumes,
                                       array_1d<double,3>  signs,
                                       Matrix Ngauss,
                                       const double Area)
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
				rDampingMatrix(FirstRow,FirstCol) += Weight * ( FourThirds * rShapeDeriv(i,0) * rShapeDeriv(j,0) + rShapeDeriv(i,1) * rShapeDeriv(j,1) );
				rDampingMatrix(FirstRow,FirstCol+1) += Weight * ( nTwoThirds * rShapeDeriv(i,0) * rShapeDeriv(j,1) + rShapeDeriv(i,1) * rShapeDeriv(j,0) );

				// Second Row
				rDampingMatrix(FirstRow+1,FirstCol) += Weight * ( nTwoThirds * rShapeDeriv(i,1) * rShapeDeriv(j,0) + rShapeDeriv(i,0) * rShapeDeriv(j,1) );
				rDampingMatrix(FirstRow+1,FirstCol+1) += Weight * ( FourThirds * rShapeDeriv(i,1) * rShapeDeriv(j,1) + rShapeDeriv(i,0) * rShapeDeriv(j,0) );

				// Update Counter
				FirstRow += 2;
			}
			FirstRow = 0;
			FirstCol += 2;
		}
		*/
		
		//boost::numeric::ublas::bounded_matrix<double, 3, 2> weighted_ShapeDeriv = ZeroMatrix(3,2);
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
		rDampingMatrix = prod(B_matrix, temp_matrix );
	}
	
	//using bounded_matrix, for the explict step:
	void PFEM22D::AddViscousTerm(boost::numeric::ublas::bounded_matrix<double, (2-1)*6, (2-1)*6 >& rDampingMatrix,
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
		rDampingMatrix = prod(B_matrix, temp_matrix );
		
	}
	
	//************************************************************************************

} // Namespace Kratos
