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

	    double delta_t = rCurrentProcessInfo[DELTA_TIME];	

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

			//saving everything
			previous_vel_in_mesh(iii*2) = velocity[0];
			previous_vel_in_mesh(iii*2+1) = velocity[1];

			previous_vel_and_press(iii*3) = velocity[0];
			previous_vel_and_press(iii*3+1) = velocity[1];
			previous_vel_and_press(iii*3+2) = GetGeometry()[iii].FastGetSolutionStepValue(PRESSURE);
		

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
		double density=this->GetValue(DENSITY);
		const double lambda = this->GetValue(YOUNG_MODULUS) * this->GetValue(POISSON_RATIO) / ( (1.0+this->GetValue(POISSON_RATIO))*(1.0-2.0*this->GetValue(POISSON_RATIO)) );	 
		const double mu = this->GetValue(YOUNG_MODULUS)/ (2.0*(1.0+this->GetValue(POISSON_RATIO))); 
		const double mu_integral = mu*solid_area;
		const double Bulk_modulus = (2.0/3.0 * mu + lambda);	
		
		this->AddElasticityTerm(rLeftHandSideMatrix,DN_DX, (mu_integral*delta_t));

		
		for (unsigned int i=0; i<TNumNodes ; i++)
		{

			Mass_matrix (i*3,i*3) = density*Area*one_third; //nodal_mass_fluid(i)+nodal_mass_solid(i);
			Mass_matrix (i*3+1,i*3+1) = density*Area*one_third; //nodal_mass_fluid(i)+nodal_mass_solid(i);

			Lumped_Pressure_Mass_matrix(i*3+2,i*3+2) = one_third * solid_area /Bulk_modulus;

		}
			
		//this->GetValue(DENSITY)=density;
		
		for (unsigned int i = 0; i < 3; i++)
		{
			for (unsigned int j = 0; j < 3; j++)
			{
				
				rLeftHandSideMatrix(i*3+2, j*3+0 ) -= (D_matrix(i,j*2) )*Area;     
				rLeftHandSideMatrix(i*3+2, j*3+1 ) -= (D_matrix(i,j*2+1) )*Area;    
				
				rLeftHandSideMatrix(j*3+0, i*3+2 ) -=  D_matrix(i,j*2)*Area;     
				rLeftHandSideMatrix(j*3+1, i*3+2 ) -= D_matrix(i,j*2+1)*Area; 
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


		noalias(rRightHandSideVector) += prod((Mass_matrix),((1.0)*previous_vel_and_press/delta_t));
		noalias(rLeftHandSideMatrix) += (1.0)*Mass_matrix/delta_t;   
				
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
					GeometryUtils::CalculateGeometryData(geom, DN_DX, N, Area);
					
					array_1d<double,(2+1)> nodal_masses = ZeroVector(2+1);

					array_1d<double,(2+1)> pressure;
					for (unsigned int i=0; i<(2+1);i++)
						pressure(i) = geom[i].FastGetSolutionStepValue(PRESSURE);
						
					
					boost::numeric::ublas::bounded_matrix<double, (2+1), (2-1)*6 > G_matrix; //(gradient)
					noalias(G_matrix) = ZeroMatrix((2+1), (2-1)*6);	
					
					boost::numeric::ublas::bounded_matrix<double, (2+1), (2-1)*6 > G_matrix_no_ro; //(gradient)
					noalias(G_matrix_no_ro) = ZeroMatrix((2+1),(2-1)*6);		
					
					
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
				
				
				Vector & element_stresses= this->GetValue(ELEMENT_MEAN_STRESS);
				
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
