//   
//   Project Name:        Kratos       
//   Last modified by:    $Author: it's me! $
//   Date:                $Date: 2008-08-08 23:58:38 $
//   Revision:            $Revision: 1.0 $
//
//

// Project includes 
#include "includes/define.h"
#include "custom_elements/fsi_3d.h"
#include "custom_utilities/pfem_particle.h"
#include "pfem_2_application.h"
#include "utilities/math_utils.h"
#include "utilities/geometry_utilities.h" 
#include "utilities/enrichment_utilities.h"
#include "utilities/discont_utils.h"

//#include "utilities/enrich_2d_2dofs.h"

#include "processes/process.h"
#include "includes/node.h"
//#include "includes/element.h"
#include "includes/model_part.h"



namespace Kratos
{
	


	//************************************************************************************
	//************************************************************************************
	FsiPFEM23D::FsiPFEM23D(IndexType NewId, GeometryType::Pointer pGeometry)
		: Element(NewId, pGeometry)
	{		
		//DO NOT ADD DOFS HERE!!!
	}

	//************************************************************************************
	//************************************************************************************
	FsiPFEM23D::FsiPFEM23D(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
		: Element(NewId, pGeometry, pProperties)
	{

	}

	Element::Pointer FsiPFEM23D::Create(IndexType NewId, NodesArrayType const& ThisNodes,  PropertiesType::Pointer pProperties) const
	{
		return Element::Pointer(new FsiPFEM23D(NewId, GetGeometry().Create(ThisNodes), pProperties));
	}

	FsiPFEM23D::~FsiPFEM23D()
	{
	}
	
	//************************************************************************************
	//************************************************************************************
	
	
	void FsiPFEM23D::AddExplicitContribution(ProcessInfo& rCurrentProcessInfo)
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
	
	void FsiPFEM23D::CalculateLocalSystem(MatrixType& rLeftHandSideMatrix, VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
	{
		KRATOS_TRY
		array_1d<double,3> gravity= rCurrentProcessInfo[GRAVITY];

		const double pressure_increasing_factor=1.0e0;

		//WE MUST CALCULATE Mass and G(or D) matrixes

	    double delta_t = rCurrentProcessInfo[DELTA_TIME];	
	    const int offset = rCurrentProcessInfo[PARTICLE_POINTERS_OFFSET];	
		//const bool split_element = (this->GetValue(SPLIT_ELEMENT));

		double theta=1.0;
		
		//element mean stress, not used for calculations but rather only to reseed particles if needed.
		Vector & stresses = this->GetValue(ELEMENT_MEAN_STRESS);
		if (stresses.size()!=6)
		{
			stresses.resize(6);
			stresses=ZeroVector(6);
		}

		const unsigned int LocalSize = GetGeometry().size()*4; //3 vel + PRESSURE!
		unsigned int TNumNodes = GetGeometry().size();
		
        // Check sizes and initialize
        if (rRightHandSideVector.size() != LocalSize)
            rRightHandSideVector.resize(LocalSize, false);
       
        if(rLeftHandSideMatrix.size1() != LocalSize)
			rLeftHandSideMatrix.resize(LocalSize,LocalSize,false);
	    // Calculate this element's geometric parameters
		double Area;
		array_1d<double, 4> N; //dimension = number of nodes
		boost::numeric::ublas::bounded_matrix<double, 4, 3> DN_DX;
		GeometryUtils::CalculateGeometryData(this->GetGeometry(), DN_DX, N, Area);

		noalias(rRightHandSideVector) = ZeroVector(LocalSize);	
		noalias(rLeftHandSideMatrix) = ZeroMatrix(LocalSize,LocalSize);	
		
		int number_of_particles_in_elem = this->GetValue(NUMBER_OF_PARTICLES);
		
		if( (number_of_particles_in_elem>0) )  //ok! active element!
		//if ((this->GetGeometry()[0].FastGetSolutionStepValue(DISTANCE))<0.0 || (this->GetGeometry()[1].FastGetSolutionStepValue(DISTANCE))<0.0 || (this->GetGeometry()[2].FastGetSolutionStepValue(DISTANCE))<0.0)
		{
			//we start by restting the LHS and RHS

			
			array_1d<double,12>  previous_vel_in_mesh = ZeroVector(12); //to calculate the deformation Gradient F. Dimension = velocity dofs
			array_1d<double,12>  previous_accel_in_mesh= ZeroVector(12); //to improve the computation of the displacements calculation, (used in the computation of the deformation gradient F too) Dimension = velocity dofs
			array_1d<double,4>  total_N_from_particles = ZeroVector(4); //to weight the contribution of the particles to each of the nodes. Dimension = number of nodes.
			array_1d<double,4>  total_N_from_solid_particles = ZeroVector(4); //to weight the contribution of the particles to each of the nodes. Dimension = number of nodes.
			array_1d<double,4>  solid_pressure_from_particles = ZeroVector(4); //to weight the contribution of the particles to each of the nodes. Dimension = number of nodes.

			array_1d<double,16>  previous_vel_and_solid_press = ZeroVector(16); //to add to the righthandside the unkwnonws from the previous time step. Dimension = total dofs
			array_1d<double,16>  previous_vel_and_press = ZeroVector(16); //to add to the righthandside the unkwnonws from the previous time step. Dimension = total dofs
			array_1d<double,16>  previous_accel= ZeroVector(16);          //same as above, but pressure spaces are left blank (unused). Dimension = total dofs
			boost::numeric::ublas::bounded_matrix<double,4, 3 > coords; //coordinates of the nodes
			array_1d<double,4> distances;
			bool has_negative_node=false;
			bool has_positive_node=false;
			
			for(unsigned int iii = 0; iii<4; iii++)
			{
				array_1d<double,3>& velocity =  GetGeometry()[iii].FastGetSolutionStepValue(VELOCITY);  //reference
				array_1d<double,3>& accel =  GetGeometry()[iii].FastGetSolutionStepValue(ACCELERATION); //reference
				//array_1d<double,3> accel = (GetGeometry()[iii].FastGetSolutionStepValue(VELOCITY,1)- GetGeometry()[iii].GetSolutionStepValue(VELOCITY,2) )/delta_t;
				//velocity -= gravity*delta_t; //it was added in the particles, but since we need it here on the rhs instead of inside the velocity, we subtract it in order to add it correctly in the rhs.
				
				//saving everything
				previous_vel_in_mesh(iii*3) = velocity[0];
				previous_vel_in_mesh(iii*3+1) = velocity[1];
				previous_vel_in_mesh(iii*3+2) = velocity[2];
				previous_accel_in_mesh(iii*3) = accel[0];
				previous_accel_in_mesh(iii*3+1) = accel[1];
				previous_accel_in_mesh(iii*3+2) = accel[2];
				//
				previous_vel_and_press(iii*4) = velocity[0];
				previous_vel_and_press(iii*4+1) = velocity[1];
				previous_vel_and_press(iii*4+2) = velocity[2];
				previous_vel_and_press(iii*4+3) = GetGeometry()[iii].FastGetSolutionStepValue(PRESSURE);
				//
				previous_vel_and_solid_press(iii*4) = velocity[0];
				previous_vel_and_solid_press(iii*4+1) = velocity[1];
				previous_vel_and_solid_press(iii*4+2) = velocity[2];
				previous_vel_and_solid_press(iii*4+3) = GetGeometry()[iii].FastGetSolutionStepValue(SOLID_PRESSURE);
				//previous_vel_and_press(iii*3+3) = GetGeometry()[iii].FastGetSolutionStepValue(SOLID_PRESSURE);
				//
				previous_accel(iii*4+0) = accel[0];
				previous_accel(iii*4+1) = accel[1]; 
				previous_accel(iii*4+2) = accel[2]; 
				//				
				const array_1d<double, 3 > & xyz = this->GetGeometry()[iii].Coordinates();
				for (unsigned int j = 0; j < 3; j++)
					coords(iii, j) = xyz[j];
				//
				distances[iii] = this->GetGeometry()[iii].FastGetSolutionStepValue(DISTANCE);
				if (distances[iii]<0.0)
					has_negative_node=true;
				else
					has_positive_node=true;		
			}
			
			bool split_element=false;
			if (has_positive_node && has_negative_node)
				split_element=true;
			
			const double one_quarter = 1.0/4.0; //in 3d we would need one_quarter
        
			boost::numeric::ublas::bounded_matrix<double, 4, 4> Laplacian_matrix; //pressure dof^2
			noalias(Laplacian_matrix) = ZeroMatrix(4,4);	
			boost::numeric::ublas::bounded_matrix<double, 4, 12 > D_matrix; //(divergence)
			noalias(D_matrix) = ZeroMatrix(4,12);	
			boost::numeric::ublas::bounded_matrix<double, 4, 12 > G_matrix; //(gradient)
			noalias(G_matrix) = ZeroMatrix(4,12);	
			boost::numeric::ublas::bounded_matrix<double, 16, 16 > Mass_matrix; //2 vel + 1 pressure per node
			noalias(Mass_matrix) = ZeroMatrix(16,16);	
			boost::numeric::ublas::bounded_matrix<double, 16, 16 > Pressure_Mass_matrix; //2 vel + 1 pressure per node
			noalias(Pressure_Mass_matrix) = ZeroMatrix(16,16);	
			boost::numeric::ublas::bounded_matrix<double, 16, 16 > Lumped_Pressure_Mass_matrix; //2 vel + 1 pressure per node
			noalias(Lumped_Pressure_Mass_matrix) = ZeroMatrix(16,16);	
			
			//we start by calculating the non-enriched functions of divergence and gradient. as i write this, gradient integrated by parts.
			for (unsigned int i = 0; i < 4; i++)
			{
				for (unsigned int j = 0; j < 4; j++)
				{
					D_matrix(i, (j*3  ) ) =  DN_DX(j,0)*one_quarter;     
					D_matrix(i, (j*3+1) ) = DN_DX(j,1)*one_quarter;
					D_matrix(i, (j*3+2) ) = DN_DX(j,2)*one_quarter;
					G_matrix(i, (j*3  ) ) = - DN_DX(i,0)*one_quarter;     
					G_matrix(i, (j*3+1) ) = - DN_DX(i,1)*one_quarter;
					G_matrix(i, (j*3+2) ) = - DN_DX(i,2)*one_quarter;
				}
			}
			
			
			//area/volume integrals using the particles  in the elements:
			
			array_1d<double,6> total_integral_stresses = ZeroVector(6);
			double total_integral_mean_solid_pressure = 0.0;
			double total_integral_mean_fluid_pressure = 0.0;
			bool has_solid_particle=false;
			bool has_fluid_particle=false;
			double fluid_area=Area*0.00000000000000000000000001; //to avoid problems when we have 2 elements 
			double solid_area=Area*0.00000000000000000000000001;
			double density_solid_integral=0.01*solid_area;
			double density_fluid_integral=0.01*fluid_area;
			double viscosity_integral = 1.0000*fluid_area;
			double mu_integral = 1000000.0*solid_area;
			double Bulk_modulus_integral_solid =  1000000.0*solid_area;
			double Bulk_modulus_integral_fluid =  1000000.0*fluid_area;
			array_1d<double,4> nodal_density = ZeroVector(4);

			int non_updated_solid_particles=0;
			
			ParticlePointerVector&  element_particle_pointers = (this->GetValue(PARTICLE_POINTERS)); //self-explainatory name
			//std::cout << "elem " << ii << " with " << (unsigned int)number_of_particles_in_elem << " particles" << std::endl;
			
			int number_of_useful_particles_in_elem=0;
			
			for (unsigned int iii=0; iii<number_of_particles_in_elem ; iii++ )
			{
				PFEM_Particle & pparticle = element_particle_pointers[offset+iii];	
				if (pparticle.GetEraseFlag()==false)
				{
					number_of_useful_particles_in_elem++;
				}
			}
			if 	(number_of_useful_particles_in_elem != number_of_particles_in_elem)
				KRATOS_WATCH("artt");
			
			double particle_area = Area/double(number_of_particles_in_elem); //giving to each particle equal fractions of area/volume.

			
			for (unsigned int iii=0; iii<number_of_particles_in_elem ; iii++ )
			{
				PFEM_Particle & pparticle = element_particle_pointers[offset+iii];	
				if (pparticle.GetEraseFlag()==false)
				{
					array_1d<double,3> particle_position = pparticle.Coordinates(); //it's a copy!! 
					array_1d<double,4> particle_N;
					bool found_it= CalculatePosition(coords, particle_position(0), particle_position(1), particle_position(2), particle_N);
					
					if ( pparticle.GetDistance()<0.0 ) //solid particle;
					{
						if( pparticle.HasUpdatedStresses()==true)
						{
							has_solid_particle=true;
							solid_area +=particle_area;	
							array_1d<double,6> particle_stresses = pparticle.GetSigma(); //it's a copy!! i guess it is faster than creating a reference, then we copy it back to the particle again;
							double particle_pressure = pparticle.GetPressure();        //copy pressure
							const double particle_mu = pparticle.GetShearModulus();    //copy mu
							//this->UpdateStressesToNewConfiguration(particle_stresses,particle_pressure,DN_DX,previous_vel_in_mesh+previous_accel_in_mesh*delta_t, particle_mu , delta_t); //updating stresses
							total_integral_stresses += particle_area*particle_stresses;               //ading particle contribution
							total_integral_mean_solid_pressure += particle_area*particle_pressure;    //adding particle contribution
							mu_integral += pparticle.GetShearModulus()*particle_area;				  //adding particle contribution
							Bulk_modulus_integral_solid += pparticle.GetBulkModulus()*particle_area;  //adding particle contribution
							density_solid_integral += pparticle.GetDensity()*particle_area;				  //adding particle contribution
													
							pparticle.GetSigma() = particle_stresses;		// updating particle information							
							pparticle.GetPressure() = particle_pressure;	// updating particle information					
							
							//////////(LINES XX)  DOES NOT WORK: NODAL VALUES PROVIDE BETTER RESULTS THAN USING THE PARTICLES TU CONSTRUCT THE PREVIOUS VELOCITY, ACCEL AND PRESSURE
							//array_1d<double,3> particle_velocity = pparticle.GetVelocity(); //it's a copy!!
							//array_1d<double,3> particle_acceleration = pparticle.GetAcceleration(); //it's a copy!! 
							
							for (unsigned int i=0; i<TNumNodes ; i++)
							{
								//for (unsigned int j=0; j<2 ; j++) //x,y
								//{
									//previous_vel_and_press(i*3+j) += particle_velocity(j)*particle_N(i);
									//previous_accel(i*3+j) += particle_acceleration(j)*particle_N(i);
								//}
								//previous_vel_and_press(i*3+2) += particle_pressure*particle_N(i);
								solid_pressure_from_particles(i) += particle_pressure*particle_N(i);
								total_N_from_solid_particles(i) += particle_N(i);
							}
						}
						else
							non_updated_solid_particles++;
							

						
					}
					else //fluid particle 
					{
						has_fluid_particle=true;
						fluid_area +=particle_area;	
						double particle_pressure = pparticle.GetPressure();
						//total_integral_stresses += 0.0*particle_stresses; //no stresses!
						total_integral_mean_fluid_pressure += particle_area*particle_pressure; //we keep the pressure, but unused if incompressible fluid is used.
							
						viscosity_integral += pparticle.GetShearModulus()*particle_area;
						Bulk_modulus_integral_fluid += pparticle.GetBulkModulus()*particle_area;
						density_fluid_integral += pparticle.GetDensity()*particle_area;
					}
					
					for (unsigned int i=0; i<TNumNodes ; i++)
					{
						nodal_density(i)+=pparticle.GetDensity()*particle_N(i);
						total_N_from_particles(i) += particle_N(i);
					}
				}

			}
			
			for (unsigned int i=0; i<TNumNodes ; i++)
				nodal_density(i) /=  total_N_from_particles(i);
			
			stresses=total_integral_stresses/solid_area;
			if(non_updated_solid_particles>0)
			{
				for (unsigned int iii=0; iii<number_of_particles_in_elem ; iii++ )
				{
					PFEM_Particle & pparticle = element_particle_pointers[offset+iii];	
					if (pparticle.GetEraseFlag()==false)
					{
						if ( pparticle.GetDistance()<0.0 ) //solid particle;
						{
							if( pparticle.HasUpdatedStresses()==false)
							{
								pparticle.GetSigma() = stresses;		// updating particle information							
								pparticle.GetPressure() =total_integral_mean_solid_pressure/solid_area;	// updating particle information	
								pparticle.HasUpdatedStresses()=true;			
							}
						}
					}
				}
			}
			
			//now we artificially enlarge the ammount of particles that have information
			const double ratio = (solid_area+double(non_updated_solid_particles)*particle_area)/solid_area;
			solid_area*=ratio;
			
			total_integral_stresses *= ratio;               //ading particle contribution
			total_integral_mean_solid_pressure *= ratio;    //adding particle contribution
			mu_integral *= ratio;				  //adding particle contribution
			Bulk_modulus_integral_solid *= ratio;  //adding particle contribution
			density_solid_integral *= ratio;				  //adding particle contribution
			
			
			
			
			this->GetValue(ELEMENT_MEAN_STRESS)=stresses;
			
			//we keep the integral of the variables in the whole element, it doesn't make any sense dividinf it to multiply it later when calculating elastic or viscosity matrices
			//mu /= solid_area;      
			//viscosity /= fluid_area;
			//const double Bulk_modulus= Bulk_modulus_integral_solid/solid_area; //+Bulk_modulus_integral_fluid/fluid_area; //unique for solid and fluid. might have to change this!
			//const double Bulk_modulus= Bulk_modulus_integral_solid/solid_area+Bulk_modulus_integral_fluid/fluid_area; //unique for solid and fluid. might have to change this!
			const double Bulk_modulus_solid = Bulk_modulus_integral_solid/solid_area; //solid
			const double Bulk_modulus_fluid = Bulk_modulus_integral_fluid/fluid_area; //fluid
			const double mu = mu_integral/solid_area;
			const double solid_pressure = total_integral_mean_solid_pressure/solid_area;
			const double fluid_pressure = total_integral_mean_fluid_pressure/fluid_area;
			const double density_solid = density_solid_integral/solid_area;
			const double density_fluid = density_fluid_integral/fluid_area;
			/*
			//printing stuff 
			KRATOS_WATCH(viscosity_integral/Area)
			KRATOS_WATCH(mu_integral/Area)
			KRATOS_WATCH(total_integral_mean_fluid_pressure)
			KRATOS_WATCH(total_integral_mean_solid_pressure)
			KRATOS_WATCH(density_integral/Area)
			KRATOS_WATCH(Bulk_modulus)
			*/
			
			//////////USEFUL WITH (LINES XX) SOME LINES ABOVE
			for (unsigned int i=0; i<TNumNodes ; i++)
			{
				//for (unsigned int j=0; j<2 ; j++) //x,y
				//{
					//previous_vel_and_press(i*3+j) /= total_N_from_particles(i);
					//previous_accel(i*3+j) /= total_N_from_particles(i);
				//}
				//previous_vel_and_press(i*3+2) /= total_N_from_particles(i);
				solid_pressure_from_particles(i) /=total_N_from_solid_particles(i);
			}
			
			
			////to replace the particle stresses with the constant one of the elment 
			//array_1d<double,3> particle_stresses = stresses;
			//this->AddLinearDeltaStresses(particle_stresses,DN_DX,previous_step_vel, mu , Poisson,delta_t);
			//stresses=particle_stresses;
			//total_integral_stresses=stresses*Area;
			//for (unsigned int iii=0; iii<number_of_particles_in_elem ; iii++ )
			//{
				//PFEM_Particle & pparticle = element_particle_pointers[offset+iii];	
				//if (pparticle.GetEraseFlag()==false && pparticle.GetDistance()<0.0 ) 
					//pparticle.GetSigma() = stresses;
			//}

			//first approximation: continuous velocity field: both elastic and pressure term remain constant and based only on the ammount of solid and fluid the element has.	
			

			double add_accel_stab_term=0.0; //(true or false)
			
			//if (split_element && (fabs(distances(0))>0.02 && fabs(distances(1))>0.02 && fabs(distances(2))>0.02 && fabs(distances(3))>0.02))
			if (split_element && has_solid_particle && has_fluid_particle)
			//if (false==true)
			//if (split_element)
			{
				add_accel_stab_term=0.0; //we turn it off
				
				
				//TOOLS NEEDED TO FIND THE ENRICHMENT SHAPE FUNCTIONS
				//get position of the cut surface
				
				//TOOLS NEEDED TO FIND THE ENRICHMENT SHAPE FUNCTIONS
				//get position of the cut surface
				array_1d<double,6>  densities(0);
				array_1d<double,6>  inv_densities(0);
				array_1d<double,6>  bulks(0);
				
				boost::numeric::ublas::bounded_matrix<double,6, 2> Nenriched;
				array_1d<double,6>  volumes(0);
				boost::numeric::ublas::bounded_matrix<double,4, 3 > coords;
				boost::numeric::ublas::bounded_matrix<double,6, 4 > Ngauss;
				array_1d<double,6>  signs(0);
				std::vector< Matrix > gauss_gradients(6);
				//fill coordinates
			   
				//unsigned int single_triangle_node = 1;
				for (unsigned int i = 0; i < TNumNodes; i++)
				{
					const array_1d<double, 3 > & xyz = this->GetGeometry()[i].Coordinates();
					//KRATOS_WATCH(distances(i));
					for (unsigned int j = 0; j < 3; j++)
						coords(i, j) = xyz[j];
				}
				for (unsigned int i = 0; i < 6; i++)
					gauss_gradients[i].resize(2, 3, false);  //2 values for the 2 shape functions, and derivates in (xy) direction).
				unsigned int ndivisions = EnrichmentUtilities::CalculateEnrichedShapeFuncions(coords, DN_DX, distances, volumes, Ngauss, signs, gauss_gradients, Nenriched);
				
				
				double distance_based_solid_area=0.0;
				double distance_based_fluid_area=0.0;

				for (unsigned int i = 0; i < ndivisions; i++)  //partition
				{
					if (signs(i)<0.0)
					{
						densities(i) = density_solid;
						inv_densities(i) = 1.0/density_solid;
						bulks(i) = Bulk_modulus_solid;
						distance_based_solid_area+=volumes(i);
					}
					else
					{
						densities(i) = density_fluid;
						inv_densities(i) = 1.0/density_fluid;
						bulks(i) = Bulk_modulus_solid;             //warning!, setting the water bulk modulus just like the solid one to have a "nicer" matrix
						distance_based_fluid_area+=volumes(i);
					}
				}
				
				array_1d<double,4>  lumped_mass = ZeroVector(4);
				array_1d<double,4>  lumped_pressure_mass = ZeroVector(4);
				
				double Density=0.0;
				for (unsigned int i = 0; i < 4; i++)  // local node (or shape function)
				{
					//KRATOS_WATCH(enrich_lhs(i));
					for (unsigned int j = 0; j < ndivisions; j++) //partitions
					{
						lumped_mass(i) +=  volumes(j)*Ngauss(j,i)*densities(j);	 
						//lumped_pressure_mass(i) += volumes(j)*Ngauss(j,i)*1.0/bulks(j);	
					}
					

					Mass_matrix(4*i,4*i) += lumped_mass(i);
					Mass_matrix(4*i+1,4*i+1) += lumped_mass(i);
					Mass_matrix(4*i+2,4*i+2) += lumped_mass(i);
					
					//Lumped_Pressure_Mass_matrix(i*4+3,i*4+3) = one_quarter * solid_area /Bulk_modulus_solid;//lumped_pressure_mass(i); //one_third*distance_based_solid_area/Bulk_modulus_solid;//
					//Pressure_Mass_matrix(i*4+3,i*4+3) = lumped_pressure_mass(i);  //one_third*distance_based_solid_area/Bulk_modulus_solid;////we only have 1 gp per partition, so in splitted elems we also lump the consistent mass matrix
					Density+=densities(i)*volumes(i);

				}
				Density*=1.0/Area;
				
				const double viscosity = viscosity_integral/fluid_area;
				const double density = (density_fluid_integral+density_solid_integral)/Area;
				//this->AddElasticityTerm(rLeftHandSideMatrix,DN_DX, (mu_integral*delta_t));
				this->AddElasticityTerm(rLeftHandSideMatrix,DN_DX, (mu*distance_based_solid_area*delta_t));
				//and for the viscous term:
				//
				//rLeftHandSideMatrix*=0.5;
				
				//noalias(rRightHandSideVector) += prod(rLeftHandSideMatrix,previous_vel_and_press);
				
				this->AddViscousTerm(rLeftHandSideMatrix,DN_DX, (viscosity*distance_based_fluid_area));				
				
				//PRESSURE TERMS:
				/////////////////// calculating tau
				double AdvVelNorm = sqrt(pow(this->GetGeometry()[0].FastGetSolutionStepValue(VELOCITY_X),0)+pow(this->GetGeometry()[0].FastGetSolutionStepValue(VELOCITY_Y),0)+pow(this->GetGeometry()[0].FastGetSolutionStepValue(VELOCITY_Z),0));
				//const double TimeFactor = rCurrentProcessInfo.GetValue(DYNAMIC_TAU);
				double mElemSize;
				array_1d<double,3> Edge(3,0.0);
				Edge = this->GetGeometry()[1].Coordinates() - this->GetGeometry()[0].Coordinates();
				mElemSize = Edge[0]*Edge[0];
				for (SizeType d = 1; d < 3; d++)
					mElemSize += Edge[d]*Edge[d];
				for (SizeType i = 2; i < 4; i++)
					for(SizeType j = 0; j < i; j++)
					{
						Edge = this->GetGeometry()[i].Coordinates() - this->GetGeometry()[j].Coordinates();
						double Length = Edge[0]*Edge[0];
						for (SizeType d = 1; d < 3; d++)
							Length += Edge[d]*Edge[d];
						if (Length < mElemSize) mElemSize = Length;
					}
					
				mElemSize = sqrt(mElemSize);
				
				//const double   TauOne = (2.0 / ( ( 1.0 / delta_t + 4.0 * viscosity / (density * mElemSize * mElemSize) + 2.0 * AdvVelNorm / mElemSize) ) );
				//changed definition to account for the "viscosity" of the elastic part, which is G*delta_t
				//this will reduce the TauOne, adding less artificial diffussion and improving mass conservation.
				//const double apparent_viscosity = (viscosity * distance_based_fluid_area +  mu*distance_based_solid_area*delta_t)/Area;
				const double apparent_viscosity = (viscosity * fluid_area +  mu*solid_area*delta_t)/Area;

				const double TauOne = 0.0*(1.0 / ( ( 1.0 / delta_t + 4.0 * apparent_viscosity / (density * mElemSize * mElemSize) + 2.0 * AdvVelNorm / mElemSize) ) );

				const double   TauOne_solid = 0.0 / ( 10.0 * mu*delta_t / (density * mElemSize * mElemSize) ) ;
				////////////////////
				//TauOne=1.0;
				array_1d<double,6> TAUs(0);
				for (unsigned int i = 0; i < ndivisions; i++)  //partition
				{
					if (signs(i)<0.0)
						TAUs(i)=TauOne*0.0;
					else
						TAUs(i)=TauOne;
				}

				

				//and now the rest of the things needed
				
				//boost::numeric::ublas::bounded_matrix<double, 2, 2> DN_DX_enrich;
				//boost::numeric::ublas::bounded_matrix<double, 2, 2> DN_DX_enrich_negative;
				const int offset = 0;
				const int Pdof = 4;
				
				boost::numeric::ublas::bounded_matrix<double, Pdof, Pdof >   condensed_pressure_mass;
				noalias(condensed_pressure_mass) = ZeroMatrix(Pdof,Pdof);	
				boost::numeric::ublas::bounded_matrix<double, Pdof, Pdof > Laplacian_enrich;
				noalias(Laplacian_enrich) = ZeroMatrix(Pdof,Pdof);	
				boost::numeric::ublas::bounded_matrix<double, Pdof, 4 > mixed_Laplacian;
				noalias(mixed_Laplacian) = ZeroMatrix(Pdof,4);	
				//To calculate D* =int (N*DN_DX_enrich).  we must use the N at the gauss point of the partition (returned by the enrichment util) and multiply by the area. 
				//notice that this first row is only the upper part of the mixed D. we also need the lower one for the jump 
				boost::numeric::ublas::bounded_matrix<double, Pdof, 12 > G_matrix_mixed;
				noalias(G_matrix_mixed) = ZeroMatrix(Pdof,12);	
				boost::numeric::ublas::bounded_matrix<double, Pdof, 12 > D_matrix_mixed;
				noalias(D_matrix_mixed) = ZeroMatrix(Pdof,12);	
				
				array_1d<double,12> mass_stabilization_terms=ZeroVector(12);
				
				//double rhs_enrich = 0.0; //the rhs term of the enriched dof 
				array_1d<double,Pdof> rhs_enrich;
				noalias(rhs_enrich) = ZeroVector(Pdof);	
				
				noalias(D_matrix) = ZeroMatrix(4,12);	
				noalias(G_matrix) = ZeroMatrix(4,12);
				//D_matrix*=Area;	
				noalias(Laplacian_matrix) = ZeroMatrix(4,4);
				
				for (unsigned int i = 0; i < ndivisions; i++)  //we go through the 3 partitions of the domain
				{
					
					for (unsigned int j = 0; j < Pdof; j++) //we go through the 3 pressure shape functions, they can be either enriched or not depending on the sign(partition)*sign(shape_function)
					{
						if(signs(i)*distances(j)>0.0)
						{
							Lumped_Pressure_Mass_matrix(j*4+3,j*4+3)+=volumes(i)*Ngauss(i,j)*1.0/bulks(i);	
							
							for (unsigned int k = 0; k < Pdof; k++) //we go through the 3 velocity shape functions, they can be either enriched or not depending on the sign(partition)*sign(shape_function)
							{
									D_matrix(j,k*3) += DN_DX(k,0) *volumes(i)*Ngauss(i,j);
									D_matrix(j,k*3+1) += DN_DX(k,1) *volumes(i)*Ngauss(i,j);	
									D_matrix(j,k*3+2) += DN_DX(k,2) *volumes(i)*Ngauss(i,j);	
									
									G_matrix(j,k*3) -= DN_DX(j,0) *volumes(i)*Ngauss(i,k);
									G_matrix(j,k*3+1) -= DN_DX(j,1) *volumes(i)*Ngauss(i,k);	
									G_matrix(j,k*3+2) -= DN_DX(j,2) *volumes(i)*Ngauss(i,k);	
							}
						}
						//if(signs(i)*distances(j)<0.0)
						else
						{
							condensed_pressure_mass(j,j)+=volumes(i)*Ngauss(i,j)*1.0/bulks(i);
							
							for (unsigned int k = 0; k < Pdof; k++) //we go through the 3 shape functions, they can be either enriched or not depending on the sign(partition)*sign(shape_function)
							{
									D_matrix_mixed(j,k*3) += DN_DX(k,0) *volumes(i)*Ngauss(i,j);
									D_matrix_mixed(j,k*3+1) += DN_DX(k,1) *volumes(i)*Ngauss(i,j);		
									D_matrix_mixed(j,k*3+2) += DN_DX(k,2) *volumes(i)*Ngauss(i,j);		
									
									G_matrix_mixed(j,k*3) -= DN_DX(j,0) *volumes(i)*Ngauss(i,k);
									G_matrix_mixed(j,k*3+1) -= DN_DX(j,1) *volumes(i)*Ngauss(i,k);		
									G_matrix_mixed(j,k*3+2) -= DN_DX(j,2) *volumes(i)*Ngauss(i,k);		
							}
						}
					}
				}

				double number_of_positive_nodes=0.0;
				double number_of_negative_nodes=0.0;
				double solid_pressure_from_nodes=0.0;
				double fluid_pressure_from_nodes=0.0;
				for (unsigned int k = 0; k < 4; k++)
				{
					if (distances(k)>0.0)
					{
						number_of_positive_nodes+=1.0;
						fluid_pressure_from_nodes+=GetGeometry()[k].FastGetSolutionStepValue(PRESSURE);
					}
					else
					{
						number_of_negative_nodes+=1.0;
						solid_pressure_from_nodes+=GetGeometry()[k].FastGetSolutionStepValue(PRESSURE);
					}
				}
				solid_pressure_from_nodes/=number_of_negative_nodes;
				fluid_pressure_from_nodes/=number_of_positive_nodes;
				
				
				for (unsigned int k = 0; k < Pdof; k++) //we go through the 4 enrichments
				{
					//Lumped_Pressure_Mass_matrix(k*3+2,k*3+2)+=one_third*Area*1.0/bulks(k);	
					//condensed_pressure_mass(k,k)+=one_third*Area*1.0/bulks(k);
					
					
					if (distances(k)<0.0)
						rhs_enrich(k) -= condensed_pressure_mass(k,k)/delta_t * solid_pressure_from_nodes; //solid_pressure;
					else
						rhs_enrich(k) -= condensed_pressure_mass(k,k)/delta_t * fluid_pressure_from_nodes; //fluid_pressure:
					
				}	
				
				
				boost::numeric::ublas::bounded_matrix<double, Pdof, 16 > condensed_rows=ZeroMatrix(Pdof,16); //Vx1,Vy1,p1,Vx2,...
				boost::numeric::ublas::bounded_matrix<double, Pdof, 16 > condensed_columns=ZeroMatrix(Pdof,16); //Vx1,Vy1,p1,Vx2,...
				for (unsigned int i = 0; i < 4; i++)
				{
					//condensed_row(0,i*3+0)=-D_matrix_mixed(0,i*2+0);
					//condensed_row(0,i*3+1)=-D_matrix_mixed(0,i*2+1);
					for (unsigned int k = 0; k < Pdof; k++) //we go through the 4 enrichments
					{
						condensed_rows(k,i*4+0)= - D_matrix_mixed(k,i*3+0);//+mass_stabilization_terms(i*2+0);
						condensed_rows(k,i*4+1)= - D_matrix_mixed(k,i*3+1);//+mass_stabilization_terms(i*2+1);
						condensed_rows(k,i*4+2)= - D_matrix_mixed(k,i*3+2);//+mass_stabilization_terms(i*2+1);
						//condensed_rows(k,i*4+3)=mixed_Laplacian(k,i);
						
						condensed_columns(k,i*4+0)= - D_matrix_mixed(k,i*3+0);
						condensed_columns(k,i*4+1)= - D_matrix_mixed(k,i*3+1);
						condensed_columns(k,i*4+2)= - D_matrix_mixed(k,i*3+2);
						//condensed_columns(k,i*4+3)=mixed_Laplacian(k,i);
					}
				}
				
				boost::numeric::ublas::bounded_matrix<double, Pdof , Pdof  > condensed_block;
				condensed_block = - condensed_pressure_mass/delta_t;
				boost::numeric::ublas::bounded_matrix<double, Pdof , Pdof  > inverse_block = ZeroMatrix(Pdof,Pdof);
				
				
				this->InvertMatrix( condensed_block,  inverse_block);
				//inverse_block(0,0) = 1.0/(condensed_block(0,0));
				
				
				//rLeftHandSideMatrix=ZeroMatrix(9,9);
				boost::numeric::ublas::bounded_matrix<double, 16 , Pdof  > temp_matrix;
				temp_matrix = prod(trans(condensed_columns),inverse_block);
				
				rLeftHandSideMatrix -=  prod(temp_matrix,condensed_rows);
				noalias(rRightHandSideVector) -= prod(temp_matrix,rhs_enrich);
					


				for (unsigned int i = 0; i < 4; i++)
				{
					for (unsigned int j = 0; j < 4; j++)
					{
						rLeftHandSideMatrix(i*4+3, j*4+0 ) -= D_matrix(i,j*3) ;     
						rLeftHandSideMatrix(i*4+3, j*4+1 ) -= D_matrix(i,j*3+1) ;     
						rLeftHandSideMatrix(i*4+3, j*4+2 ) -= D_matrix(i,j*3+2) ;
						
						rLeftHandSideMatrix(j*4+0, i*4+3 ) -= D_matrix(i,j*3);     
						rLeftHandSideMatrix(j*4+1, i*4+3 ) -= D_matrix(i,j*3+1);
						rLeftHandSideMatrix(j*4+2, i*4+3 ) -= D_matrix(i,j*3+2);
						//if (has_fluid_particle)
						//rLeftHandSideMatrix(i*4+3, j*4+3 ) -= Laplacian_matrix(i,j);

					}
				}
				
				
				//and finally the RHS of the velocity plus the RHS on the pressure due to the stabilzation terms:
				array_1d<double,3> mean_velocity = ZeroVector(3);
				//double divergence_n = 0.0;
				//double mean_pressure= 0.0;
				for (unsigned int i = 0; i < 4; i++)
				{
					mean_velocity += GetGeometry()[i].FastGetSolutionStepValue(VELOCITY)*one_quarter;
					//divergence_n += one_third*Area*(DN_DX(i,0)*GetGeometry()[i].FastGetSolutionStepValue(VELOCITY_X)+DN_DX(i,1)*GetGeometry()[i].FastGetSolutionStepValue(VELOCITY_Y));
					//mean_pressure +=one_third*GetGeometry()[i].FastGetSolutionStepValue(PRESSURE);
				}	

				const double total_mass = lumped_mass(0)+lumped_mass(1)+lumped_mass(2)+lumped_mass(3);

				//double x_solid_force;
				//double y_solid_force;
				for (unsigned int i = 0; i < 4; i++)
				{
					rRightHandSideVector(i*4+0) += lumped_mass(i)*gravity(0);
					rRightHandSideVector(i*4+1) += lumped_mass(i)*gravity(1);
					rRightHandSideVector(i*4+2) += lumped_mass(i)*gravity(2);
					//if (distances(i)<0.0)
					//{
						for (unsigned int j = 0; j < ndivisions; j++) //partitions
						{
							//double  x_force = - ( gauss_gradients[j](i,0) * (stresses(0)-solid_pressure) + gauss_gradients[j](i,1)*stresses(2)) * volumes(j);
							//double y_force = - ( gauss_gradients[j](i,1) * (stresses(1)-solid_pressure) + gauss_gradients[j](i,0)*stresses(2)) * volumes(j);
							double x_force = 0.0;
							double y_force = 0.0;
							double z_force = 0.0;
							double fluid_stab=0.0;
							
							if (signs(j)<0.0)
							{
								x_force = - ( DN_DX(i,0) * (stresses(0)) + DN_DX(i,1)*stresses(3) + DN_DX(i,2)*stresses(5) ) * volumes(j);
								y_force = - ( DN_DX(i,1) * (stresses(1)) + DN_DX(i,0)*stresses(3) + DN_DX(i,2)*stresses(4) ) * volumes(j);
								z_force = - ( DN_DX(i,2) * (stresses(2)) + DN_DX(i,1)*stresses(4) + DN_DX(i,0)*stresses(5) ) * volumes(j);
								//if (distances(i)<0.0)
								{
									//x_force  += (gauss_gradients[j](i,0) * (solid_pressure) ) * volumes(j);
									//y_force  += (gauss_gradients[j](i,1) * (solid_pressure) ) * volumes(j);
									//x_force  += (DN_DX(i,0) * (solid_pressure) ) * volumes(j);
									//y_force  += (DN_DX(i,1) * (solid_pressure) ) * volumes(j);
								}
							}
							else
							{
								//if ( distances(i)>0.0)
								//	fluid_stab -= TauOne*(gauss_gradients[j](i,0)*gravity(0) + gauss_gradients[j](i,1)*gravity(1)  ) *volumes(j);
								//if ( distances(i)>0.0)
								//	fluid_stab -= TauOne*(DN_DX(i,0)*gravity(0)+DN_DX(i,1)*gravity(1))*volumes(j);
							}
							//fluid_stab -= TauOne*fluid_area*(DN_DX(i,0)*gravity(0)+DN_DX(i,1)*gravity(1)+DN_DX(i,2)*gravity(2));
							
							
							rRightHandSideVector(i*4+0) += x_force;
							rRightHandSideVector(i*4+1) += y_force;
							rRightHandSideVector(i*4+2) += z_force;
							rRightHandSideVector(i*4+3) += fluid_stab;

							//rRightHandSideVector(i*3+2) -= TauOneNode * ( gauss_gradients[j](i,0)*x_force + gauss_gradients[j](i,1)*y_force );
							
							//rRightHandSideVector(i*3+2) -= TauOneNode * ( gauss_gradients[j](i,0) * gravity(0) + gauss_gradients[j](i,1) * gravity(1) )  * (total_mass/Area) *volumes(i);
						}
					//}
					
					
				}
				
				
				array_1d<double,4> total_forces =ZeroVector(4);
				for (unsigned int i = 0; i < 4; i++)
				{
					total_forces(0) += rRightHandSideVector(i*4+0);
					total_forces(1) += rRightHandSideVector(i*4+1);
					total_forces(2) += rRightHandSideVector(i*4+2);
					//rRightHandSideVector(i*3+2) -= TauOne*distance_based_fluid_area*(DN_DX(i,0)*gravity(0)+DN_DX(i,1)*gravity(1));
				}
				
				//noalias(rRightHandSideVector) -= prod((Lumped_Pressure_Mass_matrix),(previous_vel_and_solid_press/delta_t));
				noalias(rLeftHandSideMatrix) -= 1.0*Lumped_Pressure_Mass_matrix/delta_t;  
				noalias(rRightHandSideVector) -= prod((Lumped_Pressure_Mass_matrix),(previous_vel_and_solid_press/delta_t));

				//noalias(rRightHandSideVector) -= prod((Lumped_Pressure_Mass_matrix),(previous_vel_and_press/delta_t));
				//KRATOS_WATCH(previous_vel_and_press);
				//KRATOS_WATCH(previous_vel_and_solid_press);
				//noalias(rLeftHandSideMatrix) -= Lumped_Pressure_Mass_matrix/delta_t;  
				
				const double gamma=1.0;//;765; //must be between 0.5 and 1. the higher the more dissipative
				noalias(rRightHandSideVector) += prod((Mass_matrix),((1.0/gamma)*previous_vel_and_press/delta_t+(1.0-gamma)/gamma*previous_accel));
				noalias(rLeftHandSideMatrix) += (1.0/gamma)*Mass_matrix/delta_t;   

				//KRATOS_WATCH(rLeftHandSideMatrix)
				//rLeftHandSideMatrix += 1.0*(Area*mu*delta_t/(mElemSize*mElemSize))*prod(trans(penalty_vector),penalty_vector);
				this->GetValue(TEMPERATURE)=0.0;
				/*
				KRATOS_WATCH(ndivisions)
				KRATOS_WATCH(rLeftHandSideMatrix);
				KRATOS_WATCH(rRightHandSideVector);
				for (unsigned int i = 0; i < ndivisions; i++)
					KRATOS_WATCH(gauss_gradients[i]);
				KRATOS_WATCH(condensed_block);
				KRATOS_WATCH(inverse_block);
				KRATOS_WATCH(Laplacian_enrich);
				*/
				//KRATOS_WATCH(inverse_block)
				//KRATOS_WATCH(rRightHandSideVector)

			}
			else // ( meaning if split_element==false) NON SPLIT FLUID ELEMENT
			{
				double Weight;
				double density;
				double Bulk_modulus;
				double element_distance= one_quarter*(distances(0)+distances(1)+distances(2)+distances(3));
				const double viscosity = viscosity_integral/fluid_area;
				bool is_solid_element=false;
				
				density = (density_solid_integral+density_fluid_integral)/Area;
				//density = density_fluid;
				//fluid_area=Area;
				Weight = one_quarter*density*Area;
				Bulk_modulus= Bulk_modulus_integral_solid/solid_area;
				this->AddElasticityTerm(rLeftHandSideMatrix,DN_DX, (mu_integral*delta_t));
				//rLeftHandSideMatrix*=0.5;
				this->AddViscousTerm(rLeftHandSideMatrix,DN_DX, (viscosity_integral));
				//noalias(rRightHandSideVector) -= 0.25*prod(rLeftHandSideMatrix,previous_vel_and_press);
				//rLeftHandSideMatrix*=0.75;
				
				
				if (has_solid_particle==true)
					is_solid_element=true;
					
				
				this->GetValue(TEMPERATURE)= (-solid_area+fluid_area)/Area;
				
				
				for (unsigned int i=0; i<TNumNodes ; i++)
				{
					//Mass_matrix (i*3,i*3) = Weight;
					//Mass_matrix (i*3+1,i*3+1) = Weight;
					//if (element_distance>1.0)
					//{
						Mass_matrix (i*4,i*4) = nodal_density(i)*Area*one_quarter; //nodal_mass_fluid(i)+nodal_mass_solid(i);
						Mass_matrix (i*4+1,i*4+1) = nodal_density(i)*Area*one_quarter; //nodal_mass_fluid(i)+nodal_mass_solid(i);
						Mass_matrix (i*4+2,i*4+2) = nodal_density(i)*Area*one_quarter; //nodal_mass_fluid(i)+nodal_mass_solid(i);
						
					//}	
					//Mass_matrix (i*3,i*3) = Weight;
					//Mass_matrix (i*3+1,i*3+1) = Weight;
					Lumped_Pressure_Mass_matrix(i*4+3,i*4+3) = one_quarter * solid_area /Bulk_modulus;
					for (unsigned int j=0; j<TNumNodes ; j++)
						Pressure_Mass_matrix(i*4+3,j*4+3) = one_quarter * 0.25 * (Area /Bulk_modulus );
					Pressure_Mass_matrix(i*4+3,i*4+3) = one_quarter * 0.5 * (Area /Bulk_modulus );
				}	
				
				/////////////////// calculating tau
				double AdvVelNorm = sqrt(pow(this->GetGeometry()[0].FastGetSolutionStepValue(VELOCITY_X),2)+pow(this->GetGeometry()[0].FastGetSolutionStepValue(VELOCITY_Y),2)+pow(this->GetGeometry()[0].FastGetSolutionStepValue(VELOCITY_Z),2));
				//const double TimeFactor = rCurrentProcessInfo.GetValue(DYNAMIC_TAU);
				double mElemSize;
				array_1d<double,3> Edge(3,0.0);
				Edge = this->GetGeometry()[1].Coordinates() - this->GetGeometry()[0].Coordinates();
				mElemSize = Edge[0]*Edge[0];
				for (SizeType d = 1; d < 3; d++)
					mElemSize += Edge[d]*Edge[d];
				for (SizeType i = 2; i < 4; i++)
					for(SizeType j = 0; j < i; j++)
					{
						Edge = this->GetGeometry()[i].Coordinates() - this->GetGeometry()[j].Coordinates();
						double Length = Edge[0]*Edge[0];
						for (SizeType d = 1; d < 3; d++)
							Length += Edge[d]*Edge[d];
						if (Length < mElemSize) mElemSize = Length;
					}
					
				mElemSize = sqrt(mElemSize);
				
				/*
				double TauOne = 0.0;
				if (is_solid_element==false) //it is a fluid element. otherwise we will add no stabilization
					TauOne = (1.0 / ( ( 1.0 / delta_t + 4.0 * viscosity / (density * mElemSize * mElemSize) + 2.0 * AdvVelNorm / mElemSize) ));
				else
					TauOne = 1.0 /  ( 10.0 * mu * delta_t / (density * mElemSize * mElemSize) + 2.0 * AdvVelNorm / mElemSize );
				*/
				//double TauOneFluid = 2.0 / ( ( 1.0 / delta_t + 4.0 * viscosity / (density * mElemSize * mElemSize)));// + 2.0 * AdvVelNorm / mElemSize) ));
				//if (density<10.0)
				//	TauOneFluid*=1.0;
				//changed definition to account for the "viscosity" of the elastic part, which is G*delta_t
				//this will reduce the TauOne, adding less artificial diffussion and improving mass conservation.
				const double apparent_viscosity = (viscosity * fluid_area +  100.0* mu*solid_area*delta_t)/Area;
				const double TauOneFluid = (1.0 / ( ( 1.0 / delta_t + 4.0 * apparent_viscosity / (density * mElemSize * mElemSize) + 2.0 * AdvVelNorm / mElemSize) ) );
				
				double TauOneSolid = 0.0 /  (1.0 / delta_t +  10.0 * mu * delta_t / (density * mElemSize * mElemSize) + 2.0 * AdvVelNorm / mElemSize );
				double TauOne = (solid_area*TauOneSolid+fluid_area*TauOneFluid)/Area;
				//double TauOne = TauOneFluid;
				//KRATOS_WATCH(TauOne);
				//if ((fluid_area+solid_area)<(Area*0.99))
				//	KRATOS_WATCH("wrong!");
					
				this->GetValue(TAU)=TauOne;
				this->GetValue(DENSITY)=density;
				this->GetValue(VISCOSITY)=viscosity;
				
				
				//if (element_distance>1.0 && element_distance<2.0)
				//	TauOne*= 1.0 - (1.0-element_distance)*(element_distance-2.0)*4.0*0.9;
				//if (fluid_area>(particle_area*0.5))
				//	fluid_area=Area;
				
				//if (has_solid_particle==false)
					Laplacian_matrix = prod(DN_DX,trans(DN_DX))/density;  // B B^T  (standard laplacian, now we must add the contribution of the condensed new nodes.
					//Laplacian_matrix += prod(DN_DX,trans(DN_DX))*solid_area*TauOne/density;  // B B^T  (standard laplacian, now we must add the contribution of the condensed new nodes.
				for (unsigned int i = 0; i < 4; i++)
				{
					for (unsigned int j = 0; j < 4; j++)
					{
						
						rLeftHandSideMatrix(i*4+3, j*4+0 ) -= (D_matrix(i,j*3) + add_accel_stab_term *TauOne * DN_DX(i,0)*one_quarter/delta_t)*Area;     
						rLeftHandSideMatrix(i*4+3, j*4+1 ) -= (D_matrix(i,j*3+1) + add_accel_stab_term *TauOne * DN_DX(i,1)*one_quarter*Area/delta_t)*Area;    
						rLeftHandSideMatrix(i*4+3, j*4+2 ) -= (D_matrix(i,j*3+2) + add_accel_stab_term *TauOne * DN_DX(i,2)*one_quarter*Area/delta_t)*Area;    
						
						rLeftHandSideMatrix(j*4+0, i*4+3 ) -=  D_matrix(i,j*3)*Area;     
						rLeftHandSideMatrix(j*4+1, i*4+3 ) -= D_matrix(i,j*3+1)*Area; 
						rLeftHandSideMatrix(j*4+2, i*4+3 ) -= D_matrix(i,j*3+2)*Area; 
						
						rLeftHandSideMatrix(i*4+3, j*4+3 ) -= Laplacian_matrix(i,j)*(fluid_area*TauOneFluid+solid_area*TauOneSolid);
						/*
						if (has_fluid_particle)
						{
							rLeftHandSideMatrix(i*4+2, j*4+0 ) -= (D_matrix(i,j*2) + add_accel_stab_term *TauOne * DN_DX(i,0)*one_third/delta_t)*Area;     
							rLeftHandSideMatrix(i*4+2, j*4+1 ) -= (D_matrix(i,j*2+1) + add_accel_stab_term *TauOne * DN_DX(i,1)*one_third*Area/delta_t)*Area;    
							
							rLeftHandSideMatrix(j*4+0, i*4+2 ) -=  G_matrix(i,j*2)*fluid_area;     
							rLeftHandSideMatrix(j*4+1, i*4+2 ) -= G_matrix(i,j*2+1)*fluid_area; 
							rLeftHandSideMatrix(i*4+2, j*4+2 ) -= Laplacian_matrix(i,j)*Area*TauOneFluid;
						}
												
						if (has_solid_particle)
						{
							rLeftHandSideMatrix(i*4+3, j*4+0 ) -= (D_matrix(i,j*2) + add_accel_stab_term *TauOne * DN_DX(i,0)*one_third/delta_t)*Area;     
							rLeftHandSideMatrix(i*4+3, j*4+1 ) -= (D_matrix(i,j*2+1) + add_accel_stab_term *TauOne * DN_DX(i,1)*one_third/delta_t)*Area;     
							
							rLeftHandSideMatrix(j*4+0, i*4+3 ) -=  G_matrix(i,j*2)*solid_area;     
							rLeftHandSideMatrix(j*4+1, i*4+3 ) -= G_matrix(i,j*2+1)*solid_area;
						
							rLeftHandSideMatrix(i*4+3, j*4+3 ) -= Laplacian_matrix(i,j)*Area*TauOneSolid;
							*/ 
						}
					}
					/*
					if (has_fluid_particle==false && ((GetGeometry()[i].FastGetSolutionStepValue(YP))-(GetGeometry()[i].FastGetSolutionStepValue(SOLID_YP)))<0.001 )
						rLeftHandSideMatrix(i*4+2, i*4+2 ) -= 1.0/(mElemSize)*Area*1.0;  // 1.0/mElemSize*Area is supposed to be in the order of D, this way we ensure an almost zero diagonal for a better incompressibility
					if (has_solid_particle==false  && GetGeometry()[i].FastGetSolutionStepValue(SOLID_YP)<0.001 )
						rLeftHandSideMatrix(i*4+3, i*4+3 ) -= 1.0/(mElemSize)*Area*1.0;	
					*/
				
				
				
				//and finally the RHS of the velocity plus the RHS on the pressure due to the stabilzation terms:
				array_1d<double,3> mean_velocity = ZeroVector(3);
				//double divergence_n = 0.0;
				double mean_pressure= 0.0;
				for (unsigned int i = 0; i < 4; i++)
				{
					mean_velocity += GetGeometry()[i].FastGetSolutionStepValue(VELOCITY)*one_quarter;
					//divergence_n += one_third*Area*(DN_DX(i,0)*GetGeometry()[i].FastGetSolutionStepValue(VELOCITY_X)+DN_DX(i,1)*GetGeometry()[i].FastGetSolutionStepValue(VELOCITY_Y));
					mean_pressure +=one_quarter*GetGeometry()[i].FastGetSolutionStepValue(PRESSURE);
				}	

				//array_1d<double,3> extra_solid_forces =ZeroVector(3);
				for (unsigned int i = 0; i < 4; i++)
				{
					rRightHandSideVector(i*4+0) += gravity(0)*Mass_matrix (i*4,i*4);
					rRightHandSideVector(i*4+1) += gravity(1)*Mass_matrix (i*4+1,i*4+1);
					rRightHandSideVector(i*4+2) += gravity(2)*Mass_matrix (i*4+2,i*4+2);
					if(is_solid_element)
					{
						
						rRightHandSideVector(i*4+0) -= (DN_DX(i,0)*total_integral_stresses(0)+DN_DX(i,1)*total_integral_stresses(3)+DN_DX(i,2)*total_integral_stresses(5));
						rRightHandSideVector(i*4+1) -= (DN_DX(i,1)*total_integral_stresses(1)+DN_DX(i,0)*total_integral_stresses(3)+DN_DX(i,2)*total_integral_stresses(4));
						rRightHandSideVector(i*4+2) -= (DN_DX(i,2)*total_integral_stresses(2)+DN_DX(i,1)*total_integral_stresses(4)+DN_DX(i,0)*total_integral_stresses(5));

						//rRightHandSideVector(i*3+0) -= (DN_DX(i,0)*stresses(0)+DN_DX(i,1)*stresses(2))*Area;
						//rRightHandSideVector(i*3+1) -= (DN_DX(i,1)*stresses(1)+DN_DX(i,0)*stresses(2))*Area;
						//rRightHandSideVector(i*3+0) += DN_DX(i,0) * (total_integral_mean_solid_pressure);
						//rRightHandSideVector(i*3+1) += DN_DX(i,1) * (total_integral_mean_solid_pressure);
						//rRightHandSideVector(i*3+0) += DN_DX(i,0) * (mean_pressure)*Area;
						//rRightHandSideVector(i*3+1) += DN_DX(i,1) * (mean_pressure)*Area;
						
						//extra_solid_forces(0) += - (DN_DX(i,0)*total_integral_stresses(0)+DN_DX(i,1)*total_integral_stresses(2)) +  DN_DX(i,0) * (total_integral_mean_solid_pressure);
						//extra_solid_forces(1) += - (DN_DX(i,1)*total_integral_stresses(1)+DN_DX(i,0)*total_integral_stresses(2)) +  DN_DX(i,1) * (total_integral_mean_solid_pressure);
						
						//rRightHandSideVector(i*3+0) += DN_DX(i,0) * (solid_pressure)*Area;
						//rRightHandSideVector(i*3+1) += DN_DX(i,1) * (solid_pressure)*Area;
						
						//rRightHandSideVector(i*3+2) -= Laplacian_matrix(i,i)*(solid_pressure);
					}
					//if(has_solid_particle)
					//{
					//	rRightHandSideVector(i*3+0) += (DN_DX(i,0)*total_integral_mean_solid_pressure)*Area/solid_area;
					//	rRightHandSideVector(i*3+1) += (DN_DX(i,1)*total_integral_mean_solid_pressure)*Area/solid_area;
					//}

					//if((this->GetGeometry()[0].X())>9.9999)
					//	rRightHandSideVector(i*3+1) += 0.001/34.0;
					//if (has_solid_particle==false)
					//{
					//rRightHandSideVector(i*3+2) -= TauOne*(DN_DX(i,0)*(gravity(0)-mean_velocity(0))+DN_DX(i,1)*(gravity(1)-mean_velocity(1)))*Area;
					//rRightHandSideVector(i*3+2) -= TauOne*(DN_DX(i,0)*gravity(0)+DN_DX(i,1)*gravity(1))*Area;
					//rRightHandSideVector(i*3+2) -= 3.0*TauOne*( DN_DX(i,0)*rRightHandSideVector(i*3+0)+DN_DX(i,1)*rRightHandSideVector(i*3+1)) ;

					//}
				}
				//having all the forces, we calculate the mean x and y forces:
				
				array_1d<double, 3 > vel_gauss = ZeroVector(3);
				//Vector pressure_rhs_stab(3);
				for (unsigned int i = 0; i < 4; i++)
				{
					array_1d<double,3>& node_press_proj = GetGeometry()[i].FastGetSolutionStepValue(PRESS_PROJ);
					vel_gauss[0] += node_press_proj(0)*one_quarter;
					vel_gauss[1] += node_press_proj(1)*one_quarter;
					vel_gauss[2] += node_press_proj(2)*one_quarter;
				}
				//noalias(pressure_rhs_stab) = prod(DN_DX, vel_gauss);
				//*/
				
				const bool use_press_proj=true;
				
				for (unsigned int i = 0; i < 4; i++) 
				{
					//if (is_solid_element==true)
					//	rRightHandSideVector(i*3+2) -= pressure_rhs_stab(i)*TauOne*Area;
					//else
					//rRightHandSideVector(i*3+2) -= TauOne*(DN_DX(i,0)*(Area*density*gravity(0)+add_accel_stab_term*(mean_velocity(0)*Area*density/delta_t))+DN_DX(i,1)*(Area*density*gravity(1)+add_accel_stab_term*mean_velocity(1)*Area*density/delta_t)) /density;
					if (use_press_proj)
					{
						rRightHandSideVector(i*4+3) -= TauOneFluid*(DN_DX(i,0)*vel_gauss(0)+DN_DX(i,1)*vel_gauss(1)+DN_DX(i,2)*vel_gauss(2))*fluid_area;
					}
					else
					{
						rRightHandSideVector(i*4+3) -= TauOneFluid*fluid_area*(DN_DX(i,0)*gravity(0)+DN_DX(i,1)*gravity(1)+DN_DX(i,2)*gravity(2));
						//rRightHandSideVector(i*3+2) -= TauOneSolid*(DN_DX(i,0)*(gravity(0)*solid_area+(extra_solid_forces(0))/density)+DN_DX(i,1)*(gravity(1)*solid_area+(extra_solid_forces(1))/density));
					}
				}
				

				double gamma=1.0; //must be between 0.5 and 1. the higher the more dissipative
				if (is_solid_element && has_positive_node==false)
					gamma=0.5;
					
				noalias(rRightHandSideVector) += prod((Mass_matrix),((1.0/gamma)*previous_vel_and_press/delta_t+(1.0-gamma)/gamma*previous_accel));
				noalias(rLeftHandSideMatrix) += (1.0/gamma)*Mass_matrix/delta_t;   
				//if (has_solid_particle)
				//{
					
				Lumped_Pressure_Mass_matrix *= 1.0+TauOneSolid;	
				if (is_solid_element)
				{	
					noalias(rRightHandSideVector) -= prod((Lumped_Pressure_Mass_matrix),(previous_vel_and_solid_press/delta_t));
					noalias(rLeftHandSideMatrix) -= Lumped_Pressure_Mass_matrix/delta_t;  
					//for (unsigned int i = 0; i < 3; i++) 
					//{
					//	rLeftHandSideMatrix(i*4+3,i*4+3) -= solid_pressure_from_particles(i)*Lumped_Pressure_Mass_matrix(i*4+3,i*4+3)/delta_t;
					//	rRightHandSideVector(i*4+3) -= solid_pressure_from_particles(i)*Lumped_Pressure_Mass_matrix(i*4+3,i*4+3)/delta_t;
					//}
					
					//KRATOS_WATCH(Bulk_modulus)
				}
				//KRATOS_WATCH(rRightHandSideVector)
				//}
				
			} //finished pure solid or fluid element.
			/*
			for (unsigned int i = 0; i < 3; i++)
			{
				for (unsigned int j = 0; j < 3; j++)
				{
					//rLeftHandSideMatrix(i*3+2,j*3+0) *=pressure_increasing_factor;
					//rLeftHandSideMatrix(i*3+2,j*3+1) *=pressure_increasing_factor;
				 	
					rLeftHandSideMatrix(i*4+0,j*3+2) *=pressure_increasing_factor;
					rLeftHandSideMatrix(i*4+1,j*3+2) *=pressure_increasing_factor;
					
					rLeftHandSideMatrix(i*4+3,j*4+3) *=pressure_increasing_factor;	
				}
				//rRightHandSideVector(i*3+2) *=pressure_increasing_factor;
			}
			*/
			
			noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix,previous_vel_and_press);
			
		} 
		else
		{	//nothing in the element. let's create 'air'
						//viscous term:
			KRATOS_WATCH("ADDING STOKES AIR BECAUSE NO PARTICLES WERE FOUND IN THIS ELEMENT")
			double viscosity = 0.0;
			double density = 1.0;
			const double one_quarter=1.0/3.0;
				
			const double Weight = one_quarter * Area * density;
			
			rLeftHandSideMatrix = (Weight/delta_t)*identity_matrix<double> (16);
			rRightHandSideVector = ZeroVector(16);
			
			boost::numeric::ublas::bounded_matrix<double, 4, 9 > D_matrix; //(divergence)
			noalias(D_matrix) = ZeroMatrix(4,9);	
			boost::numeric::ublas::bounded_matrix<double, 4, 4 > Laplacian_matrix; //(gradient)
			noalias(Laplacian_matrix) = ZeroMatrix(4,4);	

			
			//we start by calculating the non-enriched functions of divergence and gradient. as i write this, gradient integrated by parts.
			for (unsigned int i = 0; i < 4; i++)
			{
				for (unsigned int j = 0; j < 4; j++)
				{
					D_matrix(i, (j*3) ) =  DN_DX(j,0)*Area*one_quarter;     
					D_matrix(i, (j*3+1) ) = DN_DX(j,1)*Area*one_quarter;
					D_matrix(i, (j*3+2) ) = DN_DX(j,2)*Area*one_quarter;
				}
			}
			
			//rLeftHandSideMatrix += Mass_matrix;
			
			//PRESSURE TERMS:

				/////////////////// calculating tau
				double AdvVelNorm = sqrt(pow(this->GetGeometry()[0].FastGetSolutionStepValue(VELOCITY_X),2)+pow(this->GetGeometry()[0].FastGetSolutionStepValue(VELOCITY_Y),2)+pow(this->GetGeometry()[0].FastGetSolutionStepValue(VELOCITY_Z),2));
				//const double TimeFactor = rCurrentProcessInfo.GetValue(DYNAMIC_TAU);
				double mElemSize;
				array_1d<double,3> Edge(3,0.0);
				Edge = this->GetGeometry()[1].Coordinates() - this->GetGeometry()[0].Coordinates();
				mElemSize = Edge[0]*Edge[0];
				for (SizeType d = 1; d < 3; d++)
					mElemSize += Edge[d]*Edge[d];
				for (SizeType i = 2; i < 4; i++)
					for(SizeType j = 0; j < i; j++)
					{
						Edge = this->GetGeometry()[i].Coordinates() - this->GetGeometry()[j].Coordinates();
						double Length = Edge[0]*Edge[0];
						for (SizeType d = 1; d < 3; d++)
							Length += Edge[d]*Edge[d];
						if (Length < mElemSize) mElemSize = Length;
					}
					
				mElemSize = sqrt(mElemSize);
				
			const double   TauOne = 1.0 / ( ( 0.01 / delta_t + 4.0 * viscosity / (density * mElemSize * mElemSize) + 2.0 * AdvVelNorm / mElemSize) );

			
			Laplacian_matrix = prod(DN_DX,trans(DN_DX))*Area*TauOne/density;  // B B^T  (standard laplacian, now we must add the contribution of the condensed new nodes.

			for (unsigned int i = 0; i < 4; i++)
			{
				for (unsigned int j = 0; j < 4; j++)
				{

					rLeftHandSideMatrix(i*4+3, j*4+0 ) =- D_matrix(i,j*3);//  + DN_DX(i,0)*TauOne*one_third*Area/delta_t;     
					rLeftHandSideMatrix(i*4+3, j*4+1 ) =- D_matrix(i,j*3+1);// + DN_DX(i,1)*TauOne*one_third*Area/delta_t;
					rLeftHandSideMatrix(i*4+3, j*4+2 ) =- D_matrix(i,j*3+2);// + DN_DX(i,1)*TauOne*one_third*Area/delta_t;

					
					rLeftHandSideMatrix(j*4+0, i*4+3 ) =- D_matrix(i,j*3);     
					rLeftHandSideMatrix(j*4+1, i*4+3 ) =- D_matrix(i,j*3+1);
					rLeftHandSideMatrix(j*4+2, i*4+3 ) =- D_matrix(i,j*3+2);

					
					rLeftHandSideMatrix(i*4+3, j*4+3 ) =- Laplacian_matrix(i,j);
					//rLeftHandSideMatrix(i*4+3, j*4+3 ) =- Laplacian_matrix(i,j);
					//rLeftHandSideMatrix(i*3+2, j*3+2 ) -= one_third * Area / (100.0 * delta_t);
					
					
					
					
				}
				/*
				if (GetGeometry()[i].FastGetSolutionStepValue(SOLID_YP)<0.001)
						rLeftHandSideMatrix(i*4+3,i*4+3) = -1.0 * TauOne * Area /(mElemSize*mElemSize)/density;
				else
					rLeftHandSideMatrix(i*4+3,i*4+3) = 0.0;
					*/ 
				
			}
			
			
			//and finally the RHS of the velocity plus the RHS on the pressure due to the stabilzation terms:
			array_1d<double,16>  previous_vel_and_press = ZeroVector(16);
			for(unsigned int iii = 0; iii<4; iii++)
			{
				array_1d<double,3>& velocity =  GetGeometry()[iii].FastGetSolutionStepValue(VELOCITY);  //reference
				previous_vel_and_press(iii*4) = velocity[0];
				previous_vel_and_press(iii*4+1) = velocity[1];
				previous_vel_and_press(iii*4+2) = velocity[2];

				previous_vel_and_press(iii*4+3) = GetGeometry()[iii].FastGetSolutionStepValue(PRESSURE);
				//previous_vel_and_press(iii*4+3) = GetGeometry()[iii].FastGetSolutionStepValue(SOLID_PRESSURE);
			}
			
			
				
			for (unsigned int i = 0; i < 4; i++)
			{
				rRightHandSideVector(i*4+0) = one_quarter*Area*gravity(0)*density + previous_vel_and_press(i*4)*Weight /delta_t;
				rRightHandSideVector(i*4+1) = one_quarter*Area*gravity(1)*density + previous_vel_and_press(i*4+1)*Weight /delta_t;
				rRightHandSideVector(i*4+2) = one_quarter*Area*gravity(2)*density + previous_vel_and_press(i*4+2)*Weight /delta_t;

				//rRightHandSideVector(i*3+2) += divergence_n;
				//rRightHandSideVector(i*3+2) -= delta_t*DN_DX(i,1)*gravity(1)*Area*one_third*0.5;
				//rRightHandSideVector(i*3+2) -= TauOne*(DN_DX(i,0)*(gravity(0)+previous_vel_and_press(3*i+0)/delta_t)+DN_DX(i,1)*(gravity(1)+previous_vel_and_press(3*i+1)/delta_t))*Area;
				//rRightHandSideVector(i*3+2) -= TauOne*(DN_DX(i,0)*(gravity(0)-mean_velocity(0))+DN_DX(i,1)*(gravity(1)-mean_velocity(1)))*Area;
				rRightHandSideVector(i*4+3) =- TauOne*(DN_DX(i,0)*(gravity(0))+DN_DX(i,1)*(gravity(1))+DN_DX(i,2)*(gravity(2)))*Area;
				//rRightHandSideVector(i*4+3) =- TauOne*(DN_DX(i,0)*(gravity(0))+DN_DX(i,1)*(gravity(1)))*Area;
				//rRightHandSideVector(i*3+2) -= one_third * Area / (100.0 * delta_t) * GetGeometry()[i].FastGetSolutionStepValue(PRESSURE);

			}

			
			noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix,previous_vel_and_press);
			
			//************  choice 2: we just add garbage to the diagaonal of the matrix (small garbage)
			rLeftHandSideMatrix = ZeroMatrix(16,16);
			for (unsigned int i = 0; i < 4; i++)
			{
				//if (GetGeometry()[i].FastGetSolutionStepValue(YP)<0.001)
				{
					rLeftHandSideMatrix(i*4+0,i*4+0) = 1.0 * Area/delta_t ;
					rLeftHandSideMatrix(i*4+1,i*4+1) = 1.0 * Area/delta_t ;
					rLeftHandSideMatrix(i*4+2,i*4+2) = 1.0 * Area/delta_t ;

					rLeftHandSideMatrix(i*4+3,i*4+3) = -1.0 * Area /(mElemSize*mElemSize)*TauOne;
					//rLeftHandSideMatrix(i*4+3,i*4+3) = -1.0 * Area /(mElemSize*mElemSize)*TauOne;
				}
				//if (GetGeometry()[i].FastGetSolutionStepValue(SOLID_YP)<0.001)
				//	rLeftHandSideMatrix(i*3+3,i*4+3) = -1.0 * Area /(mElemSize*mElemSize)*TauOne;
				//for (unsigned int j = 0; j < 3; j++)
				//	rRightHandSideVector(i*3+j)=rLeftHandSideMatrix(i*3+j,i*3+j);
			}
			rLeftHandSideMatrix *= 0.0001;
			noalias(rRightHandSideVector) = ZeroVector(16);
			//previous_vel_and_press = ZeroVector(9);
			
			noalias(rRightHandSideVector) -= prod(rLeftHandSideMatrix,previous_vel_and_press);
			
			this->GetValue(TEMPERATURE)=2.0; 
		}
		

		
		KRATOS_CATCH("");
	}
	
	
	
	//************************************************************************************
	//************************************************************************************
	// this subroutine calculates the nodal contributions for the explicit steps of the 
	// fractional step procedure
	void FsiPFEM23D::InitializeSolutionStep(ProcessInfo& CurrentProcessInfo)
	{
		KRATOS_TRY
		KRATOS_CATCH("");
	}
	
	//************************************************************************************
	//************************************************************************************

	void FsiPFEM23D::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& CurrentProcessInfo)
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
	void FsiPFEM23D::GetDofList(DofsVectorType& ElementalDofList,ProcessInfo& CurrentProcessInfo)
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
	
	//non partitioned elements using MatrixType
	void FsiPFEM23D::AddElasticityTerm(MatrixType& OutputMatrix,
                                       const boost::numeric::ublas::bounded_matrix<double, 4, 3>& rShapeDeriv,
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
		
		boost::numeric::ublas::bounded_matrix<double, 6, 6 > C_matrix = ZeroMatrix(6,6);
		
		C_matrix(0,0)=4.0/3.0;
		C_matrix(1,1)=4.0/3.0;
		C_matrix(2,2)=4.0/3.0;
		C_matrix(3,3)=1.0;
		C_matrix(4,4)=1.0;
		C_matrix(5,5)=1.0;
		
		C_matrix(0,1)=-2.0/3.0;
		C_matrix(0,2)=-2.0/3.0;
		C_matrix(1,0)=-2.0/3.0;
		C_matrix(1,2)=-2.0/3.0;
		C_matrix(2,0)=-2.0/3.0;
		C_matrix(2,1)=-2.0/3.0;
		
		C_matrix *= Weight;
		
		boost::numeric::ublas::bounded_matrix<double, 6 , 12  > temp_matrix = prod(C_matrix,trans(B_matrix));
		boost::numeric::ublas::bounded_matrix<double, 12 , 12  > viscosity_matrix = prod(B_matrix, temp_matrix );
		for (unsigned int i=0; i!=4; i++) //i node
		{
			for (unsigned int j=0; j!=4; j++) //j neighbour
			{
				for (unsigned int k=0; k!=3; k++) //xyz
					for (unsigned int l=0; l!=3; l++) //xyz
						OutputMatrix(i*4+k,j*4+l)+=viscosity_matrix(i*3+k,j*3+l); 
			}
		}
		
		//OutputMatrix += viscosity_matrix;
	}
	

	
		//non partitioned elements using MatrixType
	void FsiPFEM23D::AddViscousTerm(MatrixType& rDampMatrix,
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
	}	
	
		//non partitioned elements using MatrixType
	void FsiPFEM23D::UpdateStressesToNewConfiguration(array_1d<double,6>& OutputVector,
									   double& pressure,
                                       const boost::numeric::ublas::bounded_matrix<double, 4, 3>& rShapeDeriv,
                                       const array_1d<double,12>&  previous_step_vel,
                                       const double ShearModulus,
                                       const double delta_t)
	{
		boost::numeric::ublas::bounded_matrix<double, 12, 1 > velocities;
		for (unsigned int i=0;i!=12;i++)
			velocities(i,0)=previous_step_vel(i);	 //hopefully this is something like V(n+1)
		
		
		//having added the new component, we do sigma_new_config = F.S.F where F is the DeformationGradient
		boost::numeric::ublas::bounded_matrix<double, 4 , 3 > displacements;
		for (unsigned int i=0;i!=4;i++)
		{
			displacements(i,0)=velocities(i*3+0,0)*delta_t;
			displacements(i,1)=velocities(i*3+1,0)*delta_t;
			displacements(i,2)=velocities(i*3+2,0)*delta_t;
		}	
		
		boost::numeric::ublas::bounded_matrix<double, 3 , 3 > DeformationGradient = identity_matrix<double> (3) ; //ZeroMatrix(2,2);// = prod(trans(displacements),rShapeDeriv);
		
		DeformationGradient+=prod(trans(displacements),rShapeDeriv);

        double J = MathUtils<double>::Det3( DeformationGradient );
		
		if (J<=0.8 )
		{
			double amplifying_factor=1.0;
			double reducing_factor;
			while (J<=0.8 )
			{
				reducing_factor=1.0/(1.0+0.1*amplifying_factor);
				DeformationGradient = identity_matrix<double> (3) + reducing_factor*prod(trans(displacements),rShapeDeriv);
				J = MathUtils<double>::Det3( DeformationGradient );
				amplifying_factor++;
			}
			std::cout << "elem " << this->Id() << ", reduced timestep by " << reducing_factor << std::endl;
			
		}
		boost::numeric::ublas::bounded_matrix<double, 3 , 3 > tensorial_stresses;
		tensorial_stresses(0,0)=OutputVector(0);
		tensorial_stresses(1,1)=OutputVector(1);
		tensorial_stresses(2,2)=OutputVector(2);
		tensorial_stresses(1,0)=OutputVector(3);
		tensorial_stresses(0,1)=OutputVector(3);
		tensorial_stresses(1,2)=OutputVector(4);
		tensorial_stresses(2,1)=OutputVector(4);
		tensorial_stresses(2,0)=OutputVector(5);
		tensorial_stresses(0,2)=OutputVector(5);
		
		boost::numeric::ublas::bounded_matrix<double, 3 , 3 > temp_matrix = prod(tensorial_stresses, trans(DeformationGradient));
		tensorial_stresses = prod((DeformationGradient),temp_matrix);
		tensorial_stresses /=J;
		OutputVector(0)=tensorial_stresses(0,0);
		OutputVector(1)=tensorial_stresses(1,1);
		OutputVector(2)=tensorial_stresses(2,2);
		OutputVector(3)=tensorial_stresses(1,0);
		OutputVector(4)=tensorial_stresses(1,2);
		OutputVector(5)=tensorial_stresses(2,0);
	}
	
	
	
	
	
	void FsiPFEM23D::CalculatePressureProjection(ProcessInfo& CurrentProcessInfo)
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
					array_1d<double,4>  pressures = ZeroVector(4); //to calculate the deformation Gradient F. Dimension = velocity dofs
					boost::numeric::ublas::bounded_matrix<double,4, 3 > coords; //coordinates of the nodes
					bool has_negative_node=false;
					bool has_positive_node=false;
					
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
					}
		
					bool split_element=false;
					if (has_positive_node && has_negative_node)
						split_element=true;
			
					if (has_negative_node==false) //we only calculate if we have a pure fluid element
					{

						
						boost::numeric::ublas::bounded_matrix<double, (3+1), (3-1)*6 > G_matrix; //(gradient)
						noalias(G_matrix) = ZeroMatrix((3+1), (3-1)*6);	
						
						double density = this->GetValue(DENSITY);

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
	
	
	
	void FsiPFEM23D::Calculate(const Variable<Vector> &rVariable,
                                     Vector &rOutput,
                                     const ProcessInfo &rCurrentProcessInfo)
	{
		
		
		if (rVariable == ELEMENT_MEAN_STRESS)
		{
				const bool use_failure_criteria=false;
			
				double delta_t = rCurrentProcessInfo[DELTA_TIME];	
				const int offset = rCurrentProcessInfo[PARTICLE_POINTERS_OFFSET];	
	
			//we will add the delta stresses to each of the particles inside the elements.
				Geometry<Node<3> >& geom = this->GetGeometry(); 
				
				ParticlePointerVector&  element_particle_pointers =  (this->GetValue(PARTICLE_POINTERS));
				int & number_of_particles_in_elem=this->GetValue(NUMBER_OF_PARTICLES);
				//std::cout << "elem " << this->Id() << " with " << (unsigned int)number_of_particles_in_elem << " particles" << std::endl;
				////////////////////////
				boost::numeric::ublas::bounded_matrix<double, 12, 6 > B_matrix = ZeroMatrix(12,6);
				boost::numeric::ublas::bounded_matrix<double, 12, 1 >velocities = ZeroMatrix(12,1);
				double Area;
				array_1d<double, 4> N;
				array_1d<double, 4> distances;
				boost::numeric::ublas::bounded_matrix<double, 4, 3> rShapeDeriv;
				GeometryUtils::CalculateGeometryData(geom, rShapeDeriv, N, Area);
				

				for (unsigned int i=0; i!=4; i++) //i node
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
						
						array_1d<double, 3 >velocity =  ( geom[i].FastGetSolutionStepValue(VELOCITY)) * 1.0 ;//+geom[i].FastGetSolutionStepValue(MESH_VELOCITY) ) * 0.5;
						velocities(i*3+0,0)=velocity(0);
						velocities(i*3+1,0)=velocity(1);
						velocities(i*3+2,0)=velocity(2);
						
						distances(i)=geom[i].FastGetSolutionStepValue(DISTANCE);
				}
				
				double reduced_delta_t=delta_t;
				//CalculateReducedTimeStepForPositiveJacobian(rShapeDeriv,velocities,delta_t,reduced_delta_t);
				
				
				double number_of_solid_particles=1.e-10; //almost zero
				Vector & element_stresses= this->GetValue(ELEMENT_MEAN_STRESS);
				array_1d<double,6> sum_particle_stresses=ZeroVector(6);
				
				for (unsigned int iii=0; iii<number_of_particles_in_elem ; iii++ )
				{
					//KRATOS_WATCH(iii)
					//if (iii>mmaximum_number_of_particles) //it means we are out of our portion of the array, abort loop!
					//	break; 

					PFEM_Particle & pparticle = element_particle_pointers[offset+iii];
					
					
					bool erase_flag= pparticle.GetEraseFlag();
					if (erase_flag==false)
					{
							
						if ( (pparticle.GetDistance())<0.0)
						{
							UpdateParticlePressure(pparticle,geom);
							
							UpdateParticleStresses(pparticle,geom,velocities,distances,B_matrix,reduced_delta_t);

							if(use_failure_criteria)
								TestParticleWithDruckerPrager(pparticle);
							
							number_of_solid_particles++;
							sum_particle_stresses += pparticle.GetSigma();
						}	
							
					}
					
					
				}
				element_stresses = sum_particle_stresses/number_of_solid_particles;
				
			
			
		}
	}
	
		
	void FsiPFEM23D::UpdateParticlePressure(
						 PFEM_Particle & pparticle,
						 Geometry< Node<3> >& geom)
	{
		//bool is_found;

		//array_1d<double,3> position;
		array_1d<double,4> N;

		//we start with the first position, then it will enter the loop.
		array_1d<double,3> coords = pparticle.Coordinates();// + (pparticle)->FastGetSolutionStepValue(DISPLACEMENT); //initial coordinates

		bool is_found = CalculatePosition(geom,coords[0],coords[1],coords[2],N);
		if(is_found == false)
		{
			KRATOS_WATCH(N)
			for (int j=0 ; j!=(4); j++)
								if (N[j]<0.0 && N[j]> -1e-5)
									N[j]=1e-10;
			//KRATOS_THROW_ERROR(std::logic_error, "PARTICLE IN WRONG ELEMENT!", "");
		}
		
		double pressure_change = 0.0;
		double total_N = 1e-6;
		//bool have_water_node=false;
		for(unsigned int j=0; j<(4); j++)
		{
			if(geom[j].FastGetSolutionStepValue(DISTANCE)<0.0)
			{
				total_N += N[j];
				pressure_change += (geom[j].FastGetSolutionStepValue(PRESSURE)-geom[j].FastGetSolutionStepValue(PRESSUREAUX))*N[j];
			}
		}
		
		pressure_change/=total_N;
		
		pparticle.GetPressure() = pparticle.GetPressure()+pressure_change;

	}
	

	void FsiPFEM23D::UpdateParticleStresses(
						 PFEM_Particle & pparticle,
						 Geometry< Node<3> >& geom,
						 const boost::numeric::ublas::bounded_matrix<double, 12, 1 >& velocities,
						 const array_1d<double,4> & distances,
						 const boost::numeric::ublas::bounded_matrix<double, 12, 6 >& B_matrix,
						 const double delta_t
						 )
	{
		//bool is_found;

		//array_1d<double,3> position;
		array_1d<double,4> N;

		
		//we start with the first position, then it will enter the loop.
		array_1d<double,3> coords = pparticle.Coordinates();// + (pparticle)->FastGetSolutionStepValue(DISPLACEMENT); //initial coordinates

		bool is_found = CalculatePosition(geom,coords[0],coords[1],coords[2],N);
		if(is_found == false)
		{
			KRATOS_WATCH(N)
			for (int j=0 ; j!=(3); j++)
								if (N[j]<0.0 && N[j]> -1e-5)
									N[j]=1e-10;
			//KRATOS_THROW_ERROR(std::logic_error, "PARTICLE IN WRONG ELEMENT!", "");
		}
		
		double mesh_distance = 0.0;
		for (unsigned int i=0;i!=4;i++)
			mesh_distance+=N(i)*distances(i);
			
		//if(mesh_distance<0.0)	
		const int TDim=3;
		if(true)
		{
			boost::numeric::ublas::bounded_matrix<double, (TDim-1)*3, (TDim-1)*3 > C_matrix = ZeroMatrix((TDim-1)*3,(TDim-1)*3);
			
			C_matrix(0,0)=4.0/3.0;
			C_matrix(1,1)=4.0/3.0;
			C_matrix(2,2)=4.0/3.0;
			C_matrix(3,3)=1.0;
			C_matrix(4,4)=1.0;
			C_matrix(5,5)=1.0;
			
			C_matrix(0,1)=-2.0/3.0;
			C_matrix(0,2)=-2.0/3.0;
			C_matrix(1,0)=-2.0/3.0;
			C_matrix(1,2)=-2.0/3.0;
			C_matrix(2,0)=-2.0/3.0;
			C_matrix(2,1)=-2.0/3.0;
			
			C_matrix *= pparticle.GetShearModulus();
		
			boost::numeric::ublas::bounded_matrix<double, 6, 1 > delta_strains = prod(trans(B_matrix),velocities)*delta_t;
			
			boost::numeric::ublas::bounded_matrix<double, 6, 1 > delta_stresses = prod(C_matrix,delta_strains); //for the moment is DELTA_stresses
			
			
			//const double delta_stresses_pressure = - (delta_stresses(0,0)+delta_stresses(1,0))*0.5;
			//if (coords[0]<0.05 && coords[1]>0.95)
			//{
			//	KRATOS_WATCH(delta_stresses_pressure);
			//	KRATOS_WATCH(delta_stresses(0,0));
			//	KRATOS_WATCH(delta_stresses(1,0));
			//}
			
			/*
			for (unsigned int i=0;i!=2;i++)
			{
				delta_stresses(i,0) +=delta_stresses_pressure;
			}
			*/
			array_1d<double,6>& particle_stress = pparticle.GetSigma();
			for (unsigned int i=0;i<6;i++)
				//particle_stress(i)=(pelement->GetValue(ELEMENT_MEAN_STRESS))(i)+delta_stresses(i,0);
				particle_stress(i) += delta_stresses(i,0); 
				
			pparticle.HasUpdatedStresses()=true;
		}
		else
			pparticle.HasUpdatedStresses()=false;
			

	}
	
	void FsiPFEM23D::TestParticleWithDruckerPrager(
						 PFEM_Particle & pparticle
						 )
	{
			array_1d<double,6>& particle_stress = pparticle.GetSigma();
			const double& particle_pressure = pparticle.GetPressure();
			const double von_misses_stress =  sqrt(1.0/3.0* (0.25 * pow((particle_stress(0)-particle_stress(1)),2) + pow((particle_stress(2)),2) ) ) ;

		    const double  theta =  pparticle.GetTheta(); //radians
		    double& cohesion =  pparticle.GetCohesion();
		    //if (pparticle.GetDistance()> -0.9);
			//	cohesion = 0.0;
			double A = 6.0 * cohesion *sin(theta) / ( 3.0 - sin(theta));
			double B = 2.0 *sin(theta) / ( 3.0 - sin(theta));
			
			double failure = von_misses_stress - A - B * particle_pressure;
			if (failure>0.0)
			{
				pparticle.GetDistance()= 1.0; //now it's liquid!
				pparticle.GetShearModulus() = 100.0;				
			}
				
			

	}
	
	void FsiPFEM23D::CalculateReducedTimeStepForPositiveJacobian(const boost::numeric::ublas::bounded_matrix<double, 4, 3>& rShapeDeriv,const boost::numeric::ublas::bounded_matrix<double, 12, 1 >& velocities, const double& delta_t, double& reduced_delta_t)
    {
		boost::numeric::ublas::bounded_matrix<double, 4 , 3 > displacements;
		for (unsigned int i=0;i!=4;i++)
		{
			displacements(i,0)=velocities(i*3+0,0)*delta_t;
			displacements(i,1)=velocities(i*3+1,0)*delta_t;
			displacements(i,2)=velocities(i*3+2,0)*delta_t;
		}	
		
		boost::numeric::ublas::bounded_matrix<double, 3 , 3 > DeformationGradient = identity_matrix<double> (3) ; //ZeroMatrix(2,2);// = prod(trans(displacements),rShapeDeriv);
		
		DeformationGradient+=prod(trans(displacements),rShapeDeriv);

        double J = MathUtils<double>::Det3( DeformationGradient );
		
		double reducing_factor=1.0;
		if (J<=0.8 )
		{
			double amplifying_factor=1.0;
			
			while (J<=0.8 )
			{
				reducing_factor=1.0/(1.0+0.1*amplifying_factor);
				DeformationGradient = identity_matrix<double> (3) + reducing_factor*prod(trans(displacements),rShapeDeriv);
				J = MathUtils<double>::Det3( DeformationGradient );
				amplifying_factor++;
			}
			//std::cout << "elem " << this->Id() << ", reduced timestep by " << reducing_factor << std::endl;
			
		}
		
		reduced_delta_t = reducing_factor*delta_t;
		
	}

	inline bool FsiPFEM23D::CalculatePosition(const bounded_matrix<double, 4, 3 > & coordinates,
                const double xc, const double yc, const double zc,
                array_1d<double, 4 > & N
                )
        {

            double x0 = coordinates(0,0);
			double y0 = coordinates(0,1);
			double z0 = coordinates(0,2);
			double x1 = coordinates(1,0);
			double y1 = coordinates(1,1);
			double z1 = coordinates(1,2);
			double x2 = coordinates(2,0);
			double y2 = coordinates(2,1);
			double z2 = coordinates(2,2);
			double x3 = coordinates(3,0);
			double y3 = coordinates(3,1);
			double z3 = coordinates(3,2);

            double vol = CalculateVol(x0, y0, z0, x1, y1, z1, x2, y2, z2, x3, y3, z3);

            double inv_vol = 0.0;
            if (vol <= 0.0)
            {
				KRATOS_WATCH(coordinates)
				KRATOS_WATCH(vol)
                KRATOS_THROW_ERROR(std::logic_error, "element with zero vol found", "");
            } else
            {
                inv_vol = 1.0 / vol;
            }

            N[0] = CalculateVol(x1, y1, z1, x3, y3, z3, x2, y2, z2, xc, yc, zc) * inv_vol;
            N[1] = CalculateVol(x0, y0, z0, x1, y1, z1, x2, y2, z2, xc, yc, zc) * inv_vol;
            N[2] = CalculateVol(x3, y3, z3, x1, y1, z1, x0, y0, z0, xc, yc, zc) * inv_vol;
            N[3] = CalculateVol(x3, y3, z3, x0, y0, z0, x2, y2, z2, xc, yc, zc) * inv_vol;


            if (N[0] >= 0.0 && N[1] >= 0.0 && N[2] >= 0.0 && N[3] >= 0.0 &&
                    N[0] <= 1.0 && N[1] <= 1.0 && N[2] <= 1.0 && N[3] <= 1.0)
                //if the xc yc zc is inside the tetrahedron return true
                return true;

            return false;
        }
        
        inline bool FsiPFEM23D::CalculatePosition(Geometry<Node < 3 > >&geom,
                const double xc, const double yc, const double zc,
                array_1d<double, 4 > & N
                )
        {
            double x0 = geom[0].X();
            double y0 = geom[0].Y();
            double z0 = geom[0].Z();
            double x1 = geom[1].X();
            double y1 = geom[1].Y();
            double z1 = geom[1].Z();
            double x2 = geom[2].X();
            double y2 = geom[2].Y();
            double z2 = geom[2].Z();
            double x3 = geom[3].X();
            double y3 = geom[3].Y();
            double z3 = geom[3].Z();

            double vol = CalculateVol(x0, y0, z0, x1, y1, z1, x2, y2, z2, x3, y3, z3);

            double inv_vol = 0.0;
            if (vol <= 0.0)
            {
                KRATOS_THROW_ERROR(std::logic_error, "element with zero vol found", "");
            } else
            {
                inv_vol = 1.0 / vol;
            }

            N[0] = CalculateVol(x1, y1, z1, x3, y3, z3, x2, y2, z2, xc, yc, zc) * inv_vol;
            N[1] = CalculateVol(x0, y0, z0, x1, y1, z1, x2, y2, z2, xc, yc, zc) * inv_vol;
            N[2] = CalculateVol(x3, y3, z3, x1, y1, z1, x0, y0, z0, xc, yc, zc) * inv_vol;
            N[3] = CalculateVol(x3, y3, z3, x0, y0, z0, x2, y2, z2, xc, yc, zc) * inv_vol;


            if (N[0] >= 0.0 && N[1] >= 0.0 && N[2] >= 0.0 && N[3] >= 0.0 &&
                    N[0] <= 1.0 && N[1] <= 1.0 && N[2] <= 1.0 && N[3] <= 1.0)
                //if the xc yc zc is inside the tetrahedron return true
                return true;

            return false;
        }


    inline double FsiPFEM23D::CalculateVol(const double x0, const double y0, const double z0,
                const double x1, const double y1, const double z1,
                const double x2, const double y2, const double z2,
                const double x3, const double y3, const double z3
                )
        {
            double x10 = x1 - x0;
            double y10 = y1 - y0;
            double z10 = z1 - z0;

            double x20 = x2 - x0;
            double y20 = y2 - y0;
            double z20 = z2 - z0;

            double x30 = x3 - x0;
            double y30 = y3 - y0;
            double z30 = z3 - z0;

            double detJ = x10 * y20 * z30 - x10 * y30 * z20 + y10 * z20 * x30 - y10 * x20 * z30 + z10 * x20 * y30 - z10 * y20 * x30;
            return detJ * 0.1666666666666666666667;
        }
    
    void FsiPFEM23D::GetValueOnIntegrationPoints(const Variable<double>& rVariable,
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

	template<class T>
	bool FsiPFEM23D::InvertMatrix(const T& input, T& inverse)
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
