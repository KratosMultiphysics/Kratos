/*
==============================================================================
KratosStructuralApplication
A library based on:
Kratos
A General Purpose Software for Multi-Physics Finite Element Analysis
Version 1.0 (Released on march 05, 2007).

Copyright 2007
Pooyan Dadvand, Riccardo Rossi, Janosch Stascheit, Felix Nagel
pooyan@cimne.upc.edu
rrossi@cimne.upc.edu
janosch.stascheit@rub.de
nagel@sd.rub.de
- CIMNE (International Center for Numerical Methods in Engineering),
Gran Capita' s/n, 08034 Barcelona, Spain
- Ruhr-University Bochum, Institute for Structural Mechanics, Germany


Permission is hereby granted, free  of charge, to any person obtaining
a  copy  of this  software  and  associated  documentation files  (the
"Software"), to  deal in  the Software without  restriction, including
without limitation  the rights to  use, copy, modify,  merge, publish,
distribute,  sublicense and/or  sell copies  of the  Software,  and to
permit persons to whom the Software  is furnished to do so, subject to
the following condition:

Distribution of this code for  any  commercial purpose  is permissible
ONLY BY DIRECT ARRANGEMENT WITH THE COPYRIGHT OWNERS.

The  above  copyright  notice  and  this permission  notice  shall  be
included in all copies or substantial portions of the Software.

THE  SOFTWARE IS  PROVIDED  "AS  IS", WITHOUT  WARRANTY  OF ANY  KIND,
EXPRESS OR  IMPLIED, INCLUDING  BUT NOT LIMITED  TO THE  WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT  SHALL THE AUTHORS OR COPYRIGHT HOLDERS  BE LIABLE FOR ANY
CLAIM, DAMAGES OR  OTHER LIABILITY, WHETHER IN AN  ACTION OF CONTRACT,
TORT  OR OTHERWISE, ARISING  FROM, OUT  OF OR  IN CONNECTION  WITH THE
SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

==============================================================================
 */

//
//   Project Name:        Kratos
//   Last modified by:    $Author: M.Santasusana $
//   Date:                $Date: 2014-03-01 14:39:59 $
//   Revision:            $Revision: 1.27 $
//
//


// System includes

// External includes


// Project includes
#include "custom_elements/DEM_FEM_Particle.h"
#include "includes/define.h"
#include "custom_utilities/GeometryFunctions.h"
#include "includes/constitutive_law.h"
#include "DEM_application.h"

//#include <omp.h>

namespace Kratos
{
    using namespace GeometryFunctions;

    DEM_FEM_Particle::DEM_FEM_Particle( IndexType NewId, GeometryType::Pointer pGeometry )
            : SphericParticle( NewId, pGeometry )
    {
        //DO NOT ADD DOFS HERE!!!
    }

    //************************************************************************************
    //************************************************************************************

    DEM_FEM_Particle::DEM_FEM_Particle( IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties )
            : SphericParticle( NewId, pGeometry, pProperties )
    {
        //         const unsigned int dim = GetGeometry().WorkingSpaceDimension();

    }

    Element::Pointer DEM_FEM_Particle::Create( IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties ) const
    {
        return Element::Pointer( new DEM_FEM_Particle( NewId, GetGeometry().Create( ThisNodes ), pProperties ) );
    }

    DEM_FEM_Particle::~DEM_FEM_Particle()
    {
    }


    void DEM_FEM_Particle::Initialize()
    {

        SphericParticle::Initialize();
		
		SetInitialBallNeighbor();
		SetInitialRigidFaceNeighbor();

        //Extra functions;
		mfcohesion = this->GetGeometry()(0)->FastGetSolutionStepValue(PARTICLE_COHESION);
		mftension  = this->GetGeometry()(0)->FastGetSolutionStepValue(PARTICLE_TENSION );
		
    }

    void DEM_FEM_Particle::SetInitialBallNeighbor()
	{
		
        ParticleWeakVectorType& rNeighbours  = this->GetValue(NEIGHBOUR_ELEMENTS);
        unsigned int new_size                = rNeighbours.size();
	   
	    maInitialBallNeighborID.resize(new_size);
		maInitialBallNeighborFailureType.resize(new_size);
		
	   unsigned int neighbour_counter = 0;
       
       for (ParticleWeakIteratorType i = rNeighbours.begin(); i != rNeighbours.end(); i++){

           maInitialBallNeighborID[neighbour_counter] = static_cast<int>(i->Id());
           maInitialBallNeighborFailureType[neighbour_counter] = 0;

		   neighbour_counter++;
        }
		
	}
	
	
	
     //////RigidFace
	void DEM_FEM_Particle::SetInitialRigidFaceNeighbor()
     {

       ConditionWeakVectorType& rNeighbours  = this->GetValue(NEIGHBOUR_RIGID_FACES);
       unsigned int new_size                = rNeighbours.size();
	   
		maInitialRigidFaceNeighborID.resize(new_size);
		maInitialRigidFaceNeighborFailureType.resize(new_size);
		
		
		unsigned int neighbour_counter = 0;
       
       for (ConditionWeakIteratorType i = rNeighbours.begin(); i != rNeighbours.end(); i++)
		{

           maInitialRigidFaceNeighborID[neighbour_counter] = static_cast<int>(i->Id());
           maInitialRigidFaceNeighborFailureType[neighbour_counter] = 0;

		   neighbour_counter++;
        }
      }
	
	
	
	
	


    ////************************************************************************************
    ////************************************************************************************

    void DEM_FEM_Particle::InitializeSolutionStep( ProcessInfo& CurrentProcessInfo )
    {
        SphericParticle::InitializeSolutionStep(CurrentProcessInfo);
		
		if (CurrentProcessInfo[VIRTUAL_MASS_OPTION])
		{
	
          //double& mass              = GetGeometry()(0)->FastGetSolutionStepValue(NODAL_MASS);
          double& moment_of_inertia = GetGeometry()(0)->FastGetSolutionStepValue(PARTICLE_MOMENT_OF_INERTIA);
		  
		  double mass  = mYoung * M_PI * mRadius; 
                  GetGeometry()(0)->FastGetSolutionStepValue(SQRT_OF_MASS) = sqrt(mass);
		  
		  if(mRotationOption)
		  {
			  mass = mass * 2.5;
		  }
		  moment_of_inertia = 0.4 * mass * mRadius * mRadius;  
		  
		  mMomentOfInertia = moment_of_inertia;
		}

    }


    ////************************************************************************************
    ////************************************************************************************

    void DEM_FEM_Particle::FinalizeSolutionStep( ProcessInfo& CurrentProcessInfo )
    {
		SphericParticle::FinalizeSolutionStep(CurrentProcessInfo);
    }


  void DEM_FEM_Particle::ComputeBallToBallContactForce(array_1d<double, 3>& rContactForce,
                                                          array_1d<double, 3>& rContactMoment,
                                                          array_1d<double, 3>& rElasticForce,
                                                          array_1d<double, 3>& rInitialRotaMoment,
                                                          ProcessInfo& rCurrentProcessInfo)
      {
		  
		  
          KRATOS_TRY

          ParticleWeakVectorType& rNeighbours    = this->GetValue(NEIGHBOUR_ELEMENTS);

          // KINEMATICS

          double dt = rCurrentProcessInfo[DELTA_TIME];
          double dt_i = 1 / dt;

          const array_1d<double, 3>& vel         = this->GetGeometry()(0)->FastGetSolutionStepValue(VELOCITY);
          const array_1d<double, 3>& delta_displ = this->GetGeometry()(0)->FastGetSolutionStepValue(DELTA_DISPLACEMENT);
          const array_1d<double, 3>& ang_vel     = this->GetGeometry()(0)->FastGetSolutionStepValue(ANGULAR_VELOCITY);
          double RotaAcc[3]                      = {0.0};
          double InitialRotaMoment[3]            = {0.0};

          if (mRotationOption){
              RotaAcc[0]                         = ang_vel[0] * dt_i;
              RotaAcc[1]                         = ang_vel[1] * dt_i;
              RotaAcc[2]                         = ang_vel[2] * dt_i;

              InitialRotaMoment[0]               = RotaAcc[0] * mMomentOfInertia;
              InitialRotaMoment[1]               = RotaAcc[1] * mMomentOfInertia;
              InitialRotaMoment[2]               = RotaAcc[2] * mMomentOfInertia;
          }
		  

          //KRATOS_WATCH("Ball2ball entered") //SALVA
          //LOOP OVER NEIGHBOURS BEGINS
		  
          size_t i_neighbour_count = 0;

          for (ParticleWeakIteratorType neighbour_iterator = rNeighbours.begin(); neighbour_iterator != rNeighbours.end(); neighbour_iterator++)
		 {
              
              // BASIC CALCULATIONS              
              array_1d<double, 3> other_to_me_vect    = this->GetGeometry()(0)->Coordinates() - neighbour_iterator->GetGeometry()(0)->Coordinates();
              double distance                       = sqrt(other_to_me_vect[0] * other_to_me_vect[0] +
                                                          other_to_me_vect[1] * other_to_me_vect[1] +
                                                          other_to_me_vect[2] * other_to_me_vect[2]);
              const double &other_radius              = neighbour_iterator->GetGeometry()(0)->FastGetSolutionStepValue(RADIUS);
 
              //double distance                         = sqrt(other_to_me_vect[0] * other_to_me_vect[0] + other_to_me_vect[1] * other_to_me_vect[1] + other_to_me_vect[2] * other_to_me_vect[2]);
              double radius_sum                       = mRadius + other_radius;
			  
			  
              double equiv_radius                     = radius_sum * 0.5;
             // double indentation                      = radius_sum - distance;
	      double kn;
	      double kt;
              
              double equiv_area                       =  M_PI * equiv_radius * equiv_radius;
			  
			  
	      double other_cohesion = neighbour_iterator->GetGeometry()(0)->FastGetSolutionStepValue(PARTICLE_COHESION);
	      double other_tension  = neighbour_iterator->GetGeometry()(0)->FastGetSolutionStepValue(PARTICLE_TENSION );
		  
	      double equiv_cohesion = (mfcohesion + other_cohesion) * 0.5;
	      double equiv_tension  = (mftension  + other_tension ) * 0.5;

              double equiv_visco_damp_coeff_normal;
              double equiv_visco_damp_coeff_tangential;
              double equiv_tg_of_fri_ang;

              CalculateEquivalentConstitutiveParameters(other_to_me_vect, other_radius, radius_sum, kn, kt, equiv_visco_damp_coeff_normal, equiv_visco_damp_coeff_tangential, equiv_tg_of_fri_ang, neighbour_iterator);
              
              double DeltDisp[3]                       = {0.0};
              double LocalDeltDisp[3]                  = {0.0};
              double RelVel[3]                         = {0.0};
              double LocalRelVel[3]                    = {0.0};
              //double NormalDir[3]                      = {0.0};
              //double OldNormalDir[3]                   = {0.0};
              double LocalCoordSystem[3][3]            = {{0.0}, {0.0}, {0.0}};
              double OldLocalCoordSystem[3][3]         = {{0.0}, {0.0}, {0.0}};
			  
			  
	      const double &other_young       = neighbour_iterator->GetGeometry()(0)->FastGetSolutionStepValue(YOUNG_MODULUS);
	      const double &other_poisson     = neighbour_iterator->GetGeometry()(0)->FastGetSolutionStepValue(POISSON_RATIO);
	      double equiv_young              = (mYoung + other_young) * 0.5;
	      double equiv_poisson            = (mPoisson + other_poisson) * 0.5;
	      double equiv_shearM             = equiv_young / (2.0 * (1.0 + equiv_poisson));
	      kn                              = equiv_young  * equiv_area / radius_sum;
	      kt                              = equiv_shearM * equiv_area / radius_sum;
			  

              EvaluateDeltaDisplacement(DeltDisp, RelVel, /*NormalDir, OldNormalDir,*/ LocalCoordSystem, OldLocalCoordSystem, other_to_me_vect, vel, delta_displ, neighbour_iterator, distance);

              if (mRotationOption)
			  {
                  DisplacementDueToRotation(DeltDisp, /*NormalDir,*/ LocalCoordSystem, other_radius, dt, ang_vel, neighbour_iterator);
			  }


              double LocalElasticContactForce[3]       = {0.0};
              double GlobalElasticContactForce[3]      = {0.0};


              GlobalElasticContactForce[0]             = mOldNeighbourContactForces[i_neighbour_count][0];
              GlobalElasticContactForce[1]             = mOldNeighbourContactForces[i_neighbour_count][1];
              GlobalElasticContactForce[2]             = mOldNeighbourContactForces[i_neighbour_count][2];

              GeometryFunctions::VectorGlobal2Local(LocalCoordSystem, GlobalElasticContactForce, LocalElasticContactForce); // Here we recover the old local forces projected in the new coordinates in the way they were in the old ones; Now they will be increased if its the necessary
              GeometryFunctions::VectorGlobal2Local(LocalCoordSystem, DeltDisp, LocalDeltDisp);
              GeometryFunctions::VectorGlobal2Local(LocalCoordSystem, RelVel, LocalRelVel);



            LocalElasticContactForce[0] +=  - kt * LocalDeltDisp[0];
            LocalElasticContactForce[1] +=  - kt * LocalDeltDisp[1];
            LocalElasticContactForce[2] +=  - kn * LocalDeltDisp[2];
			
			
			
			
			////Check if_initial_neighbor
			unsigned int inb;
			for(inb = 0; inb < maInitialBallNeighborID.size(); inb++)
			{
				if( maInitialBallNeighborID[inb] == static_cast<int>(neighbour_iterator->Id()) )
				{
					if(maInitialBallNeighborFailureType[inb] > 0)
					{
					  equiv_cohesion = 0.0;
					  equiv_tension  = 0.0;
					}
					break;
				}
				  
			}
			// Not the initial contact 
			if(inb == maInitialBallNeighborID.size())
			{
				equiv_cohesion = 0.0;
				equiv_tension  = 0.0;		
			}
			

          //  bool If_sliding = false;
			
			int failure_type= 0;
			
            if (-LocalElasticContactForce[2] > equiv_area * equiv_tension)
            {
                LocalElasticContactForce[0] = 0.0;
                LocalElasticContactForce[1]  = 0.0;
                LocalElasticContactForce[2]  = 0.0;
				
				failure_type = 1;
            }
            else
            {

                double ShearForceMax = LocalElasticContactForce[2] * equiv_tg_of_fri_ang + equiv_area * equiv_cohesion;
				
                double ShearForceNow = sqrt(LocalElasticContactForce[0] * LocalElasticContactForce[0]
                                     +      LocalElasticContactForce[1] * LocalElasticContactForce[1]);


                //Cfeng: for shear failure
                if(ShearForceMax <= 0.0)
                {
                    LocalElasticContactForce[0] = 0.0;
                    LocalElasticContactForce[1] = 0.0;
                }
                else if(ShearForceNow > ShearForceMax)
                {
                    LocalElasticContactForce[0] = ShearForceMax / ShearForceNow * LocalElasticContactForce[0];
                    LocalElasticContactForce[1] = ShearForceMax / ShearForceNow * LocalElasticContactForce[1];
					
				//	If_sliding = true;
					
					failure_type = 2;
                }
            }
			
			

            GeometryFunctions::VectorLocal2Global(LocalCoordSystem, LocalElasticContactForce, GlobalElasticContactForce);

            mOldNeighbourContactForces[i_neighbour_count][0] = GlobalElasticContactForce[0];
            mOldNeighbourContactForces[i_neighbour_count][1] = GlobalElasticContactForce[1];
            mOldNeighbourContactForces[i_neighbour_count][2] = GlobalElasticContactForce[2];
			
			

            rContactForce[0] += GlobalElasticContactForce[0];
            rContactForce[1] += GlobalElasticContactForce[1];
            rContactForce[2] += GlobalElasticContactForce[2];
			
			
			
			/////////////////////////Damping Force Calculation//////////////////////////7
		/*	if(LocalElasticContactForce[2] > 0.0) // Compression , use visc damp
			{			
				double ViscoDampingLocalContactForce[3] = {0.0};
			
				CalculateViscoDamping(LocalRelVel, ViscoDampingLocalContactForce,indentation, equiv_visco_damp_coeff_normal,
									  equiv_visco_damp_coeff_tangential, If_sliding);
				
				double GlobalDampForce[3]   = {0.0};
				GeometryFunctions::VectorLocal2Global(LocalCoordSystem, ViscoDampingLocalContactForce, GlobalDampForce);	
				
				rContactForce[0] += GlobalDampForce[0];
				rContactForce[1] += GlobalDampForce[1];
				rContactForce[2] += GlobalDampForce[2];
			}
*/


              if (mRotationOption)
			  {
					
				double MA[3] = {0.0};
				GeometryFunctions::CrossProduct(LocalCoordSystem[2], GlobalElasticContactForce, MA);
				rContactMoment[0] -= MA[0] * mRadius;
				rContactMoment[1] -= MA[1] * mRadius;
				rContactMoment[2] -= MA[2] * mRadius;
		      }

              i_neighbour_count++;
			  
			  
			  
			  if(inb < maInitialBallNeighborID.size())
			  {
				  if(maInitialBallNeighborFailureType[inb] == 0)
				  {
					  maInitialBallNeighborFailureType[inb] = failure_type;
				  }
			  }

          }// for each neighbour
		  

          rInitialRotaMoment[0] = InitialRotaMoment[0];
          rInitialRotaMoment[1] = InitialRotaMoment[1];
          rInitialRotaMoment[2] = InitialRotaMoment[2];

          KRATOS_CATCH("")
      }// ComputeBallToBallContactForce
	  
	  
	  
	  
	  
  void DEM_FEM_Particle::ComputeBallToRigidFaceContactForce(array_1d<double, 3>& rContactForce,
                                                          array_1d<double, 3>& rContactMoment,
                                                          array_1d<double, 3>& rElasticForce,
                                                          array_1d<double, 3>& rInitialRotaMoment,
                                                          ProcessInfo& rCurrentProcessInfo)
	{
		  
		  KRATOS_TRY
		  

          ConditionWeakVectorType& rNeighbours    = this->GetValue(NEIGHBOUR_RIGID_FACES);
          
		
		double mTimeStep    = rCurrentProcessInfo[DELTA_TIME];
		/////int CalRotateOption = rCurrentProcessInfo[RIGID_FACE_FLAG];

		double Friction       = mTgOfFrictionAngle;
		double young          = mYoung;
		double poisson        = mPoisson;
		double radius         = GetGeometry()(0)->FastGetSolutionStepValue(RADIUS);
		double area           = M_PI * radius * radius;
		double kn             = young * area / (2.0 * radius);
		double ks             = kn / (2.0 * (1.0 + poisson));
		
		double equiv_cohesion = mfcohesion;
	    double equiv_tension  = mftension;
		
        std::size_t iRigidFaceNeighbour = 0;

        for(ConditionWeakIteratorType ineighbour = rNeighbours.begin(); ineighbour != rNeighbours.end(); ineighbour++)
        {
			
            double LocalContactForce[3]  = {0.0};
            double GlobalContactForce[3] = {0.0};
            double GlobalContactForceOld[3] = {0.0};
			

            array_1d<double, 3 > vel = GetGeometry()(0)->FastGetSolutionStepValue(VELOCITY);

            array_1d<double, 3 > other_to_me_vel;
            noalias(other_to_me_vel) = ZeroVector(3);

            double LocalCoordSystem[3][3] = {{0.0}, {0.0}, {0.0}};
			
			double DistPToB = 0.0;


			ComputeRigidFaceToMeVelocity(ineighbour, iRigidFaceNeighbour, LocalCoordSystem, DistPToB, other_to_me_vel);			
			

            double DeltDisp[3] = {0.0};
            double DeltVel [3] = {0.0};

            DeltVel[0] = (vel[0] - other_to_me_vel[0]);
            DeltVel[1] = (vel[1] - other_to_me_vel[1]);
            DeltVel[2] = (vel[2] - other_to_me_vel[2]);

            // For translation movement delt displacement
            DeltDisp[0] = DeltVel[0] * mTimeStep;
            DeltDisp[1] = DeltVel[1] * mTimeStep;
            DeltDisp[2] = DeltVel[2] * mTimeStep;


            if (mRotationOption)
            {
                double velA[3]   = {0.0};
                double dRotaDisp[3] = {0.0};

                array_1d<double, 3 > AngularVel= GetGeometry()(0)->FastGetSolutionStepValue(ANGULAR_VELOCITY);
                double Vel_Temp[3] = { AngularVel[0], AngularVel[1], AngularVel[2]};
                GeometryFunctions::CrossProduct(Vel_Temp, LocalCoordSystem[2], velA);

                dRotaDisp[0] = -velA[0] * radius;
                dRotaDisp[1] = -velA[1] * radius;
                dRotaDisp[2] = -velA[2] * radius;

                //////contribution of the rotation vel
                DeltDisp[0] += dRotaDisp[0] * mTimeStep;
                DeltDisp[1] += dRotaDisp[1] * mTimeStep;
                DeltDisp[2] += dRotaDisp[2] * mTimeStep;
            }


            double LocalDeltDisp[3] = {0.0};
            GeometryFunctions::VectorGlobal2Local(LocalCoordSystem, DeltDisp, LocalDeltDisp);
			
			

            //////120323,for global storage
			
		    GlobalContactForceOld[0] = mFemOldNeighbourContactForces[iRigidFaceNeighbour][0];
		    GlobalContactForceOld[1] = mFemOldNeighbourContactForces[iRigidFaceNeighbour][1];
		    GlobalContactForceOld[2] = mFemOldNeighbourContactForces[iRigidFaceNeighbour][2];
			
			
            GeometryFunctions::VectorGlobal2Local(LocalCoordSystem, GlobalContactForceOld, LocalContactForce);
            LocalContactForce[0] +=  - ks * LocalDeltDisp[0];
            LocalContactForce[1] +=  - ks * LocalDeltDisp[1];
            LocalContactForce[2] +=  - kn * LocalDeltDisp[2];
			
		
		    ////Check if_initial_neighbor
			unsigned int inb;
			for(inb = 0; inb < maInitialRigidFaceNeighborID.size(); inb++)
			{
				if( maInitialRigidFaceNeighborID[inb] == static_cast<int>(ineighbour->Id()) )
				{
					if(maInitialRigidFaceNeighborFailureType[inb] > 0)
					{
					  equiv_cohesion = 0.0;
					  equiv_tension  = 0.0;
					}
					break;
				}
				  
			}
			// Not the initial contact 
			if(inb == maInitialRigidFaceNeighborID.size())
			{
				equiv_cohesion = 0.0;
				equiv_tension  = 0.0;		
			}
			
			
			
			int failure_type= 0;
			
            if (-LocalContactForce[2] > area * equiv_tension)
            {
                LocalContactForce[0] = 0.0;
                LocalContactForce[1]  = 0.0;
                LocalContactForce[2]  = 0.0;
				
				failure_type = 1;
            }
            else
            {

                double ShearForceMax = LocalContactForce[2] * Friction + area * equiv_cohesion;
				
                double ShearForceNow = sqrt(LocalContactForce[0] * LocalContactForce[0]
                                     +      LocalContactForce[1] * LocalContactForce[1]);


                //Cfeng: for shear failure
                if(ShearForceMax <= 0.0)
                {
                    LocalContactForce[0] = 0.0;
                    LocalContactForce[1] = 0.0;
                }
                else if(ShearForceNow > ShearForceMax)
                {
                    LocalContactForce[0] = ShearForceMax / ShearForceNow * LocalContactForce[0];
                    LocalContactForce[1] = ShearForceMax / ShearForceNow * LocalContactForce[1];
				
					
					failure_type = 2;
                }
            }
			
			
			
			

            GeometryFunctions::VectorLocal2Global(LocalCoordSystem, LocalContactForce, GlobalContactForce);

            mFemOldNeighbourContactForces[iRigidFaceNeighbour][0] = GlobalContactForce[0];
            mFemOldNeighbourContactForces[iRigidFaceNeighbour][1] = GlobalContactForce[1];
            mFemOldNeighbourContactForces[iRigidFaceNeighbour][2] = GlobalContactForce[2];
			
            ///Global stored contact force between rigid face and particle, used by fem elements
            
            Vector& neighbour_rigid_faces_elastic_contact_force = this->GetValue(NEIGHBOUR_RIGID_FACES_ELASTIC_CONTACT_FORCE);
            neighbour_rigid_faces_elastic_contact_force[3 * iRigidFaceNeighbour + 0] = GlobalContactForce[0];
            neighbour_rigid_faces_elastic_contact_force[3 * iRigidFaceNeighbour + 1] = GlobalContactForce[1];
            neighbour_rigid_faces_elastic_contact_force[3 * iRigidFaceNeighbour + 2] = GlobalContactForce[2];
			

            rContactForce[0] += GlobalContactForce[0];
            rContactForce[1] += GlobalContactForce[1];
            rContactForce[2] += GlobalContactForce[2];
		
			
			/////////////////////////////////////////////////////////////////////////////

            if ( mRotationOption)
            {
                double MA[3] = {0.0};
                GeometryFunctions::CrossProduct(LocalCoordSystem[2], GlobalContactForce, MA);
                rContactMoment[0] -= MA[0] * radius;
                rContactMoment[1] -= MA[1] * radius;
                rContactMoment[2] -= MA[2] * radius;
            }


            iRigidFaceNeighbour++;
			
			
			  if(inb < maInitialRigidFaceNeighborID.size())
			  {
				  if(maInitialRigidFaceNeighborFailureType[inb] == 0)
				  {
					  maInitialRigidFaceNeighborFailureType[inb] = failure_type;
				  }
			  }

        }
		  
          KRATOS_CATCH("")
		  
      }// ComputeBallToRigidFaceContactForce
	  
	  
	  
	  
	  



void DEM_FEM_Particle::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
      {
          KRATOS_TRY

          SphericParticle::CalculateRightHandSide(rRightHandSideVector, rCurrentProcessInfo);
		  
		  
		  
		  double local_damp_ratio   = this->GetGeometry()(0)->FastGetSolutionStepValue(LOCAL_DAMP_RATIO);
		  array_1d<double,3> vel      = this->GetGeometry()(0)->FastGetSolutionStepValue(VELOCITY);
		  array_1d<double,3> rota_vel = this->GetGeometry()(0)->FastGetSolutionStepValue(ANGULAR_VELOCITY);
		  
		  
		 // KRATOS_WATCH(rRightHandSideVector[1])
		  
		  for(unsigned int i = 0; i < 3; i++)
		  {
			  if(vel[i] > 0.0)
			  {
				  rRightHandSideVector[i] -= local_damp_ratio * fabs(rRightHandSideVector[i]);
			  }
			  else
			  {
				  rRightHandSideVector[i] += local_damp_ratio * fabs(rRightHandSideVector[i]);				  
			  }
			  
			  if (mRotationOption)
			  {
				  if(rota_vel[i] > 0.0)
				  {
					  rRightHandSideVector[i + 3] -= local_damp_ratio * fabs(rRightHandSideVector[i + 3]);
				  }
				  else
				  {
					  rRightHandSideVector[i + 3] += local_damp_ratio * fabs(rRightHandSideVector[i + 3]);				  
				  }
			  }
			  
		  }
		 // KRATOS_WATCH(rRightHandSideVector[1])
		 // KRATOS_WATCH(local_damp_ratio)

          KRATOS_CATCH( "" )
      }


	
	  void DEM_FEM_Particle::Calculate(const Variable<double>& rVariable, double& Output, const ProcessInfo& rCurrentProcessInfo)
	  {
		  SphericParticle::Calculate(rVariable, Output, rCurrentProcessInfo);
	  }
      void DEM_FEM_Particle::Calculate(const Variable<array_1d<double, 3 > >& rVariable, array_1d<double, 3 > & Output, const ProcessInfo& rCurrentProcessInfo)
	  {
		  SphericParticle::Calculate(rVariable, Output, rCurrentProcessInfo);		  
	  }
      void DEM_FEM_Particle::Calculate(const Variable<Vector >& rVariable, Vector& Output, const ProcessInfo& rCurrentProcessInfo)
	  {
		  SphericParticle::Calculate(rVariable, Output, rCurrentProcessInfo);		  
	  }
      void DEM_FEM_Particle::Calculate(const Variable<Matrix >& rVariable, Matrix& Output, const ProcessInfo& rCurrentProcessInfo)
	  {
		  SphericParticle::Calculate(rVariable, Output, rCurrentProcessInfo);		  
	  }

    //************************************************************************************
    //************************************************************************************
  
    void DEM_FEM_Particle::save( Serializer& rSerializer ) const
    {
//  std::cout << "Saving the Particle #" << Id() << std::endl;
        rSerializer.save( "Name", "Particle" );
        KRATOS_SERIALIZE_SAVE_BASE_CLASS( rSerializer, Element );
    }

    void DEM_FEM_Particle::load( Serializer& rSerializer )
    {
        KRATOS_SERIALIZE_LOAD_BASE_CLASS( rSerializer, Element );
    }



} // Namespace Kratos


