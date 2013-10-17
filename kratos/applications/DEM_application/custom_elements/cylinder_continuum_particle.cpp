//
//   Project Name:        Kratos
//   Last Modified by:    $Author: M.Santasusana $
//   Date:                $Date: 2006-11-27 16:07:33 $
//   Revision:            $Revision: 1.1.1.1 $
//
//


// System includes
#include <string>
#include <iostream>

// External includes


// Project includes
#include "includes/define.h"
#include "cylinder_continuum_particle.h"
#include "custom_utilities/GeometryFunctions.h"
#include "custom_utilities/AuxiliaryFunctions.h"
#include "DEM_application.h"



namespace Kratos
{
     // using namespace GeometryFunctions;

      CylinderContinuumParticle::CylinderContinuumParticle()
      : SphericContinuumParticle(){mInitializedVariablesFlag = 0;}

      CylinderContinuumParticle::CylinderContinuumParticle(IndexType NewId, GeometryType::Pointer pGeometry)
      : SphericContinuumParticle(NewId, pGeometry){mInitializedVariablesFlag = 0;}

      CylinderContinuumParticle::CylinderContinuumParticle(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
      : SphericContinuumParticle(NewId, pGeometry, pProperties){mInitializedVariablesFlag = 0;}

      CylinderContinuumParticle::CylinderContinuumParticle(IndexType NewId, NodesArrayType const& ThisNodes)
      : SphericContinuumParticle(NewId, ThisNodes){mInitializedVariablesFlag = 0;}

      Element::Pointer CylinderContinuumParticle::Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const
      {
           return SphericContinuumParticle::Pointer(new CylinderContinuumParticle(NewId, GetGeometry().Create(ThisNodes), pProperties));

      }

      /// Destructor.
      CylinderContinuumParticle::~CylinderContinuumParticle(){}

      void CylinderContinuumParticle::Initialize()
      {
          KRATOS_TRY

          mDimension                = 2;
          mRadius                   = GetGeometry()(0)->FastGetSolutionStepValue(RADIUS);
          mYoung                    = GetGeometry()(0)->FastGetSolutionStepValue(YOUNG_MODULUS);         
          mPoisson                  = GetGeometry()(0)->FastGetSolutionStepValue(POISSON_RATIO);
          mTgOfFrictionAngle        = GetGeometry()(0)->FastGetSolutionStepValue(PARTICLE_FRICTION);
          mLnOfRestitCoeff          = GetGeometry()(0)->FastGetSolutionStepValue(LN_OF_RESTITUTION_COEFF);
          double& density           = GetGeometry()(0)->FastGetSolutionStepValue(PARTICLE_DENSITY);
          double& mass              = GetGeometry()(0)->FastGetSolutionStepValue(NODAL_MASS);
          double& sqrt_of_mass      = GetGeometry()(0)->FastGetSolutionStepValue(SQRT_OF_MASS);
          double& moment_of_inertia = GetGeometry()(0)->FastGetSolutionStepValue(PARTICLE_MOMENT_OF_INERTIA);

          mass                      = M_PI * density * mRadius * mRadius * 1.0;
          sqrt_of_mass              = sqrt(mass);
          moment_of_inertia         = 0.5 * mass * mRadius * mRadius;
          mRealMass                 = mass;          
          mSqrtOfRealMass           = sqrt_of_mass;
          mMomentOfInertia          = moment_of_inertia;

          CustomInitialize();

          KRATOS_CATCH( "" )
      }

      //**************************************************************************************************************************************************
      //**************************************************************************************************************************************************

      void CylinderContinuumParticle::CalculateRightHandSide(VectorType& rRightHandSideVector, ProcessInfo& rCurrentProcessInfo)
      {
          KRATOS_TRY

          array_1d<double, 3> contact_force;
          array_1d<double, 3> contact_moment;
          array_1d<double, 3> additionally_applied_force;
          array_1d<double, 3> additionally_applied_moment;
          array_1d<double, 3> initial_rotation_moment;     
          array_1d<double, 3>& elastic_force = this->GetGeometry()[0].GetSolutionStepValue(ELASTIC_FORCES);

          contact_force.clear();
          contact_moment.clear();
          additionally_applied_force.clear();
          additionally_applied_moment.clear();
          initial_rotation_moment.clear();
          elastic_force.clear();

          ComputeNewNeighboursHistoricalData();

          ComputeBallToBallContactForce(contact_force, contact_moment, elastic_force, initial_rotation_moment, rCurrentProcessInfo);

          if (mLimitSurfaceOption > 0){

              for (int surface_num = 0; surface_num < mLimitSurfaceOption; surface_num++){
                  ComputeBallToSurfaceContactForce(contact_force, contact_moment, initial_rotation_moment, surface_num, rCurrentProcessInfo);
              }
          }
          
          if (mLimitCylinderOption > 0){

              for (int cylinder_num = 0; cylinder_num < mLimitCylinderOption; cylinder_num++){
                  ComputeBallToCylinderContactForce(contact_force, contact_moment, initial_rotation_moment, cylinder_num, rCurrentProcessInfo);
              }
          }
          
          ComputeAdditionalForces(contact_force, contact_moment, additionally_applied_force, additionally_applied_moment, rCurrentProcessInfo);

          rRightHandSideVector[0] = contact_force[0]  + additionally_applied_force[0];
          rRightHandSideVector[1] = contact_force[1]  + additionally_applied_force[1];
          rRightHandSideVector[2] = contact_force[2]  + additionally_applied_force[2];
          //rRightHandSideVector[2] = contact_force[2]  + additionally_applied_force[2];
          rRightHandSideVector[3] = contact_moment[0] + additionally_applied_moment[0];
          rRightHandSideVector[4] = contact_moment[1] + additionally_applied_moment[0];
          rRightHandSideVector[5] = contact_moment[2] + additionally_applied_moment[0];
          //rRightHandSideVector[5] = contact_moment[2] + additionally_applied_moment[0];

          KRATOS_CATCH( "" )
      }
      
      
       void SphericContinuumParticle::ContactAreaWeighting2D(const ProcessInfo& rCurrentProcessInfo) //MISMI 10: POOYAN this could be done by calculating on the bars. not looking at the neighbous of my neighbours.
      { 

        double alpha = 1.0;
        double sphere_perimeter = 2*M_PI*mRadius;  
        
        double total_equiv_perimeter = 0.0;

        ParticleWeakVectorType r_continuum_ini_neighbours    = this->GetValue(CONTINUUM_INI_NEIGHBOUR_ELEMENTS);
        int cont_ini_neighbours_size                         = r_continuum_ini_neighbours.size();
        
        mcont_ini_neigh_area.resize(cont_ini_neighbours_size);  //NOTE: Here we use "mcont_ini_neigh_area" just becouse in the general 3D particle this is the name used.
        
        //computing the total equivalent area
        
        size_t index = 0;
        
        for(ParticleWeakIteratorType ini_cont_neighbour_iterator = r_continuum_ini_neighbours.begin();     // MSIMSI 99:Could this loop be done during the bar creation in the strategy and so avoid another repetition?
            ini_cont_neighbour_iterator != r_continuum_ini_neighbours.end(); ini_cont_neighbour_iterator++)
        {   
            double other_radius     = ini_cont_neighbour_iterator->GetGeometry()(0)->GetSolutionStepValue(RADIUS);
            double equiv_radius     = 2*mRadius * other_radius / (mRadius + other_radius);        
            //double equiv_area       = (0.25)*M_PI * equiv_radius * equiv_radius; //we now take 1/2 of the efective mRadius.
            total_equiv_perimeter  += equiv_radius;
        
            mcont_ini_neigh_area[index] = equiv_radius; //*  //this is consistent since in 2D, we work with cylinders of depth unit 1.0.
            index++; //*
            
            
        } //for every neighbour
      
        if (cont_ini_neighbours_size >= 3)
        {
            if(!*mSkinSphere)
            {
            
              
              AuxiliaryFunctions::CalculateAlphaFactor2D(cont_ini_neighbours_size, sphere_perimeter, total_equiv_perimeter, alpha); 
              
              size_t not_skin_index = 0;
          
              for(ParticleWeakIteratorType ini_cont_neighbour_iterator = r_continuum_ini_neighbours.begin();
                  ini_cont_neighbour_iterator != r_continuum_ini_neighbours.end(); ini_cont_neighbour_iterator++)
                  
                  {      
                      mcont_ini_neigh_area[not_skin_index] = alpha*mcont_ini_neigh_area[not_skin_index];
                      not_skin_index++;  
                      
                  } //for every neighbour

            }
            
            else //skin sphere 
            {
 
                size_t skin_index = 0; 
                
                for(ParticleWeakIteratorType ini_cont_neighbour_iterator = r_continuum_ini_neighbours.begin();
                  ini_cont_neighbour_iterator != r_continuum_ini_neighbours.end(); ini_cont_neighbour_iterator++)
                  
                  {

                          alpha            = 2.0*(1.10266)*(sphere_perimeter/total_equiv_perimeter)*((double(cont_ini_neighbours_size))/6); // 6 is mean coordination number.
                          mcont_ini_neigh_area[skin_index] = alpha*mcont_ini_neigh_area[skin_index];
     
                    
                  skin_index++;
                  
                  }//loop on cont neighs       
                  
            }//skin particles.
            
        }//if 3 neighbours or more.
   
      } //Contact Area Weighting

      //**************************************************************************************************************************************************
      //**************************************************************************************************************************************************
      void CylinderContinuumParticle::Calculate(const Variable<double>& rVariable, double& Output, const ProcessInfo& rCurrentProcessInfo){}
      void CylinderContinuumParticle::Calculate(const Variable<array_1d<double, 3> >& rVariable, array_1d<double, 3>& Output, const ProcessInfo& rCurrentProcessInfo){}
      void CylinderContinuumParticle::Calculate(const Variable<Vector >& rVariable, Vector& Output, const ProcessInfo& rCurrentProcessInfo){}
      void CylinderContinuumParticle::Calculate(const Variable<Matrix >& rVariable, Matrix& Output, const ProcessInfo& rCurrentProcessInfo){}

      //**************************************************************************************************************************************************
      //**************************************************************************************************************************************************

}  // namespace Kratos.

