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
      : SphericContinuumParticle(){/*mInitializedVariablesFlag = 0;*/}

      CylinderContinuumParticle::CylinderContinuumParticle(IndexType NewId, GeometryType::Pointer pGeometry)
      : SphericContinuumParticle(NewId, pGeometry){/*mInitializedVariablesFlag = 0;*/}

      CylinderContinuumParticle::CylinderContinuumParticle(IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
      : SphericContinuumParticle(NewId, pGeometry, pProperties){/*mInitializedVariablesFlag = 0;*/}

      CylinderContinuumParticle::CylinderContinuumParticle(IndexType NewId, NodesArrayType const& ThisNodes)
      : SphericContinuumParticle(NewId, ThisNodes){/*mInitializedVariablesFlag = 0;*/}

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
          mRadius                   = GetGeometry()[0].FastGetSolutionStepValue(RADIUS);
          double density            = GetDensity();

          double& mass              = GetGeometry()[0].FastGetSolutionStepValue(NODAL_MASS);
          mass                      = KRATOS_M_PI_3 * density * mRadius * mRadius * 1.0;
          mRealMass                 = mass;

          if (this->Is(DEMFlags::HAS_ROTATION) ){
            double moment_of_inertia = 0.5 * mass * mRadius * mRadius;   
            GetGeometry()[0].FastGetSolutionStepValue(PARTICLE_MOMENT_OF_INERTIA) = moment_of_inertia;
          }                                                                        

          CustomInitialize();

          KRATOS_CATCH( "" )
      }

      //**************************************************************************************************************************************************
      //**************************************************************************************************************************************************

           
      
       void CylinderContinuumParticle::ContactAreaWeighting2D() //MISMI 10: POOYAN this could be done by calculating on the bars. not looking at the neighbous of my neighbours.
      { 

        double alpha = 1.0;
        double sphere_perimeter = 2*KRATOS_M_PI*mRadius;  
        
        double total_equiv_perimeter = 0.0;

        ParticleWeakVectorType r_continuum_ini_neighbours    = this->GetValue(CONTINUUM_INI_NEIGHBOUR_ELEMENTS);
        int cont_ini_neighbours_size                         = r_continuum_ini_neighbours.size();
        
        mcont_ini_neigh_area.resize(cont_ini_neighbours_size);  //NOTE: Here we use "mcont_ini_neigh_area" just because in the general 3D particle this is the name used.ContactAreaWeighting2D
        
        //computing the total equivalent area
        
        size_t index = 0;
        
        for(ParticleWeakIteratorType ini_cont_neighbour_iterator = r_continuum_ini_neighbours.begin();     // MSIMSI 99:Could this loop be done during the bar creation in the strategy and so avoid another repetition?
            ini_cont_neighbour_iterator != r_continuum_ini_neighbours.end(); ini_cont_neighbour_iterator++)
        {   
            double other_radius     = ini_cont_neighbour_iterator->GetGeometry()[0].FastGetSolutionStepValue(RADIUS);
            double equiv_radius     = 2*mRadius * other_radius / (mRadius + other_radius);        
            //double equiv_area       = (0.25)*KRATOS_M_PI * equiv_radius * equiv_radius; //we now take 1/2 of the efective mRadius.
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

                          alpha            = 1.30*(1.10266)*(sphere_perimeter/total_equiv_perimeter)*((double(cont_ini_neighbours_size))/6); // 6 is mean coordination number.
                          mcont_ini_neigh_area[skin_index] = alpha*mcont_ini_neigh_area[skin_index];
     
                    
                  skin_index++;
                  
                  }//loop on cont neighs       
                  
            }//skin particles.
            
        }//if 3 neighbours or more.
   
      } //Contact Area Weighting
      
      
      void CylinderContinuumParticle::AddNeighbourContributionToStressTensor(double GlobalElasticContactForce[3],
                                                                            array_1d<double,3> &other_to_me_vect,
                                                                            const double &distance,
                                                                            const double &radius_sum,
                                                                            const double &calculation_area,
                                                                            ParticleWeakIteratorType neighbour_iterator, 
                                                                            ProcessInfo& rCurrentProcessInfo, 
                                                                            double &rRepresentative_Volume)
      {

        double gap                  = distance - radius_sum;
      
        array_1d<double,3> normal_vector_on_contact =  -1 * other_to_me_vect; //outwards
      
        double Dummy_Dummy = 0.0;
        GeometryFunctions::normalize(normal_vector_on_contact,Dummy_Dummy); // Normalize to unitary module

        array_1d<double,3> x_centroid      = (mRadius + 0.5*gap) * normal_vector_on_contact;
      
        array_1d<double,3> surface_baricenter = x_centroid;
        
        double result_product = GeometryFunctions::DotProduct(surface_baricenter,normal_vector_on_contact); 
        
        //Aproximation with error: surface_baricenter should be the baricenter of each surface, which can no be calculated because the surfaces are imaginary.
      
        rRepresentative_Volume = rRepresentative_Volume + 0.5 * (result_product * calculation_area);
        
        for (int i=0; i<3; i++)
        {
            for (int j=0; j<3; j++)
            {   
                mStressTensor[i][j] += (x_centroid[j]) * GlobalElasticContactForce[i]; //ref: Katalin Bagi 1995 Mean stress tensor           
            
              
            }
        }

          
      } //AddNeighbourContributionToStressTensor

      //**************************************************************************************************************************************************
      //**************************************************************************************************************************************************
//       void CylinderContinuumParticle::Calculate(const Variable<double>& rVariable, double& Output, const ProcessInfo& rCurrentProcessInfo){}
//       void CylinderContinuumParticle::Calculate(const Variable<array_1d<double, 3> >& rVariable, array_1d<double, 3>& Output, const ProcessInfo& rCurrentProcessInfo){}
//       void CylinderContinuumParticle::Calculate(const Variable<Vector >& rVariable, Vector& Output, const ProcessInfo& rCurrentProcessInfo){}
//       void CylinderContinuumParticle::Calculate(const Variable<Matrix >& rVariable, Matrix& Output, const ProcessInfo& rCurrentProcessInfo){}

      //**************************************************************************************************************************************************
      //**************************************************************************************************************************************************

}  // namespace Kratos.

