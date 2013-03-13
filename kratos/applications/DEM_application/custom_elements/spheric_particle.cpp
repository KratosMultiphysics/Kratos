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
#include <iomanip>

// External includes


// Project includes
#include "includes/define.h"
#include "spheric_particle.h"
#include "custom_utilities/GeometryFunctions.h"
#include "custom_utilities/AuxiliaryFunctions.h"
#include "DEM_application.h"


namespace Kratos
{
     // using namespace GeometryFunctions;

      SphericParticle::SphericParticle() : DiscreteElement(){}

      SphericParticle::SphericParticle( IndexType NewId, GeometryType::Pointer pGeometry) : DiscreteElement(NewId, pGeometry){}

      SphericParticle::SphericParticle( IndexType NewId, GeometryType::Pointer pGeometry,  PropertiesType::Pointer pProperties)
      : DiscreteElement(NewId, pGeometry, pProperties){}

      SphericParticle::SphericParticle(IndexType NewId, NodesArrayType const& ThisNodes)
      : DiscreteElement(NewId, ThisNodes){}

      Element::Pointer SphericParticle::Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const
      {
         return DiscreteElement::Pointer(new SphericParticle(NewId, GetGeometry().Create( ThisNodes ), pProperties) );
      }

      /// Destructor.
      SphericParticle::~SphericParticle(){}


      void SphericParticle::Initialize(){

        KRATOS_TRY

        mpFailureId = &(this->GetValue(PARTICLE_FAILURE_ID));

        //mDimension = this->GetGeometry().WorkingSpaceDimension();
        mDimension = 3; ///WARNING: I: to be revised.
        
        double density      = GetGeometry()(0)->FastGetSolutionStepValue(PARTICLE_DENSITY);
        double radius       = GetGeometry()(0)->FastGetSolutionStepValue(RADIUS);
        double& mass        = GetGeometry()(0)->FastGetSolutionStepValue(NODAL_MASS);

        double & Inertia         = GetGeometry()(0)->FastGetSolutionStepValue(PARTICLE_INERTIA);
        double & MomentOfInertia = GetGeometry()(0)->FastGetSolutionStepValue(PARTICLE_MOMENT_OF_INERTIA);
        
        double& Representative_Volume = this->GetGeometry()[0].GetSolutionStepValue(REPRESENTATIVE_VOLUME);
        Representative_Volume = 0.0;
        
        mContinuumGroup     = this->GetGeometry()[0].GetSolutionStepValue(PARTICLE_CONTINUUM);

        //TO BE IMPROVED: i would like to work with *mpFailureId as integer. the problem is that it has to be exported to GID to be plotted.

        //provisional way is:
        if(mContinuumGroup==0)  {*mpFailureId=1;}
        else                    {*mpFailureId=0;}

        //easiest way is:   
        //*mpFailureId          = !(mContinuumGroup);
           
        if(mDimension ==2)
        {
            mass     = M_PI * radius * radius * density;

            mRealMass = mass;

            Inertia = 0.25 * M_PI * radius * radius * radius  * radius ;

            MomentOfInertia = 0.5 * radius * radius * mass;

        }
        else
        {
            mass     = 4.0 / 3.0 * M_PI * radius * radius * radius * density;

            mRealMass = mass;

            Inertia = 0.25 * M_PI * radius * radius * radius  * radius ;

            MomentOfInertia = 0.4 * radius * radius * mass;

        }
        

        KRATOS_CATCH( "" )

      }

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************


      void SphericParticle::CalculateRightHandSide(VectorType& rRightHandSideVector,ProcessInfo& rCurrentProcessInfo){

        ComputeParticleContactForce(rCurrentProcessInfo);

        if( (rCurrentProcessInfo[ROTATION_OPTION] != 0) && (rCurrentProcessInfo[ROTATION_SPRING_OPTION] != 0) )
        {
            ComputeParticleRotationSpring(rCurrentProcessInfo);
        }

        //CharacteristicParticleFailureId(rCurrentProcessInfo);

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

      }
      void SphericParticle::EquationIdVector(EquationIdVectorType& rResult, ProcessInfo& rCurrentProcessInfo){}

      void SphericParticle::MassMatrix(MatrixType& rMassMatrix, ProcessInfo& rCurrentProcessInfo)

      {

          double radius = GetGeometry()(0)->GetSolutionStepValue(RADIUS);
          double volume =   1.333333333333333*M_PI*radius*radius*radius;
          double density = GetGeometry()(0)->GetSolutionStepValue(PARTICLE_DENSITY);
          rMassMatrix.resize(1,1);
          rMassMatrix(0,0) = volume*density;

      }

        void SphericParticle::SetInitialContacts(int case_opt, ProcessInfo& rCurrentProcessInfo  ) //vull ficar que sigui zero si no son veins cohesius.
        {   
            /*
             * 
             * ELEMENT_NEIGHBOURS / NEIGHBOURS_IDS: the ones important for calculating forces.
             * INI_NEIGHBOURS_IDS: the ones to be treated specially due to initial delta or continuum case.
             * INI_CONTINUUM_NEIGHBOURS_IDS: only the ones that are continuum at 0 step and we should treat the possible detachment.
             * 
             * These 3 classes do NOT coincide at t=0!
 
            */
                                
            bool delta_OPTION;
            bool continuum_simulation_OPTION;
            int contact_mesh_OPTION = rCurrentProcessInfo[CONTACT_MESH_OPTION];
            
            AuxiliaryFunctions::SwitchCase(case_opt, delta_OPTION, continuum_simulation_OPTION);
            
            // DEFINING THE REFERENCES TO THE MAIN PARAMETERS

            ParticleWeakVectorType& r_neighbours                = this->GetValue(NEIGHBOUR_ELEMENTS);
            
            ParticleWeakVectorType& r_continuum_ini_neighbours    = this->GetValue(CONTINUUM_INI_NEIGHBOUR_ELEMENTS);
            
            r_continuum_ini_neighbours.clear();
            
            size_t ini_size = 0;
            size_t continuum_ini_size =0;

            unsigned int i=0;

            //default
            *mpFailureId=1;

            //SAVING THE INICIAL NEIGHBOURS, THE DELTAS AND THE FAILURE ID

            
            for(ParticleWeakIteratorType_ptr ineighbour = r_neighbours.ptr_begin();  //loop over the neighbours and store into a initial_neighbours vector.
            ineighbour != r_neighbours.ptr_end(); ineighbour++)
            {
           
                array_1d<double,3> other_to_me_vect = this->GetGeometry()(0)->Coordinates() - ((*ineighbour).lock())->GetGeometry()(0)->Coordinates();
                double distance                     = sqrt(other_to_me_vect[0] * other_to_me_vect[0] +
                                                     other_to_me_vect[1] * other_to_me_vect[1] +
                                                     other_to_me_vect[2] * other_to_me_vect[2]);

                double radius_sum                   = this->GetGeometry()(0)->GetSolutionStepValue(RADIUS) + ((*ineighbour).lock())->GetGeometry()(0)->GetSolutionStepValue(RADIUS);
                double initial_delta                = radius_sum - distance;

                int r_other_continuum_group = ((*ineighbour).lock())->GetGeometry()(0)->GetSolutionStepValue(PARTICLE_CONTINUUM);

                if( ( (r_other_continuum_group == mContinuumGroup) && (mContinuumGroup != 0) ) || ( fabs(initial_delta)>1.0e-6 ) ) // which would be the appropiate tolerance?
                //THESE ARE THE CASES THAT NEED TO STORE THE INITIAL NEIGHBOURS
                {
                    //Number of contacts accounting.
              
                    int& total_number_of_contacts = rCurrentProcessInfo[TOTAL_CONTACTS];

                    total_number_of_contacts++;

                    //initial neighbours

                    ini_size++;

                    (this->GetValue(INI_NEIGHBOURS_IDS)).resize(ini_size);
                    (this->GetValue(PARTICLE_INITIAL_DELTA)).resize(ini_size);
                    (this->GetValue(PARTICLE_INITIAL_FAILURE_ID)).resize(ini_size);
                                   
                 
                    this->GetValue(PARTICLE_INITIAL_DELTA)[ini_size - 1]      =  0.0;
                    this->GetValue(PARTICLE_INITIAL_FAILURE_ID)[ini_size - 1] = 1;
                    this->GetValue(INI_NEIGHBOURS_IDS)[ini_size - 1] = ((*ineighbour).lock())->Id();
                   
                  
                    if (delta_OPTION == true)
                    {
            
                        this->GetValue(PARTICLE_INITIAL_DELTA)[ini_size - 1]    =   initial_delta;
       
                        this->GetValue(PARTICLE_CONTACT_DELTA)[i]               =   initial_delta; //these variables are different and need to be kept.

                    }

                    if (continuum_simulation_OPTION == true)
                    {
                        
                        if ( (r_other_continuum_group == mContinuumGroup) && (mContinuumGroup != 0) )
                        {
                            this->GetValue(PARTICLE_CONTACT_FAILURE_ID)[i]=0;
                            this->GetValue(PARTICLE_INITIAL_FAILURE_ID)[ini_size - 1]=0;
                            *mpFailureId=0; // if a cohesive contact exist, the FailureId becomes 0.          
                            continuum_ini_size++;
                                
                            r_continuum_ini_neighbours.push_back(*ineighbour);
                            (this->GetValue(CONTINUUM_INI_NEIGHBOURS_IDS)).resize(continuum_ini_size);
                            (this->GetValue(CONTINUUM_PARTICLE_INITIAL_FAILURE_ID)).resize(continuum_ini_size);
                            
                              this->GetValue(CONTINUUM_INI_NEIGHBOURS_IDS)[continuum_ini_size - 1] = ((*ineighbour).lock())->Id();
                             this->GetValue(CONTINUUM_PARTICLE_INITIAL_FAILURE_ID)[continuum_ini_size - 1] = 0;
                             
                            if(contact_mesh_OPTION)
                            {

                                (this->GetGeometry()(0))->GetValue(NODE_TO_NEIGH_ELEMENT_POINTER).resize(continuum_ini_size);
                   
                            }                          
                        }

                    }//for continuum_simulation_OPTION      //hi havia: r_VectorContactFailureId[i]=1; //generally detached    //diferent group

                } // FOR THE CASES THAT NEED STORING INITIAL NEIGHBOURS

                i++;
              
            } //end for: ParticleWeakIteratorType ineighbour
            
            if (continuum_simulation_OPTION == true)
            {
             
              ContactAreaWeighting(rCurrentProcessInfo);
                  
            }
            
       
        }//SetInitialContacts

       void SphericParticle::ContactAreaWeighting(const ProcessInfo& rCurrentProcessInfo ) //only for the continuum_case

       { 
     
          int skin_sphere     = this->GetValue(SKIN_SPHERE);

          double alpha = 1.0;
          double radius = this->GetGeometry()(0)->GetSolutionStepValue(RADIUS);
          double external_sphere_area = 4*M_PI*radius*radius;  
          
          mtotal_equiv_area = 0.0;
              
          int cont_ini_neighbours_size = this->GetValue(CONTINUUM_INI_NEIGHBOUR_ELEMENTS).size();
              
          ParticleWeakVectorType r_continuum_ini_neighbours    = this->GetValue(CONTINUUM_INI_NEIGHBOUR_ELEMENTS);

          mcont_ini_neigh_area.resize(cont_ini_neighbours_size);
          
          //computing the total equivalent area
          
          for(ParticleWeakIteratorType ini_cont_neighbour_iterator = r_continuum_ini_neighbours.begin();
              ini_cont_neighbour_iterator != r_continuum_ini_neighbours.end(); ini_cont_neighbour_iterator++)
          {   
              double other_radius     = ini_cont_neighbour_iterator->GetGeometry()(0)->GetSolutionStepValue(RADIUS);
              double equiv_radius    = 2*radius * other_radius / (radius + other_radius);        
              double equiv_area      = (0.25)*M_PI * equiv_radius * equiv_radius; //we now take 1/2 of the efective radius.
              mtotal_equiv_area   += equiv_area;
          
          } //for every neighbour
       
       
          if(skin_sphere != 1)
          {
          
            AuxiliaryFunctions::CalculateAlphaFactor(cont_ini_neighbours_size, external_sphere_area, mtotal_equiv_area, alpha); 
          
          }
          else //skin sphere //AIXO HAURIA DE CANVIAR PER L'AREA COPIADA BONA DEL VEI KE NO ES 
          {
              //alpha            = 1.40727*4*M_PI*radius*radius*n_neighbours/(11*total_equiv_area);
              //alpha            = (1.40727)*(external_sphere_area/mtotal_equiv_area)*((double(cont_ini_neighbours_size))/11);
            
              alpha            = 1.0*(1.40727)*(external_sphere_area/mtotal_equiv_area)*((double(cont_ini_neighbours_size))/11);
                
          }
          
              
          size_t index = 0;
        
          for(ParticleWeakIteratorType ini_cont_neighbour_iterator = r_continuum_ini_neighbours.begin();
              ini_cont_neighbour_iterator != r_continuum_ini_neighbours.end(); ini_cont_neighbour_iterator++)
          {   
              double other_radius     = ini_cont_neighbour_iterator->GetGeometry()(0)->GetSolutionStepValue(RADIUS);
              double equiv_radius    = 2*radius * other_radius / (radius + other_radius);        
              double equiv_area      = (0.25)*M_PI * equiv_radius * equiv_radius;
              double corrected_area  = alpha*equiv_area;
              
              mcont_ini_neigh_area[index] = corrected_area;
              
              index++;    
          } //for every neighbour
          
          
          
          
          
      } //Contact Area Weighting
        
        
      void SphericParticle::ComputeParticleContactForce(ProcessInfo& rCurrentProcessInfo )
      {
          KRATOS_TRY

          ParticleWeakVectorType& r_neighbours                = this->GetValue(NEIGHBOUR_ELEMENTS);
          //ParticleWeakVectorType& r_continuum_ini_neighbours  = this->GetValue(CONTINUUM_INI_NEIGHBOUR_ELEMENTS);

          vector<double>& r_VectorContactInitialDelta         = this->GetValue(PARTICLE_CONTACT_DELTA);

          
          
          //optimitzation:
          
          VectorArray3Double&  GlobalContactForceMatrix = this->GetValue(PARTICLE_CONTACT_FORCES);
   
          
          // PROCESS INFO

          double magic_factor                 = rCurrentProcessInfo[DEM_MAGIC_FACTOR];
          
          const array_1d<double,3>& gravity   = rCurrentProcessInfo[GRAVITY];

          double dt                           = rCurrentProcessInfo[DELTA_TIME];
          int damp_id                         = rCurrentProcessInfo[DAMP_TYPE];
          int force_calculation_type_id       = rCurrentProcessInfo[FORCE_CALCULATION_TYPE];
          int rotation_OPTION                 = rCurrentProcessInfo[ROTATION_OPTION]; //M:  it's 1/0, should be a boolean
          int global_variables_OPTION         = rCurrentProcessInfo[GLOBAL_VARIABLES_OPTION]; //M:  it's 1/0, should be a boolean

          int case_OPTION                     = rCurrentProcessInfo[CASE_OPTION];
          bool delta_OPTION;
          bool continuum_simulation_OPTION;

          AuxiliaryFunctions::SwitchCase(case_OPTION, delta_OPTION, continuum_simulation_OPTION);

          
          // GETTING PARTICLE PROPERTIES
             
          double radius               = this->GetGeometry()[0].GetSolutionStepValue(RADIUS);
          double mass                 = mRealMass;

          double young                = this->GetGeometry()[0].GetSolutionStepValue(YOUNG_MODULUS);
          double poisson              = this->GetGeometry()[0].GetSolutionStepValue(POISSON_RATIO);
          double FriAngle             = this->GetGeometry()[0].GetSolutionStepValue(PARTICLE_FRICTION);
          
          double restitution_coeff    = this->GetGeometry()[0].GetSolutionStepValue(RESTITUTION_COEFF);
      
          array_1d<double,3>& rhs             = this->GetGeometry()[0].GetSolutionStepValue(RHS);
          array_1d<double,3>& total_forces    = this->GetGeometry()[0].GetSolutionStepValue(TOTAL_FORCES);
          array_1d<double,3>& damp_forces     = this->GetGeometry()[0].GetSolutionStepValue(DAMP_FORCES);

          array_1d<double, 3 > vel            = this->GetGeometry()(0)->GetSolutionStepValue(VELOCITY);
          array_1d<double, 3 > delta_displ            = this->GetGeometry()(0)->GetSolutionStepValue(DELTA_DISPLACEMENT);
          
          //Aplied Force for pressure:
          
          array_1d<double,3> external_total_applied_force;
          
          double& Representative_Volume = this->GetGeometry()[0].GetSolutionStepValue(REPRESENTATIVE_VOLUME);
          
          
          external_total_applied_force[1] = 0.0;
          external_total_applied_force[1] = 0.0;
          external_total_applied_force[1] = 0.0;
          
          array_1d<double,3>& applied_force = this->GetGeometry()[0].GetSolutionStepValue(APPLIED_FORCE);
          
          if (rCurrentProcessInfo[INT_DUMMY_2]==1) //activated external force
          {
            
            external_total_applied_force = this->GetGeometry()[0].GetSolutionStepValue(EXTERNAL_APPLIED_FORCE);
          
            double initial_time     = rCurrentProcessInfo[INITIAL_PRESSURE_TIME];
            double final_time   = 0.01*rCurrentProcessInfo[TIME_INCREASING_RATIO] * rCurrentProcessInfo[FINAL_SIMULATION_TIME]; 
            double current_time     = rCurrentProcessInfo[TIME];
          
            int& dummy_switch = rCurrentProcessInfo[INT_DUMMY_1];
      
            applied_force = AuxiliaryFunctions::LinearTimeIncreasingFunction(external_total_applied_force,initial_time,current_time,final_time,dummy_switch);

          }
          
              
          //temporaly modulation of the applied force:
                   
          rhs  = mass*gravity + applied_force;

          total_forces = rhs;
                  
          size_t iContactForce = 0;

          array_1d<double, 3 > AngularVel = this->GetGeometry()(0)->FastGetSolutionStepValue(ANGULAR_VELOCITY);
          double momentofinertia = GetGeometry()(0)->FastGetSolutionStepValue(PARTICLE_MOMENT_OF_INERTIA);
          
                        
          double RotaAcc[3] = {0.0};
          double Initial_Rota_Moment[3] = {0.0};
          double Max_Rota_Moment[3] = {0.0};
                
          if ( rotation_OPTION == 1 )
          {
              RotaAcc[0] = AngularVel[0] / dt;
              RotaAcc[1] = AngularVel[1] / dt;
              RotaAcc[2] = AngularVel[2] / dt;
      
              Initial_Rota_Moment[0] = RotaAcc[0] * momentofinertia;
              Initial_Rota_Moment[1] = RotaAcc[1] * momentofinertia;
              Initial_Rota_Moment[2] = RotaAcc[2] * momentofinertia;
                
              Max_Rota_Moment[0] = Initial_Rota_Moment[0];
              Max_Rota_Moment[1] = Initial_Rota_Moment[1];
              Max_Rota_Moment[2] = Initial_Rota_Moment[2];
            }  
                      
             
          for(ParticleWeakIteratorType neighbour_iterator = r_neighbours.begin();
              neighbour_iterator != r_neighbours.end(); neighbour_iterator++)
          {
            

              // GETTING NEIGHBOUR PROPERTIES

              
              double other_radius                 = neighbour_iterator->GetGeometry()(0)->GetSolutionStepValue(RADIUS);
              double other_density                = neighbour_iterator->GetGeometry()(0)->GetSolutionStepValue(PARTICLE_DENSITY);
              double other_mass                   = 4.0 / 3.0 * M_PI * other_radius * other_radius * other_radius * other_density;
              double equiv_mass                   = sqrt(mass*other_mass);//(mass*other_mass*(mass+other_mass)) / ((mass+other_mass)*(mass+other_mass)); //I: calculated by Roberto Flores
              
              double other_restitution_coeff      = neighbour_iterator->GetGeometry()(0)->GetSolutionStepValue(RESTITUTION_COEFF);
              double equiv_visco_damp_coeff_normal;       //= (visco_damp_coeff + other_visco_damp_coeff) / 2.0;   //M: is it correct to be a simple mean.
              double equiv_visco_damp_coeff_tangential;
              double equiv_restitution_coeff      = sqrt(restitution_coeff * other_restitution_coeff); //I: we assume this.
            
              double other_young                  = neighbour_iterator->GetGeometry()[0].GetSolutionStepValue(YOUNG_MODULUS);
              double other_poisson                = neighbour_iterator->GetGeometry()[0].GetSolutionStepValue(POISSON_RATIO);
              //double other_tension                = neighbour_iterator->GetGeometry()[0].GetSolutionStepValue(PARTICLE_TENSION);
              //double other_cohesion               = neighbour_iterator->GetGeometry()[0].GetSolutionStepValue(PARTICLE_COHESION);
              double other_FriAngle               = neighbour_iterator->GetGeometry()[0].GetSolutionStepValue(PARTICLE_FRICTION);
              double equiv_FriAngle               = (FriAngle + other_FriAngle) * 0.5; 
              
              
              // CONTINUUM SIMULATING PARAMETERS:

              
              double initial_delta = 0.0;
              //double CTension = 0.0;
              //double CCohesion = 0.0;

              
              //SLIDING

              
              bool sliding = false;

              if( delta_OPTION && (iContactForce < r_VectorContactInitialDelta.size()) )
              {
                  initial_delta = r_VectorContactInitialDelta[iContactForce];
              }

              
              // BASIC CALCULATIONS

              
              array_1d<double,3> other_to_me_vect = this->GetGeometry()(0)->Coordinates() - neighbour_iterator->GetGeometry()(0)->Coordinates();
              double distance                     = sqrt(other_to_me_vect[0] * other_to_me_vect[0] +
                                                    other_to_me_vect[1] * other_to_me_vect[1] +
                                                    other_to_me_vect[2] * other_to_me_vect[2]);

              double radius_sum                   = radius + other_radius;
              double indentation                  = radius_sum - distance - initial_delta; //M: Here, Initial_delta is expected to be positive if it is embeding and negative if it's separation.

              double equiv_radius     = 2*radius * other_radius / (radius + other_radius);
              //we now take 1/2 of the efective radius.
              double equiv_area       = (0.25)*M_PI * equiv_radius * equiv_radius; // 0.25 is becouse we take only the half of the equivalent radius, corresponding to the case of one sphere with radius Requivalent and other = radius 0.
              double equiv_poisson    = 2* poisson * other_poisson / (poisson + other_poisson);
              double equiv_young      = 2 * young * other_young / (young + other_young);
              
              double corrected_area = equiv_area;
              double this_poisson_contribution = 0.0;
              double neigh_poisson_contribution = 0.0;
              
              int size_ini_cont_neigh = this->GetValue(CONTINUUM_INI_NEIGHBOURS_IDS).size();

              
              //for the updating steps...... //teporarily set as only the first 
              
              Element::Pointer lock_p_weak;
              
              if(continuum_simulation_OPTION)
              {

                  for (int index_area=0; index_area<size_ini_cont_neigh; index_area++)
                  {

                    //MPI_CARLOS_MIQUEL DESCOMENTAR EL ELSE I EL IF, QUE SEMPRE FACI EL IF...............
                    
                      if ( this->GetValue(CONTINUUM_INI_NEIGHBOURS_IDS)[index_area] == int(neighbour_iterator->Id()) ) 
                      {
                           if(rCurrentProcessInfo[CONTACT_MESH_OPTION]==1)
                           {
                           lock_p_weak = (this->GetGeometry()[0].GetValue(NODE_TO_NEIGH_ELEMENT_POINTER)(index_area)).lock();
                           }
                          
                       // if(rCurrentProcessInfo[TIME_STEPS]==0 || rCurrentProcessInfo[CONTACT_MESH_OPTION]==0 )//(1<2)//rCurrentProcessInfo[TIME_STEPS]==0) //MIQUEL: NO BARRES:
                          {

                              corrected_area = mcont_ini_neigh_area[index_area];
                             
                             this_poisson_contribution = 0.0;
                             neigh_poisson_contribution = 0.0;
                             
                              break;
                              
                          } //for the updating steps //THESE STEPS SHOULD BE DONE OUTSIDE THE CALCULATION BECOUSE THEY WOULD HAVE DIFFERENT FORCES.
                                                  
                        
                        /*
                        else 
                          {

                              Element::Pointer lock_p_weak = (this->GetGeometry()[0].GetValue(NODE_TO_NEIGH_ELEMENT_POINTER)(index_area)).lock();

                              corrected_area = lock_p_weak->GetValue(MEAN_CONTACT_AREA);
                              
                              
                              if (this->Id() < neighbour_iterator->Id())
                              {
                                
                         
                                this_poisson_contribution   = lock_p_weak->GetValue(LOW_POISSON_FORCE); 
                                neigh_poisson_contribution  = lock_p_weak->GetValue(HIGH_POISSON_FORCE); 
                                
                              }
                              
                              else
                              {
                         
                                this_poisson_contribution   = lock_p_weak->GetValue(HIGH_POISSON_FORCE) ; 
                                neigh_poisson_contribution  = lock_p_weak->GetValue(LOW_POISSON_FORCE) ;
   
                                if(rCurrentProcessInfo[TIME_STEPS] ==500)
                                {
                                KRATOS_WATCH(this->Id())
                                KRATOS_WATCH(lock_p_weak->GetValue(HIGH_POISSON_FORCE))
                                KRATOS_WATCH(lock_p_weak->GetValue(LOW_POISSON_FORCE))
                                }
                              }
                              
                             break;
 
                          }//for the known steps....
                       */
                       //MPI_CARLOS_MIQUEL DESCOMENTAR EL ELSE I EL IF, QUE SEMPRE FACI EL IF...............
                       
                       
                      }// if ( this->GetValue(CONTINUUM_INI_NEIGHBOURS_IDS)[index_area] == int(neighbour_iterator->Id()) ) 
                      
                
                    
                  }//for every neighbour      

              }//if(continuum_simulation_OPTION)

              
         

         
              //MACRO PARAMETERS

              double kn               = magic_factor*equiv_young*corrected_area/(radius + other_radius); //M_PI * 0.5 * equiv_young * equiv_radius; //M: CANET FORMULA    
              double ks               = kn / (2.0 * (1.0 + equiv_poisson));
           

              if (global_variables_OPTION == 1) //globally defined parameters       // ha de ser canviat aixo pk ara rn i rt no entren al calcul
              {

                // KRATOS_WATCH("-------GLOBAL VARIABLES----------")

                  kn       = rCurrentProcessInfo[GLOBAL_KN];
                  ks       = rCurrentProcessInfo[GLOBAL_KT];
                  //RN       = rCurrentProcessInfo[GLOBAL_RN];
                  //RT_base  = rCurrentProcessInfo[GLOBAL_RT];
                  //Friction = tan( (rCurrentProcessInfo[GLOBAL_FRI_ANG])*(M_PI/180) );

              }

              //historical minimun K for the critical time:
              if (rCurrentProcessInfo[CRITICAL_TIME_OPTION]==1)

              {
                  double historic = rCurrentProcessInfo[HISTORICAL_MIN_K];
                  if( (kn<historic) || (ks<historic))
                  {
                      historic = GeometryFunctions::min(kn,ks);
                  }

              }

              if(equiv_restitution_coeff>0)
              {
                  equiv_visco_damp_coeff_normal  = -( (2*log(equiv_restitution_coeff)*sqrt(equiv_mass*kn)) / (sqrt( (log(equiv_restitution_coeff)*log(equiv_restitution_coeff)) + (M_PI*M_PI) )) );
                  equiv_visco_damp_coeff_tangential  = -( (2*log(equiv_restitution_coeff)*sqrt(equiv_mass*ks)) / (sqrt( (log(equiv_restitution_coeff)*log(equiv_restitution_coeff)) + (M_PI*M_PI) )) );
              }
              else 
              {
                  equiv_visco_damp_coeff_normal  = ( 2*sqrt(equiv_mass*kn) );
                  equiv_visco_damp_coeff_tangential  = ( 2*sqrt(equiv_mass*ks) );
              }
              
                
              // FORMING LOCAL CORDINATES
              
              //Notes: Since we will normally inherit the mesh from GiD, we respect the global system X,Y,Z [0],[1],[2]
              //In the local coordinates we will define the normal direction of the contact as the [2] component!!!!!
              //the way the normal direction is defined (other_to_me_vect) compression is positive

              double NormalDir[3]           = {0.0};
              double LocalCoordSystem[3][3] = {{0.0}, {0.0}, {0.0}};
              NormalDir[0] = other_to_me_vect[0];  
              NormalDir[1] = other_to_me_vect[1];
              NormalDir[2] = other_to_me_vect[2];
              
              
              GeometryFunctions::ComputeContactLocalCoordSystem(NormalDir, LocalCoordSystem); //new Local Coord System
        
              // FORMING OLD LOCAL CORDINATES
            
             array_1d<double,3> old_coord_target     = this->GetGeometry()(0)->GetInitialPosition() + this->GetGeometry()(0)->GetSolutionStepValue(DISPLACEMENT,1);
             array_1d<double,3> old_coord_neigh      = neighbour_iterator->GetGeometry()(0)->GetInitialPosition()+neighbour_iterator->GetGeometry()(0)->GetSolutionStepValue(DISPLACEMENT,1);
             array_1d<double,3> Old_other_to_me_vect = old_coord_target - old_coord_neigh; 
             //
             double OldNormalDir[3]           = {0.0};
              double OldLocalCoordSystem[3][3] = {{0.0}, {0.0}, {0.0}};
              OldNormalDir[0] = Old_other_to_me_vect[0];   // M. this way the compresion is positive.
              OldNormalDir[1] = Old_other_to_me_vect[1];
              OldNormalDir[2] = Old_other_to_me_vect[2];
              GeometryFunctions::ComputeContactLocalCoordSystem(OldNormalDir, OldLocalCoordSystem); //Old Local Coord System
              
              
              // VELOCITIES AND DISPLACEMENTS

              
              array_1d<double, 3 > other_vel      = neighbour_iterator->GetGeometry()(0)->GetSolutionStepValue(VELOCITY);            
              array_1d<double, 3 > other_delta_displ      = neighbour_iterator->GetGeometry()(0)->GetSolutionStepValue(DELTA_DISPLACEMENT);

             
              
              double DeltDisp[3] = {0.0};
              double RelVel [3] = {0.0};

              RelVel[0] = (vel[0] - other_vel[0]);
              RelVel[1] = (vel[1] - other_vel[1]);
              RelVel[2] = (vel[2] - other_vel[2]);

              //DeltDisp in global cordinates

              DeltDisp[0] = (delta_displ[0] - other_delta_displ[0]);
              DeltDisp[1] = (delta_displ[1] - other_delta_displ[1]);
              DeltDisp[2] = (delta_displ[2] - other_delta_displ[2]);


              if ( rotation_OPTION == 1 )
              {

                  double velA[3]   = {0.0};
                  double velB[3]   = {0.0};
                  double dRotaDisp[3] = {0.0};

                  array_1d<double, 3 > AngularVel       = this->GetGeometry()(0)->FastGetSolutionStepValue(ANGULAR_VELOCITY);
                  array_1d<double, 3 > Other_AngularVel = neighbour_iterator->GetGeometry()(0)->FastGetSolutionStepValue(ANGULAR_VELOCITY);

                  double Vel_Temp[3]       = {      AngularVel[0],       AngularVel[1],       AngularVel[2]};
                  double Other_Vel_Temp[3] = {Other_AngularVel[0], Other_AngularVel[1], Other_AngularVel[2]};
                  GeometryFunctions::CrossProduct(Vel_Temp,             OldLocalCoordSystem[2], velA); //it was Local Coordinate system, now we do OLD. 
                  GeometryFunctions::CrossProduct(Other_Vel_Temp, OldLocalCoordSystem[2], velB);

                  dRotaDisp[0] = -velA[0] * radius - velB[0] * other_radius;
                  dRotaDisp[1] = -velA[1] * radius - velB[1] * other_radius;
                  dRotaDisp[2] = -velA[2] * radius - velB[2] * other_radius;
                  //////contribution of the rotation vel
                  DeltDisp[0] += dRotaDisp[0] * dt;
                  DeltDisp[1] += dRotaDisp[1] * dt;
                  DeltDisp[2] += dRotaDisp[2] * dt;

              }//if rotation_OPTION

              double LocalDeltDisp[3] = {0.0};
              double LocalContactForce[3]  = {0.0};
              double GlobalContactForce[3] = {0.0};
              double LocalRelVel[3] = {0.0};
              //double GlobalContactForceOld[3] = {0.0};


              GlobalContactForce[0] = GlobalContactForceMatrix[iContactForce][0];   //M:aqui tenim guardades les del neighbour calculator.
              GlobalContactForce[1] = GlobalContactForceMatrix[iContactForce][1];
              GlobalContactForce[2] = GlobalContactForceMatrix[iContactForce][2];

              /*
              
              if(this->Id() == 4368 && neighbour_iterator->Id() == 4619)
              {   
                  //std::cout << "12 IND " << indentation - (1 - sqrt(other_to_me_vect[0] * other_to_me_vect[0] + other_to_me_vect[1] * other_to_me_vect[1] + other_to_me_vect[2] * other_to_me_vect[2])) << " " << GlobalContactForce[1] << std::endl;
                  std::cout << "4368 MY_COORD " << std::setprecision(13) << this->GetGeometry()(0)->Coordinates()[1] << std::endl;    
                  std::cout << "4368 NEIGH_CO0RD " << std::setprecision(13) << neighbour_iterator->GetGeometry()(0)->Coordinates()[1] << std::endl; 
              }
              if(this->Id() == 4619 && neighbour_iterator->Id() ==4368)
              {   
                  //std::cout << "11 IND " << indentation - (1 - sqrt(other_to_me_vect[0] * other_to_me_vect[0] + other_to_me_vect[1] * other_to_me_vect[1] + other_to_me_vect[2] * other_to_me_vect[2])) << " " << GlobalContactForce[1] << std::endl;
                  std::cout << "4619 MY_COORD " << std::setprecision(13) << this->GetGeometry()(0)->Coordinates()[1] << std::endl;    
                  std::cout << "4619 NEIGH_CO0RD " << std::setprecision(13) << neighbour_iterator->GetGeometry()(0)->Coordinates()[1] << std::endl; 
                  
                  KRATOS_WATCH("____________________________________________________")
                  
              }
              
              */              
                            
              /*
              if(this->Id() == 12 && neighbour_iterator->Id() == 11)
              {   
                  //std::cout << "12 IND " << indentation - (1 - sqrt(other_to_me_vect[0] * other_to_me_vect[0] + other_to_me_vect[1] * other_to_me_vect[1] + other_to_me_vect[2] * other_to_me_vect[2])) << " " << GlobalContactForce[1] << std::endl;
                  std::cout << "12 MY_COORD " << std::setprecision(13) << this->GetGeometry()(0)->Coordinates()[1] << std::endl;    
                  std::cout << "12 NEIGH_CO0RD " << std::setprecision(13) << neighbour_iterator->GetGeometry()(0)->Coordinates()[1] << std::endl; 
              }
              if(this->Id() == 11 && neighbour_iterator->Id() == 12)
              {   
                  //std::cout << "11 IND " << indentation - (1 - sqrt(other_to_me_vect[0] * other_to_me_vect[0] + other_to_me_vect[1] * other_to_me_vect[1] + other_to_me_vect[2] * other_to_me_vect[2])) << " " << GlobalContactForce[1] << std::endl;
                  std::cout << "11 MY_COORD " << std::setprecision(13) << this->GetGeometry()(0)->Coordinates()[1] << std::endl;    
                  std::cout << "11 NEIGH_CO0RD " << std::setprecision(13) << neighbour_iterator->GetGeometry()(0)->Coordinates()[1] << std::endl; 
              }
            
*/
             //GeometryFunctions::VectorGlobal2Local(OldLocalCoordSystem, GlobalContactForce, LocalContactForce); //we recover this way the old local forces projected in the new coordinates in the way they were in the old ones; Now they will be increased if its the necessary
             
           
             
             GeometryFunctions::VectorGlobal2Local(OldLocalCoordSystem, GlobalContactForce, LocalContactForce); //we recover this way the old local forces projected in the new coordinates in the way they were in the old ones; Now they will be increased if its the necessary
             GeometryFunctions::VectorGlobal2Local(OldLocalCoordSystem, DeltDisp, LocalDeltDisp);
             GeometryFunctions::VectorGlobal2Local(OldLocalCoordSystem, RelVel, LocalRelVel);
             
              // FORCES
           
              if ( (indentation > 0.0) || (this->GetValue(PARTICLE_CONTACT_FAILURE_ID)[iContactForce] == 0) )   // for detached particles we enter only if the indentation is > 0.
              {                                                                                                 // for attached particles we enter only if the particle is still attached.
            
                  // NORMAL FORCE
                  
                  switch (force_calculation_type_id) //  0---linear comp & tension ; 1 --- Hertzian (no linear comp, linear tension)
                  {
                        case 0:
                            
                            LocalContactForce[2]= kn * indentation;
                            
                            if(rCurrentProcessInfo[NON_LINEAR_OPTION])
                            {
                             
                              
                                if (this->Id() == 1 && rCurrentProcessInfo[TIME_STEPS] == 4 ) KRATOS_WATCH( "MUST BE IMPROVED, THE CALCULATION OF STRESS IN THE NON_LINEAR_OPTION ONLY TAKES INTO ACCOUNT THE INDENTATION, IS IT OKAY??" )
                                  
                                double kn_b = kn/rCurrentProcessInfo[SLOPE_FRACTION_N1];
                                double kn_c = kn/rCurrentProcessInfo[SLOPE_FRACTION_N2];
                                
                                
                                double compression_limit_1 = rCurrentProcessInfo[SLOPE_LIMIT_COEFF_C1]*rCurrentProcessInfo[CONTACT_SIGMA_MAX]*1e6;
                                double compression_limit_2 = rCurrentProcessInfo[SLOPE_LIMIT_COEFF_C2]*rCurrentProcessInfo[CONTACT_SIGMA_MAX]*1e6;
                                
                                double sigma_a = (kn * indentation)/(corrected_area);
//                                 if(this->Id()==3927)
//                                 {
//                                    KRATOS_WATCH(rCurrentProcessInfo[SLOPE_LIMIT_COEFF_C1])
//                                   KRATOS_WATCH(rCurrentProcessInfo[SLOPE_FRACTION_N1])
//                                 KRATOS_WATCH(sigma_a)
//                                 KRATOS_WATCH(compression_limit_1)
//                                 KRATOS_WATCH(" ")
//                                 }
                             
                                double sigma_b = compression_limit_1 + kn_b*(indentation/corrected_area - compression_limit_1/kn);
                                
                                
                              
                                if( (indentation >= 0.0) && (sigma_a < compression_limit_1) ) 
                                {
                                    LocalContactForce[2]= kn * indentation;
                                    
                                    if(rCurrentProcessInfo[CONTACT_MESH_OPTION])
                                    
                                    {
                                      lock_p_weak->GetValue(NON_ELASTIC_STAGE) = 1.0;
                                    }
                                    
                                }
                                
                                else if( (indentation >= 0.0) && (sigma_a >= compression_limit_1) && ( sigma_b < compression_limit_2 ) ) 
                                {
                                    LocalContactForce[2]= compression_limit_1*corrected_area + kn_b*(indentation - compression_limit_1*corrected_area/kn);
                                     
                                    if(rCurrentProcessInfo[CONTACT_MESH_OPTION])
                                    
                                    {
                                    lock_p_weak->GetValue(NON_ELASTIC_STAGE) = 2.0;
                                    
                                    }
                                }
                                
                                else if ( indentation >= 0.0 ) 
                                {
                                    
                                    LocalContactForce[2]= compression_limit_2*corrected_area + kn_c*(indentation - corrected_area*(compression_limit_1/kn + compression_limit_2/kn_b - compression_limit_1/kn_b));
  
                                }
                                
                                else {LocalContactForce[2]= kn * indentation; }
                            }
                            break;
                                        
                        case 1:

                            if(indentation >= 0.0) {LocalContactForce[2]= kn * pow(indentation, 1.5); }
                            else { LocalContactForce[2] = kn * indentation; }
                            break;

                        default:
                            if(indentation >= 0.0) {LocalContactForce[2]= kn * indentation; }
                            else {LocalContactForce[2]= kn * indentation; }
                            break;

                  }
                    
                  // TANGENTIAL FORCE: incremental calculation. YADE develops a complicated "absolute method"

                  LocalContactForce[0] += - ks * LocalDeltDisp[0];  // 0: first tangential
                  LocalContactForce[1] += - ks * LocalDeltDisp[1];  // 1: second tangential
              }


              if ( (indentation <= 0.0) && (this->GetValue(PARTICLE_CONTACT_FAILURE_ID)[iContactForce] != 0) )
              {                      
                  LocalContactForce[0] = 0.0;  // 0: first tangential
                  LocalContactForce[1] = 0.0;  // 1: second tangential
                  LocalContactForce[2] = 0.0;  // 2: normal force
              }
                

              //EVALUATING THE POSSIBLE FAILURE FOR THE CONTINUUM CONTACTS
              
              double contact_tau = 0.0;
              double contact_sigma = 0.0;
              
              double failure_criterion_state = 0.0;
                        
              double DYN_FRI_ANG = equiv_FriAngle*M_PI/180; //entrar el valor des de el material de gid
              
              //double compression_limit        = rCurrentProcessInfo[CONTACT_SIGMA_MAX]*1e6;
              double tension_limit            = rCurrentProcessInfo[CONTACT_SIGMA_MIN]*1e6;
                                      
              double ShearForceNow = sqrt(LocalContactForce[0] * LocalContactForce[0]
                                    +      LocalContactForce[1] * LocalContactForce[1]); 
              
              double Frictional_ShearForceMax = tan(DYN_FRI_ANG) * LocalContactForce[2];
              
              
              if(Frictional_ShearForceMax < 0.0){Frictional_ShearForceMax = 0.0;}
              
              
              if ( this->GetValue(PARTICLE_CONTACT_FAILURE_ID)[iContactForce] != 0 )
              {
                  failure_criterion_state = 1.0;
                                                        
                  if( (ShearForceNow >  Frictional_ShearForceMax) && (ShearForceNow != 0.0) ) 
                  {
                      
                      LocalContactForce[0] = (Frictional_ShearForceMax / ShearForceNow) * LocalContactForce[0];
                      LocalContactForce[1] = (Frictional_ShearForceMax / ShearForceNow )* LocalContactForce[1];
                      sliding = true ;
                                  
                  }
                  
              }

  
              if(corrected_area < 1e-09) {KRATOS_WATCH(corrected_area) KRATOS_WATCH(this->Id())}
              
              if(this->GetValue(PARTICLE_CONTACT_FAILURE_ID)[iContactForce] == 0)         
              {  
                  int failure_criterion_OPTION = rCurrentProcessInfo[FAILURE_CRITERION_OPTION];
            
                  ///(1) MOHR-COULOMB FAILURE: (we don't consider rotational spring!!!!! here) need to be thought.
          
                  if (failure_criterion_OPTION==1)  //MOHR-COULOMB
                  {   
                      contact_tau = ShearForceNow/(corrected_area);
                      contact_sigma = LocalContactForce[2]/(corrected_area);

                      double sigma_max, sigma_min;

                      if (LocalContactForce[2]>=0)
                      {
                          sigma_max = contact_sigma;
                          sigma_min = 0;
                      }
                      else 
                      {
                          sigma_max = 0;
                          sigma_min = contact_sigma;
                      }

                      //change into principal stresses

                      double centre = 0.5*(sigma_max + sigma_min);
                      double radius = sqrt( (sigma_max - centre)*(sigma_max - centre) + contact_tau*contact_tau   ) ;

                      double sigma_I = centre + radius;
                      double sigma_II = centre - radius;
   
                      // Check:

                      //double tau_zero = 0.5*sqrt(compression_limit*tension_limit); 
                      double tau_zero = rCurrentProcessInfo[CONTACT_TAU_ZERO]*1e6;
            
                      //double Failure_FriAngle =  atan((compression_limit-tension_limit)/(2*sqrt(compression_limit*tension_limit)));
                      double Failure_FriAngle =  rCurrentProcessInfo[CONTACT_INTERNAL_FRICC]*M_PI/180;
          
                      double distance_to_failure = ( tau_zero/(tan(Failure_FriAngle)) + centre )*sin(Failure_FriAngle);
                  
                      failure_criterion_state = radius/distance_to_failure;
            
                     if ( sigma_I - sigma_II >= 2*tau_zero*cos(Failure_FriAngle) + (sigma_I + sigma_II)*sin(Failure_FriAngle) )
                      {

                          //breaks

                          this->GetValue(PARTICLE_CONTACT_FAILURE_ID)[iContactForce] = 5; //mohr coulomb
                            
                          //tangential mapping, divide 2 tangent
            
                          if((this->GetValue(PARTICLE_CONTACT_FAILURE_ID)[iContactForce])!=0)
                          {
                              //KRATOS_WATCH(this->GetValue(PARTICLE_CONTACT_FAILURE_ID)[iContactForce])
                              //KRATOS_WATCH(lock_p_weak->GetValue(CONTACT_FAILURE_LOW))
                          }
                
                          
                          sliding = true ;

                      }
                      else
                      {
                          // doesn't brake
                      }
                  } //MOHR-COULOMB
            
                  ///(2) UNCOUPLED FRACTURE
            
                  ///vam decidir amb miguel angel de no fer el mapping de les shear fins al pas seguent.. esta correcte? afecta quan trenca?
            
                  if (failure_criterion_OPTION==2)//UNCOUPLED FRACTURE
                  {    
                    
                      contact_tau = ShearForceNow/(corrected_area);
                      contact_sigma = LocalContactForce[2]/(corrected_area);

                      //double tau_zero = 0.5*sqrt(compression_limit*tension_limit); 
                      double tau_zero = rCurrentProcessInfo[CONTACT_TAU_ZERO]*1e6;
            
                      //double Failure_FriAngle =  atan((compression_limit-tension_limit)/(2*sqrt(compression_limit*tension_limit)));
                      double Failure_FriAngle =  rCurrentProcessInfo[CONTACT_INTERNAL_FRICC]*M_PI/180;
            
                      if (LocalContactForce[2]>=0)
                      {
                          double tau_strength = tau_zero+tan(Failure_FriAngle)*contact_sigma;
                          
                          failure_criterion_state = contact_tau/tau_strength;
                                                
                          if(contact_tau>tau_strength)
                          {
                              this->GetValue(PARTICLE_CONTACT_FAILURE_ID)[iContactForce] = 2;
                              sliding = true;
                          }
                      } //positive sigmas
                      else //negative sigmas
                      {
                          double tau_strength = tau_zero;
                
                          failure_criterion_state = GeometryFunctions::max(contact_tau/tau_strength, -contact_sigma/tension_limit) ;
                
                          if(contact_tau > tau_strength)
                          {
                              this->GetValue(PARTICLE_CONTACT_FAILURE_ID)[iContactForce] = 3;  //shear failure
                              sliding = true;
                  
                              if(contact_sigma<-tension_limit)
                              {
                                  this->GetValue(PARTICLE_CONTACT_FAILURE_ID)[iContactForce] = 12;
                              } //both shear and tension
                  
                          }
                          else if(contact_sigma<-tension_limit)
                          {
                              this->GetValue(PARTICLE_CONTACT_FAILURE_ID)[iContactForce] = 4; //tension failure
                              sliding = true;
                          }
                          
                      } //negative values of sigma              
              
                  } //UNCOUPLED FRACTURE
                    
              }// if (this->GetValue(PARTICLE_CONTACT_FAILURE_ID)[iContactForce] == 0)
                
              if(this->GetValue(PARTICLE_CONTACT_FAILURE_ID)[iContactForce] != 0 && rCurrentProcessInfo[ACTIVATE_SEARCH]==0)
              {
                  rCurrentProcessInfo.SetValue(ACTIVATE_SEARCH, 1);

                  KRATOS_WATCH(" ")
                  KRATOS_WATCH("-------->From now on, searching neighbours, some contacs have failed<-------")
                  KRATOS_WATCH("Time step:")
                  KRATOS_WATCH(rCurrentProcessInfo[TIME_STEPS])
                  KRATOS_WATCH("Particle_1")
                  KRATOS_WATCH(this->Id())
                  KRATOS_WATCH("Particle_2")
                  KRATOS_WATCH(neighbour_iterator->Id())      
                  KRATOS_WATCH(" ")
              }


          // VISCODAMPING (applyied locally)

                    //*** the compbrobation is component-wise since localContactForce and RelVel have in principle no relationship.
                    // the visco force can be higher than the contact force only if they go to the same direction. (in my opinion)
                    // but in oposite direction the visco damping can't overpass the force...


              double ViscoDampingLocalContactForce[3]    = {0.0};
        
              if ( (damp_id > 0  ) && ( (indentation > 0.0) || (this->GetValue(PARTICLE_CONTACT_FAILURE_ID)[iContactForce] == 0) ) )
              {
                
                  if (damp_id == 11 || damp_id == 10)
                  {
                      ViscoDampingLocalContactForce[2] = - equiv_visco_damp_coeff_normal * LocalRelVel[2];
                  }
                  
                  for (unsigned int index = 0; index < 2; index++)
                  {

                      if(sliding == false && ( damp_id == 1 || damp_id ==11) ) //only applied when no sliding to help to the regularized friccion law or the spring convergence
                      {
                          ViscoDampingLocalContactForce[index] = - equiv_visco_damp_coeff_tangential * LocalRelVel[index];  
                      }
                  }
              }

              // TRANSFORMING TO GLOBAL FORCES AND ADDING UP


              double LocalResultantContactForce[3] ={0.0};
              double ViscoDampingGlobalContactForce[3] = {0.0};
              double GlobalResultantContactForce[3] = {0.0};

              for (unsigned int index = 0; index < 3; index++)
              {
                  LocalResultantContactForce[index] = LocalContactForce[index]  + ViscoDampingLocalContactForce[index];
              }
              
              double PoissonContactForce[3]        = {0.0};
              double GlobalContactPoissonForce[3]  = {0.0};
              PoissonContactForce [2] = neigh_poisson_contribution;

              GeometryFunctions::VectorLocal2Global(LocalCoordSystem, LocalContactForce, GlobalContactForce);
              GeometryFunctions::VectorLocal2Global(LocalCoordSystem, ViscoDampingLocalContactForce, ViscoDampingGlobalContactForce);
              GeometryFunctions::VectorLocal2Global(LocalCoordSystem, LocalResultantContactForce, GlobalResultantContactForce);
              GeometryFunctions::VectorLocal2Global(LocalCoordSystem, PoissonContactForce, GlobalContactPoissonForce);

              rhs[0] += GlobalContactForce[0]; //RHS
              rhs[1] += GlobalContactForce[1];
              rhs[2] += GlobalContactForce[2];

              total_forces[0] += GlobalResultantContactForce[0] + GlobalContactPoissonForce[0];
              total_forces[1] += GlobalResultantContactForce[1] + GlobalContactPoissonForce[1];
              total_forces[2] += GlobalResultantContactForce[2] + GlobalContactPoissonForce[2];

              damp_forces[0] += ViscoDampingGlobalContactForce[0];
              damp_forces[1] += ViscoDampingGlobalContactForce[1];
              damp_forces[2] += ViscoDampingGlobalContactForce[2];

                          
              // SAVING CONTACT FORCES FOR NEXT STEPS

              this->GetValue(PARTICLE_CONTACT_FORCES)[iContactForce][0] = GlobalContactForce[0];
              this->GetValue(PARTICLE_CONTACT_FORCES)[iContactForce][1] = GlobalContactForce[1];
              this->GetValue(PARTICLE_CONTACT_FORCES)[iContactForce][2] = GlobalContactForce[2];
          
              if ( rotation_OPTION == 1 )
              {
                
                     array_1d<double, 3 > & mRota_Moment = this->GetGeometry()[0].GetSolutionStepValue(PARTICLE_MOMENT);
                    double Rota_Moment[3] = {0.0};
                        
                    Rota_Moment[0] = mRota_Moment[0];
                    Rota_Moment[1] = mRota_Moment[1];
                    Rota_Moment[2] = mRota_Moment[2];
                       
                    double MA[3] = {0.0};
                        
                    //GeometryFunctions::CrossProduct(LocalCoordSystem[2], GlobalResultantContactForce, MA);
                    GeometryFunctions::CrossProduct(LocalCoordSystem[2], GlobalContactForce, MA);
                                          
                    Rota_Moment[0] -= MA[0] * radius;
                    Rota_Moment[1] -= MA[1] * radius;
                    Rota_Moment[2] -= MA[2] * radius;
                   
                    if(rCurrentProcessInfo[ROTA_DAMP_TYPE]==2)  //Rolling friccion type   
                    {
                       double RollingFriction              = this->GetGeometry()[0].GetSolutionStepValue(ROLLING_FRICTION);
                       double RollingFrictionCoeff         = RollingFriction * radius;  
                       double other_RollingFriction        = neighbour_iterator->GetGeometry()[0].GetSolutionStepValue(ROLLING_FRICTION);
                       double other_RollingFrictionCoeff   = other_RollingFriction * other_radius;
                       double equiv_RollingFrictionCoeff   = 2*RollingFrictionCoeff * other_RollingFrictionCoeff / (RollingFrictionCoeff + other_RollingFrictionCoeff); 
                    
                       if (equiv_RollingFrictionCoeff != 0.0)
                       {
                          Max_Rota_Moment[0] += Rota_Moment[0];
                          Max_Rota_Moment[1] += Rota_Moment[1];
                          Max_Rota_Moment[2] += Rota_Moment[2];                     
                          
                          double CoordSystemMoment[3] = {0.0};
                          double MR[3] = {0.0};
                          
                          double NormalForce[3] = {0.0};      
                          
                          //NormalForce[0] = LocalCoordSystem[2][0] * fabs(LocalResultantContactForce[2]);
                          //NormalForce[1] = LocalCoordSystem[2][1] * fabs(LocalResultantContactForce[2]);
                          //NormalForce[2] = LocalCoordSystem[2][2] * fabs(LocalResultantContactForce[2]);
                          
                          NormalForce[0] = LocalCoordSystem[2][0] * fabs(LocalContactForce[2]);
                          NormalForce[1] = LocalCoordSystem[2][1] * fabs(LocalContactForce[2]);
                          NormalForce[2] = LocalCoordSystem[2][2] * fabs(LocalContactForce[2]);
                          
                          GeometryFunctions::CrossProduct(LocalCoordSystem[2], Max_Rota_Moment, CoordSystemMoment);
                          
                          double DetCoordSystemMoment = sqrt(CoordSystemMoment[0] * CoordSystemMoment[0] + CoordSystemMoment[1] * CoordSystemMoment[1] + CoordSystemMoment[2] * CoordSystemMoment[2]);
                          
                          CoordSystemMoment[0] = CoordSystemMoment[0] / DetCoordSystemMoment;
                          CoordSystemMoment[1] = CoordSystemMoment[1] / DetCoordSystemMoment;
                          CoordSystemMoment[2] = CoordSystemMoment[2] / DetCoordSystemMoment;                            
                                                    
                          GeometryFunctions::CrossProduct(NormalForce, CoordSystemMoment, MR);

                          double DetMR = sqrt( MR[0] * MR[0] + MR[1] * MR[1] + MR[2] * MR[2] );
                          double MR_now = DetMR * equiv_RollingFrictionCoeff;
                          double MR_max = sqrt( Max_Rota_Moment[0] * Max_Rota_Moment[0] + Max_Rota_Moment[1] * Max_Rota_Moment[1] + Max_Rota_Moment[2] * Max_Rota_Moment[2] );
                          
                          if ( MR_max > MR_now )
                          {
                              Rota_Moment[0] += MR[0] * equiv_RollingFrictionCoeff;
                              Rota_Moment[1] += MR[1] * equiv_RollingFrictionCoeff;
                              Rota_Moment[2] += MR[2] * equiv_RollingFrictionCoeff;
                              
                              Max_Rota_Moment[0] += MR[0] * equiv_RollingFrictionCoeff;
                              Max_Rota_Moment[1] += MR[1] * equiv_RollingFrictionCoeff;
                              Max_Rota_Moment[2] += MR[2] * equiv_RollingFrictionCoeff;                            
                          }
                        
                          else
                          {
                              Rota_Moment[0] = -Initial_Rota_Moment[0];
                              Rota_Moment[1] = -Initial_Rota_Moment[1];
                              Rota_Moment[2] = -Initial_Rota_Moment[2];
                          }                        
                       } // if (equiv_RollingFrictionCoeff != 0.0)
                    } //  if(rCurrentProcessInfo[ROTA_DAMP_TYPE]==2)
              
                    mRota_Moment[0] = Rota_Moment[0];
                    mRota_Moment[1] = Rota_Moment[1];
                    mRota_Moment[2] = Rota_Moment[2];               
               
                
          
                    
                    
              } //if ( rotation_OPTION == 1 )     
              
        
            //COMPUTE THE MEAN STRESS TENSOR:
            if(rCurrentProcessInfo[INT_DUMMY_9]==1)    //IF stress_strain_options ON
            {
              
              double gap                  = distance - radius_sum;
              
              array_1d<double,3> normal_vector_on_contact =  -1 * other_to_me_vect; //outwards
              
              double Dummy_Dummy = 0.0;
              GeometryFunctions::norm(normal_vector_on_contact,Dummy_Dummy); // Normalize to unitary module

              array_1d<double,3> coord_target    = this->GetGeometry()(0)->Coordinates();
              array_1d<double,3> coord_neigh     = neighbour_iterator->GetGeometry()(0)->Coordinates();
              //double neigh_radius                = neighbour_iterator->GetGeometry()(0)->GetSolutionStepValue(RADIUS);
              double target_radius               = this->GetGeometry()(0)->GetSolutionStepValue(RADIUS);
              array_1d<double,3> x_centroid      = (target_radius + 0.5*gap) * normal_vector_on_contact;
              
              //KRATOS_WATCH(kn)
              double result_product = GeometryFunctions::DotProduct(x_centroid,normal_vector_on_contact);
              
              Representative_Volume = Representative_Volume + 0.33333333333333 * (result_product * corrected_area);
              
                
                for (int i=0; i<3; i++)
                {
                  
                  for (int j=0; j<3; j++)
                  {
                    
                    mStressTensor[i][j] += (x_centroid[j]) * GlobalContactForce[i]; //ref: Katalin Bagi 1995 Mean stress tensor           
                    
                  }
                  
                }
            
            }   
  
        //CONTACT ELEMENT
              
              
              // Transfer values to the contact element
                    
              //obtenir el punter a la barra 
                   
              if(rCurrentProcessInfo[CONTACT_MESH_OPTION]==1)    
              {
                  int size_ini_cont_neigh = this->GetValue(CONTINUUM_INI_NEIGHBOURS_IDS).size();

                  for (int iii=0; iii<size_ini_cont_neigh; iii++)
                  {
                      if ( this->GetValue(CONTINUUM_INI_NEIGHBOURS_IDS)[iii] == int(neighbour_iterator->Id()) ) 
                      {
                          //obtaining pointer to contact element.
                                
                          Element::Pointer lock_p_weak = (this->GetGeometry()[0].GetValue(NODE_TO_NEIGH_ELEMENT_POINTER)(iii)).lock();
                          
                          if( this->Id() < neighbour_iterator->Id() )  // if id pequeña
                          {
                              //COPY VARIABLES LOW
                                                          
                              //storing values:
                                 
                              //HIGH-LOW variables
                                    
                              lock_p_weak->GetValue(LOCAL_CONTACT_FORCE_LOW)[0] = LocalContactForce[0];
                              lock_p_weak->GetValue(LOCAL_CONTACT_FORCE_LOW)[1] = LocalContactForce[1];
                              lock_p_weak->GetValue(LOCAL_CONTACT_FORCE_LOW)[2] = LocalContactForce[2];
                              
                              if(rCurrentProcessInfo[TIME_STEPS]==0)
                              {
                              lock_p_weak->GetValue(LOCAL_CONTACT_AREA_LOW) = corrected_area;
                              
                              }
                              
                              //COMBINED MEAN          
                    
                              lock_p_weak->GetValue(CONTACT_SIGMA) += 0.5*contact_sigma;
                              lock_p_weak->GetValue(CONTACT_TAU)   += 0.5*contact_tau;
                                                                                       
                              //UNIQUE VALUES
                              
                              //1) failure
                                  lock_p_weak->GetValue(CONTACT_FAILURE) = (this->GetValue(PARTICLE_CONTACT_FAILURE_ID)[iContactForce]);                                        
                                        
                                  if(failure_criterion_state<=1.0)
                                  {
                                      lock_p_weak->GetValue(FAILURE_CRITERION_STATE) = failure_criterion_state; 
                                  }
                                  
                                  else
                                  {
                                      //KRATOS_WATCH (failure_criterion_state )
                                  }   
                                                 
                          } // if Target Id < Neigh Id
                          else   
                          {
                              //COPY VARIABLES HIGH 
                                    
                              //HIGH-LOW variables
                                    
                              lock_p_weak->GetValue(LOCAL_CONTACT_FORCE_HIGH)[0] = LocalContactForce[0];
                              lock_p_weak->GetValue(LOCAL_CONTACT_FORCE_HIGH)[1] = LocalContactForce[1];
                              lock_p_weak->GetValue(LOCAL_CONTACT_FORCE_HIGH)[2] = LocalContactForce[2];
                              
                              if(rCurrentProcessInfo[TIME_STEPS]==0)
                              {

                              lock_p_weak->GetValue(LOCAL_CONTACT_AREA_HIGH) = corrected_area;

                              }
                              
                                                                    
                              //COMBINED MEAN       
                    
                              lock_p_weak->GetValue(CONTACT_SIGMA)                += 0.5*contact_sigma;
                              lock_p_weak->GetValue(CONTACT_TAU)                  += 0.5*contact_tau;
                              
                                                    
                          }
                          
                          //CONTACT AREA
                          
                            if ( ( rCurrentProcessInfo[TIME_STEPS]==0 ) && ( this->GetValue(SKIN_SPHERE)==0 ) && ( neighbour_iterator->GetValue(SKIN_SPHERE)==0 ) )
                               {
                                                             
                                 lock_p_weak->GetValue(MEAN_CONTACT_AREA)   += 0.5*corrected_area;

                               }
                               
                               else if ( ( rCurrentProcessInfo[TIME_STEPS]==0 ) && ( this->GetValue(SKIN_SPHERE)==1 ) && ( neighbour_iterator->GetValue(SKIN_SPHERE)==1 ) )
                               {
                            
                                 lock_p_weak->GetValue(MEAN_CONTACT_AREA)   += 0.5*corrected_area;
                                                             
                               }
                               
                                        
                               
                               else if ( ( rCurrentProcessInfo[TIME_STEPS]==0 ) && ( this->GetValue(SKIN_SPHERE)==0 ) && ( neighbour_iterator->GetValue(SKIN_SPHERE)==1 ) )
                               {
                               
                                lock_p_weak->GetValue(MEAN_CONTACT_AREA)   = corrected_area;
                                 
                               }
                              
                               

                      } //copying the data only to the initial neighbours.
                            
                      /***************************************************************************/         
                      //QUESTION!::: M.S.I.     
                      //what happens with the initial continuum contacs that now are not found becouse they are broken....
                      //should be assured that they become 0 when they break and this value keeps.
                      /************************************************************************************/      

                  } //for continuum initial neighbours
                        
              } // if(rCurrentProcessInfo[CONTACT_MESH_OPTION]==1)

              iContactForce++;

          }//for each neighbour

        if(rCurrentProcessInfo[INT_DUMMY_9] == 1) // if stress_strain_options ON 
        {
            if ( ( Representative_Volume <= 0.0 ))// && ( this->GetValue(SKIN_SPHERE) == 0 ) )
            {
        
              this->GetGeometry()(0)->GetSolutionStepValue(GROUP_ID) = 15;
              KRATOS_WATCH(this->Id())
              KRATOS_WATCH("Negatiu volume")
              KRATOS_WATCH(rCurrentProcessInfo[TIME_STEPS])
              
              
            }
                
            else
            {
                  for (int i=0; i<3; i++)
                  {
              
                          for (int j=0; j<3; j++)
                          {
                            //KRATOS_WATCH(Representative_Volume)
                            //KRATOS_WATCH(mStressTensor[i][j])
                                mStressTensor[i][j] = (1/Representative_Volume)*0.5*(mStressTensor [i][j] + mStressTensor[j][i]);  //THIS WAY THE TENSOR BECOMES SYMMETRIC
                          }
                  }       
            
            }
            
        }//if stress_strain_options
        

         
         KRATOS_CATCH("")
      }//ComputeParticleContactForce

      void SphericParticle::ApplyLocalMomentsDamping(const ProcessInfo& rCurrentProcessInfo )
      {

          KRATOS_TRY
      
          array_1d<double, 3 > & RotaMoment       = this->GetGeometry()[0].GetSolutionStepValue(PARTICLE_MOMENT);
          double RotaDampRatio                    = this->GetGeometry()[0].GetSolutionStepValue(PARTICLE_ROTATION_DAMP_RATIO);

          // LOCAL DAMPING OPTION FOR THE UNBALANCED FORCES (IN GLOBAL CORDINATES).
        
          for (int iDof = 0; iDof < 3; iDof++)
          {
              if (this->GetGeometry()(0)->GetSolutionStepValue(ANGULAR_VELOCITY)[iDof] > 0.0)
              {
                  RotaMoment[iDof] = RotaMoment[iDof] - RotaDampRatio * fabs(RotaMoment[iDof]);
              }
              else
              {
                  RotaMoment[iDof] = RotaMoment[iDof] + RotaDampRatio * fabs(RotaMoment[iDof]);
              }
          }

          KRATOS_CATCH("")

      } //ApplyLocalMomentsDamping


      void SphericParticle::CharacteristicParticleFailureId(const ProcessInfo& rCurrentProcessInfo )
      {  

          KRATOS_TRY

          // SETTING THE CARACHTERISTIC FAILURE TYPE OF THE PARTICLE (THE DOMINANT ONE) FOR VISUALITZATION.

            /* ContactFailureId = 2 (partially detached) it's never considerer in a contact, but it can be a desired output on the visualitzation.
             *
             * If the particle has no neighbour, but it is not from ContinuumGroup = 0 its becouse it has lost all the neighbours,
             * we are not interested in its global failure type but in its neighbours failure type. See: //// DETECTING A COHESIVE PARTICLE THAT HAS BEEN COMPLETELY DETACHED. in AfterForceCalculation.
             *
             * If one of the failure types completely dominates over the others this is the FailureId choosen for the particle
             *
             * If its partially detached by some type(s) of failure but the dominant contact state is "still atached", so failureId =0, we chose FailureId = 2.
             *
             * If the particle is not simulating the continuum, the failure type is 1 (default set from the particle whose ContinuumGroup is 0); It won't be printed as output.
             *
             */

            /*
             *   *mpFailureId values:
             *      0 := all neighbours attached
             *      1 := General detachment (all neighbours lost or no initial continuum case)
             *      2 := Partially detached and no dominance
             *      3 := tensile dominating on the contact detachment
             *      4 := shear dominating on the contact detachment
             */


          int tempType[5] = {0,0,0,0,0};

          if (*mpFailureId != 1)  // for *mpFailureId == 1 there's no failure to represent, the particle is not a continuum-simulating one or has been completelly detached already.
          {
              vector<int>& r_initial_neighbours_id          = this->GetValue(INI_NEIGHBOURS_IDS);

              int ini_neighbours_size = r_initial_neighbours_id.size();
     
              for(int index = 0; index < ini_neighbours_size; index++)    //loop over neighbours, just to run the index.
              {    
                  if( this->GetValue(PARTICLE_INITIAL_FAILURE_ID)[index] == 0)
                  {
                      tempType[0]++;
                  }
                  else if( this->GetValue(PARTICLE_INITIAL_FAILURE_ID)[index] == 1)
                  {
                      tempType[1]++;
                  }

                  // mContactFailureId == 2 intentionally skipped!

                  else if( this->GetValue(PARTICLE_INITIAL_FAILURE_ID)[index] == 3)
                  {
                      tempType[3]++;
                  }
                  else if( this->GetValue(PARTICLE_INITIAL_FAILURE_ID)[index] == 4)
                  {
                      tempType[4]++;
                  }
              }

              if ( tempType[0] == 0)  //no neighbour is attached
              {
                  *mpFailureId = 1;
              }   // no one neighbour is attached (but maybe still contacting).
              else if( (tempType[3] > tempType[4]) ) //some neighbour attached but failure 3 dominates over 4.
              {
                  *mpFailureId = 3;
              }
              else if( (tempType[4] > tempType[3]) ) // the same but 4 dominates over 3.
              {
                  *mpFailureId = 4;
              }
              else if ( (tempType[4] > 0) || (tempType[3] > 0) ) // no 3 neither 4 dominates but one of them may exist.
              {
                  *mpFailureId = 2;  // Partially detached / mix case.
              }
              else
              {
                  *mpFailureId = 0;  // last option: no one detached.
              }

          }// if (*mpFailureId != 1)

          KRATOS_CATCH("")

      } //CharacteristicParticleFailureId


      void SphericParticle::ComputeParticleRotationSpring(const ProcessInfo& rCurrentProcessInfo) // shan de corregir areees etc etc etc...
      {
        //double dt                           = rCurrentProcessInfo[DELTA_TIME]; //C.F.: neew
        /*
                    c=objecte_contacte(particula,vei)

            força=(c.RADI)*3;  //M: idea: to create a class contact, create objects of contacts with all the paramaters. easier...
                                /no puc amb MPI oi? pk hauria de passar punters...
          */
          double magic_factor =  rCurrentProcessInfo[DEM_MAGIC_FACTOR];
          
          double Tension        = this->GetGeometry()[0].GetSolutionStepValue(PARTICLE_TENSION);
          double Cohesion       = this->GetGeometry()[0].GetSolutionStepValue(PARTICLE_COHESION);
          double young          = this->GetGeometry()[0].GetSolutionStepValue(YOUNG_MODULUS);
          double poisson        = this->GetGeometry()[0].GetSolutionStepValue(POISSON_RATIO);
          double radius         = this->GetGeometry()[0].GetSolutionStepValue(RADIUS);
          double inertia        = this->GetGeometry()[0].GetSolutionStepValue(PARTICLE_INERTIA);

          array_1d<double, 3 > & mRota_Moment = GetGeometry()(0)->FastGetSolutionStepValue(PARTICLE_MOMENT);

          ParticleWeakVectorType& rE             = this->GetValue(NEIGHBOUR_ELEMENTS);

          Vector & mRotaSpringFailureType  = this->GetValue(PARTICLE_ROTATE_SPRING_FAILURE_TYPE);

          size_t iContactForce = 0;

          for(ParticleWeakIteratorType ineighbour = rE.begin(); ineighbour != rE.end(); ineighbour++)
          {

              //if(mIfInitalContact[iContactForce] == 1 && mRotaSpringFailureType[iContactForce] == 0) ///M.S:NEWWWW, IF THE SPRING BRAKES... NO MORE CONTRIBUION.
              //if( mRotaSpringFailureType[iContactForce] == 0) //M.S: CAL FICAR A INITIALIZE QUE SIGUI 1 I DESPRES INITIAL CONTACTS POSAR 0 SI NECESITEN, IGUAL QUE FAILURE NORMAL.
              //mmm.. what about the other failure types? if a contact is broken due to shear or tensile, it cant be a bending
              {
                
                
                  array_1d<double, 3 > & mRotaSpringMoment  = this->GetValue(PARTICLE_ROTATE_SPRING_MOMENT)[ iContactForce ];

                  double other_radius    = ineighbour->GetGeometry()(0)->FastGetSolutionStepValue(RADIUS);
                  double other_young     = ineighbour->GetGeometry()[0].GetSolutionStepValue(YOUNG_MODULUS);
                  double other_poisson   = ineighbour->GetGeometry()[0].GetSolutionStepValue(POISSON_RATIO);
                  double other_tension   = ineighbour->GetGeometry()[0].GetSolutionStepValue(PARTICLE_TENSION);
                  double other_cohesion  = ineighbour->GetGeometry()[0].GetSolutionStepValue(PARTICLE_COHESION);
                  double other_inertia   = ineighbour->GetGeometry()(0)->FastGetSolutionStepValue(PARTICLE_INERTIA);

                  Tension  = (Tension  + other_tension ) * 0.5;
                  Cohesion = (Cohesion + other_cohesion) * 0.5;

                  double equiv_radius     = (radius + other_radius) * 0.5 ;
                  double equiv_area       = M_PI * equiv_radius * equiv_radius;
                  double equiv_poisson    = (poisson + other_poisson) * 0.5 ;
                  double equiv_young      = (young  + other_young)  * 0.5;

                  double kn               = magic_factor*equiv_young * equiv_area / (2.0 * equiv_radius);
                  double ks               = kn / (2.0 * (1.0 + equiv_poisson));

                  array_1d<double,3> other_to_me_vect = GetGeometry()(0)->Coordinates() - ineighbour->GetGeometry()(0)->Coordinates();

                  /////Cfeng: Forming the Local Contact Coordinate system
                  double NormalDir[3]           = {0.0};
                  double LocalCoordSystem[3][3] = {{0.0}, {0.0}, {0.0}};
                  NormalDir[0] = other_to_me_vect[0];
                  NormalDir[1] = other_to_me_vect[1];
                  NormalDir[2] = other_to_me_vect[2];
                  GeometryFunctions::ComputeContactLocalCoordSystem(NormalDir, LocalCoordSystem);

                  double LocalRotaSpringMoment[3]     = {0.0};
                  double GlobalRotaSpringMoment[3]    = {0.0};
                  double GlobalRotaSpringMomentOld[3] = {0.0};

                  double DeltRotaDisp[3] = {0.0};
                  //double DeltRotaDisp2[3] = {0.0};

                  double LocalDeltRotaDisp[3] = {0.0};

                  double TargetDeltRotaDist[3] = {0.0};
                  double NeighbourDeltRotaDist[3] = {0.0};
                
                  for (int i=0;i<3;i++)
                  {
                      TargetDeltRotaDist[i] = this->GetGeometry()(0)->FastGetSolutionStepValue(DELTA_ROTA_DISPLACEMENT)[i];
                      NeighbourDeltRotaDist[i] = ineighbour->GetGeometry()(0)->FastGetSolutionStepValue(DELTA_ROTA_DISPLACEMENT)[i];
                      DeltRotaDisp[i] =  - ( TargetDeltRotaDist[i] - NeighbourDeltRotaDist[i] );
                  }
              
                  GeometryFunctions::VectorGlobal2Local(LocalCoordSystem, DeltRotaDisp, LocalDeltRotaDisp);

                  GlobalRotaSpringMomentOld[0] = mRotaSpringMoment[ 0 ];

                  GlobalRotaSpringMomentOld[1] = mRotaSpringMoment[ 1 ];

                  GlobalRotaSpringMomentOld[2] = mRotaSpringMoment[ 2 ];

                  GeometryFunctions::VectorGlobal2Local(LocalCoordSystem, GlobalRotaSpringMomentOld, LocalRotaSpringMoment);


                  double Inertia_I = (inertia + other_inertia) * 0.5;

                  double Inertia_J = Inertia_I * 2.0;


                  LocalRotaSpringMoment[0] +=  - Inertia_I * LocalDeltRotaDisp[0] * kn / equiv_area;

                  LocalRotaSpringMoment[1] +=  - Inertia_I * LocalDeltRotaDisp[1] * kn / equiv_area;

                  LocalRotaSpringMoment[2] +=  - Inertia_J * LocalDeltRotaDisp[2] * ks / equiv_area;


                  ////Judge if the rotate spring is broken or not
                  double GlobalContactForce[3]  = {0.0};
                  double LocalContactForce [3]  = {0.0};

                  GlobalContactForce[0] = this->GetValue(PARTICLE_CONTACT_FORCES)[ iContactForce ][ 0 ];
                  GlobalContactForce[1] = this->GetValue(PARTICLE_CONTACT_FORCES)[ iContactForce ][ 1 ];
                  GlobalContactForce[2] = this->GetValue(PARTICLE_CONTACT_FORCES)[ iContactForce ][ 2 ]; //GlobalContactForce[2] = mContactForces[3 * iContactForce  + 2 ];
                  GeometryFunctions::VectorGlobal2Local(LocalCoordSystem, GlobalContactForce, LocalContactForce);

                  double ForceN  = LocalContactForce[2];
                  double ForceS  = sqrt( LocalContactForce[0] * LocalContactForce[0] + LocalContactForce[1] * LocalContactForce[1]);
                  double MomentS = sqrt(LocalRotaSpringMoment[0] * LocalRotaSpringMoment[0] + LocalRotaSpringMoment[1] * LocalRotaSpringMoment[1]);
                  double MomentN = LocalRotaSpringMoment[2];

                  //////bending stress and axial stress add together, use edge of the bar will failure first
                  double TensiMax = -ForceN / equiv_area + MomentS        / Inertia_I * equiv_radius;
                  double ShearMax = ForceS  / equiv_area + fabs(MomentN)  / Inertia_J * equiv_radius;

                  if(TensiMax > Tension || ShearMax > Cohesion)
                  {
                      mRotaSpringFailureType[iContactForce] = 1;

                      LocalRotaSpringMoment[0] = 0.0;
                      LocalRotaSpringMoment[1] = 0.0;
                      LocalRotaSpringMoment[2] = 0.0;
                  }

                  GeometryFunctions::VectorLocal2Global(LocalCoordSystem, LocalRotaSpringMoment, GlobalRotaSpringMoment);

                  mRotaSpringMoment[ 0 ] = GlobalRotaSpringMoment[0];
                  mRotaSpringMoment[ 1 ] = GlobalRotaSpringMoment[1];
                  mRotaSpringMoment[ 2 ] = GlobalRotaSpringMoment[2];

                  ////feedback, contact moment----induce by rotation spring
                  mRota_Moment[0] -= GlobalRotaSpringMoment[0];
                  mRota_Moment[1] -= GlobalRotaSpringMoment[1];
                  mRota_Moment[2] -= GlobalRotaSpringMoment[2];
              }

              iContactForce++;
          }
      }//ComputeParticleRotationSpring

     
     void SphericParticle::SymmetrizeTensor(const ProcessInfo& rCurrentProcessInfo)
     {

        for (int i=0; i<3; i++)
        {
            for (int j=0; j<3; j++)
            {
                mSymmStressTensor[i][j] = 0.5*(mStressTensor[i][j]+mStressTensor[j][i]);
                
            }

        }

       //composició 
          double radius                             = this->GetGeometry()(0)->GetSolutionStepValue(RADIUS);
          ParticleWeakVectorType& r_neighbours      = this->GetValue(NEIGHBOUR_ELEMENTS);
          double poisson                            = this->GetGeometry()[0].GetSolutionStepValue(POISSON_RATIO);
          
          for(ParticleWeakIteratorType neighbour_iterator = r_neighbours.begin();
              neighbour_iterator != r_neighbours.end(); neighbour_iterator++)
          {            
              // GETTING NEIGHBOUR PROPERTIES

              //searching for the area
              double other_radius     = neighbour_iterator->GetGeometry()(0)->GetSolutionStepValue(RADIUS);
              double other_poisson    = neighbour_iterator->GetGeometry()[0].GetSolutionStepValue(POISSON_RATIO);
              double equiv_radius     = 2*radius * other_radius / (radius + other_radius);
              int size_ini_cont_neigh = this->GetValue(CONTINUUM_INI_NEIGHBOURS_IDS).size();
              double equiv_area       = (0.25)*M_PI * equiv_radius * equiv_radius; // 0.25 is becouse we take only the half of the equivalent radius, corresponding to the case of one sphere with radius Requivalent and other = radius 0.
              double equiv_poisson    = 2* poisson * other_poisson / (poisson + other_poisson);
              //double equiv_young      = 2 * young * other_young / (young + other_young);
              //bool is_continuum       = false;
               double corrected_area = equiv_area;
               
               bool found = false;
               
               for (int index_area=0; index_area<size_ini_cont_neigh; index_area++)
               {

                      if ( this->GetValue(CONTINUUM_INI_NEIGHBOURS_IDS)[index_area] == int(neighbour_iterator->Id()) ) 
                      {
                         
                        //MPI_CARLOS_ en MPI no pot funcionar lo de symmetrize de moment

                              Element::Pointer lock_p_weak = (this->GetGeometry()[0].GetValue(NODE_TO_NEIGH_ELEMENT_POINTER)(index_area)).lock();

                              corrected_area = lock_p_weak->GetValue(MEAN_CONTACT_AREA);

                               array_1d<double,3> other_to_me_vect = this->GetGeometry()(0)->Coordinates() - neighbour_iterator->GetGeometry()(0)->Coordinates();
              
                              double NormalDir[3]           = {0.0};
                              double LocalCoordSystem[3][3] = {{0.0}, {0.0}, {0.0}};
                              NormalDir[0] = other_to_me_vect[0];  
                              NormalDir[1] = other_to_me_vect[1];
                              NormalDir[2] = other_to_me_vect[2];
                              
                              
                              double Auxiliar[3][2] = {{0.0}, {0.0}, {0.0}};
                              double stress_projection[2] = {0.0};
                              GeometryFunctions::ComputeContactLocalCoordSystem(NormalDir, LocalCoordSystem);  
                              
                                  
                              //assumció de que la suma de sigmaX i sigmaZ en un pla sempre dona el mateix valor.
                                            
                              //vector ortogonal 1 = LocalCoordSystem[0]
                              //vector ortogonal 2 = LocalCoordSystem[1]
                              
                              //NOTE:puc fer un clean si un valor es massa petit
                              
                              for (int i=0;i<2;i++)//only for 0 and 1, dos auxiliars
                              { 
                                          
                                  for (int j=0;j<3;j++)//for 0,1,2. Component dels auxiliars
                                  {
                                    for (int u=0;u<3;u++)
                                    {
                    
                                    
                                    Auxiliar[j][i] += mSymmStressTensor[j][u]*LocalCoordSystem[i][u];
                                  
                                    }

                                      for (int k=0;k<3;k++)//for 0,1,2.
                                      {
                                                      
                                          stress_projection[i] += Auxiliar[k][i]*LocalCoordSystem[i][k];
                                    
                                      } 
                                  }  
                                  
                              }
                              
                              double Normal_Contact_Contribution = -1.0*corrected_area*equiv_poisson*(stress_projection[0]+stress_projection[1]); 
                              
                              //storing the value in the bar and doing the mean
                            
                              if(this->Id() < neighbour_iterator->Id())
                              {
                                
                                lock_p_weak->GetValue(LOW_POISSON_FORCE) = Normal_Contact_Contribution; //crec que ja té la direcció cap a on ha danar el veí.
                                 
                                 if(rCurrentProcessInfo[TIME_STEPS]==500)
                                 {
                                 KRATOS_WATCH(lock_p_weak->GetValue(LOW_POISSON_FORCE))             
                                 }
                                
                                 
                              }
                              else
                              {
                                
                                lock_p_weak->GetValue(HIGH_POISSON_FORCE) = Normal_Contact_Contribution;
                                 if(rCurrentProcessInfo[TIME_STEPS]==500)
                                 {
                                 KRATOS_WATCH(lock_p_weak->GetValue(HIGH_POISSON_FORCE))   
                                 
                                 }            
                              }
                                              
                              found = true;
                              
                              break;
                              
             
                        
                      }// if ( this->GetValue(CONTINUUM_INI_NEIGHBOURS_IDS)[index_area] == int(neighbour_iterator->Id()) ) 
                                         
                  }//for every ini_neighbour      
              
              
              if(found == false)
              {
                
                KRATOS_WATCH("ERROR!!!!!!!!!!!!!   THIS IS ONLY VALID FOR KNOWN INITIAL NEIGHBOURS, NO NEW ONES")KRATOS_WATCH(rCurrentProcessInfo[TIME_STEPS])
                
              }
              
              //NOTE: this area is provisionally obtained diferent from each side of the contact. bars are not yet ready for mpi.
              
             
              
               
               //this->GetValue(POISSON_NORMAL_FORCE) += 0.5*Normal_Contact_Contribution*LocalCoordSystem[2]; //seria el mateix que other to me normalitzat oi?
               //this->GetValue(POISSON_NORMAL_FORCE) += 0.5*Normal_Contact_Contribution*LocalCoordSystem[2];
               
               //neighbour_iterator->GetValue(POISSON_NORMAL_FORCE) += 0.5*Normal_Contact_Contribution;
                
                
         
              
          } // for every neighbour    
     }
     
     
      void SphericParticle::DampMatrix(MatrixType& rDampMatrix, ProcessInfo& rCurrentProcessInfo){}

      void SphericParticle::GetDofList(DofsVectorType& ElementalDofList, ProcessInfo& CurrentProcessInfo)
      {
          ElementalDofList.resize( 0 );

          for ( unsigned int i = 0; i < GetGeometry().size(); i++ )
          {
              ElementalDofList.push_back( GetGeometry()[i].pGetDof( DISPLACEMENT_X ) );
              ElementalDofList.push_back( GetGeometry()[i].pGetDof( DISPLACEMENT_Y ) );

              if ( GetGeometry().WorkingSpaceDimension() == 3 )
              {
                  ElementalDofList.push_back( GetGeometry()[i].pGetDof( DISPLACEMENT_Z ) );
              }
          }
      }

      void SphericParticle::InitializeSolutionStep(ProcessInfo& rCurrentProcessInfo)
      {
        
        
          int case_opt                  = rCurrentProcessInfo[CASE_OPTION];        
          int neighbours_initialized    = rCurrentProcessInfo[NEIGH_INITIALIZED];

          if( (case_opt!=0) && (neighbours_initialized == 0) )
          {
              SetInitialContacts(case_opt, rCurrentProcessInfo);
          }

          array_1d<double,3>& rhs               = this->GetGeometry()[0].GetSolutionStepValue(RHS);//RHS forces, we reset to 0. and we calculate again.
          noalias(rhs)                          = ZeroVector(3);
          array_1d<double,3>& total_forces      = this->GetGeometry()[0].GetSolutionStepValue(TOTAL_FORCES);
          noalias(total_forces)                 = ZeroVector(3);
          array_1d<double,3>& damp_forces       = this->GetGeometry()[0].GetSolutionStepValue(DAMP_FORCES);
          noalias(damp_forces)  = ZeroVector(3);
          
          if(rCurrentProcessInfo[ROTATION_OPTION]==1)
          {
            array_1d<double, 3 > & mRota_Moment   = this->GetGeometry()[0].GetSolutionStepValue(PARTICLE_MOMENT);
            noalias(mRota_Moment) = ZeroVector(3);
          }
          
          double& Representative_Volume = this->GetGeometry()[0].GetSolutionStepValue(REPRESENTATIVE_VOLUME);
         
          Representative_Volume = 0.0;
          
          for (int i=0; i<3; i++)
          {
          
              for (int j=0; j<3; j++)
              {
                   mStressTensor[i][j] = 0.0;
                   mSymmStressTensor[i][j] = 0.0;
              }
          
           }
    
           
           //DEBUG MEDICIÓ
           /*
                  double& area_vertical_tapa = rCurrentProcessInfo[AREA_VERTICAL_TAPA];
                  double& area_vertical_centre = rCurrentProcessInfo[AREA_VERTICAL_CENTRE];
                  
                  if(rCurrentProcessInfo[TIME_STEPS]==2)
                      
                  {
                     
                     if( (this->GetGeometry()[0].GetSolutionStepValue(GROUP_ID)==1) || (this->GetGeometry()[0].GetSolutionStepValue(GROUP_ID)==5))
                            
                      {
                      
                       ParticleWeakVectorType& r_neighbours                = this->GetValue(NEIGHBOUR_ELEMENTS);
           
                         for(ParticleWeakIteratorType_ptr ineighbour = r_neighbours.ptr_begin();  //loop over the neighbours and store into a initial_neighbours vector.
                             ineighbour != r_neighbours.ptr_end(); ineighbour++)
                  
                           {
                             
                             if( (((*ineighbour).lock())->GetGeometry()[0].GetSolutionStepValue(GROUP_ID)!=1) && (((*ineighbour).lock())->GetGeometry()[0].GetSolutionStepValue(GROUP_ID)!=5) ) 
                   
                             {
                                 
                                 int size_ini_cont_neigh = this->GetValue(CONTINUUM_INI_NEIGHBOURS_IDS).size();
                                 
                                  for (int index_area=0; index_area<size_ini_cont_neigh; index_area++)
                               {

                                 if ( this->GetValue(CONTINUUM_INI_NEIGHBOURS_IDS)[index_area] == (int)( ((*ineighbour).lock())->Id() ) )
                                {     
                                     Element::Pointer lock_p_weak = (this->GetGeometry()[0].GetValue(NODE_TO_NEIGH_ELEMENT_POINTER)(index_area)).lock();
                                                 double corrected_area = lock_p_weak->GetValue(MEAN_CONTACT_AREA);

                                                  array_1d<double,3> other_to_me_vect = this->GetGeometry()(0)->Coordinates() - ((*ineighbour).lock())->GetGeometry()(0)->Coordinates();
                                      array_1d<double,3> normal_vector_on_contact =  -1 * other_to_me_vect; //outwards     
                                      double Dummy_Dummy = 0.0;
                                      GeometryFunctions::norm(normal_vector_on_contact,Dummy_Dummy); // Normalize to unitary module
                                                 
                                                   if ((this->GetGeometry()[0].GetSolutionStepValue(GROUP_ID)==1)){
                                                 area_vertical_tapa += corrected_area*fabs(normal_vector_on_contact[1]);
                                                   
                                                  }
                                                  else if ( this->GetGeometry()[0].GetSolutionStepValue(GROUP_ID)==5 ) 
                                                   {
                                                      area_vertical_centre += 0.5*corrected_area*fabs(normal_vector_on_contact[1]);
                                                  }
                                                  
                                                  break;
                                       }// if ( this->GetValue(CONTINUUM_INI_NEIGHBOURS_IDS)[index_area] == int(neighbour_iterator->Id()) ) 
                      
                      
                                  }//for every neighbour      
                                 
                                 
                                 
                             }
                             
                         }
                       
                     }
                  }
           */
               
        
      }
        
      void SphericParticle::FinalizeSolutionStep(ProcessInfo& rCurrentProcessInfo) 
      {
          
    
          //this->GetGeometry()[0].FastGetSolutionStepValue(EXPORT_PARTICLE_FAILURE_ID) = double(this->GetValue(PARTICLE_FAILURE_ID)); //temporarily unused
          if(rCurrentProcessInfo[INT_DUMMY_3]==1)
          {
            this->GetGeometry()[0].FastGetSolutionStepValue(EXPORT_ID) = double(this->Id());
          }
          if(rCurrentProcessInfo[INT_DUMMY_4]==1)
          {
            
            this->GetGeometry()[0].FastGetSolutionStepValue(EXPORT_SKIN_SPHERE) = double(this->GetValue(SKIN_SPHERE));  
            
          }
          
          this->GetGeometry()[0].FastGetSolutionStepValue(NUM_OF_NEIGH) = this->GetValue(NEIGHBOUR_ELEMENTS).size();
         
          if( rCurrentProcessInfo[CONTACT_MESH_OPTION]==1 && rCurrentProcessInfo[INT_DUMMY_9] )
          {
      
          this->GetGeometry()[0].FastGetSolutionStepValue(DEM_STRESS_XX) =  mStressTensor[0][0];
          this->GetGeometry()[0].FastGetSolutionStepValue(DEM_STRESS_XY) =  mStressTensor[0][1];
          this->GetGeometry()[0].FastGetSolutionStepValue(DEM_STRESS_XZ) =  mStressTensor[0][2];
          this->GetGeometry()[0].FastGetSolutionStepValue(DEM_STRESS_YX) =  mStressTensor[1][0];
          this->GetGeometry()[0].FastGetSolutionStepValue(DEM_STRESS_YY) =  mStressTensor[1][1];
          this->GetGeometry()[0].FastGetSolutionStepValue(DEM_STRESS_YZ) =  mStressTensor[1][2];
          this->GetGeometry()[0].FastGetSolutionStepValue(DEM_STRESS_ZX) =  mStressTensor[2][0];
          this->GetGeometry()[0].FastGetSolutionStepValue(DEM_STRESS_ZY) =  mStressTensor[2][1];
          this->GetGeometry()[0].FastGetSolutionStepValue(DEM_STRESS_ZZ) =  mStressTensor[2][2];
          
          }
          
           if(rCurrentProcessInfo[INT_DUMMY_8]==1)
          {
            
            double X = this->GetGeometry()[0].GetSolutionStepValue(DISPLACEMENT_X);
            double Z = this->GetGeometry()[0].GetSolutionStepValue(DISPLACEMENT_Z);
            
            this->GetGeometry()[0].GetSolutionStepValue(RADIAL_DISPLACEMENT) = sqrt(X*X+Z*Z);
            
          }
          
          

           // the elemental variable is copied to a nodal variable in order to export the results onto GiD Post. Also a casting to double is necessary for GiD interpretation.
      }

      void SphericParticle::Calculate(const Variable<double>& rVariable, double& Output, const ProcessInfo& rCurrentProcessInfo)
      {
          //KRATOS_WATCH(rVariable)

          if (rVariable == DELTA_TIME)
          {
              double mass  = mRealMass;
              double coeff = rCurrentProcessInfo[NODAL_MASS_COEFF];

              if (coeff>1.0)
              {
                  KRATOS_ERROR(std::runtime_error,"The coefficient assigned for vitual mass is larger than one, virtual_mass_coeff= ",coeff)
              }

              else if ((coeff==1.0) && (rCurrentProcessInfo[VIRTUAL_MASS_OPTION]))
              {
                  Output = 9.0E09;
              }

              else
              {
                  if (rCurrentProcessInfo[VIRTUAL_MASS_OPTION])
                  {
                      mass = mass/(1-coeff);
                  }

                  double E = this->GetGeometry()(0)->FastGetSolutionStepValue(YOUNG_MODULUS);
                  double K = E * M_PI * this->GetGeometry()(0)->FastGetSolutionStepValue(RADIUS); //M. Error, should be the same that the local definition.

                  if (rCurrentProcessInfo[GLOBAL_VARIABLES_OPTION]==1)
                      K = rCurrentProcessInfo[GLOBAL_KN];

                  Output = 0.34 * sqrt( mass / K);

                  if(rCurrentProcessInfo[ROTATION_OPTION] == 1)
                  {
                      Output = Output * 0.5; //factor for critical time step when rotation is allowed.
                  }
              }
          }//CRITICAL DELTA CALCULATION

          if (rVariable == PARTICLE_ROTATION_DAMP_RATIO)
          {
              ApplyLocalMomentsDamping( rCurrentProcessInfo );
          } //DAMPING

          if (rVariable == DEM_STRESS_XX)  //operations with the stress_strain tensors
          {
            
              SymmetrizeTensor( rCurrentProcessInfo );
             // SymmetrizeTensor( rCurrentProcessInfo );
             // SymmetrizeTensor( rCurrentProcessInfo );
            
          } //EULER_ANGLES

          
          
          

      }//calculate

      void SphericParticle::Calculate(const Variable<array_1d<double, 3 > >& rVariable, array_1d<double, 3 > & Output, const ProcessInfo& rCurrentProcessInfo){}
      void SphericParticle::Calculate(const Variable<Vector >& rVariable, Vector& Output, const ProcessInfo& rCurrentProcessInfo){}
      void SphericParticle::Calculate(const Variable<Matrix >& rVariable, Matrix& Output, const ProcessInfo& rCurrentProcessInfo){}

}  // namespace Kratos.
