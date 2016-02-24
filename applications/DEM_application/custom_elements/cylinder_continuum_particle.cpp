//
// Author: Miquel Santasusana msantasusana@cimne.upc.edu
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

          mRadius                   = GetGeometry()[0].FastGetSolutionStepValue(RADIUS);
          double density            = GetDensity();

          double& mass              = GetGeometry()[0].FastGetSolutionStepValue(NODAL_MASS);
          mass                      = density * GetVolume();
          mRealMass                 = mass;

          if (this->Is(DEMFlags::HAS_ROTATION) ){
            double moment_of_inertia = 0.5 * mass * GetRadius() * GetRadius();   
            GetGeometry()[0].FastGetSolutionStepValue(PARTICLE_MOMENT_OF_INERTIA) = moment_of_inertia;
          }                                                                        

          CustomInitialize();

          KRATOS_CATCH( "" )
      }

      //**************************************************************************************************************************************************
      //**************************************************************************************************************************************************

           
      
        void CylinderContinuumParticle::ContactAreaWeighting() //MISMI 10: POOYAN this could be done by calculating on the bars. not looking at the neighbous of my neighbours.
        { 

            double alpha = 1.0;
            double sphere_perimeter = 2*KRATOS_M_PI * GetRadius();       
            double total_equiv_perimeter = 0.0;
            //int cont_ini_neighbours_size = mContinuumIniNeighbourElements.size();
            unsigned int continuous_initial_neighbours_size = mContinuumInitialNeighborsSize;
        
            for (unsigned int i = 0; i < continuous_initial_neighbours_size; i++) {
                //SphericParticle* ini_cont_neighbour_iterator = mContinuumIniNeighbourElements[i];
                SphericParticle* ini_cont_neighbour_iterator = mNeighbourElements[i];
                double other_radius     = ini_cont_neighbour_iterator->GetGeometry()[0].FastGetSolutionStepValue(RADIUS);
                double area = mContinuumConstitutiveLawArray[i]->CalculateContactArea(GetRadius(), other_radius, mContIniNeighArea); //This call fills the vector of areas only if the Constitutive Law wants.         
                total_equiv_perimeter += area;
            } //for every neighbour
      
            if (continuous_initial_neighbours_size >= 3) {
            
                if (!*mSkinSphere) {
                    AuxiliaryFunctions::CalculateAlphaFactor2D(continuous_initial_neighbours_size, sphere_perimeter, total_equiv_perimeter, alpha); 
                    for (unsigned int i = 0; i < mContIniNeighArea.size(); i++) {
                        mContIniNeighArea[i] = alpha*mContIniNeighArea[i];                      
                    } //for every neighbour
                }
                else { //skin sphere                              
                    for (unsigned int i = 0; i < mContIniNeighArea.size(); i++) {
                        alpha            = 1.30*(1.10266)*(sphere_perimeter/total_equiv_perimeter)*((double(continuous_initial_neighbours_size))/6); // 6 is mean coordination number.
                        mContIniNeighArea[i] = alpha*mContIniNeighArea[i];
                    }     //loop on cont neighs
                }           //skin particles.
            }               //if 3 neighbours or more.
        }                 //Contact Area Weighting
      
      


      void CylinderContinuumParticle::ComputeBallToBallContactForce(array_1d<double, 3>& rElasticForce,
              array_1d<double, 3 > & rContactForce,
              array_1d<double, 3>& rInitialRotaMoment,
              ProcessInfo& rCurrentProcessInfo,
              double dt,
              const bool multi_stage_RHS) {
          KRATOS_TRY

          const double dt_i = 1 / dt;
          const int time_steps = rCurrentProcessInfo[TIME_STEPS];

          int& search_control = rCurrentProcessInfo[SEARCH_CONTROL];
          vector<int>& search_control_vector = rCurrentProcessInfo[SEARCH_CONTROL_VECTOR];

          /* Initializations */

          const array_1d<double, 3>& vel         = this->GetGeometry()[0].FastGetSolutionStepValue(VELOCITY);
          const array_1d<double, 3>& delta_displ = this->GetGeometry()[0].FastGetSolutionStepValue(DELTA_DISPLACEMENT);
          const array_1d<double, 3>& ang_vel     = this->GetGeometry()[0].FastGetSolutionStepValue(ANGULAR_VELOCITY);
          const double moment_of_inertia         = this->GetGeometry()[0].FastGetSolutionStepValue(PARTICLE_MOMENT_OF_INERTIA);
          double RotaAcc[3] = {0.0};

          if (this->Is(DEMFlags::HAS_ROTATION)) {
              RotaAcc[0] = 0;
              RotaAcc[1] = 0;
              RotaAcc[2] = ang_vel[2] * dt_i;

              rInitialRotaMoment[0] = 0;
              rInitialRotaMoment[1] = 0;
              rInitialRotaMoment[2] = RotaAcc[2] * moment_of_inertia;
          }

          for (unsigned int i_neighbour_count = 0; i_neighbour_count < mNeighbourElements.size(); i_neighbour_count++) {

              SphericContinuumParticle* neighbour_iterator = dynamic_cast<SphericContinuumParticle*>(mNeighbourElements[i_neighbour_count]);

              unsigned int neighbour_iterator_id = neighbour_iterator->Id();

              array_1d<double, 3> other_to_me_vect = this->GetGeometry()[0].Coordinates() - neighbour_iterator->GetGeometry()[0].Coordinates();
              const double &other_radius = neighbour_iterator->GetRadius();
              double distance = sqrt(other_to_me_vect[0] * other_to_me_vect[0] +
                      other_to_me_vect[1] * other_to_me_vect[1]);
              double radius_sum = GetRadius() + other_radius;
              double initial_delta = mNeighbourDelta[i_neighbour_count]; //*
              double initial_dist = (radius_sum - initial_delta);
              double indentation = initial_dist - distance; //#1
              double myYoung = GetYoung();
              double myPoisson = GetPoisson();

              double kn_el;
              double kt_el;
              double DeltDisp[3] = {0.0};
              double RelVel[3] = {0.0};
              double LocalCoordSystem[3][3]         = {{0.0}, {0.0}, {0.0}};
              double OldLocalCoordSystem[3][3]      = {{0.0}, {0.0}, {0.0}};
              bool sliding = false;

              //const int mapping_new_cont = mMappingNewCont[i_neighbour_count];

              double contact_tau = 0.0;
              double contact_sigma = 0.0;
              double failure_criterion_state = 0.0;
              double acumulated_damage = 0.0;

              // Getting neighbor properties
              double other_young = neighbour_iterator->GetYoung();
              double other_poisson = neighbour_iterator->GetPoisson();
              double equiv_poisson;
              if ((myPoisson + other_poisson) != 0.0) {
                  equiv_poisson = 2.0 * myPoisson * other_poisson / (myPoisson + other_poisson);
              } else {
                  equiv_poisson = 0.0;
              }

              double equiv_young = 2.0 * myYoung * other_young / (myYoung + other_young);
              double calculation_area = 0.0;

              //if (mapping_new_cont != -1) {
              if (i_neighbour_count < mContinuumInitialNeighborsSize) {
                  mContinuumConstitutiveLawArray[i_neighbour_count]-> CalculateContactArea(GetRadius(), other_radius, calculation_area);
                  mContinuumConstitutiveLawArray[i_neighbour_count]-> CalculateElasticConstants(kn_el, kt_el, initial_dist, equiv_young, equiv_poisson, calculation_area);
              } else {
                  mDiscontinuumConstitutiveLaw -> CalculateContactArea(GetRadius(), other_radius, calculation_area);
                  mDiscontinuumConstitutiveLaw -> CalculateElasticConstants(kn_el, kt_el, initial_dist, equiv_young, equiv_poisson, calculation_area);
              }

              EvaluateDeltaDisplacement(DeltDisp, RelVel, LocalCoordSystem, OldLocalCoordSystem, other_to_me_vect, vel, delta_displ, neighbour_iterator, distance);

              if (this->Is(DEMFlags::HAS_ROTATION)) {
                  DisplacementDueToRotationMatrix(DeltDisp, RelVel, OldLocalCoordSystem, other_radius, dt, ang_vel, neighbour_iterator);
              }

              double LocalDeltDisp[3] = {0.0};
              double LocalElasticContactForce[3] = {0.0}; // 0: first tangential, // 1: second tangential, // 2: normal force
              double GlobalElasticContactForce[3] = {0.0};

              GlobalElasticContactForce[0] = mNeighbourElasticContactForces[i_neighbour_count][0];
              GlobalElasticContactForce[1] = mNeighbourElasticContactForces[i_neighbour_count][1];
              GlobalElasticContactForce[2] = 0;

              GeometryFunctions::VectorGlobal2Local(OldLocalCoordSystem, GlobalElasticContactForce, LocalElasticContactForce);
              //recover old local forces projected over new coordinates
              GeometryFunctions::VectorGlobal2Local(OldLocalCoordSystem, DeltDisp, LocalDeltDisp);
              
              int failure_id = mIniNeighbourFailureId[i_neighbour_count];
              
              /* Translational Forces */

              if (indentation > 0.0 || (mIniNeighbourFailureId[i_neighbour_count] == 0))       //#3
              {
                  //if (mapping_new_cont != -1) {                                             //Normal Forces
                  if (i_neighbour_count < mContinuumInitialNeighborsSize) {
                      mContinuumConstitutiveLawArray[i_neighbour_count]-> CalculateForces(
                              rCurrentProcessInfo,
                              LocalElasticContactForce,
                              LocalDeltDisp,
                              kn_el,
                              kt_el,
                              contact_sigma,
                              contact_tau,
                              failure_criterion_state,
                              equiv_young,
                              indentation,
                              calculation_area,
                              acumulated_damage,
                              this,
                              neighbour_iterator,
                              i_neighbour_count,
                              rCurrentProcessInfo[TIME_STEPS],
                              sliding,
                              search_control,
                              search_control_vector);
                  } else {
                      mDiscontinuumConstitutiveLaw -> CalculateForces(
                              rCurrentProcessInfo,
                              LocalElasticContactForce,
                              LocalDeltDisp,
                              kn_el,
                              kt_el,
                              indentation,
                              failure_criterion_state,
                              sliding,
                              this,
                              neighbour_iterator,
                              mIniNeighbourFailureId[i_neighbour_count]/*, mapping_new_cont*/);
                  }
              } // compression or cohesive contact

              //  Viscodamping (applied locally)
              double ViscoDampingLocalContactForce[3] = {0.0};
              double equiv_visco_damp_coeff_normal;
              double equiv_visco_damp_coeff_tangential;

              if (indentation > 0.0 || (mIniNeighbourFailureId[i_neighbour_count] == 0)) {

                  double LocalRelVel[3] = {0.0};
                  GeometryFunctions::VectorGlobal2Local(LocalCoordSystem, RelVel, LocalRelVel);
                  //if (mapping_new_cont != -1) {
                  if (i_neighbour_count < mContinuumInitialNeighborsSize) {
                      mContinuumConstitutiveLawArray[i_neighbour_count]->CalculateViscoDampingCoeff(equiv_visco_damp_coeff_normal,
                              equiv_visco_damp_coeff_tangential,
                              this,
                              neighbour_iterator,
                              kn_el,
                              kt_el);

                      mContinuumConstitutiveLawArray[i_neighbour_count]->CalculateViscoDamping(LocalRelVel,
                              ViscoDampingLocalContactForce,
                              indentation,
                              equiv_visco_damp_coeff_normal,
                              equiv_visco_damp_coeff_tangential,
                              sliding,
                              failure_id);
                  } else {
                      mDiscontinuumConstitutiveLaw -> CalculateViscoDampingCoeff(equiv_visco_damp_coeff_normal,
                              equiv_visco_damp_coeff_tangential,
                              this,
                              neighbour_iterator,
                              kn_el,
                              kt_el);

                      mDiscontinuumConstitutiveLaw -> CalculateViscoDamping(LocalRelVel,
                              ViscoDampingLocalContactForce,
                              indentation,
                              equiv_visco_damp_coeff_normal,
                              equiv_visco_damp_coeff_tangential,
                              sliding);
                  }
              }

              // Transforming to global forces and adding up
              double LocalContactForce[3] = {0.0};
              double GlobalContactForce[3] = {0.0};

              //if (rCurrentProcessInfo[STRESS_STRAIN_OPTION] && mapping_new_cont != -1) {AddPoissonContribution(equiv_poisson,
              if (rCurrentProcessInfo[STRESS_STRAIN_OPTION] && (i_neighbour_count < mContinuumInitialNeighborsSize)) {AddPoissonContribution(equiv_poisson,
                                                                                          LocalCoordSystem,
                                                                                          LocalElasticContactForce[2],
                                                                                          calculation_area);}

              AddUpForcesAndProject(OldLocalCoordSystem, LocalCoordSystem, LocalContactForce, LocalElasticContactForce, GlobalContactForce,
                      GlobalElasticContactForce, ViscoDampingLocalContactForce, 0.0, rElasticForce, rContactForce, i_neighbour_count);

              array_1d<double, 3> temp_force = ZeroVector(3);

              temp_force[0] = GlobalContactForce[0];
              temp_force[1] = GlobalContactForce[1];
              temp_force[2] = 0;

              if (this->Is(DEMFlags::HAS_ROTATION)) {
                  ComputeMoments(LocalElasticContactForce[2], temp_force, rInitialRotaMoment, LocalCoordSystem[2], neighbour_iterator, indentation);
              }

              if (rCurrentProcessInfo[CONTACT_MESH_OPTION] == 1 && (i_neighbour_count < mContinuumInitialNeighborsSize) && this->Id() < neighbour_iterator_id) {

                  CalculateOnContactElements(i_neighbour_count,
                                             LocalElasticContactForce,
                                             contact_sigma,
                                             contact_tau,
                                             failure_criterion_state,
                                             acumulated_damage,
                                             time_steps);
              }

              if (rCurrentProcessInfo[STRESS_STRAIN_OPTION] && (i_neighbour_count < mContinuumInitialNeighborsSize)) {
                  AddNeighbourContributionToStressTensor(GlobalElasticContactForce,
                                                          LocalCoordSystem[2],distance, radius_sum);}

              AddContributionToRepresentativeVolume(distance, radius_sum, calculation_area);



          } //for each neighbor

          KRATOS_CATCH("")

      } //ComputeBallToBallContactForce
      
      double CylinderContinuumParticle::GetVolume(){
          return KRATOS_M_PI * GetRadius() * GetRadius();
      }



//      void CylinderContinuumParticle::AddNeighbourContributionToStressTensor(double GlobalElasticContactForce[3],
//                                                                            array_1d<double,3> &other_to_me_vect,
//                                                                            const double &distance,
//                                                                            const double &radius_sum,
//                                                                            const double &calculation_area,
//                                                                            ParticleWeakIteratorType neighbour_iterator,
//                                                                            ProcessInfo& rCurrentProcessInfo,
//                                                                            double &rRepresentative_Volume){
//        double gap  = distance - radius_sum;
//        array_1d<double,3> normal_vector_on_contact =  -1 * other_to_me_vect; //outwards
//        double value = 0.0;
//        GeometryFunctions::normalize(normal_vector_on_contact,value); // Normalize to unitary module
//        array_1d<double,3> x_centroid      = (GetRadius() + 0.5*gap) * normal_vector_on_contact;
//        array_1d<double,3> surface_baricenter = x_centroid;
//        double result_product = GeometryFunctions::DotProduct(surface_baricenter,normal_vector_on_contact);
//        //Aproximation with error: surface_baricenter should be the baricenter of each surface, which can no be calculated because the surfaces are imaginary.
//        rRepresentative_Volume = rRepresentative_Volume + 0.5 * (result_product * calculation_area);
        
//        for (int i=0; i<3; i++){
//            for (int j=0; j<3; j++){
//                (*mStressTensor)(i,j) += (x_centroid[j]) * GlobalElasticContactForce[i]; //ref: Katalin Bagi 1995 Mean stress tensor
//            }
//        }
//      } //AddNeighbourContributionToStressTensor




      //**************************************************************************************************************************************************
      //**************************************************************************************************************************************************
//       void CylinderContinuumParticle::Calculate(const Variable<double>& rVariable, double& Output, const ProcessInfo& rCurrentProcessInfo){}
//       void CylinderContinuumParticle::Calculate(const Variable<array_1d<double, 3> >& rVariable, array_1d<double, 3>& Output, const ProcessInfo& rCurrentProcessInfo){}
//       void CylinderContinuumParticle::Calculate(const Variable<Vector >& rVariable, Vector& Output, const ProcessInfo& rCurrentProcessInfo){}
//       void CylinderContinuumParticle::Calculate(const Variable<Matrix >& rVariable, Matrix& Output, const ProcessInfo& rCurrentProcessInfo){}

      //**************************************************************************************************************************************************
      //**************************************************************************************************************************************************

}  // namespace Kratos.

