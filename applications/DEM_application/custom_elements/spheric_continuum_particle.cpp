//
// Authors: Miquel Santasusana msantasusana@cimne.upc.edu
//          Miguel Angel Celigueta maceli@cimne.upc.edu
//

// System includes
#include <string>
#include <iostream>
#include <iomanip> // to improve std::cout precision

// Project includes
#include "includes/define.h"
#include "spheric_particle.h"
#include "custom_utilities/GeometryFunctions.h"
#include "custom_utilities/AuxiliaryFunctions.h"
#include "DEM_application.h"
#include "utilities/openmp_utils.h"
#include "utilities/timer.h"

#define CUSTOMTIMER 0  // ACTIVATES AND DISABLES ::TIMER:::::

#ifdef CUSTOMTIMER
#define KRATOS_TIMER_START(t) Timer::Start(t);
#define KRATOS_TIMER_STOP(t) Timer::Stop(t);
#else
#define KRATOS_TIMER_START(t)
#define KRATOS_TIMER_STOP(t)
#endif

namespace Kratos {

    SphericContinuumParticle::SphericContinuumParticle():SphericParticle() {}

    SphericContinuumParticle::SphericContinuumParticle(IndexType NewId, GeometryType::Pointer pGeometry) : SphericParticle(NewId, pGeometry){}

    SphericContinuumParticle::SphericContinuumParticle(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
    : SphericParticle(NewId, pGeometry, pProperties){}

    SphericContinuumParticle::SphericContinuumParticle(IndexType NewId, NodesArrayType const& ThisNodes)
    : SphericParticle(NewId, ThisNodes){}

    Element::Pointer SphericContinuumParticle::Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const
    {
        return SphericParticle::Pointer(new SphericContinuumParticle(NewId, GetGeometry().Create(ThisNodes), pProperties));
    }

    /// Destructor

    SphericContinuumParticle::~SphericContinuumParticle() {}

    
    void SphericContinuumParticle::SearchSkinParticles(ProcessInfo& r_process_info) {
    
        array_1d<double, 3> other_to_me_vect_total = ZeroVector(3);
        
        for (unsigned int i = 0; i < mContinuumInitialNeighborsSize; i++) {
            
            if (mNeighbourElements[i] == NULL) continue;
            
            SphericContinuumParticle* neighbour_iterator = dynamic_cast<SphericContinuumParticle*>(mNeighbourElements[i]);
            
            array_1d<double, 3> other_to_me_vect;
            noalias(other_to_me_vect) = this->GetGeometry()[0].Coordinates() - neighbour_iterator->GetGeometry()[0].Coordinates();

            GeometryFunctions::normalize(other_to_me_vect);
                    
            DEM_ADD_SECOND_TO_FIRST(other_to_me_vect_total, other_to_me_vect)
        }
        
        double other_to_me_vect_total_modulus = DEM_MODULUS_3(other_to_me_vect_total);
        
        if ((other_to_me_vect_total_modulus < 0.85) || (mContinuumInitialNeighborsSize > 9)) {
            this->GetGeometry()[0].FastGetSolutionStepValue(SKIN_SPHERE) = 0.0;
        }
    } // SearchSkinParticles
    
    void SphericContinuumParticle::SetInitialSphereContacts(ProcessInfo& r_process_info) {
         
        std::vector<SphericContinuumParticle*> ContinuumInitialNeighborsElements;
        std::vector<SphericContinuumParticle*> DiscontinuumInitialNeighborsElements;
        std::vector<int> DiscontinuumInitialNeighborsIds;
        std::vector<double> DiscontinuumInitialNeighborsDeltas;
        mIniNeighbourFailureId.clear(); // We will have to build this vector, we still don't know its size, it applies only to continuum particles
        size_t continuum_ini_size    = 0;
        size_t discontinuum_ini_size = 0;
        unsigned int neighbours_size = mNeighbourElements.size();
        mIniNeighbourIds.resize(neighbours_size);
        mIniNeighbourDelta.resize(neighbours_size);

        for (unsigned int i = 0; i < mNeighbourElements.size(); i++) {
            
            SphericContinuumParticle* neighbour_iterator = dynamic_cast<SphericContinuumParticle*>(mNeighbourElements[i]);
            array_1d<double, 3> other_to_me_vect;
            noalias(other_to_me_vect) = this->GetGeometry()[0].Coordinates() - neighbour_iterator->GetGeometry()[0].Coordinates();

            double distance = DEM_MODULUS_3(other_to_me_vect);
            double radius_sum = GetRadius() + neighbour_iterator->GetRadius();
            double initial_delta = radius_sum - distance;
            int r_other_continuum_group = neighbour_iterator->mContinuumGroup; // finding out neighbor's Continuum Group Id
            
            if ((r_other_continuum_group  == this->mContinuumGroup) && (this->mContinuumGroup != 0)) {
                
                mIniNeighbourIds[continuum_ini_size]   = neighbour_iterator->Id();
                mIniNeighbourDelta[continuum_ini_size] = initial_delta;
                mIniNeighbourFailureId.push_back(0);
                ContinuumInitialNeighborsElements.push_back(neighbour_iterator);
                continuum_ini_size++;

            } else {
                
                DiscontinuumInitialNeighborsIds.push_back(neighbour_iterator->Id());
                DiscontinuumInitialNeighborsDeltas.push_back(initial_delta);
                DiscontinuumInitialNeighborsElements.push_back(neighbour_iterator);
                discontinuum_ini_size++;
            }            
        }

        mContinuumInitialNeighborsSize = continuum_ini_size;
        mInitialNeighborsSize = neighbours_size;
        
        for (unsigned int j = 0; j < continuum_ini_size; j++) {
            mNeighbourElements[j] = ContinuumInitialNeighborsElements[j];
        }
            
        for (unsigned int k = 0; k < discontinuum_ini_size; k++) {
            
            mIniNeighbourIds[continuum_ini_size + k]   = DiscontinuumInitialNeighborsIds[k];
            mIniNeighbourDelta[continuum_ini_size + k] = DiscontinuumInitialNeighborsDeltas[k];
            mNeighbourElements[continuum_ini_size + k] = DiscontinuumInitialNeighborsElements[k];
        }
    }//SetInitialSphereContacts    
    
    void SphericContinuumParticle::CreateContinuumConstitutiveLaws() {
        
        unsigned int continuous_neighbor_size = mContinuumInitialNeighborsSize;
        mContinuumConstitutiveLawArray.resize(continuous_neighbor_size);

        for (unsigned int i = 0; i < continuous_neighbor_size; i++) {
            DEMContinuumConstitutiveLaw::Pointer NewContinuumConstitutiveLaw = GetProperties()[DEM_CONTINUUM_CONSTITUTIVE_LAW_POINTER]-> Clone();
            mContinuumConstitutiveLawArray[i] = NewContinuumConstitutiveLaw;
            mContinuumConstitutiveLawArray[i]->Initialize();
        }
    }
    
    void SphericContinuumParticle::SetInitialFemContacts() {
        
        std::vector<DEMWall*>& rFemNeighbours = this->mNeighbourRigidFaces;

        unsigned int fem_neighbours_size = rFemNeighbours.size();

        mFemIniNeighbourIds.resize(fem_neighbours_size);
        mFemIniNeighbourDelta.resize(fem_neighbours_size);
           
        for (unsigned int i = 0; i < rFemNeighbours.size(); i++) {
            
            double LocalCoordSystem[3][3]            = {{0.0}, {0.0}, {0.0}};
            array_1d<double, 3> wall_delta_disp_at_contact_point = ZeroVector(3);
            array_1d<double, 3> wall_velocity_at_contact_point = ZeroVector(3);    
            double DistPToB = 0.0;
            int ContactType = -1;
            array_1d<double, 4> Weight = ZeroVector(4);                
        
            ComputeConditionRelativeData(rFemNeighbours[i], LocalCoordSystem, DistPToB, Weight, wall_delta_disp_at_contact_point, wall_velocity_at_contact_point, ContactType);
          
            double initial_delta = -(DistPToB - GetRadius());

            mFemIniNeighbourIds[i] = rFemNeighbours[i]->Id();
            mFemIniNeighbourDelta[i] = initial_delta;
        }
    }//SetInitialFemContacts              
    
    void SphericContinuumParticle::ContactAreaWeighting() //MISMI 10: POOYAN this could be done by calculating on the bars. not looking at the neighbors of my neighbors.
    {
        double alpha = 1.0;
        double external_sphere_area = 4 * KRATOS_M_PI * GetRadius()*GetRadius();
        double total_equiv_area = 0.0;

        int cont_ini_neighbours_size = mContinuumInitialNeighborsSize;

        for (int i = 0; i < cont_ini_neighbours_size; i++) {
            SphericParticle* ini_cont_neighbour_iterator = mNeighbourElements[i];
            double other_radius = ini_cont_neighbour_iterator->GetRadius();
            double area = mContinuumConstitutiveLawArray[i]->CalculateContactArea(GetRadius(), other_radius, mContIniNeighArea); //This call fills the vector of areas only if the Constitutive Law wants.
            total_equiv_area += area;
        } 

        if (cont_ini_neighbours_size >= 4) { //more than 3 neighbors
            if (!*mSkinSphere) {
                AuxiliaryFunctions::CalculateAlphaFactor3D(cont_ini_neighbours_size, external_sphere_area, total_equiv_area, alpha);
                for (unsigned int i = 0; i < mContIniNeighArea.size(); i++) {
                    mContIniNeighArea[i] = alpha * mContIniNeighArea[i];                   
                } //for every neighbor
            }

            else {//skin sphere             
                for (unsigned int i = 0; i < mContIniNeighArea.size(); i++) {
                    alpha = 1.00 * (1.40727)*(external_sphere_area / total_equiv_area)*((double(cont_ini_neighbours_size)) / 11.0);
                    mContIniNeighArea[i] = alpha * mContIniNeighArea[i];
                } //for every neighbor
            }
        }//if more than 3 neighbors
    } 
    
    void SphericContinuumParticle::ComputeBallToBallContactForce(array_1d<double, 3>& rElasticForce,
                                                                 array_1d<double, 3>& rContactForce,
                                                                 array_1d<double, 3>& rInitialRotaMoment,
                                                                 ProcessInfo& r_process_info,
                                                                 const double dt,
                                                                 const bool multi_stage_RHS) {
        KRATOS_TRY

        const int time_steps = r_process_info[TIME_STEPS];
        const int& search_control = r_process_info[SEARCH_CONTROL];
        vector<int>& search_control_vector = r_process_info[SEARCH_CONTROL_VECTOR];

        const array_1d<double, 3>& vel         = this->GetGeometry()[0].FastGetSolutionStepValue(VELOCITY);
        const array_1d<double, 3>& delta_displ = this->GetGeometry()[0].FastGetSolutionStepValue(DELTA_DISPLACEMENT);
        const array_1d<double, 3>& ang_vel     = this->GetGeometry()[0].FastGetSolutionStepValue(ANGULAR_VELOCITY);
                        
        for (unsigned int i = 0; i < mNeighbourElements.size(); i++) {
            
            if (mNeighbourElements[i] == NULL) continue;
            if (this->Is(NEW_ENTITY) && mNeighbourElements[i]->Is(NEW_ENTITY)) continue;  
                                
            SphericContinuumParticle* neighbour_iterator = dynamic_cast<SphericContinuumParticle*>(mNeighbourElements[i]);
            
            unsigned int neighbour_iterator_id = neighbour_iterator->Id();
            
            array_1d<double, 3> other_to_me_vect;
            noalias(other_to_me_vect) = this->GetGeometry()[0].Coordinates() - neighbour_iterator->GetGeometry()[0].Coordinates();
                            
            const double& other_radius = neighbour_iterator->GetRadius();
 
            double distance = DEM_MODULUS_3(other_to_me_vect);
            double radius_sum = GetRadius() + other_radius;
                        
            double initial_delta = GetInitialDelta(i);
            
            double initial_dist = radius_sum - initial_delta;
            double indentation = initial_dist - distance;
            double myYoung = GetYoung();
            double myPoisson = GetPoisson();

            double kn_el;
            double kt_el;            
            double DeltDisp[3] = {0.0};
            double RelVel[3] = {0.0};
            double LocalCoordSystem[3][3]    = {{0.0}, {0.0}, {0.0}};
            double OldLocalCoordSystem[3][3] = {{0.0}, {0.0}, {0.0}};
            bool sliding = false;

            double contact_tau = 0.0;
            double contact_sigma = 0.0;
            double failure_criterion_state = 0.0;
            double acumulated_damage = 0.0;

            // Getting neighbour properties
            double other_young = neighbour_iterator->GetYoung();
            double other_poisson = neighbour_iterator->GetPoisson();
            double equiv_poisson;
            if ((myPoisson + other_poisson) != 0.0) { equiv_poisson = 2.0 * myPoisson * other_poisson / (myPoisson + other_poisson); }
            else { equiv_poisson = 0.0; }

            double equiv_young = 2.0 * myYoung * other_young / (myYoung + other_young);
            double calculation_area = 0.0;
            
            if (i < mContinuumInitialNeighborsSize) {
                mContinuumConstitutiveLawArray[i]->GetContactArea(GetRadius(), other_radius, mContIniNeighArea, i, calculation_area);
                mContinuumConstitutiveLawArray[i]->CalculateElasticConstants(kn_el, kt_el, initial_dist, equiv_young, equiv_poisson, calculation_area);
            } 
             
            EvaluateDeltaDisplacement(DeltDisp, RelVel, LocalCoordSystem, OldLocalCoordSystem, other_to_me_vect, vel, delta_displ, neighbour_iterator, distance);

            if (this->Is(DEMFlags::HAS_ROTATION)) {
                DisplacementDueToRotationMatrix(DeltDisp, RelVel, OldLocalCoordSystem, other_radius, dt, ang_vel, neighbour_iterator);
            }

            double LocalDeltDisp[3] = {0.0};
            double LocalElasticContactForce[3] = {0.0}; // 0: first tangential, // 1: second tangential, // 2: normal force
            double GlobalElasticContactForce[3] = {0.0};
            double OldLocalElasticContactForce[3] = {0.0};
            
            GeometryFunctions::VectorGlobal2Local(LocalCoordSystem, DeltDisp, LocalDeltDisp);
            
            RotateOldContactForces(OldLocalCoordSystem, LocalCoordSystem, mNeighbourElasticContactForces[i]);
            
            // Here we recover the old local forces projected in the new coordinates in the way they were in the old ones; Now they will be increased if necessary
            GeometryFunctions::VectorGlobal2Local(LocalCoordSystem, mNeighbourElasticContactForces[i], OldLocalElasticContactForce);
        
            GlobalElasticContactForce[0] = mNeighbourElasticContactForces[i][0];
            GlobalElasticContactForce[1] = mNeighbourElasticContactForces[i][1];
            GlobalElasticContactForce[2] = mNeighbourElasticContactForces[i][2];
                        
            GeometryFunctions::VectorGlobal2Local(LocalCoordSystem, GlobalElasticContactForce, LocalElasticContactForce);
            //we recover this way the old local forces projected in the new coordinates in the way they were in the old ones; Now they will be increased if its the necessary
                        
            double ViscoDampingLocalContactForce[3] = {0.0};
            double equiv_visco_damp_coeff_normal;
            double equiv_visco_damp_coeff_tangential;
            double ElasticLocalRotationalMoment[3] = {0.0};
            double ViscoLocalRotationalMoment[3] = {0.0};
                
            double LocalRelVel[3] = {0.0};
            GeometryFunctions::VectorGlobal2Local(LocalCoordSystem, RelVel, LocalRelVel);

            if (i < mContinuumInitialNeighborsSize) {
                int failure_id = mIniNeighbourFailureId[i];
                mContinuumConstitutiveLawArray[i]->CalculateForces(r_process_info, LocalElasticContactForce, LocalDeltDisp, kn_el, kt_el,
                                                                                    contact_sigma, contact_tau, failure_criterion_state, equiv_young, indentation,
                                                                                    calculation_area, acumulated_damage, this, neighbour_iterator, i,
                                                                                    r_process_info[TIME_STEPS], sliding, search_control, search_control_vector);

                mContinuumConstitutiveLawArray[i]->CalculateViscoDampingCoeff(equiv_visco_damp_coeff_normal, equiv_visco_damp_coeff_tangential, this, neighbour_iterator, kn_el, kt_el);
                mContinuumConstitutiveLawArray[i]->CalculateViscoDamping(LocalRelVel, ViscoDampingLocalContactForce, indentation, equiv_visco_damp_coeff_normal, equiv_visco_damp_coeff_tangential, sliding, failure_id);

            } else if (indentation > 0.0) { 
                double cohesive_force =  0.0;                
                const double previous_indentation = indentation + LocalDeltDisp[2];
                mDiscontinuumConstitutiveLaw->CalculateForces(r_process_info, OldLocalElasticContactForce, LocalElasticContactForce, LocalDeltDisp, LocalRelVel, indentation, previous_indentation, ViscoDampingLocalContactForce, cohesive_force, this, neighbour_iterator, sliding);                      
            }
        
            // Transforming to global forces and adding up
            double LocalContactForce[3] = {0.0};
            double GlobalContactForce[3] = {0.0};

            if (this->Is(DEMFlags::HAS_STRESS_TENSOR) && (i < mContinuumInitialNeighborsSize)) { // We leave apart the discontinuum neighbors (the same for the walls). The neighbor would not be able to do the same if we activate it. 
                mContinuumConstitutiveLawArray[i]->AddPoissonContribution(equiv_poisson, LocalCoordSystem, LocalElasticContactForce[2], calculation_area, mSymmStressTensor, this, neighbour_iterator);
            }

            AddUpForcesAndProject(OldLocalCoordSystem, LocalCoordSystem, LocalContactForce, LocalElasticContactForce, GlobalContactForce,
                                  GlobalElasticContactForce, ViscoDampingLocalContactForce, 0.0, rElasticForce, rContactForce, i);
            
            if (this->Is(DEMFlags::HAS_ROTATION)) {
                if (this->Is(DEMFlags::HAS_ROLLING_FRICTION) && !multi_stage_RHS) {
                    const double coeff_acc      = this->GetGeometry()[0].FastGetSolutionStepValue(PARTICLE_MOMENT_OF_INERTIA) / dt;
                    noalias(rInitialRotaMoment) = coeff_acc * ang_vel; // the moment needed to stop the spin in a single time step
                }   
                ComputeMoments(LocalElasticContactForce[2], mNeighbourTotalContactForces[i], rInitialRotaMoment, LocalCoordSystem[2], neighbour_iterator, indentation);
                if (i < mContinuumInitialNeighborsSize) {       
                    if (mIniNeighbourFailureId[i] == 0) {
                        mContinuumConstitutiveLawArray[i]->ComputeParticleRotationalMoments(this, neighbour_iterator, equiv_young, distance, calculation_area,
                                                                                            LocalCoordSystem, ElasticLocalRotationalMoment, ViscoLocalRotationalMoment);
                    }
                }
            }                            
            
            AddUpMomentsAndProject(LocalCoordSystem, ElasticLocalRotationalMoment, ViscoLocalRotationalMoment);
            
            if (r_process_info[CONTACT_MESH_OPTION] == 1 && (i < mContinuumInitialNeighborsSize) && this->Id() < neighbour_iterator_id) {
                CalculateOnContactElements(i, LocalElasticContactForce, contact_sigma, contact_tau, failure_criterion_state, acumulated_damage, time_steps);
            }

            if (this->Is(DEMFlags::HAS_STRESS_TENSOR) && (i < mContinuumInitialNeighborsSize)) {
                AddNeighbourContributionToStressTensor(GlobalElasticContactForce, LocalCoordSystem[2], distance, radius_sum);
            }

            AddContributionToRepresentativeVolume(distance, radius_sum, calculation_area);

        } // for each neighbor

            KRATOS_CATCH("")
    } //  ComputeBallToBallContactForce   

    void SphericContinuumParticle::AddContributionToRepresentativeVolume(const double distance,
                                                                         const double radius_sum,
                                                                         const double contact_area) {

        KRATOS_TRY

        double gap = distance - radius_sum;
        double real_distance = GetInteractionRadius() + 0.5 * gap;
        double& rRepresentative_Volume = this->GetGeometry()[0].FastGetSolutionStepValue(REPRESENTATIVE_VOLUME);
        rRepresentative_Volume += 0.33333333333333 * (real_distance * contact_area);

        KRATOS_CATCH("")
    }

    void SphericContinuumParticle::FinalizeSolutionStep(ProcessInfo& r_process_info) {
        KRATOS_TRY
        
        SphericParticle::FinalizeSolutionStep(r_process_info);
                
        //Update sphere mass and inertia taking into account the real volume of the represented volume:
        SetMass(this->GetGeometry()[0].FastGetSolutionStepValue(REPRESENTATIVE_VOLUME) * GetDensity());
        if (this->Is(DEMFlags::HAS_ROTATION) ){
            GetGeometry()[0].FastGetSolutionStepValue(PARTICLE_MOMENT_OF_INERTIA) = CalculateMomentOfInertia();
        }

        // the elemental variable is copied to a nodal variable in order to export the results onto GiD Post. Also a casting to double is necessary for GiD interpretation.
        KRATOS_CATCH("")
    }

    void SphericContinuumParticle::ReorderAndRecoverInitialPositionsAndFilter(std::vector<SphericParticle*>& TempNeighbourElements) {
        
        KRATOS_TRY

        unsigned int current_neighbors_size = mNeighbourElements.size();
        unsigned int initial_neighbors_size = mIniNeighbourIds.size();
        TempNeighbourElements.resize(initial_neighbors_size);

        for (unsigned int i = 0; i < initial_neighbors_size; i++) {
            TempNeighbourElements[i] = NULL;
        }
            
        // Loop over current neighbors
        for (unsigned int i = 0; i < current_neighbors_size; i++) {            
            SphericParticle* i_neighbour = mNeighbourElements[i];
            bool found = false;
            // Loop over initial neighbors
            for (unsigned int k = 0; k < initial_neighbors_size; k++) {                
    
                if (static_cast<int>(i_neighbour->Id()) == mIniNeighbourIds[k]) {                    
                    TempNeighbourElements[k] = i_neighbour;
                    found = true;
                    break;                 
                }
            }
            
            if (!found) {
                double other_radius = i_neighbour->GetInteractionRadius();
                double radius_sum = GetInteractionRadius() + other_radius;
                array_1d<double, 3> other_to_me_vect;
                noalias(other_to_me_vect) = this->GetGeometry()[0].Coordinates() - i_neighbour->GetGeometry()[0].Coordinates();
                double distance = DEM_MODULUS_3(other_to_me_vect);
                double indentation = radius_sum - distance;

                if (indentation > 0.0) {
                    TempNeighbourElements.push_back(i_neighbour);
                }
            }
        }
                                  
        mNeighbourElements.swap(TempNeighbourElements);
        
        if (mBondElements.size()) {
            for (unsigned int i = 0; i < mContinuumInitialNeighborsSize; i++) { 
                if (mNeighbourElements[i] == NULL) {
                    mBondElements[i] = NULL;
                }
                if ((mNeighbourElements[i] == NULL) && (mIniNeighbourFailureId[i] == 0)) {
                    mIniNeighbourFailureId[i] = 6;
                }
            }
        }

        KRATOS_CATCH("")
    }
    
    void SphericContinuumParticle::ReorderFEMneighbours() {        
        KRATOS_TRY

        unsigned int current_neighbors_size = mNeighbourRigidFaces.size();
        unsigned int initial_neighbors_size = mFemIniNeighbourIds.size();
        
        std::vector<DEMWall*> TempNeighbourElements;
        TempNeighbourElements.resize(initial_neighbors_size);

        for (unsigned int i = 0; i < initial_neighbors_size; i++) { TempNeighbourElements[i] = NULL; }
            
        // Loop over current neighbors
        for (unsigned int i = 0; i < current_neighbors_size; i++) {            
            DEMWall* i_neighbour = mNeighbourRigidFaces[i];
            bool found = false;
            // Loop over initial neighbors
            for (unsigned int k = 0; k < initial_neighbors_size; k++) {                
                if (static_cast<int>(i_neighbour->Id()) == mFemIniNeighbourIds[k]) {                    
                    TempNeighbourElements[k] = i_neighbour;
                    found = true;
                    break;                 
                }
            }
            
            if (!found) { TempNeighbourElements.push_back(i_neighbour); }
        }
                                  
        mNeighbourRigidFaces.swap(TempNeighbourElements);

        KRATOS_CATCH("")        
    }

    void SphericContinuumParticle::CalculateMeanContactArea(const bool has_mpi, const ProcessInfo& r_process_info) {
        KRATOS_TRY
        
	  if(!mContIniNeighArea.size()) return; //TODO: ugly fix //This means that for this case the areas are not being saved (because the consitutive law is not filling this vector)

        for (unsigned int i = 0; i < mContinuumInitialNeighborsSize; i++) {           
            SphericContinuumParticle* r_continuum_ini_neighbour = dynamic_cast<SphericContinuumParticle*>(mNeighbourElements[i]);            
            if (r_continuum_ini_neighbour == NULL) continue; //The initial neighbor was deleted at some point in time!!
            if(this->Id() > r_continuum_ini_neighbour->Id()) continue; //The Sphere with lower Id will do the job only                                   

            int index_of_the_neighbour_that_is_me = -1;
            for (unsigned int j = 0; j <  r_continuum_ini_neighbour->mContinuumInitialNeighborsSize; j++){
                if(this == r_continuum_ini_neighbour->mNeighbourElements[j]) index_of_the_neighbour_that_is_me = (int) j;
            }
            
            bool neigh_is_skin = (bool) r_continuum_ini_neighbour->mSkinSphere;
            if ((mSkinSphere && neigh_is_skin) || (!mSkinSphere && !neigh_is_skin)) { //both skin or both inner.
                double mean_area =  0.5 * (mContIniNeighArea[i] + r_continuum_ini_neighbour->mContIniNeighArea[index_of_the_neighbour_that_is_me]);
                mContIniNeighArea[i] = mean_area;
                r_continuum_ini_neighbour->mContIniNeighArea[index_of_the_neighbour_that_is_me] = mean_area;                
            }
            else if (!mSkinSphere && neigh_is_skin) {//we will store both the same only coming from the inner to the skin.                                    
                r_continuum_ini_neighbour->mContIniNeighArea[index_of_the_neighbour_that_is_me] = mContIniNeighArea[i];                    
            }
            else {
                mContIniNeighArea[i] = r_continuum_ini_neighbour->mContIniNeighArea[index_of_the_neighbour_that_is_me];
            }            
        }//loop neigh.

        return;
        KRATOS_CATCH("")
    }

    void SphericContinuumParticle::Calculate(const Variable<double>& rVariable, double& Output, const ProcessInfo& r_process_info) {
        KRATOS_TRY

        if (rVariable == DELTA_TIME) {
            double coeff = r_process_info[NODAL_MASS_COEFF];
            double mass = GetMass();

            if (coeff > 1.0) { KRATOS_THROW_ERROR(std::runtime_error, "The coefficient assigned for virtual mass is larger than one, virtual_mass_coeff= ", coeff) }
            else if ((coeff == 1.0) && (r_process_info[VIRTUAL_MASS_OPTION])) { Output = 9.0E09; }
            else {
                if (r_process_info[VIRTUAL_MASS_OPTION]) { mass /= 1 - coeff; }
                double K = GetYoung() * KRATOS_M_PI * GetRadius();
                Output = 0.34 * sqrt(mass / K);
                if (r_process_info[ROTATION_OPTION] == 1) { Output = Output * 0.5; } //factor for critical time step when rotation is allowed.                
            }
            return;
        }//CRITICAL DELTA CALCULATION

        KRATOS_CATCH("")
    }//Calculate

    void SphericContinuumParticle::CustomInitialize() {
        mSkinSphere     = &(this->GetGeometry()[0].FastGetSolutionStepValue(SKIN_SPHERE));
        this->GetGeometry()[0].FastGetSolutionStepValue(SKIN_SPHERE) = 1.0;
        mContinuumGroup = this->GetGeometry()[0].FastGetSolutionStepValue(COHESIVE_GROUP);
    }
    
    double SphericContinuumParticle::GetInitialDelta(int index) {
        if (index < (int) mIniNeighbourDelta.size()) return mIniNeighbourDelta[index];
        else  return 0.0;                 
    }

    double SphericContinuumParticle::GetInitialDeltaWithFEM(int index) {
        if(index< (int) mFemIniNeighbourDelta.size()) return mFemIniNeighbourDelta[index];
        else return 0.0;
    }

    void SphericContinuumParticle::CalculateOnContactElements(size_t i, double LocalElasticContactForce[3], double contact_sigma, double contact_tau, double failure_criterion_state, double acumulated_damage, int time_steps) {
        KRATOS_TRY

        if (!mBondElements.size()) return; // we skip this function if the vector of bonds hasn't been filled yet.
        ParticleContactElement* bond = mBondElements[i];
        if (bond == NULL) return; //This bond was never created (happens in some MPI cases, see CreateContactElements() in explicit_solve_continumm.h)

        bond->mLocalContactForce[0] = LocalElasticContactForce[0];
        bond->mLocalContactForce[1] = LocalElasticContactForce[1];
        bond->mLocalContactForce[2] = LocalElasticContactForce[2];
        bond->mContactSigma = contact_sigma;
        bond->mContactTau = contact_tau;
        bond->mContactFailure = mIniNeighbourFailureId[i];
        bond->mFailureCriterionState = failure_criterion_state;

        if ((time_steps == 0) || (acumulated_damage > bond->mUnidimendionalDamage)) {
            bond->mUnidimendionalDamage = acumulated_damage;
        }
        KRATOS_CATCH("")
    }                                  
    
} // namespace Kratos
