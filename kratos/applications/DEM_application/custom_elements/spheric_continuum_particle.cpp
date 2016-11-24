//
// Authors: Miquel Santasusana msantasusana@cimne.upc.edu
//          Miguel Angel Celigueta maceli@cimne.upc.edu
//

#include "spheric_continuum_particle.h"
#include <cmath>

namespace Kratos {

    SphericContinuumParticle::SphericContinuumParticle():SphericParticle() {
      //mOldSymmStressTensor = NULL;
    }

    SphericContinuumParticle::SphericContinuumParticle(IndexType NewId, GeometryType::Pointer pGeometry) : SphericParticle(NewId, pGeometry){
      //mOldSymmStressTensor = NULL;
    }

    SphericContinuumParticle::SphericContinuumParticle(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
    : SphericParticle(NewId, pGeometry, pProperties){
      //mOldSymmStressTensor = NULL;
    }

    SphericContinuumParticle::SphericContinuumParticle(IndexType NewId, NodesArrayType const& ThisNodes)
    : SphericParticle(NewId, ThisNodes){
      //mOldSymmStressTensor = NULL;
    }

    Element::Pointer SphericContinuumParticle::Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const {
      return SphericParticle::Pointer(new SphericContinuumParticle(NewId, GetGeometry().Create(ThisNodes), pProperties));
    }

    /// Destructor

    SphericContinuumParticle::~SphericContinuumParticle() {
        /*if (mOldSymmStressTensor!=NULL) {
            delete mOldSymmStressTensor;
            mOldSymmStressTensor = NULL;
        }*/
    }

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
                array_1d<double, 3> vector_of_zeros(3,0.0);
                //mArrayOfOldDeltaDisplacements.push_back(vector_of_zeros);
                //mArrayOfDeltaDisplacements.push_back(vector_of_zeros);
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

        mContinuumConstitutiveLawArray.resize(mContinuumInitialNeighborsSize);

        for (unsigned int i = 0; i < mContinuumInitialNeighborsSize; i++) {
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
            array_1d<double, 4>& Weight = this->mContactConditionWeights[i];

            ComputeConditionRelativeData(i,rFemNeighbours[i], LocalCoordSystem, DistPToB, Weight, wall_delta_disp_at_contact_point, wall_velocity_at_contact_point, ContactType);

            double initial_delta = -(DistPToB - GetRadius());

            mFemIniNeighbourIds[i] = rFemNeighbours[i]->Id();
            mFemIniNeighbourDelta[i] = initial_delta;
        }
    }//SetInitialFemContacts

    void SphericContinuumParticle::ContactAreaWeighting() //MISMI 10: POOYAN this could be done by calculating on the bars. not looking at the neighbors of my neighbors.
    {
        double alpha = 1.0;
        //double external_sphere_area = 4 * KRATOS_M_PI * GetRadius()*GetRadius();
        double effectiveVolumeRadius = EffectiveVolumeRadius();  //calculateEffectiveVolumeRadius
        double external_sphere_area = 4 * KRATOS_M_PI * effectiveVolumeRadius * effectiveVolumeRadius;
        double total_equiv_area = 0.0;
        double total_mContIniNeighArea = 0.0;
        int cont_ini_neighbours_size = mContinuumInitialNeighborsSize;
        Vector& cont_ini_neigh_area = GetValue(NEIGHBOURS_CONTACT_AREAS);
        bool print_debug_files = false;

        for (int i = 0; i < cont_ini_neighbours_size; i++) {
            SphericParticle* ini_cont_neighbour_iterator = mNeighbourElements[i];
            double other_radius = ini_cont_neighbour_iterator->GetRadius();
            double area = mContinuumConstitutiveLawArray[i]->CalculateContactArea(GetRadius(), other_radius, cont_ini_neigh_area); //This call fills the vector of areas only if the Constitutive Law wants.
            total_equiv_area += area;

        }
        if (print_debug_files) {
            std::ofstream outputfile("external_sphere_area-total_equiv_area.txt", std::ios_base::out | std::ios_base::app);
            outputfile << external_sphere_area << "  " << total_equiv_area << "\n";
            outputfile.close();
        }
        if (cont_ini_neighbours_size >= 6) { 
            if (!IsSkin()) {
                AuxiliaryFunctions::CalculateAlphaFactor3D(cont_ini_neighbours_size, external_sphere_area, total_equiv_area, alpha);
                if (print_debug_files) {
                    std::ofstream outputfile("alpha.txt", std::ios_base::out | std::ios_base::app);
                    outputfile << alpha << "\n";
                    outputfile.close();
                }
                for (unsigned int i = 0; i < cont_ini_neigh_area.size(); i++) {
                    cont_ini_neigh_area[i] = alpha * cont_ini_neigh_area[i];
                    total_mContIniNeighArea += cont_ini_neigh_area[i];
                } //for every neighbor

                if (print_debug_files) {
                    std::ofstream outputfile2("total_mContIniNeighArea-total_equiv_area.txt", std::ios_base::out | std::ios_base::app);
                    outputfile2 << total_mContIniNeighArea << "  " << total_equiv_area << "\n";
                    outputfile2.close();
                }
            }

            else {//skin sphere
                for (unsigned int i = 0; i < cont_ini_neigh_area.size(); i++) {
                    alpha = 1.00 * (1.40727)*(external_sphere_area / total_equiv_area)*((double(cont_ini_neighbours_size)) / 11.0);

                    cont_ini_neigh_area[i] = alpha * cont_ini_neigh_area[i];
                } //for every neighbor
            }
        }//if more than 3 neighbors
    }

    double SphericContinuumParticle::EffectiveVolumeRadius() {

        int cont_ini_neighbours_size = mContinuumInitialNeighborsSize;
        
        double effectiveVolumeRadiusSum = 0.0;
        
        for (int i = 0; i < cont_ini_neighbours_size; i++) {

            SphericContinuumParticle* neighbour_iterator = dynamic_cast<SphericContinuumParticle*>(mNeighbourElements[i]);
            double other_radius = neighbour_iterator->GetRadius();
            array_1d<double, 3> other_to_me_vect;
            noalias(other_to_me_vect) = this->GetGeometry()[0].Coordinates() - neighbour_iterator->GetGeometry()[0].Coordinates();
            double distance = DEM_MODULUS_3(other_to_me_vect);
            effectiveVolumeRadiusSum += 0.5 * (distance + GetRadius() - other_radius);
        }
        
        double effectiveVolumeRadius = effectiveVolumeRadiusSum / cont_ini_neighbours_size;

        return effectiveVolumeRadius;
    }


    /*void SphericContinuumParticle::InitializeSolutionStep(ProcessInfo& r_process_info) {
    
    KRATOS_TRY
        
    SphericParticle::InitializeSolutionStep(r_process_info);
    
    for (unsigned int i = 0; i < mContinuumInitialNeighborsSize; i++) {
        DEM_COPY_SECOND_TO_FIRST_3(mArrayOfOldDeltaDisplacements[i], mArrayOfDeltaDisplacements[i]);
    }
    
    KRATOS_CATCH("")
    }*/

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
        Vector& cont_ini_neigh_area              = this->GetValue(NEIGHBOURS_CONTACT_AREAS);
        int NeighbourSize = mNeighbourElements.size();
        GetGeometry()[0].GetSolutionStepValue(NEIGHBOUR_SIZE) = NeighbourSize;

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

            double kn_el = 0.0;
            double kt_el = 0.0;
            double DeltDisp[3] = {0.0};
            double RelVel[3] = {0.0};
            double LocalCoordSystem[3][3]    = {{0.0}, {0.0}, {0.0}};
            double OldLocalCoordSystem[3][3] = {{0.0}, {0.0}, {0.0}};
            bool sliding = false;

            double contact_tau = 0.0;
            double contact_sigma = 0.0;
            double failure_criterion_state = 0.0;
            double acumulated_damage = 0.0;

            // Getting neighbor properties
            double other_young = neighbour_iterator->GetYoung();
            double other_poisson = neighbour_iterator->GetPoisson();
            double equiv_poisson;
            if ((myPoisson + other_poisson) != 0.0) { equiv_poisson = 2.0 * myPoisson * other_poisson / (myPoisson + other_poisson); }
            else { equiv_poisson = 0.0; }

            double equiv_young = 2.0 * myYoung * other_young / (myYoung + other_young);
            double calculation_area = 0.0;
            const double equiv_shear = equiv_young / (2.0 * (1 + equiv_poisson));

            if (i < mContinuumInitialNeighborsSize) {
                mContinuumConstitutiveLawArray[i]->GetContactArea(GetRadius(), other_radius, cont_ini_neigh_area, i, calculation_area); //some Constitutive Laws get a value, some others calculate the value.
                mContinuumConstitutiveLawArray[i]->CalculateElasticConstants(kn_el, kt_el, distance, equiv_young, equiv_poisson, calculation_area, this, neighbour_iterator);
            }

            EvaluateDeltaDisplacement(DeltDisp, RelVel, LocalCoordSystem, OldLocalCoordSystem, other_to_me_vect, vel, delta_displ, neighbour_iterator, distance);

            if (this->Is(DEMFlags::HAS_ROTATION)) {
                RelativeDisplacementAndVelocityOfContactPointDueToRotationMatrix(DeltDisp, RelVel, OldLocalCoordSystem, other_radius, dt, ang_vel, neighbour_iterator);                
            }

            RelativeDisplacementAndVelocityOfContactPointDueToOtherReasons(r_process_info, DeltDisp, RelVel, OldLocalCoordSystem, LocalCoordSystem, neighbour_iterator);

            double LocalDeltDisp[3] = {0.0};
            double LocalElasticContactForce[3] = {0.0}; // 0: first tangential, // 1: second tangential, // 2: normal force
            double LocalElasticExtraContactForce[3] = {0.0};
            double GlobalElasticContactForce[3] = {0.0};
            double GlobalElasticExtraContactForce[3] = {0.0};
            double TotalGlobalElasticContactForce[3] = {0.0};
            double OldLocalElasticContactForce[3] = {0.0};

            FilterNonSignificantDisplacements(DeltDisp, RelVel, indentation);

            GeometryFunctions::VectorGlobal2Local(LocalCoordSystem, DeltDisp, LocalDeltDisp);

            RotateOldContactForces(OldLocalCoordSystem, LocalCoordSystem, mNeighbourElasticContactForces[i]);
            RotateOldContactForces(OldLocalCoordSystem, LocalCoordSystem, mNeighbourElasticExtraContactForces[i]);

            GeometryFunctions::VectorGlobal2Local(LocalCoordSystem, mNeighbourElasticContactForces[i], OldLocalElasticContactForce);

            GlobalElasticContactForce[0] = mNeighbourElasticContactForces[i][0];
            GlobalElasticContactForce[1] = mNeighbourElasticContactForces[i][1];
            GlobalElasticContactForce[2] = mNeighbourElasticContactForces[i][2];

            GeometryFunctions::VectorGlobal2Local(LocalCoordSystem, GlobalElasticContactForce, LocalElasticContactForce); //TODO: can we remove this? We should overwrite LocalElasticContactForce afterwards

            double ViscoDampingLocalContactForce[3] = {0.0};
            double equiv_visco_damp_coeff_normal;
            double equiv_visco_damp_coeff_tangential;
            double ElasticLocalRotationalMoment[3] = {0.0};
            double ViscoLocalRotationalMoment[3] = {0.0};

            double LocalRelVel[3] = {0.0};
            GeometryFunctions::VectorGlobal2Local(LocalCoordSystem, RelVel, LocalRelVel);

            if (i < mContinuumInitialNeighborsSize) {

                mContinuumConstitutiveLawArray[i]->CheckFailure(i, this, neighbour_iterator);
                
                mContinuumConstitutiveLawArray[i]->CalculateForces(r_process_info,
                                                                   OldLocalElasticContactForce,
                                                                   LocalElasticContactForce,
                                                                   LocalElasticExtraContactForce,
                                                                   LocalCoordSystem,
                                                                   LocalDeltDisp,
                                                                   kn_el,
                                                                   kt_el,
                                                                   contact_sigma,
                                                                   contact_tau,
                                                                   failure_criterion_state,
                                                                   equiv_young,
                                                                   equiv_shear,
                                                                   indentation,
                                                                   calculation_area,
                                                                   acumulated_damage,
                                                                   this,
                                                                   neighbour_iterator,
                                                                   i,
                                                                   r_process_info[TIME_STEPS],
                                                                   sliding,
                                                                   search_control,
                                                                   search_control_vector,
                                                                   equiv_visco_damp_coeff_normal,
                                                                   equiv_visco_damp_coeff_tangential,
                                                                   LocalRelVel,
                                                                   ViscoDampingLocalContactForce);

            } else if (indentation > 0.0) {
                double cohesive_force =  0.0;
                const double previous_indentation = indentation + LocalDeltDisp[2];
                mDiscontinuumConstitutiveLaw->CalculateForces(r_process_info, OldLocalElasticContactForce, LocalElasticContactForce, LocalDeltDisp, LocalRelVel, indentation, previous_indentation, ViscoDampingLocalContactForce, cohesive_force, this, neighbour_iterator, sliding);
            }

            // Transforming to global forces and adding up
            double LocalContactForce[3] = {0.0};
            double GlobalContactForce[3] = {0.0};
            
            if (this->Is(DEMFlags::HAS_STRESS_TENSOR) && (i < mContinuumInitialNeighborsSize)) { // We leave apart the discontinuum neighbors (the same for the walls). The neighbor would not be able to do the same if we activate it.
                mContinuumConstitutiveLawArray[i]->AddPoissonContribution(equiv_poisson, LocalCoordSystem, LocalElasticContactForce[2], calculation_area, mSymmStressTensor, this, neighbour_iterator, r_process_info, i, indentation);
            }

            array_1d<double, 3> other_ball_to_ball_forces(3,0.0);
            ComputeOtherBallToBallForces(other_ball_to_ball_forces);

            AddUpForcesAndProject(OldLocalCoordSystem, LocalCoordSystem, LocalContactForce, LocalElasticContactForce, LocalElasticExtraContactForce, GlobalContactForce,
                                  GlobalElasticContactForce, GlobalElasticExtraContactForce, TotalGlobalElasticContactForce,ViscoDampingLocalContactForce, 0.0, other_ball_to_ball_forces, rElasticForce, rContactForce, i, r_process_info); //TODO: replace the 0.0 with an actual cohesive force for discontinuum neighbours

            if (this->Is(DEMFlags::HAS_ROTATION)) {
                if (this->Is(DEMFlags::HAS_ROLLING_FRICTION) && !multi_stage_RHS) {
                    const double coeff_acc      = this->GetGeometry()[0].FastGetSolutionStepValue(PARTICLE_MOMENT_OF_INERTIA) / dt;
                    noalias(rInitialRotaMoment) = coeff_acc * ang_vel; // the moment needed to stop the spin in a single time step
                }
                ComputeMoments(LocalElasticContactForce[2], TotalGlobalElasticContactForce, rInitialRotaMoment, LocalCoordSystem[2], neighbour_iterator, indentation);
                if (i < mContinuumInitialNeighborsSize && mIniNeighbourFailureId[i] == 0) {
                    mContinuumConstitutiveLawArray[i]->ComputeParticleRotationalMoments(this, neighbour_iterator, equiv_young, distance, calculation_area,
                                                                                        LocalCoordSystem, ElasticLocalRotationalMoment, ViscoLocalRotationalMoment, equiv_poisson, indentation);
                }

                AddUpMomentsAndProject(LocalCoordSystem, ElasticLocalRotationalMoment, ViscoLocalRotationalMoment);
            }

            if (r_process_info[CONTACT_MESH_OPTION] == 1 && (i < mContinuumInitialNeighborsSize) && this->Id() < neighbour_iterator_id) {
                double total_local_elastic_contact_force[3] = {0.0};
                total_local_elastic_contact_force[0] = LocalElasticContactForce[0] + LocalElasticExtraContactForce[0];
                total_local_elastic_contact_force[1] = LocalElasticContactForce[1] + LocalElasticExtraContactForce[1];
                total_local_elastic_contact_force[2] = LocalElasticContactForce[2] + LocalElasticExtraContactForce[2];
                CalculateOnContactElements(i, total_local_elastic_contact_force, contact_sigma, contact_tau, failure_criterion_state, acumulated_damage, time_steps);
            }

            if (this->Is(DEMFlags::HAS_STRESS_TENSOR) && (i < mContinuumInitialNeighborsSize)) {
                AddNeighbourContributionToStressTensor(TotalGlobalElasticContactForce, LocalCoordSystem[2], distance, radius_sum, this);
            }

            AddContributionToRepresentativeVolume(distance, radius_sum, calculation_area);
            
            /*if (i < mContinuumInitialNeighborsSize) {
                DEM_COPY_SECOND_TO_FIRST_3(mArrayOfDeltaDisplacements[i], DeltDisp);
            }*/
        } // for each neighbor
                       
        KRATOS_CATCH("")
    } //  ComputeBallToBallContactForce

    void SphericContinuumParticle::FilterNonSignificantDisplacements(double DeltDisp[3], //IN GLOBAL AXES
                                                                     double RelVel[3], //IN GLOBAL AXES
                                                                     double& indentation) {

        if (fabs(DeltDisp[0]) < 1e-15){DeltDisp[0] = 0.0;}
        if (fabs(DeltDisp[1]) < 1e-15){DeltDisp[1] = 0.0;}
        if (fabs(DeltDisp[2]) < 1e-15){DeltDisp[2] = 0.0;}
        if (fabs(RelVel[0]) < 1e-15){RelVel[0] = 0.0;}
        if (fabs(RelVel[1]) < 1e-15){RelVel[1] = 0.0;}
        if (fabs(RelVel[2]) < 1e-15){RelVel[2] = 0.0;}
        if (fabs(indentation) < 1e-15){indentation = 0.0;}
    }

    void SphericContinuumParticle::AddContributionToRepresentativeVolume(const double distance,
                                                                         const double radius_sum,
                                                                         const double contact_area) {
        KRATOS_TRY

        double gap = distance - radius_sum;
        double real_distance = GetInteractionRadius() + 0.5 * gap;
        //double& rRepresentative_Volume = this->GetGeometry()[0].FastGetSolutionStepValue(REPRESENTATIVE_VOLUME);
        //KRATOS_WATCH(rRepresentative_Volume)
        //rRepresentative_Volume += 0.33333333333333 * (real_distance * contact_area);
        mPartialRepresentativeVolume += 0.33333333333333 * real_distance * contact_area;

        KRATOS_CATCH("")
    }

    /*
    void SphericContinuumParticle::CorrectRepresentativeVolume(double& rRepresentative_Volume, bool& is_smaller_than_sphere) {
    
        KRATOS_TRY
        
        SphericParticle::CorrectRepresentativeVolume(rRepresentative_Volume, is_smaller_than_sphere);
        
        //if (*mSkinSphere && is_smaller_than_sphere) rRepresentative_Volume *= 1.5; // This is the quotient between the volume of the cylinder circumscribed about a sphere and the latter
                                                                                   // So the minimum volume for a skin sphere is that of a cylinder
        //if (*mSkinSphere) rRepresentative_Volume *= 0.5;
                
        KRATOS_CATCH("")
    }
    */
        
    void SphericContinuumParticle::FinalizeSolutionStep(ProcessInfo& r_process_info) {
        KRATOS_TRY

        SphericParticle::FinalizeSolutionStep(r_process_info);

        //Update sphere mass and inertia taking into account the real volume of the represented volume:
        SetMass(this->GetGeometry()[0].FastGetSolutionStepValue(REPRESENTATIVE_VOLUME) * GetDensity());
        if (this->Is(DEMFlags::HAS_ROTATION) ){
            GetGeometry()[0].FastGetSolutionStepValue(PARTICLE_MOMENT_OF_INERTIA) = CalculateMomentOfInertia();
        }
        
        KRATOS_CATCH("")
    }
    
    void SphericContinuumParticle::GetStressTensorFromNeighbourStep1(){
        
        Set(DEMFlags::COPIED_STRESS_TENSOR, false);
        Set(DEMFlags::COPIED_STRESS_TENSOR2,false);
        if(!IsSkin()) return;
        
        for (unsigned int i=0; i<mNeighbourElements.size(); i++) { 
            if (mNeighbourElements[i] == NULL) continue;
            SphericContinuumParticle* p_neighbour = dynamic_cast<SphericContinuumParticle*>(mNeighbourElements[i]);
            if(!p_neighbour->IsSkin()) {
                *(mStressTensor) = *(p_neighbour->mStressTensor);
                *(mSymmStressTensor) = *(p_neighbour->mSymmStressTensor);
                Set(DEMFlags::COPIED_STRESS_TENSOR,true);
                break;
            }
        }
    }
    
    void SphericContinuumParticle::GetStressTensorFromNeighbourStep2(){
        
        if(!IsSkin()) return;
        if(Is(DEMFlags::COPIED_STRESS_TENSOR)) return;
        
        for (unsigned int i=0; i<mNeighbourElements.size(); i++) {
            if (mNeighbourElements[i] == NULL) continue;
            SphericContinuumParticle* p_neighbour = dynamic_cast<SphericContinuumParticle*>(mNeighbourElements[i]);
            if(p_neighbour->Is(DEMFlags::COPIED_STRESS_TENSOR)) {
                *(mStressTensor) = *(p_neighbour->mStressTensor);
                *(mSymmStressTensor) = *(p_neighbour->mSymmStressTensor);
                Set(DEMFlags::COPIED_STRESS_TENSOR2,true);
                break;
            }
        }
    }
    
    void SphericContinuumParticle::GetStressTensorFromNeighbourStep3(){
        
        if(!IsSkin()) return;
        if(Is(DEMFlags::COPIED_STRESS_TENSOR)) return;
        if(Is(DEMFlags::COPIED_STRESS_TENSOR2)) return;
        
        for (unsigned int i=0; i<mNeighbourElements.size(); i++) {
            if (mNeighbourElements[i] == NULL) continue;
            SphericContinuumParticle* p_neighbour = dynamic_cast<SphericContinuumParticle*>(mNeighbourElements[i]);
            if(p_neighbour->Is(DEMFlags::COPIED_STRESS_TENSOR2)) {
                *(mStressTensor) = *(p_neighbour->mStressTensor);
                *(mSymmStressTensor) = *(p_neighbour->mSymmStressTensor);
                break;
            }
        }
    }
    
    void SphericContinuumParticle::MarkNewSkinParticlesDueToBreakage() {
        
        KRATOS_TRY
        
        for (unsigned int i = 0; i < mContinuumInitialNeighborsSize; i++) {
            
            if (mNeighbourElements[i] == NULL) {
                *mSkinSphere = 1.0;
                break;
            }
        }

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
            for (unsigned int kk = 0; kk < initial_neighbors_size; kk++) {

                if (static_cast<int>(i_neighbour->Id()) == mIniNeighbourIds[kk]) {
                    TempNeighbourElements[kk] = i_neighbour;
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
                    mIniNeighbourFailureId[i] = 6; // Breakage due to search (neighbor not broken was not found during search)
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

    double SphericContinuumParticle::CalculateMaxSearchDistance(const bool has_mpi, const ProcessInfo& r_process_info) {
        
        KRATOS_TRY

        double max_local_search = 0.0;

        for (unsigned int i = 0; i < mContinuumInitialNeighborsSize; i++) {
	    if(mNeighbourElements[i] == NULL) continue;
            SphericContinuumParticle* r_continuum_ini_neighbour = dynamic_cast<SphericContinuumParticle*>(mNeighbourElements[i]);
            double search_dist = mContinuumConstitutiveLawArray[i]->LocalMaxSearchDistance(i, this, r_continuum_ini_neighbour);
            if (search_dist > max_local_search) max_local_search = search_dist;
        }
        
        return max_local_search;
        
        KRATOS_CATCH("")
    }

    bool SphericContinuumParticle::OverlappedParticleRemoval() {
        
        KRATOS_TRY

        for (unsigned int i = 0; i < mNeighbourElements.size(); i++) {
            if (mNeighbourElements[i] == NULL) continue;
            SphericParticle* ini_cont_neighbour_iterator = mNeighbourElements[i];
            double other_radius = ini_cont_neighbour_iterator->GetRadius();

            array_1d<double, 3> other_to_me_vect;
            noalias(other_to_me_vect) = this->GetGeometry()[0].Coordinates() - mNeighbourElements[i]->GetGeometry()[0].Coordinates();
            double distance = DEM_MODULUS_3(other_to_me_vect);

            double alpha = 1.0; // alpha = 1.0 means that the particle is completely inside another
            
            if (GetRadius() + distance < alpha * other_radius) {
                this->Set(TO_ERASE, true);
                return true;
            }
        }
        return false;
        
        KRATOS_CATCH("")
    }

    double SphericContinuumParticle::CalculateLocalMaxPeriod(const bool has_mpi, const ProcessInfo& r_process_info) {
        
        KRATOS_TRY
        
        double max_sqr_period = 0.0;

        for (unsigned int i = 0; i < mContinuumInitialNeighborsSize; i++) {
            SphericContinuumParticle* r_continuum_ini_neighbour = dynamic_cast<SphericContinuumParticle*>(mNeighbourElements[i]);
            double sqr_period = mContinuumConstitutiveLawArray[i]->LocalPeriod(i, this, r_continuum_ini_neighbour);
            if (sqr_period > max_sqr_period) { (max_sqr_period = sqr_period); }
        }
        for (unsigned int i = mContinuumInitialNeighborsSize; i < mNeighbourElements.size(); i++) {

            double sqr_period_discontinuum = mDiscontinuumConstitutiveLaw->LocalPeriod(i, this, mNeighbourElements[i]);
            if (sqr_period_discontinuum > max_sqr_period) { (max_sqr_period = sqr_period_discontinuum); }
        }
        return max_sqr_period;
        
        KRATOS_CATCH("")
    }

    void SphericContinuumParticle::CalculateMeanContactArea(const bool has_mpi, const ProcessInfo& r_process_info) {

        KRATOS_TRY
        
        Vector& cont_ini_neigh_area = GetValue(NEIGHBOURS_CONTACT_AREAS);

        if (!cont_ini_neigh_area.size()) return; //TODO: ugly fix //This means that for this case the areas are not being saved (because the constitutive law is not filling this vector)

        for (unsigned int i = 0; i < cont_ini_neigh_area.size(); i++) {

            SphericContinuumParticle* r_continuum_ini_neighbour = dynamic_cast<SphericContinuumParticle*>(mNeighbourElements[i]);
            if (r_continuum_ini_neighbour == NULL) continue; //The initial neighbor was deleted at some point in time!!

            if (this->Id() > r_continuum_ini_neighbour->Id()) continue; // The Sphere with a lower Id will do the job only                                   

            Vector& NeighbourContIniNeighArea = r_continuum_ini_neighbour->GetValue(NEIGHBOURS_CONTACT_AREAS);

            int index_of_the_neighbour_that_is_me = -1;

            for (unsigned int j = 0; j <  NeighbourContIniNeighArea.size(); j++) {
                boost::numeric::ublas::vector<int>& vector_of_ids_of_neighbours_of_this_neighbour = r_continuum_ini_neighbour->GetValue(NEIGHBOUR_IDS);
                if ((int) this->Id() == vector_of_ids_of_neighbours_of_this_neighbour[j]) index_of_the_neighbour_that_is_me = (int)j;
            }
            
            if (index_of_the_neighbour_that_is_me == -1) {               
                std::string message = "An element (Id " + std::to_string(this->Id()) + ") found a neighbor (had contact area) but the neighbor (Id " \
                                                        + std::to_string(r_continuum_ini_neighbour->Id()) + ") did not have area for that element  ";

                KRATOS_THROW_ERROR(std::runtime_error, message, 0);
            }

            bool neigh_is_skin = r_continuum_ini_neighbour->IsSkin();
            if ((IsSkin() && neigh_is_skin) || (!IsSkin() && !neigh_is_skin)) { //both skin or both inner
                double mean_area =  0.5 * (cont_ini_neigh_area[i] + NeighbourContIniNeighArea[index_of_the_neighbour_that_is_me]);
                cont_ini_neigh_area[i] = mean_area;
                NeighbourContIniNeighArea[index_of_the_neighbour_that_is_me] = mean_area;
            }
            else if (!IsSkin() && neigh_is_skin) {//we will store both the same only coming from the inner to the skin.
                NeighbourContIniNeighArea[index_of_the_neighbour_that_is_me] = cont_ini_neigh_area[i];
            }
            else {
                cont_ini_neigh_area[i] = NeighbourContIniNeighArea[index_of_the_neighbour_that_is_me];
            }
        } //loop neighbors

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

    void SphericContinuumParticle::Initialize(const ProcessInfo& r_process_info){

        SphericParticle::Initialize(r_process_info);

        /*if (this->Is(DEMFlags::HAS_STRESS_TENSOR)) {
            mOldSymmStressTensor  = new Matrix(3,3);
            *mOldSymmStressTensor = ZeroMatrix(3,3);
        }
        else {
            mOldSymmStressTensor = NULL;
        }*/
        SetValue(NEIGHBOURS_CONTACT_AREAS, Vector());

        mSkinSphere     = &(this->GetGeometry()[0].FastGetSolutionStepValue(SKIN_SPHERE));
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
