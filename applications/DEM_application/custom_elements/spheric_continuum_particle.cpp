
//   Project Name:        Kratos
//   Last Modified by:    $Author: M.Santasusana $
//   Date:                $Date: 2006-11-27 16:07:33 $
//   Revision:            $Revision: 1.1.1.1 $
//
//


// System includes
#include <string>
#include <iostream>
#include <iomanip> // to improve std::cout precision

// External includes


// Project includes
#include "includes/define.h"
#include "spheric_particle.h"
#include "custom_utilities/GeometryFunctions.h"
#include "custom_utilities/AuxiliaryFunctions.h"
#include "DEM_application.h"
#include "utilities/openmp_utils.h"

//#include "custom_constitutive/DEM_continuum_constitutive_law.h"

//TIMER....................
#include "utilities/timer.h"

#define CUSTOMTIMER 0  // ACTIVATES AND DISABLES ::TIMER:::::

#ifdef CUSTOMTIMER
#define KRATOS_TIMER_START(t) Timer::Start(t);
#define KRATOS_TIMER_STOP(t) Timer::Stop(t);
#else
#define KRATOS_TIMER_START(t)
#define KRATOS_TIMER_STOP(t)
#endif

//......................
//std::cout<<print...<<std::endl;

namespace Kratos {

      SphericContinuumParticle::SphericContinuumParticle() : SphericParticle(){}

      SphericContinuumParticle::SphericContinuumParticle( IndexType NewId, GeometryType::Pointer pGeometry) : SphericParticle(NewId, pGeometry){}

      SphericContinuumParticle::SphericContinuumParticle(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
      : SphericParticle(NewId, pGeometry, pProperties){}

      SphericContinuumParticle::SphericContinuumParticle(IndexType NewId, NodesArrayType const& ThisNodes)
      : SphericParticle(NewId, ThisNodes){}

      Element::Pointer SphericContinuumParticle::Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const
      {
        return SphericParticle::Pointer(new SphericContinuumParticle(NewId, GetGeometry().Create(ThisNodes), pProperties));
      }

    /// Destructor.

      SphericContinuumParticle::~SphericContinuumParticle(){}


    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void SphericContinuumParticle::FullInitialize(const ProcessInfo& r_process_info) {
        KRATOS_TRY        
        MemberDeclarationFirstStep(r_process_info);
        ContinuumSphereMemberDeclarationFirstStep(r_process_info);
        Initialize();
        CreateDiscontinuumConstitutiveLaws(r_process_info);
        KRATOS_CATCH("")
    }

    void SphericContinuumParticle::SetInitialSphereContacts(ProcessInfo& rCurrentProcessInfo) {

        /*
         * 
         * ELEMENT_NEIGHBOURS / NEIGHBOURS_IDS: the ones important for calculating forces.
         * INI_NEIGHBOURS_IDS: the ones to be treated specially due to initial delta or continuum case.
         * INI_CONTINUUM_NEIGHBOURS_IDS: only the ones that are continuum at 0 step and we should treat the possible detachment.
         * 
         * These 3 classes do NOT coincide at t=0!

         */
        mContinuumIniNeighbourElements.clear(); ///////////////////////

        size_t ini_size = 0;
        size_t continuum_ini_size = 0;
        size_t cont_ini_mapping_index = 0;

        unsigned int neighbours_size = mNeighbourElements.size(); //////////////////////
        mIniNeighbourIds.resize(neighbours_size);
        mIniNeighbourToIniContinuum.resize(neighbours_size);
        mIniNeighbourDelta.resize(neighbours_size);
        mIniNeighbourFailureId.resize(neighbours_size);
        mMapping_New_Ini.resize(neighbours_size);

        //SAVING THE INITIAL NEIGHBOURS, THE DELTAS AND THE FAILURE ID

        for (unsigned int i = 0; i < mNeighbourElements.size(); i++) {
            SphericContinuumParticle* neighbour_iterator = dynamic_cast<SphericContinuumParticle*> (mNeighbourElements[i]);

            array_1d<double, 3 > other_to_me_vect = this->GetGeometry()[0].Coordinates() - neighbour_iterator->GetGeometry()[0].Coordinates();

            double distance = sqrt(other_to_me_vect[0] * other_to_me_vect[0] +
                    other_to_me_vect[1] * other_to_me_vect[1] +
                    other_to_me_vect[2] * other_to_me_vect[2]);

            double radius_sum = GetRadius() + neighbour_iterator->GetRadius();
            double initial_delta = radius_sum - distance;

            int r_other_continuum_group = neighbour_iterator->mContinuumGroup;

            ini_size++;
            mIniNeighbourIds[ini_size - 1] = neighbour_iterator->Id();
            mIniNeighbourFailureId[ini_size - 1] = 1;
            mMapping_New_Ini[ini_size - 1] = ini_size - 1;
            mIniNeighbourToIniContinuum[ini_size - 1] = -1; //-1 is initial but not continuum.             
            mIniNeighbourDelta[ini_size - 1] = initial_delta;


            if ((r_other_continuum_group == mContinuumGroup) && (mContinuumGroup != 0)) {

                mIniNeighbourToIniContinuum[ini_size - 1] = cont_ini_mapping_index;

                mIniNeighbourFailureId[ini_size - 1] = 0;

                continuum_ini_size++;
                cont_ini_mapping_index++;

                mContinuumIniNeighbourElements.push_back(neighbour_iterator);

                mMapping_New_Cont.push_back(-1);


            }//if ( (r_other_continuum_group == mContinuumGroup) && (mContinuumGroup != 0) )

        } //end for: ParticleWeakIteratorType ineighbour

        ContactAreaWeighting();

    }//SetInitialSphereContacts


    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void SphericContinuumParticle::CreateContinuumConstitutiveLaws(ProcessInfo& rCurrentProcessInfo) {
        unsigned int cont_neigh_size = mContinuumIniNeighbourElements.size();

        mContinuumConstitutiveLawArray.resize(cont_neigh_size);

        for (unsigned int i = 0; i < cont_neigh_size; i++) {

            DEMContinuumConstitutiveLaw::Pointer NewContinuumConstitutiveLaw = GetProperties()[DEM_CONTINUUM_CONSTITUTIVE_LAW_POINTER]-> Clone();
            mContinuumConstitutiveLawArray[i] = NewContinuumConstitutiveLaw;
            mContinuumConstitutiveLawArray[i]->Initialize(rCurrentProcessInfo);

        }
    }

    //**************************************************************************************************************************************************
    //**************************************************************************************************************************************************

    void SphericContinuumParticle::SetInitialFemContacts() {
        const std::vector<double>& RF_Pram = this->mNeighbourRigidFacesPram;
        std::vector<DEMWall*>& rFemNeighbours = this->mNeighbourRigidFaces;

        unsigned int fem_neighbours_size = rFemNeighbours.size();

        mFemIniNeighbourIds.resize(fem_neighbours_size);
        mFemMappingNewIni.resize(fem_neighbours_size);
        mFemIniNeighbourDelta.resize(fem_neighbours_size);

        for (unsigned int i = 0; i < rFemNeighbours.size(); i++) {
            int ino1 = i * 16;
            double DistPToB = RF_Pram[ino1 + 9];
            int iNeighborID = static_cast<int> (RF_Pram[ino1 + 14]);
            double initial_delta = -(DistPToB - GetRadius());

            mFemIniNeighbourIds[i] = iNeighborID;
            mFemMappingNewIni[i] = i;
            mFemIniNeighbourDelta[i] = initial_delta;
        }

    }//SetInitialFemContacts              

    void SphericContinuumParticle::ContactAreaWeighting() //MISMI 10: POOYAN this could be done by calculating on the bars. not looking at the neighbous of my neighbours.
    {

        double alpha = 1.0;
        double external_sphere_area = 4 * KRATOS_M_PI * GetRadius()*GetRadius();
        double total_equiv_area = 0.0;

        int cont_ini_neighbours_size = mContinuumIniNeighbourElements.size();

        mcont_ini_neigh_area.resize(cont_ini_neighbours_size);

        //computing the total equivalent area

        size_t index = 0;

        for (unsigned int i = 0; i < mContinuumIniNeighbourElements.size(); i++) {
            SphericParticle* ini_cont_neighbour_iterator = mContinuumIniNeighbourElements[i];

            double other_radius = ini_cont_neighbour_iterator->GetRadius();
            double equiv_radius = 2 * GetRadius() * other_radius / (GetRadius() + other_radius);
            double equiv_area = 0.25 * KRATOS_M_PI * equiv_radius * equiv_radius; //we now take 1/2 of the effective radius.
            total_equiv_area += equiv_area;

            mcont_ini_neigh_area[index] = equiv_area; //*
            index++; //*

            //RECAREY
            //array_1d<double,3> other_to_me_vect   = this->GetGeometry()[0].Coordinates() - ini_cont_neighbour_iterator->GetGeometry()[0].Coordinates();
            //double distance                       = sqrt(other_to_me_vect[0] * other_to_me_vect[0] +
            //                                              other_to_me_vect[1] * other_to_me_vect[1] +
            //                                              other_to_me_vect[2] * other_to_me_vect[2]);
            //do this in pragma omp critical!!
            //rCurrentProcessInfo[TOTAL_CONTACT_DISTANCES] +=  0.5*distance*distance;

        } //for every neighbor

        if (cont_ini_neighbours_size >= 4) //more than 3 neighbors. 
        {
            if (!*mSkinSphere) {

                AuxiliaryFunctions::CalculateAlphaFactor3D(cont_ini_neighbours_size, external_sphere_area, total_equiv_area, alpha);

                size_t not_skin_index = 0;

                for (unsigned int i = 0; i < mContinuumIniNeighbourElements.size(); i++) {

                    mcont_ini_neigh_area[not_skin_index] = alpha * mcont_ini_neigh_area[not_skin_index];
                    not_skin_index++;

                } //for every neighbor

            }//if(!*mSkinSphere)

            else //skin sphere 
            {

                size_t skin_index = 0;

                for (unsigned int i = 0; i < mContinuumIniNeighbourElements.size(); i++) {

                    alpha = 1.00 * (1.40727)*(external_sphere_area / total_equiv_area)*((double(cont_ini_neighbours_size)) / 11);
                    mcont_ini_neigh_area[skin_index] = alpha * mcont_ini_neigh_area[skin_index];

                    skin_index++;
                } //for every neighbor

            }//skin particles.

        }//if more than 3 neighbors

    } //Contact Area Weighting           

    /**
     * Calculates all particle's ball-to-ball forces based on its neighbors
     * @param rContactForce
     * @param rContactMoment
     * @param rCurrentProcessInfo
     **/
    void SphericContinuumParticle::ComputeBallToBallContactForce(array_1d<double, 3>& rElasticForce,
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
            RotaAcc[0] = ang_vel[0] * dt_i;
            RotaAcc[1] = ang_vel[1] * dt_i;
            RotaAcc[2] = ang_vel[2] * dt_i;

            rInitialRotaMoment[0] = RotaAcc[0] * moment_of_inertia;
            rInitialRotaMoment[1] = RotaAcc[1] * moment_of_inertia;
            rInitialRotaMoment[2] = RotaAcc[2] * moment_of_inertia;
        }
               
        for (unsigned int i_neighbour_count = 0; i_neighbour_count < mNeighbourElements.size(); i_neighbour_count++) {
                                
            SphericContinuumParticle* neighbour_iterator = dynamic_cast<SphericContinuumParticle*>(mNeighbourElements[i_neighbour_count]);
            
            unsigned int neighbour_iterator_id = neighbour_iterator->Id();
            
            array_1d<double, 3> other_to_me_vect = this->GetGeometry()[0].Coordinates() - neighbour_iterator->GetGeometry()[0].Coordinates();
            const double &other_radius = neighbour_iterator->GetRadius();
            double distance = sqrt(other_to_me_vect[0] * other_to_me_vect[0] +
                    other_to_me_vect[1] * other_to_me_vect[1] +
                    other_to_me_vect[2] * other_to_me_vect[2]);
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

            //const int mapping_new_ini = mMapping_New_Ini[i_neighbour_count]; //*
            const int mapping_new_cont = mMapping_New_Cont[i_neighbour_count];

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


            //KRATOS_WATCH(mapping_new_cont)

            if (mapping_new_cont != -1) {
                mContinuumConstitutiveLawArray[mapping_new_cont]-> CalculateContactArea(GetRadius(), other_radius, calculation_area);
                mContinuumConstitutiveLawArray[mapping_new_cont]-> CalculateElasticConstants(kn_el, kt_el, initial_dist, equiv_young, equiv_poisson, calculation_area);
            } else {
                mDiscontinuumConstitutiveLaw -> CalculateContactArea(GetRadius(), other_radius, calculation_area);
                mDiscontinuumConstitutiveLaw -> CalculateElasticConstants(kn_el, kt_el, initial_dist, equiv_young, equiv_poisson, calculation_area);
            }


            //            bool equivalent_area_method = false; //  Miquel - should be moved out of SphericcontinuumParticle
            //            if (equivalent_area_method) {
            //                if (mapping_new_cont != -1) {
            //                    calculation_area = mcont_ini_neigh_area[mapping_new_cont];}
            //                kn_el = equiv_young * calculation_area * radius_sum_i; //0.9237*equiv_young * calculation_area; //MSIMSI 1: initial gap? we are only dividing by radius sum, it is not correct..
            //                kt_el = kn_el / (2.0 + equiv_poisson + equiv_poisson); //0.0*...
            //                aux_norm_to_tang = sqrt(kt_el / kn_el);}

            //if (mCriticalTimeOption){
            /*if (this->Is(DEMFlags::HAS_CRITICAL_TIME)) { //MA: I comment this out because historic or HISTORICAL_MIN_K seems to be doing nothing in the whole DEM application...
                double historic = rCurrentProcessInfo[HISTORICAL_MIN_K];

                if ((kn_el < historic) || (kt_el < historic)) {
                    historic = std::min(kn_el, kt_el);
                }
            }*/
            EvaluateDeltaDisplacement(DeltDisp, RelVel, LocalCoordSystem, OldLocalCoordSystem, other_to_me_vect, vel, delta_displ, neighbour_iterator, distance);

            if (this->Is(DEMFlags::HAS_ROTATION)) {
                //DisplacementDueToRotation(DeltDisp, RelVel, OldLocalCoordSystem, other_radius, dt, ang_vel, neighbour_iterator);
                DisplacementDueToRotationMatrix(DeltDisp, RelVel, OldLocalCoordSystem, other_radius, dt, ang_vel, neighbour_iterator);
            }

            double LocalDeltDisp[3] = {0.0};
            double LocalElasticContactForce[3] = {0.0}; // 0: first tangential, // 1: second tangential, // 2: normal force
            double GlobalElasticContactForce[3] = {0.0};

            GlobalElasticContactForce[0] = mNeighbourElasticContactForces[i_neighbour_count][0];
            GlobalElasticContactForce[1] = mNeighbourElasticContactForces[i_neighbour_count][1];
            GlobalElasticContactForce[2] = mNeighbourElasticContactForces[i_neighbour_count][2];

            GeometryFunctions::VectorGlobal2Local(OldLocalCoordSystem, GlobalElasticContactForce, LocalElasticContactForce);
            //we recover this way the old local forces projected in the new coordinates in the way they were in the old ones; Now they will be increased if its the necessary
            GeometryFunctions::VectorGlobal2Local(OldLocalCoordSystem, DeltDisp, LocalDeltDisp);

            /* Translational Forces */

            if (indentation > 0.0 || (mNeighbourFailureId[i_neighbour_count] == 0))//*  //#3
            {
                if (mapping_new_cont != -1) {//Normal Forces
                    mContinuumConstitutiveLawArray[mapping_new_cont]-> CalculateForces(
                            rCurrentProcessInfo,
                            LocalElasticContactForce,
                            LocalDeltDisp,
                            kn_el,
                            kt_el,
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
                            mNeighbourFailureId[i_neighbour_count],
                            mapping_new_cont);
                }
            } // compression or cohesive contact

            //          Viscodamping (applied locally)
            double ViscoDampingLocalContactForce[3] = {0.0};
            double equiv_visco_damp_coeff_normal;
            double equiv_visco_damp_coeff_tangential;

            if (indentation > 0.0 || (mNeighbourFailureId[i_neighbour_count] == 0)) {

                double LocalRelVel[3] = {0.0};
                GeometryFunctions::VectorGlobal2Local(LocalCoordSystem, RelVel, LocalRelVel);
                if (mapping_new_cont != -1) {
                    mContinuumConstitutiveLawArray[mapping_new_cont]->CalculateViscoDampingCoeff(equiv_visco_damp_coeff_normal,
                            equiv_visco_damp_coeff_tangential,
                            this,
                            neighbour_iterator,
                            kn_el,
                            kt_el);

                    mContinuumConstitutiveLawArray[mapping_new_cont]->CalculateViscoDamping(LocalRelVel,
                            ViscoDampingLocalContactForce,
                            indentation,
                            equiv_visco_damp_coeff_normal,
                            equiv_visco_damp_coeff_tangential,
                            sliding);
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

            if (rCurrentProcessInfo[STRESS_STRAIN_OPTION] && mapping_new_cont != -1) {AddPoissonContribution(equiv_poisson, 
                                                                                        LocalCoordSystem, 
                                                                                        LocalElasticContactForce[2], 
                                                                                        calculation_area);}
            
            AddUpForcesAndProject(OldLocalCoordSystem, LocalCoordSystem, LocalContactForce, LocalElasticContactForce, GlobalContactForce,
                    GlobalElasticContactForce, ViscoDampingLocalContactForce, 0.0, rElasticForce, rContactForce, i_neighbour_count);
                       
            array_1d<double, 3> temp_force = ZeroVector(3);
        
            temp_force[0] = GlobalContactForce[0];
            temp_force[1] = GlobalContactForce[1];
            temp_force[2] = GlobalContactForce[2];
            
            if (this->Is(DEMFlags::HAS_ROTATION)) {
                ComputeMoments(LocalElasticContactForce[2], temp_force, rInitialRotaMoment, LocalCoordSystem[2], neighbour_iterator, indentation);
            }

            if (rCurrentProcessInfo[CONTACT_MESH_OPTION] == 1 && (mapping_new_cont != -1) && this->Id() < neighbour_iterator_id) {
                CalculateOnContactElements(neighbour_iterator_id, 
                                            i_neighbour_count, 
                                            mapping_new_cont, 
                                            LocalElasticContactForce, 
                                            contact_sigma, 
                                            contact_tau, 
                                            failure_criterion_state, 
                                            acumulated_damage, 
                                            time_steps);}

            if (rCurrentProcessInfo[STRESS_STRAIN_OPTION] && mapping_new_cont != -1) {
                AddNeighbourContributionToStressTensor(GlobalElasticContactForce, 
                                                        LocalCoordSystem[2], 
                                                        distance, 
                                                        radius_sum, 
                                                        calculation_area);}            
        } //for each neighbor

        KRATOS_CATCH("")

    } //ComputeBallToBallContactForce


    void SphericContinuumParticle::ApplyLocalMomentsDamping(const ProcessInfo& rCurrentProcessInfo) {

        KRATOS_TRY

        array_1d<double, 3 > & RotaMoment = this->GetGeometry()[0].FastGetSolutionStepValue(PARTICLE_MOMENT);
        double RotaDampRatio = this->GetGeometry()[0].FastGetSolutionStepValue(PARTICLE_ROTATION_DAMP_RATIO);

        // LOCAL DAMPING OPTION FOR THE UNBALANCED FORCES (IN GLOBAL COORDINATES).        
        for (int iDof = 0; iDof < 3; iDof++) {
            if (this->GetGeometry()[0].FastGetSolutionStepValue(ANGULAR_VELOCITY)[iDof] > 0.0) {
                RotaMoment[iDof] = RotaMoment[iDof] - RotaDampRatio * fabs(RotaMoment[iDof]);
            } else {
                RotaMoment[iDof] = RotaMoment[iDof] + RotaDampRatio * fabs(RotaMoment[iDof]);
            }
        }

        KRATOS_CATCH("")

    } //ApplyLocalMomentsDamping      

    void SphericContinuumParticle::SymmetrizeTensor(const ProcessInfo& rCurrentProcessInfo) //MSIMSI10
    {

        KRATOS_TRY
        
        KRATOS_CATCH("")
    } //SymmetrizeTensor

    void SphericContinuumParticle::InitializeSolutionStep(ProcessInfo& rCurrentProcessInfo) {         
        KRATOS_TRY
        
        SphericParticle::InitializeSolutionStep(rCurrentProcessInfo);
        
        KRATOS_CATCH("")
    }//void SphericContinuumParticle::InitializeSolutionStep(ProcessInfo& r_process_info)

    void SphericContinuumParticle::FinalizeSolutionStep(ProcessInfo& rCurrentProcessInfo) {

        KRATOS_TRY
        if (rCurrentProcessInfo[PRINT_SKIN_SPHERE] == 1) {
            this->GetGeometry()[0].FastGetSolutionStepValue(EXPORT_SKIN_SPHERE) = double(*mSkinSphere);
        }        

        // the elemental variable is copied to a nodal variable in order to export the results onto GiD Post. Also a casting to double is necessary for GiD interpretation.
        KRATOS_CATCH("")
    }

    //VELL:

    void SphericContinuumParticle::ComputeNewNeighboursHistoricalData(std::vector<unsigned int>& mTempNeighboursIds, //We are passing all these temporal vectors as arguments because creating them inside the function is slower (memory allocation and deallocation)
            std::vector<array_1d<double, 3> >& mTempNeighbourElasticContactForces,
            std::vector<array_1d<double, 3> >& mTempNeighbourTotalContactForces,
            std::vector<SphericParticle*>& mTempNeighbourElements,
            std::vector<double>& mTempNeighboursDelta,
            std::vector<int>& mTempNeighboursFailureId,
            std::vector<int>& mTempNeighboursMapping,
            std::vector<int>& mTempContNeighboursMapping) {

        KRATOS_TRY
        mTempNeighbourElements.swap(mNeighbourElements);

        unsigned int temp_size = mTempNeighbourElements.size();

        mNeighbourElements.clear(); //////////////////////////////////////              

        mTempNeighboursIds.resize(temp_size);
        mTempNeighboursDelta.resize(temp_size);
        mTempNeighboursFailureId.resize(temp_size);
        mTempNeighbourElasticContactForces.resize(temp_size);
        mTempNeighbourTotalContactForces.resize(temp_size);
        mTempNeighboursMapping.resize(temp_size);
        mTempContNeighboursMapping.resize(temp_size);

        array_1d<double, 3> vector_of_zeros;
        vector_of_zeros[0] = vector_of_zeros[1] = vector_of_zeros[2] = 0.0;

        unsigned int neighbour_counter = 0; //Not increased at every iteration!! only if found as a real neighbor.

        for (unsigned int i = 0; i < mTempNeighbourElements.size(); i++) {
            SphericParticle* i_neighbour = mTempNeighbourElements[i];

            double ini_delta = 0.0;
            int failure_id = 1;
            array_1d<double, 3> neigh_elastic_forces = vector_of_zeros;
            array_1d<double, 3> neigh_total_forces = vector_of_zeros;
            int mapping_new_ini = -1;
            int mapping_new_cont = -1;

            //Loop Over Initial Neighbors        

            for (unsigned int k = 0; k != mIniNeighbourIds.size(); k++) {
                if (static_cast<int> ((i_neighbour)->Id()) == mIniNeighbourIds[k]) {
                    ini_delta = mIniNeighbourDelta[k];
                    failure_id = mIniNeighbourFailureId[k];
                    mapping_new_ini = k;
                    mapping_new_cont = mIniNeighbourToIniContinuum[k];
                    break;
                }
            }

            //Loop Over Last time-step Neighbors
            for (unsigned int j = 0; j != mOldNeighbourIds.size(); j++) {
                if ((i_neighbour)->Id() == mOldNeighbourIds[j]) {
                    neigh_elastic_forces = mNeighbourElasticContactForces[j];
                    neigh_total_forces = mNeighbourTotalContactForces[j];
                    break;
                }
            }

            //Judge if it is neighbor            
            double other_radius = i_neighbour->GetRadius();
            double radius_sum = GetRadius() + other_radius;
            array_1d<double, 3> other_to_me_vect = this->GetGeometry()[0].Coordinates() - i_neighbour->GetGeometry()[0].Coordinates();
            double distance = sqrt(other_to_me_vect[0] * other_to_me_vect[0] + other_to_me_vect[1] * other_to_me_vect[1] + other_to_me_vect[2] * other_to_me_vect[2]);
            double indentation = radius_sum - distance - ini_delta;

            if (indentation > 0.0 || failure_id == 0) //WE NEED TO SET A NUMERICAL TOLERANCE FUNCTION OF THE RADIUS.  MSIMSI 10
            {
                mNeighbourElements.push_back(i_neighbour); ///////////////////////////////*/

                mTempNeighboursIds[neighbour_counter] = (i_neighbour)->Id();
                mTempNeighboursMapping[neighbour_counter] = mapping_new_ini;
                mTempContNeighboursMapping[neighbour_counter] = mapping_new_cont;
                mTempNeighboursDelta[neighbour_counter] = ini_delta;
                mTempNeighboursFailureId[neighbour_counter] = failure_id;
                mTempNeighbourElasticContactForces[neighbour_counter] = neigh_elastic_forces;
                mTempNeighbourTotalContactForces[neighbour_counter] = neigh_total_forces;

                neighbour_counter++;

            }

        }//for ParticleWeakIteratorType i

        int final_size = mNeighbourElements.size();
        mTempNeighboursIds.resize(final_size);
        mTempNeighboursDelta.resize(final_size);
        mTempNeighboursFailureId.resize(final_size);
        mTempNeighbourElasticContactForces.resize(final_size);
        mTempNeighbourTotalContactForces.resize(final_size);
        mTempNeighboursMapping.resize(final_size);
        mTempContNeighboursMapping.resize(final_size);

        mMapping_New_Ini.swap(mTempNeighboursMapping);
        mMapping_New_Cont.swap(mTempContNeighboursMapping);
        mOldNeighbourIds.swap(mTempNeighboursIds);
        mNeighbourDelta.swap(mTempNeighboursDelta);
        mNeighbourFailureId.swap(mTempNeighboursFailureId);
        mNeighbourElasticContactForces.swap(mTempNeighbourElasticContactForces);
        mNeighbourTotalContactForces.swap(mTempNeighbourTotalContactForces);

        KRATOS_CATCH("")
    } //ComputeNewNeighboursHistoricalData

    /*
     //RIC!!!!!
     void SphericContinuumParticle::ComputeNewNeighboursHistoricalData() //NOTA: LOOP SOBRE TOTS ELS VEINS PROVISIONALS, TEN KEDERAS UNS QUANTS FENT PUSHBACK. ALS VECTORS DELTA ETC.. HI HAS DE POSAR
        //LA POSICIÓ DELS QUE SON DEFINITIUS.
        {
           KRATOS_TRY

           ParticleWeakVectorType& TempNeighbours = this->GetValue(NEIGHBOUR_ELEMENTS);

           unsigned int neighbour_counter       = 0;
           unsigned int temp_neighbour_counter  = 0;

           unsigned int temp_size = TempNeighbours.size();


           std::vector<int>&                  temp_neighbours_ids = mTempNeighboursIds;
           std::vector<double>&               temp_neighbours_delta = mTempNeighboursDelta;
           std::vector<int>&                  temp_neighbours_failure_id = mTempNeighboursFailureId;
           std::vector<array_1d<double, 3> >& temp_neighbours_contact_forces = mTempNeighbourElasticContactForces;
           std::vector<int>&                  temp_neighbours_mapping = mTempNeighboursMapping;
           std::vector<int>&                  temp_cont_neighbours_mapping = mTempContNeighboursMapping;   
        
           temp_neighbours_ids.resize(temp_size);
           temp_neighbours_delta.resize(temp_size);
           temp_neighbours_failure_id.resize(temp_size);
           temp_neighbours_contact_forces.resize(temp_size);
           temp_neighbours_mapping.resize(temp_size);
           temp_cont_neighbours_mapping.resize(temp_size);
       
           //double                ini_delta           = 0.0;
           //int                   failure_id          = 1;
           //array_1d<double, 3>   neigh_forces        (3,0.0); // **zerovector anava mes rapid
           //double                mapping_new_ini     = -1;  
           //double                mapping_new_cont    = -1;
               

           for (ParticleWeakIteratorType i = TempNeighbours.begin(); i != TempNeighbours.end(); i++)
        
           {

            double                ini_delta           = 0.0;
           int                   failure_id          = 1;
           array_1d<double, 3>   neigh_forces        (3,0.0); // **zerovector anava mes rapid
           double                mapping_new_ini     = -1;  
           double                mapping_new_cont    = -1;

             //Loop Over Initial Neighbours

             for (unsigned int k = 0; k != mIniNeighbourIds.size(); k++)
             {
                        
               if (static_cast<int>((i)->Id()) == mIniNeighbourIds[k])
               {               
              
                 ini_delta  = mIniNeighbourDelta[k];
                 failure_id = mIniNeighbourFailureId[k];
                 mapping_new_ini = k; 
                 mapping_new_cont = mIniNeighbourToIniContinuum[k];
              
                 break;
               }

             }
    
             //Judge if its neighbour
          
             double other_radius                 = i->GetGeometry()[0].FastGetSolutionStepValue(RADIUS);
             double radius_sum                   = GetRadius() + other_radius;
             array_1d<double,3> other_to_me_vect = this->GetGeometry()[0].Coordinates() - i->GetGeometry()[0].Coordinates();
             double distance                     = sqrt(other_to_me_vect[0] * other_to_me_vect[0] + other_to_me_vect[1] * other_to_me_vect[1] + other_to_me_vect[2] * other_to_me_vect[2]);
             double indentation                  = radius_sum - distance - ini_delta;
          
             if ( indentation > 0.0 || failure_id == 0 )  //WE NEED TO SET A NUMERICAL TOLERANCE FUNCTION OF THE RADIUS.  MSIMSI 10
             {
        
               //Loop Over Last time-step Neighbours
          
               for (unsigned int j = 0; j != mOldNeighbourIds.size(); j++)
               {
                 if (static_cast<int>(i->Id()) == mOldNeighbourIds[j])
                 {
                   neigh_forces = mNeighbourElasticContactForces[j];
                   break;
                 }

               }
            
               if(neighbour_counter != temp_neighbour_counter)
               {
              
                 (*(TempNeighbours.ptr_begin() + neighbour_counter)).swap( (*(TempNeighbours.ptr_begin() + temp_neighbour_counter))) ;     
                              
               }
            
               temp_neighbours_mapping[neighbour_counter]          = mapping_new_ini;
               temp_cont_neighbours_mapping[neighbour_counter]     = mapping_new_cont;
               temp_neighbours_ids[neighbour_counter]              = static_cast<int>((i)->Id());
               temp_neighbours_delta[neighbour_counter]            = ini_delta;
               temp_neighbours_failure_id[neighbour_counter]       = failure_id;
               temp_neighbours_contact_forces[neighbour_counter]   = neigh_forces;
            
               neighbour_counter++;
            
             }
            
               temp_neighbour_counter++;

           }
        
           TempNeighbours.erase(TempNeighbours.begin()+neighbour_counter, TempNeighbours.end());
     
         //  if(mMapping_New_Ini.size() != neighbour_counter)   //si en comptes de fer tot aixo fes un resize del mMapping_New_... etc... no quedaria tallada la part ke no minteressa i hagues pogut ferlo servir ampliat i ja sta. fer resize i quedarme amb lo bo.
           {
               mMapping_New_Ini.resize(neighbour_counter);
               mMapping_New_Cont.resize(neighbour_counter);
               mOldNeighbourIds.resize(neighbour_counter);
               mNeighbourDelta.resize(neighbour_counter);
               mNeighbourFailureId.resize(neighbour_counter);
               mNeighbourElasticContactForces.resize(neighbour_counter);
           }   
        
           for(unsigned int w=0; w<neighbour_counter; w++)
           {
             mMapping_New_Ini[w]           = temp_neighbours_mapping[w];
             mMapping_New_Cont[w]          = temp_cont_neighbours_mapping[w];
             mOldNeighbourIds[w]           = temp_neighbours_ids[w];
             mNeighbourDelta[w]            = temp_neighbours_delta[w];
             mNeighbourFailureId[w]        = temp_neighbours_failure_id[w];
             mNeighbourElasticContactForces[w] = temp_neighbours_contact_forces[w];
           }

         KRATOS_CATCH("")

         } //ComputeNewNeighboursHistoricalData
 
      
     */

    void SphericContinuumParticle::ComputeNewRigidFaceNeighboursHistoricalData() {

        KRATOS_TRY

        std::vector<DEMWall*> mFemTempNeighbours;
        mFemTempNeighbours.swap(mNeighbourRigidFaces);

        unsigned int fem_temp_size = mFemTempNeighbours.size();

        mNeighbourRigidFaces.clear();

        unsigned int fem_neighbour_counter = 0;

        std::vector<int> fem_temp_neighbours_ids; //these temporal vectors are very small, saving them as a member of the particle loses time (usually they consist on 1 member).
        std::vector<double> fem_temp_neighbours_delta;
        std::vector<array_1d<double, 3> > fem_temp_neighbours_contact_forces;
        std::vector<array_1d<double, 3> > fem_temp_neighbours_elastic_contact_forces;
        std::vector<int> fem_temp_neighbours_mapping;

        fem_temp_neighbours_ids.resize(fem_temp_size);
        fem_temp_neighbours_delta.resize(fem_temp_size);
        fem_temp_neighbours_contact_forces.resize(fem_temp_size);
        fem_temp_neighbours_elastic_contact_forces.resize(fem_temp_size);
        fem_temp_neighbours_mapping.resize(fem_temp_size);

        array_1d<double, 3> vector_of_zeros;
        vector_of_zeros[0] = 0.0;
        vector_of_zeros[1] = 0.0;
        vector_of_zeros[2] = 0.0;

        const std::vector<double>& RF_Pram = mNeighbourRigidFacesPram;

        for (unsigned int i = 0; i < mFemTempNeighbours.size(); i++) {

            int ino1 = i * 16;
            double DistPToB = RF_Pram[ino1 + 9];
            int iNeighborID = static_cast<int> (RF_Pram[ino1 + 14]);
            double ini_delta = 0.0;
            array_1d<double, 3> neigh_forces = vector_of_zeros;
            array_1d<double, 3> neigh_forces_elastic = vector_of_zeros;
            double mapping_new_ini = -1;

            for (unsigned int k = 0; k != mFemIniNeighbourIds.size(); k++) {
                if (iNeighborID == mFemIniNeighbourIds[k]) {
                    ini_delta = mFemIniNeighbourDelta[k];
                    mapping_new_ini = k;
                    break;
                }
            }

            for (unsigned int j = 0; j != mFemOldNeighbourIds.size(); j++) {
                if (static_cast<int> ((mFemTempNeighbours[i])->Id()) == mFemOldNeighbourIds[j]) {
                    neigh_forces_elastic = mNeighbourRigidFacesElasticContactForce[j];
                    neigh_forces = mNeighbourRigidFacesTotalContactForce[j]; 
                    break;
                }
            }

            //Judge if it is neighbor                  
            double indentation = -(DistPToB - GetRadius()) - ini_delta;

            if (indentation > 0.0) {
                mNeighbourRigidFaces.push_back(mFemTempNeighbours[i]);

                fem_temp_neighbours_ids[fem_neighbour_counter] = static_cast<int> ((mFemTempNeighbours[i])->Id());
                fem_temp_neighbours_mapping[fem_neighbour_counter] = mapping_new_ini;
                fem_temp_neighbours_delta[fem_neighbour_counter] = ini_delta;
                fem_temp_neighbours_contact_forces[fem_neighbour_counter] = neigh_forces;
                fem_temp_neighbours_elastic_contact_forces[fem_neighbour_counter] = neigh_forces_elastic;

                fem_neighbour_counter++;
            }

        }//for ConditionWeakIteratorType i

        int final_size = mNeighbourRigidFaces.size();
        fem_temp_neighbours_ids.resize(final_size);
        fem_temp_neighbours_delta.resize(final_size);
        fem_temp_neighbours_contact_forces.resize(final_size);
        fem_temp_neighbours_elastic_contact_forces.resize(final_size);
        fem_temp_neighbours_mapping.resize(final_size);

        mFemMappingNewIni.swap(fem_temp_neighbours_mapping);
        mFemOldNeighbourIds.swap(fem_temp_neighbours_ids);
        mFemNeighbourDelta.swap(fem_temp_neighbours_delta);
        mNeighbourRigidFacesElasticContactForce.swap(fem_temp_neighbours_elastic_contact_forces);
        mNeighbourRigidFacesTotalContactForce.swap(fem_temp_neighbours_contact_forces);
        mNeighbourRigidFacesPram.clear();

        KRATOS_CATCH("")
    }

    void SphericContinuumParticle::CalculateMeanContactArea(const bool has_mpi, const ProcessInfo& rCurrentProcessInfo, const bool first) {
        int my_id = this->Id();
        double my_partition_index = 0.0;
        if (has_mpi) my_partition_index = this->GetGeometry()[0].FastGetSolutionStepValue(PARTITION_INDEX);
        bool im_skin = bool(this->GetGeometry()[0].FastGetSolutionStepValue(SKIN_SPHERE));

        std::vector<SphericContinuumParticle*>& r_continuum_ini_neighbours = this->mContinuumIniNeighbourElements;

        for (unsigned int i = 0; i < r_continuum_ini_neighbours.size(); i++) {
            if (r_continuum_ini_neighbours[i] == NULL) continue; //The initial neighbor was deleted at some point in time!!
            //TODO: SHOULD WE CHECK HERE THAT BOTH BELONG TO THE SAME CONTINUUM GROUP??? It is done in other places.

            Particle_Contact_Element* bond_i = mBondElements[i];
            double other_partition_index = 0.0;
            if (has_mpi) other_partition_index = r_continuum_ini_neighbours[i]->GetGeometry()[0].FastGetSolutionStepValue(PARTITION_INDEX);

            if (first) {

                bool neigh_is_skin = bool(r_continuum_ini_neighbours[i]->GetGeometry()[0].FastGetSolutionStepValue(SKIN_SPHERE));

                int neigh_id = r_continuum_ini_neighbours[i]->Id();

                if ((im_skin && neigh_is_skin) || (!im_skin && !neigh_is_skin)) {
                    if (my_id < neigh_id) {
                        bond_i->mLocalContactAreaLow = mcont_ini_neigh_area[i];
                    }// if my id < neigh id                        
                    else {
                        if (!has_mpi) bond_i->mLocalContactAreaHigh = mcont_ini_neigh_area[i];
                        else {
                            if (other_partition_index == my_partition_index) bond_i->mLocalContactAreaHigh = mcont_ini_neigh_area[i];
                        }
                    }
                }//both skin or both inner.

                else if (!im_skin && neigh_is_skin) {//we will store both the same only coming from the inner to the skin.                    
                    if (!has_mpi) {
                        bond_i -> mLocalContactAreaHigh = mcont_ini_neigh_area[i];
                        bond_i -> mLocalContactAreaLow = mcont_ini_neigh_area[i];
                    } else {
                        if (other_partition_index == my_partition_index) {
                            bond_i -> mLocalContactAreaHigh = mcont_ini_neigh_area[i];
                            bond_i -> mLocalContactAreaLow = mcont_ini_neigh_area[i];
                        }
                    }
                } //neigh skin

            }//if(first_time)

            else {//last operation                                  
                if (!has_mpi) mcont_ini_neigh_area[i] = bond_i->mMeanContactArea;
            }

        }//loop neigh.

        return;
    }

    void SphericContinuumParticle::Calculate(const Variable<double>& rVariable, double& Output, const ProcessInfo& rCurrentProcessInfo) {

        KRATOS_TRY

        if (rVariable == DELTA_TIME) {

            double coeff = rCurrentProcessInfo[NODAL_MASS_COEFF];
            double mass = mRealMass;

            if (coeff > 1.0) {
                KRATOS_THROW_ERROR(std::runtime_error, "The coefficient assigned for virtual mass is larger than one, virtual_mass_coeff= ", coeff)
            }
            else if ((coeff == 1.0) && (rCurrentProcessInfo[VIRTUAL_MASS_OPTION])) {
                Output = 9.0E09;
            }
            else {

                if (rCurrentProcessInfo[VIRTUAL_MASS_OPTION]) {
                    mass /= 1 - coeff;
                }

                double K = GetYoung() * KRATOS_M_PI * GetRadius();

                Output = 0.34 * sqrt(mass / K);

                if (rCurrentProcessInfo[ROTATION_OPTION] == 1) {
                    Output = Output * 0.5; //factor for critical time step when rotation is allowed.
                }
            }
            return;
        }//CRITICAL DELTA CALCULATION
        ////////////////////////////////////////////////////////////////////////
        if (rVariable == PARTICLE_ROTATION_DAMP_RATIO) {
            //ApplyLocalMomentsDamping( rCurrentProcessInfo ); MSIMSI
            return;
        }
        ////////////////////////////////////////////////////////////////////////

        if (rVariable == DEM_STRESS_XX) //operations with the stress_strain tensors
        {
            SymmetrizeTensor(rCurrentProcessInfo);
            return;
        }
        ////////////////////////////////////////////////////////////////////////
        if (rVariable == DUMMY_DEBUG_DOUBLE) //Dummy variable for debugging  MSIMSI DEBUG
        {
            //CheckPairWiseBreaking();
            return;
        }
        ////////////////////////////////////////////////////////////////////////

        //        if (rVariable == CALCULATE_SET_INITIAL_DEM_CONTACTS)
        //        {
        //            SetInitialSphereContacts(rCurrentProcessInfo);

        //            CreateContinuumConstitutiveLaws();
        //            return;
        //        }

        if (rVariable == CALCULATE_SET_INITIAL_FEM_CONTACTS) {
            SetInitialFemContacts();
            return;
        }

        KRATOS_CATCH("")

    }//Calculate

    void SphericContinuumParticle::Calculate(const Variable<Vector >& rVariable, Vector& Output, const ProcessInfo& rCurrentProcessInfo) {
    }//calculate Output vector

    void SphericContinuumParticle::Calculate(const Variable<array_1d<double, 3 > >& rVariable, array_1d<double, 3 > & Output, const ProcessInfo& rCurrentProcessInfo) {
    }

    void SphericContinuumParticle::Calculate(const Variable<Matrix >& rVariable, Matrix& Output, const ProcessInfo& rCurrentProcessInfo) {
    }

    void SphericContinuumParticle::ComputeAdditionalForces(array_1d<double, 3>& additionally_applied_force,
            array_1d<double, 3>& additionally_applied_moment,
            ProcessInfo& rCurrentProcessInfo,
            const array_1d<double, 3>& gravity) {

        //const array_1d<double,3>& gravity         = rCurrentProcessInfo[GRAVITY];

        KRATOS_TRY

        if (rCurrentProcessInfo[TRIAXIAL_TEST_OPTION] && *mSkinSphere) //could be applied to selected particles.
        {
            ComputePressureForces(additionally_applied_force, rCurrentProcessInfo);
        }

        //if( mRotationOption != 0 && mRotationSpringOption != 0 )
        /*if( this->Is(DEMFlags::HAS_ROTATION) && this->Is(DEMFlags::HAS_ROTATION_SPRING)  )             
        {
            //ComputeParticleRotationSpring(); MSI: #C2
        }*/

        double mass = mRealMass;
        additionally_applied_force[0] += mass * gravity[0];
        additionally_applied_force[1] += mass * gravity[1];
        additionally_applied_force[2] += mass * gravity[2];

        KRATOS_CATCH("")
    }

    void SphericContinuumParticle::CustomInitialize() {
        distances_squared = 0.0;           
    }

    double SphericContinuumParticle::GetInitialDeltaWithFEM(int index) {
        return mFemNeighbourDelta[index];
    }

    void SphericContinuumParticle::CalculateOnContactElements(unsigned int neighbour_iterator_id,
            size_t i_neighbour_count,
            int mapping_new_cont,
            double LocalElasticContactForce[3],
            double contact_sigma,
            double contact_tau,
            double failure_criterion_state,
            double acumulated_damage,
            int time_steps) {
        KRATOS_TRY

        Particle_Contact_Element* bond = mBondElements[mapping_new_cont];
        if (bond == NULL) return; //This bond was never created (happens in some MPI cases, see CreateContactElements() in explicit_solve_continumm.h)

        bond->mLocalContactForce[0] = LocalElasticContactForce[0];
        bond->mLocalContactForce[1] = LocalElasticContactForce[1];
        bond->mLocalContactForce[2] = LocalElasticContactForce[2];
        bond->mContactSigma = contact_sigma;
        bond->mContactTau = contact_tau;
        bond->mContactFailure = (mNeighbourFailureId[i_neighbour_count]);
        bond->mFailureCriterionState = failure_criterion_state;

        if ((time_steps == 0) || (acumulated_damage > bond->mUnidimendionalDamage)) {
            bond->mUnidimendionalDamage = acumulated_damage;
        }
        // if Target Id < Neigh Id        

        KRATOS_CATCH("")

    }//CalculateOnContactElements                              
    

    void SphericContinuumParticle::AddPoissonContribution(const double equiv_poisson, double LocalCoordSystem[3][3], double& normal_force, double calculation_area) {

        double force[3];

        for (int i = 0; i < 3; i++) force[i] = (*mSymmStressTensor)(i,0) * LocalCoordSystem[0][0] + (*mSymmStressTensor)(i,1) * LocalCoordSystem[0][1] + (*mSymmStressTensor)(i,2) * LocalCoordSystem[0][2]; //StressTensor*unitaryNormal0

        double sigma_x = force[0] * LocalCoordSystem[0][0] + force[1] * LocalCoordSystem[0][1] + force[2] * LocalCoordSystem[0][2]; // projection to normal to obtain value of the normal stress

        for (int i = 0; i < 3; i++) force[i] = (*mSymmStressTensor)(i,0) * LocalCoordSystem[1][0] + (*mSymmStressTensor)(i,1) * LocalCoordSystem[1][1] + (*mSymmStressTensor)(i,2) * LocalCoordSystem[1][2]; //StressTensor*unitaryNormal1

        double sigma_y = force[0] * LocalCoordSystem[1][0] + force[1] * LocalCoordSystem[1][1] + force[2] * LocalCoordSystem[1][2]; // projection to normal to obtain value of the normal stress

        double poisson_force = calculation_area * equiv_poisson * (sigma_x + sigma_y);

        normal_force -= poisson_force;
    }

    void SphericContinuumParticle::ComputePressureForces(array_1d<double, 3>& externally_applied_force, ProcessInfo& rCurrentProcessInfo) {

        externally_applied_force = this->GetGeometry()[0].FastGetSolutionStepValue(EXTERNAL_APPLIED_FORCE);

        /*
         double time_now = rCurrentProcessInfo[TIME]; //MSIMSI 1 I tried to do a *mpTIME

         if( mFinalPressureTime <= 1e-10 )
         {
          
           //  KRATOS_WATCH("WARNING: SIMULATION TIME TO CLOSE TO ZERO")
          
         }
        
         else if ( (mFinalPressureTime > 1e-10) && (time_now < mFinalPressureTime) )
         {
  
          
           externally_applied_force = AuxiliaryFunctions::LinearTimeIncreasingFunction(total_externally_applied_force, time_now, mFinalPressureTime);

         }       
         else
         {

           externally_applied_force = total_externally_applied_force;
          
         }
       
         */

    } //SphericContinuumParticle::ComputePressureForces

    void SphericContinuumParticle::ContinuumSphereMemberDeclarationFirstStep(const ProcessInfo& rCurrentProcessInfo) {

        KRATOS_TRY

        mSkinSphere                     = &(this->GetGeometry()[0].FastGetSolutionStepValue(SKIN_SPHERE));
        mContinuumGroup                 = this->GetGeometry()[0].FastGetSolutionStepValue(COHESIVE_GROUP);
        
        KRATOS_CATCH("")

    } //ContinuumSphereMemberDeclarationFirstStep

} // namespace Kratos.


//NOTE::
/*
 * #1: Here, Initial_delta is expected to be positive if it is embedding and negative if there's a separation.
 * #2: 0.25 is because we take only the half of the equivalent radius, corresponding to the case of one sphere with radius Requivalent and other = radius 0.
 * #3: For detached particles we enter only if the indentation is > 0. For attached particles we enter only if the particle is still attached.
 * #4: we use incremental calculation. YADE develops a complicated "absolute method"
 * #5: We only store in the initial neighbors array the ones that are cohesive or the ones that have positive or negative initial indentation. In other words,
 *     the non-cohesive ones with 0 indentation (<some tolerance) don't have to be stored since we can treat it indistinctly from other new neighbours that the particle in stury would meet.
 */

//#C2: ComputeParticleRotationSpring;

//           void SphericContinuumParticle::ComputeParticleRotationSpring() // shan de corregir areees etc etc etc...
//       {
//         //double dt                           = rCurrentProcessInfo[DELTA_TIME]; //C.F.: new
//         /*
//                     c=objecte_contacte(particula,vei)
// 
//             força=(c.RADI)*3;  //M: idea: to create a class contact, create objects of contacts with all the parameters. easier...
//                                 /no puc amb MPI oi? pk hauria de passar punters...
//           */
//           ParticleWeakVectorType& mrNeighbours                = this->GetValue(NEIGHBOUR_ELEMENTS);
// 
//           array_1d<double, 3 > & mRota_Moment = GetGeometry()[0].FastGetSolutionStepValue(PARTICLE_MOMENT);
// 
// 
//           Vector & mRotaSpringFailureType  = this->GetValue(PARTICLE_ROTATE_SPRING_FAILURE_TYPE);
// 
//           size_t i_neighbour_count = 0;
// 
//           for(ParticleWeakIteratorType ineighbour = mrNeighbours.begin(); ineighbour != mrNeighbours.end(); ineighbour++)
//           {
// 
//               //if(mIfInitalContact[i_neighbour_count] == 1 && mRotaSpringFailureType[i_neighbour_count] == 0) ///M.S:NEWWWW,  if THE SPRING BRAKES... NO MORE CONTRIBUION.
//               //if( mRotaSpringFailureType[i_neighbour_count] == 0) //M.S: CAL FICAR A INITIALIZE QUE SIGUI 1 I DESPRES INITIAL CONTACTS POSAR 0 SI NECESITEN, IGUAL QUE FAILURE NORMAL.
//               //mmm.. what about the other failure types? if a contact is broken due to shear or tensile, it cant be a bending
//               {
//                 
//                 
//                   array_1d<double, 3 > & mRotaSpringMoment  = this->GetValue(PARTICLE_ROTATE_SPRING_MOMENT)[ i_neighbour_count ];
// 
//                   double other_radius    = ineighbour->GetGeometry()[0].FastGetSolutionStepValue(RADIUS);
//                   double other_young     = ineighbour->GetGeometry()[0].FastGetSolutionStepValue(YOUNG_MODULUS);
//                   double other_poisson   = ineighbour->GetGeometry()[0].FastGetSolutionStepValue(POISSON_RATIO);
//                   double other_tension   = ineighbour->GetGeometry()[0].FastGetSolutionStepValue(PARTICLE_TENSION); //MSI: these variables are eliminated.
//                   double other_cohesion  = ineighbour->GetGeometry()[0].FastGetSolutionStepValue(PARTICLE_COHESION);
//                   double other_inertia   = ineighbour->GetGeometry()[0].FastGetSolutionStepValue(PARTICLE_INERTIA);
// 
//                   double equiv_tension  = (mTension  + other_tension ) * 0.5;
//                   double equiv_cohesion = (mCohesion + other_cohesion) * 0.5;
// 
//                   double equiv_radius     = (GetRadius() + other_radius) * 0.5 ;
//                   double equiv_area       = KRATOS_M_PI * equiv_radius * equiv_radius;
//                   double equiv_poisson    = (mPoisson + other_poisson) * 0.5 ;
//                   double equiv_young      = (mYoung  + other_young)  * 0.5;
// 
//                   double kn_el               = equiv_young * equiv_area / (2.0 * equiv_radius);
//                   double ks               = kn_el / (2.0 * (1.0 + equiv_poisson));
// 
//                   array_1d<double,3> other_to_me_vect = GetGeometry()[0].Coordinates() - ineighbour->GetGeometry()[0].Coordinates();
// 
//                   /////Cfeng: Forming the Local Contact Coordinate system
//                   double NormalDir[3]           = {0.0};
//                   double LocalCoordSystem[3][3] = {{0.0}, {0.0}, {0.0}};
//                   NormalDir[0] = other_to_me_vect[0];
//                   NormalDir[1] = other_to_me_vect[1];
//                   NormalDir[2] = other_to_me_vect[2];
//                   GeometryFunctions::ComputeContactLocalCoordSystem(NormalDir, LocalCoordSystem);
// 
//                   double LocalRotaSpringMoment[3]     = {0.0};
//                   double GlobalRotaSpringMoment[3]    = {0.0};
//                   double GlobalRotaSpringMomentOld[3] = {0.0};
// 
//                   double DeltRotaDisp[3] = {0.0};
//                   //double DeltRotaDisp2[3] = {0.0};
// 
//                   double LocalDeltRotaDisp[3] = {0.0};
// 
//                   double TargetDeltRotaDist[3] = {0.0};
//                   double NeighbourDeltRotaDist[3] = {0.0};
//                 
//                   for (int i=0;i<3;i++)
//                   {
//                       TargetDeltRotaDist[i] = this->GetGeometry()[0].FastGetSolutionStepValue(DELTA_ROTA_DISPLACEMENT)[i];
//                       NeighbourDeltRotaDist[i] = ineighbour->GetGeometry()[0].FastGetSolutionStepValue(DELTA_ROTA_DISPLACEMENT)[i];
//                       DeltRotaDisp[i] =  - ( TargetDeltRotaDist[i] - NeighbourDeltRotaDist[i] );
//                   }
//               
//                   GeometryFunctions::VectorGlobal2Local(LocalCoordSystem, DeltRotaDisp, LocalDeltRotaDisp);
// 
//                   GlobalRotaSpringMomentOld[0] = mRotaSpringMoment[ 0 ];
// 
//                   GlobalRotaSpringMomentOld[1] = mRotaSpringMoment[ 1 ];
// 
//                   GlobalRotaSpringMomentOld[2] = mRotaSpringMoment[ 2 ];
// 
//                   GeometryFunctions::VectorGlobal2Local(LocalCoordSystem, GlobalRotaSpringMomentOld, LocalRotaSpringMoment);
// 
// 
//                   double Inertia_I = (mSectionalInertia + other_inertia) * 0.5;
// 
//                   double Inertia_J = Inertia_I * 2.0;
// 
// 
//                   LocalRotaSpringMoment[0] +=  - Inertia_I * LocalDeltRotaDisp[0] * kn_el / equiv_area;
// 
//                   LocalRotaSpringMoment[1] +=  - Inertia_I * LocalDeltRotaDisp[1] * kn_el / equiv_area;
// 
//                   LocalRotaSpringMoment[2] +=  - Inertia_J * LocalDeltRotaDisp[2] * ks / equiv_area;
// 
// 
//                   ////Judge if the rotate spring is broken or not
//                   double GlobalElasticContactForce[3]  = {0.0};
//                   double LocalElasticContactForce [3]  = {0.0};
// 
//                   GlobalElasticContactForce[0] = this->GetValue(PARTICLE_CONTACT_FORCES)[ i_neighbour_count ][ 0 ];
//                   GlobalElasticContactForce[1] = this->GetValue(PARTICLE_CONTACT_FORCES)[ i_neighbour_count ][ 1 ];
//                   GlobalElasticContactForce[2] = this->GetValue(PARTICLE_CONTACT_FORCES)[ i_neighbour_count ][ 2 ]; //GlobalElasticContactForce[2] = mContactForces[3 * i_neighbour_count  + 2 ];
//                   GeometryFunctions::VectorGlobal2Local(LocalCoordSystem, GlobalElasticContactForce, LocalElasticContactForce);
// 
//                   double ForceN  = LocalElasticContactForce[2];
//                   double ForceS  = sqrt( LocalElasticContactForce[0] * LocalElasticContactForce[0] + LocalElasticContactForce[1] * LocalElasticContactForce[1]);
//                   double MomentS = sqrt(LocalRotaSpringMoment[0] * LocalRotaSpringMoment[0] + LocalRotaSpringMoment[1] * LocalRotaSpringMoment[1]);
//                   double MomentN = LocalRotaSpringMoment[2];
// 
//                   //////bending stress and axial stress add together, use edge of the bar will failure first
//                   double TensiMax = -ForceN / equiv_area + MomentS        / Inertia_I * equiv_radius;
//                   double ShearMax = ForceS  / equiv_area + fabs(MomentN)  / Inertia_J * equiv_radius;
// 
//                   if(TensiMax > equiv_tension || ShearMax > equiv_cohesion)
//                   {
//                       mRotaSpringFailureType[i_neighbour_count] = 1;
// 
//                       LocalRotaSpringMoment[0] = 0.0;
//                       LocalRotaSpringMoment[1] = 0.0;
//                       LocalRotaSpringMoment[2] = 0.0;
//                   }
// 
//                   GeometryFunctions::VectorLocal2Global(LocalCoordSystem, LocalRotaSpringMoment, GlobalRotaSpringMoment);
// 
//                   mRotaSpringMoment[ 0 ] = GlobalRotaSpringMoment[0];
//                   mRotaSpringMoment[ 1 ] = GlobalRotaSpringMoment[1];
//                   mRotaSpringMoment[ 2 ] = GlobalRotaSpringMoment[2];
// 
//                   ////feedback, contact moment----induce by rotation spring
//                   mRota_Moment[0] -= GlobalRotaSpringMoment[0];
//                   mRota_Moment[1] -= GlobalRotaSpringMoment[1];
//                   mRota_Moment[2] -= GlobalRotaSpringMoment[2];
//               }
// 
//               i_neighbour_count++;
//           }
//       }//ComputeParticleRotationSpring
//      
/*         
 * 
 * 
 */

// #C3: InitializeContactElements : aquesta funcio començava abans de evaluatedeltadisplacement per poisson etc pero no crec ke faci falta ara.
