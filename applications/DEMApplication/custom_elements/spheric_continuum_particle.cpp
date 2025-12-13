//
// Authors: Miquel Santasusana msantasusana@cimne.upc.edu
//          Miguel Angel Celigueta maceli@cimne.upc.edu
//

#include "spheric_continuum_particle.h"
#include <cmath>
#include <random>

namespace Kratos {

    SphericContinuumParticle::SphericContinuumParticle():SphericParticle() {
        mContinuumInitialNeighborsSize = 0;
        mInitialNeighborsSize = 0;
    }

    SphericContinuumParticle::SphericContinuumParticle(IndexType NewId, GeometryType::Pointer pGeometry) : SphericParticle(NewId, pGeometry){
        mContinuumInitialNeighborsSize = 0;
        mInitialNeighborsSize = 0;
    }

    SphericContinuumParticle::SphericContinuumParticle(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
    : SphericParticle(NewId, pGeometry, pProperties){
        mContinuumInitialNeighborsSize = 0;
        mInitialNeighborsSize = 0;
    }

    SphericContinuumParticle::SphericContinuumParticle(IndexType NewId, NodesArrayType const& ThisNodes)
    : SphericParticle(NewId, ThisNodes){
        mContinuumInitialNeighborsSize = 0;
        mInitialNeighborsSize = 0;
    }

    Element::Pointer SphericContinuumParticle::Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const {
        return SphericParticle::Pointer(new SphericContinuumParticle(NewId, GetGeometry().Create(ThisNodes), pProperties));
    }

    /// Destructor

    SphericContinuumParticle::~SphericContinuumParticle() {
    }
    
    void SphericContinuumParticle::BeforeSetInitialSphereContacts(const ProcessInfo& r_process_info, std::set<std::pair<int, int>>& CementedContactPairsLocal) {

        for (unsigned int i = 0; i < mNeighbourElements.size(); i++) {
            SphericContinuumParticle* neighbour_iterator = dynamic_cast<SphericContinuumParticle*>(mNeighbourElements[i]);
            int r_other_continuum_group = neighbour_iterator->mContinuumGroup; // finding out neighbor's Continuum Group Id
            if ((r_other_continuum_group == this->mContinuumGroup) && (this->mContinuumGroup != 0)) {
                if (this->Id() < neighbour_iterator->Id()){
                    CementedContactPairsLocal.insert(makeKey(this->Id(), neighbour_iterator->Id()));
                }
            } 
        }
    }//BeforeSetInitialSphereContacts
    
    void SphericContinuumParticle::SetInitialSphereContacts(const ProcessInfo& r_process_info, std::map<std::pair<int, int>, double>& CementedContactAreasMapLocal) {

        std::vector<SphericContinuumParticle*> ContinuumInitialNeighborsElements;
        std::vector<SphericContinuumParticle*> DiscontinuumInitialNeighborsElements;
        std::vector<int> DiscontinuumInitialNeighborsIds;
        //std::vector<double> DiscontinuumInitialNeighborsDeltas;
        mIniNeighbourFailureId.clear(); // We will have to build this vector, we still don't know its size, it applies only to continuum particles
        size_t continuum_ini_size    = 0;
        size_t discontinuum_ini_size = 0;
        unsigned int neighbours_size = mNeighbourElements.size();
        mIniNeighbourIds.resize(neighbours_size);
        //mIniNeighbourDelta.resize(neighbours_size);

        for (unsigned int i = 0; i < mNeighbourElements.size(); i++) {

            SphericContinuumParticle* neighbour_iterator = dynamic_cast<SphericContinuumParticle*>(mNeighbourElements[i]);
            auto& central_node = GetGeometry()[0];
            auto& neighbour_node = neighbour_iterator->GetGeometry()[0];
            
            array_1d<double, 3> other_to_me_vect;
            if (!r_process_info[DOMAIN_IS_PERIODIC]){ // default infinite-domain case
                noalias(other_to_me_vect) = central_node.Coordinates() - neighbour_node.Coordinates();
            } else { // periodic domain
                double my_coors[3] = {central_node[0], central_node[1], central_node[2]};
                double other_coors[3] = {neighbour_node[0], neighbour_node[1], neighbour_node[2]};
                TransformNeighbourCoorsToClosestInPeriodicDomain(r_process_info, my_coors, other_coors);
                other_to_me_vect[0] = my_coors[0] - other_coors[0];
                other_to_me_vect[1] = my_coors[1] - other_coors[1];
                other_to_me_vect[2] = my_coors[2] - other_coors[2];
            }

            double distance = DEM_MODULUS_3(other_to_me_vect);
            double radius_sum = GetRadius() + neighbour_iterator->GetRadius();
            double initial_delta = radius_sum - distance;

            bool set_continuum_contact = false;
            if (mCementedContactPairsSetPtr->find(makeKey(this->Id(), neighbour_iterator->Id())) != mCementedContactPairsSetPtr->end()) {
                set_continuum_contact = true;
            }
            
            int r_other_continuum_group = neighbour_iterator->mContinuumGroup; // finding out neighbor's Continuum Group Id
            if ((r_other_continuum_group  == this->mContinuumGroup) && (this->mContinuumGroup != 0) && set_continuum_contact) {
                mIniNeighbourIds[continuum_ini_size] = neighbour_iterator->Id();
                //mIniNeighbourDelta[continuum_ini_size] = initial_delta;
                mIniNeighbourDelta[static_cast<int>(neighbour_iterator->Id())] = initial_delta;
                mIniNeighbourFailureId.push_back(0);
                array_1d<double, 3> vector_of_zeros(3,0.0);
                //mArrayOfOldDeltaDisplacements.push_back(vector_of_zeros);
                //mArrayOfDeltaDisplacements.push_back(vector_of_zeros);
                ContinuumInitialNeighborsElements.push_back(neighbour_iterator);
                continuum_ini_size++;

                //calculate the initial contact area
                double ini_bond_contact_area = 0.0;
                CalculateInitialBondContactArea(distance, GetSearchRadius(), 
                                                neighbour_iterator->GetSearchRadius(), 
                                                GetRadius(), 
                                                neighbour_iterator->GetRadius(), 
                                                ini_bond_contact_area,
                                                r_process_info);
                //mIniBondContactArea[static_cast<int>(neighbour_iterator->Id())] = ini_bond_contact_area;
                if (this->Id() < neighbour_iterator->Id()){
                    CementedContactAreasMapLocal[makeKey(this->Id(), neighbour_iterator->Id())] = ini_bond_contact_area;
                }
                double V_bond = 0.0;
                double R_bond = std::sqrt(ini_bond_contact_area / Globals::Pi);
                CalculateBondVolume(distance, GetRadius(), neighbour_iterator->GetRadius(), R_bond, V_bond);
                mBondVolume[static_cast<int>(neighbour_iterator->Id())] = V_bond;
            } else {
                DiscontinuumInitialNeighborsIds.push_back(neighbour_iterator->Id());
                //DiscontinuumInitialNeighborsDeltas.push_back(initial_delta);
                mIniNeighbourDelta[static_cast<int>(neighbour_iterator->Id())] = initial_delta;
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
            //mIniNeighbourDelta[continuum_ini_size + k] = DiscontinuumInitialNeighborsDeltas[k];
            mNeighbourElements[continuum_ini_size + k] = DiscontinuumInitialNeighborsElements[k];
        }

        GetGeometry()[0].FastGetSolutionStepValue(INITIAL_BOND_NUMBER) = mContinuumInitialNeighborsSize;
    }//SetInitialSphereContacts

    void SphericContinuumParticle::CreateContinuumConstitutiveLaws() {

        mContinuumConstitutiveLawArray.resize(mContinuumInitialNeighborsSize);

        for (unsigned int i = 0; i < mContinuumInitialNeighborsSize; i++) {
            Properties::Pointer properties_of_this_contact = GetProperties().pGetSubProperties(mNeighbourElements[i]->GetProperties().Id());
            mContinuumConstitutiveLawArray[i] = (*properties_of_this_contact)[DEM_CONTINUUM_CONSTITUTIVE_LAW_POINTER]-> Clone();
            SphericContinuumParticle* p_cont_neighbour_particle = dynamic_cast<SphericContinuumParticle*>(mNeighbourElements[i]);
            mContinuumConstitutiveLawArray[i]->Initialize(this, p_cont_neighbour_particle, properties_of_this_contact);
        }
    }

    void SphericContinuumParticle::SetInitialFemContacts() {

        std::vector<DEMWall*>& rFemNeighbours = this->mNeighbourRigidFaces;

        unsigned int fem_neighbours_size = rFemNeighbours.size();

        mFemIniNeighbourIds.resize(fem_neighbours_size);
        mFemIniNeighbourDelta.resize(fem_neighbours_size);
        mContactConditionWeights.resize(fem_neighbours_size);

        for (unsigned int i = 0; i < rFemNeighbours.size(); i++) {

            double LocalCoordSystem[3][3] = {{0.0}, {0.0}, {0.0}};
            array_1d<double, 3> wall_delta_disp_at_contact_point = ZeroVector(3);
            array_1d<double, 3> wall_velocity_at_contact_point = ZeroVector(3);
            double DistPToB = 0.0;
            int ContactType = -1;
            array_1d<double, 4>& Weight = this->mContactConditionWeights[i];

            rFemNeighbours[i]->ComputeConditionRelativeData(i, this, LocalCoordSystem, DistPToB, Weight, wall_delta_disp_at_contact_point, wall_velocity_at_contact_point, ContactType);

            double initial_delta = -(DistPToB - GetRadius());

            mFemIniNeighbourIds[i] = rFemNeighbours[i]->Id();
            mFemIniNeighbourDelta[i] = initial_delta;
        }
    }//SetInitialFemContacts

    void SphericContinuumParticle::ContactAreaWeighting(const ProcessInfo& r_process_info) //MISMI 10: POOYAN this could be done by calculating on the bars. not looking at the neighbors of my neighbors.
    {
        double alpha = 1.0;
        //double external_sphere_area = 4 * Globals::Pi * GetRadius()*GetRadius();
        double effectiveVolumeRadius = EffectiveVolumeRadius(r_process_info);  //calculateEffectiveVolumeRadius
        double external_sphere_area = 4 * Globals::Pi * effectiveVolumeRadius * effectiveVolumeRadius;
        double total_equiv_area = 0.0;
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
                    outputfile <<" "<<cont_ini_neighbours_size<<" "<<external_sphere_area<<" "<<total_equiv_area<<" "<< alpha << "\n";
                    outputfile.close();
                }
                for (unsigned int i = 0; i < cont_ini_neigh_area.size(); i++) {
                    cont_ini_neigh_area[i] = alpha * cont_ini_neigh_area[i];
                } //for every neighbor
            }

            else {//skin sphere
                for (unsigned int i = 0; i < cont_ini_neigh_area.size(); i++) {
                    alpha = 1.00 * (1.40727)*(external_sphere_area / total_equiv_area)*((double(cont_ini_neighbours_size)) / 11.0);

                    cont_ini_neigh_area[i] = alpha * cont_ini_neigh_area[i];
                } //for every neighbor
            }
        }//if more than 6 neighbors
    }

    void SphericContinuumParticle::CalculateInitialBondContactArea(const double distance,
                                                                const double my_search_radius,
                                                                const double other_search_radius,
                                                                const double my_radius,
                                                                const double other_radius,
                                                                double& bond_contact_area,
                                                                const ProcessInfo& r_process_info) 
    {
        const double dis_squared = distance * distance;
        const double my_search_r_squared = my_search_radius * my_search_radius;
        const double other_search_r_squared = other_search_radius * other_search_radius;
        const double squared_dis_squared = dis_squared * dis_squared;

        double numerator = 2.0 * dis_squared * (my_search_r_squared + other_search_r_squared) - std::pow(my_search_r_squared - other_search_r_squared, 2) - squared_dis_squared;
        double Rbond = 0.0;

        if (numerator > 0.0) {
            Rbond = std::sqrt(numerator) / (2.0 * distance);
            Rbond = std::min(Rbond, std::min(my_radius, other_radius));
        } else {
            Rbond = std::min(my_radius, other_radius);  // fallback if square root argument is negative
        }

        bond_contact_area = Globals::Pi * Rbond * Rbond;

        bool adjust_bond_contact_area = r_process_info[ADJUST_BOND_CONTACT_AREA_OPTION];
        std::string bond_contact_area_type = r_process_info[BOND_CONTACT_AREA_DISTRIBUTION_TYPE];
        if (adjust_bond_contact_area) {
            
            if (bond_contact_area_type == "lognormal") {
                
                const double bond_contact_area_lognormal_median = r_process_info[BOND_CONTACT_AREA_LOGNORMAL_MEDIAN];
                const double bond_contact_area_lognormal_std_dev = r_process_info[BOND_CONTACT_AREA_LOGNORMAL_STD_DEV];
                const double bond_contact_area_upper_bound = r_process_info[BOND_CONTACT_AREA_UPPER_BOUND]; // This is the upper bound for the bond contact area

                std::random_device rd;
                std::mt19937 gen(rd());

                std::lognormal_distribution<> lognorm(bond_contact_area_lognormal_median, bond_contact_area_lognormal_std_dev);
                do {
                    bond_contact_area = lognorm(gen) * 1e-12; // Convert to m^2
                } while (bond_contact_area > bond_contact_area_upper_bound);

            } else if (bond_contact_area_type == "constant") {
                const double bond_contact_area_constant_value = r_process_info[BOND_CONTACT_AREA_CONSTANT_VALUE];
                bond_contact_area = bond_contact_area_constant_value;
            }
        }
    }

    void SphericContinuumParticle::CalculateBondVolume(const double distance,
                                                    const double my_radius,
                                                    const double other_radius,
                                                    const double R_bond,
                                                    double& V_bond) 
    {

        // Step 1: h'_i, h'_j, h_bond
        double hi_prime = my_radius - std::sqrt(std::max(0.0, my_radius * my_radius - R_bond * R_bond));
        double hj_prime = other_radius - std::sqrt(std::max(0.0, other_radius * other_radius - R_bond * R_bond));
        double h_bond = distance - (my_radius + other_radius) + hi_prime + hj_prime;

        // Step 2: V_cylinder
        double V_cylinder = Globals::Pi * R_bond * R_bond * h_bond;

        // Step 3: V_ball_cylinder_intersect
        double V_ball_cylinder_intersect = Globals::Pi * hi_prime * hi_prime * (my_radius - hi_prime / 3.0)
                        + Globals::Pi * hj_prime * hj_prime * (other_radius - hj_prime / 3.0);

        // Step 4: V_ball_ball_intersect
        double V_ball_ball_intersect = 0.0;
        if (distance < (my_radius + other_radius)) {
            double distance_squared = distance * distance;
            double cos_alpha = (other_radius * other_radius + distance_squared - my_radius * my_radius) / (2.0 * other_radius * distance);
            double cos_beta  = (my_radius * my_radius + distance_squared - other_radius * other_radius) / (2.0 * my_radius * distance);
            cos_alpha = std::clamp(cos_alpha, -1.0, 1.0);
            cos_beta = std::clamp(cos_beta, -1.0, 1.0);

            double hi = other_radius * (1.0 - cos_alpha);
            double hj = my_radius * (1.0 - cos_beta);

            V_ball_ball_intersect = (Globals::Pi / 3.0) * (3 * other_radius - hi) * hi * hi + (Globals::Pi / 3.0) * (3 * my_radius - hj) * hj * hj;
        }
        
        // Step 5: Final V_bond
        V_bond = V_cylinder - V_ball_cylinder_intersect + V_ball_ball_intersect;
    }

    double SphericContinuumParticle::EffectiveVolumeRadius(const ProcessInfo& r_process_info) {

        int cont_ini_neighbours_size = mContinuumInitialNeighborsSize;

        double effectiveVolumeRadiusSum = 0.0;

        for (int i = 0; i < cont_ini_neighbours_size; i++) {

            SphericContinuumParticle* neighbour_iterator = dynamic_cast<SphericContinuumParticle*>(mNeighbourElements[i]);
            double other_radius = neighbour_iterator->GetRadius();
            auto& central_node = GetGeometry()[0];
            auto& neighbour_node = neighbour_iterator->GetGeometry()[0];
            
            array_1d<double, 3> other_to_me_vect;
            if (!r_process_info[DOMAIN_IS_PERIODIC]){ // default infinite-domain case
                noalias(other_to_me_vect) = central_node.Coordinates() - neighbour_node.Coordinates();
            } else { // periodic domain
                double my_coors[3] = {central_node[0], central_node[1], central_node[2]};
                double other_coors[3] = {neighbour_node[0], neighbour_node[1], neighbour_node[2]};
                TransformNeighbourCoorsToClosestInPeriodicDomain(r_process_info, my_coors, other_coors);
                other_to_me_vect[0] = my_coors[0] - other_coors[0];
                other_to_me_vect[1] = my_coors[1] - other_coors[1];
                other_to_me_vect[2] = my_coors[2] - other_coors[2];
            }
            double distance = DEM_MODULUS_3(other_to_me_vect);
            effectiveVolumeRadiusSum += 0.5 * (distance + GetRadius() - other_radius);
        }

        double effectiveVolumeRadius = effectiveVolumeRadiusSum / cont_ini_neighbours_size;

        return effectiveVolumeRadius;
    }

    void SphericContinuumParticle::ComputeBallToBallContactForceAndMoment(SphericParticle::ParticleDataBuffer & data_buffer,
                                                                const ProcessInfo& r_process_info,
                                                                array_1d<double, 3>& rElasticForce,
                                                                array_1d<double, 3>& rContactForce)
    {
        KRATOS_TRY

        NodeType& this_node = this->GetGeometry()[0];
        DEM_COPY_SECOND_TO_FIRST_3(data_buffer.mMyCoors, this_node)

        const int time_steps = r_process_info[TIME_STEPS];

        const array_1d<double, 3>& vel         = this->GetGeometry()[0].FastGetSolutionStepValue(VELOCITY);
        const array_1d<double, 3>& delta_displ = this->GetGeometry()[0].FastGetSolutionStepValue(DELTA_DISPLACEMENT);
        const array_1d<double, 3>& ang_vel     = this->GetGeometry()[0].FastGetSolutionStepValue(ANGULAR_VELOCITY);
        Vector& cont_ini_neigh_area            = this->GetValue(NEIGHBOURS_CONTACT_AREAS);
        int NeighbourSize = mNeighbourElements.size();
        GetGeometry()[0].GetSolutionStepValue(NEIGHBOUR_SIZE) = NeighbourSize;

        for (int i = 0; data_buffer.SetNextNeighbourOrExit(i); ++i) {

            if (mNeighbourElements[i] == NULL) continue;
            if (this->Is(NEW_ENTITY) && mNeighbourElements[i]->Is(NEW_ENTITY)) continue;
            SphericContinuumParticle* neighbour_iterator = dynamic_cast<SphericContinuumParticle*>(mNeighbourElements[i]);
            data_buffer.mpOtherParticle = neighbour_iterator;

            unsigned int neighbour_iterator_id = data_buffer.mpOtherParticle->Id();

            auto& central_node = GetGeometry()[0];
            auto& neighbour_node = neighbour_iterator->GetGeometry()[0];
            
            array_1d<double, 3> other_to_me_vect;
            if (!r_process_info[DOMAIN_IS_PERIODIC]){ // default infinite-domain case
                noalias(data_buffer.mOtherToMeVector) = central_node.Coordinates() - neighbour_node.Coordinates();
            } else { // periodic domain
                double my_coors[3] = {central_node[0], central_node[1], central_node[2]};
                double other_coors[3] = {neighbour_node[0], neighbour_node[1], neighbour_node[2]};
                TransformNeighbourCoorsToClosestInPeriodicDomain(r_process_info, my_coors, other_coors);
                data_buffer.mOtherToMeVector[0] = my_coors[0] - other_coors[0];
                data_buffer.mOtherToMeVector[1] = my_coors[1] - other_coors[1];
                data_buffer.mOtherToMeVector[2] = my_coors[2] - other_coors[2];
            }
            
            const double& other_radius = data_buffer.mpOtherParticle->GetRadius();

            data_buffer.mDistance = DEM_MODULUS_3(data_buffer.mOtherToMeVector);
            double radius_sum = GetRadius() + other_radius;

            double initial_delta = GetInitialDelta(neighbour_iterator_id);         
            
            double initial_dist = radius_sum - initial_delta;
            double indentation = initial_dist - data_buffer.mDistance;
            double indentation_particle = radius_sum - data_buffer.mDistance;
            
            if (this->mIndentationInitialOption) {

                if (data_buffer.mTime < (0.0 + 2.0 * data_buffer.mDt)){
                    if (indentation_particle > 0.0) {
                        this->mIndentationInitial[data_buffer.mpOtherParticle->Id()] = indentation_particle;
                    } 
                }

                if (this->mIndentationInitial.find(data_buffer.mpOtherParticle->Id()) != this->mIndentationInitial.end()){
                    indentation_particle -= this->mIndentationInitial[data_buffer.mpOtherParticle->Id()];
                } 

            }

            double myYoung = GetYoung();
            double myPoisson = GetPoisson();

            double kn_el = 0.0;
            double kt_el = 0.0;
            double DeltDisp[3] = {0.0};
            double RelVel[3] = {0.0};
            DEM_SET_COMPONENTS_TO_ZERO_3x3(data_buffer.mLocalCoordSystem)
            DEM_SET_COMPONENTS_TO_ZERO_3x3(data_buffer.mOldLocalCoordSystem)
            bool sliding = false;

            double contact_tau = 0.0;
            double contact_sigma = 0.0;
            double failure_criterion_state = 0.0;
            double acumulated_damage = 0.0;

            // Getting neighbor properties
            double other_young = data_buffer.mpOtherParticle->GetYoung();
            double other_poisson = data_buffer.mpOtherParticle->GetPoisson();
            double equiv_poisson;
            if ((myPoisson + other_poisson) != 0.0) { equiv_poisson = 2.0 * myPoisson * other_poisson / (myPoisson + other_poisson); }
            else { equiv_poisson = 0.0; }

            double equiv_young = 2.0 * myYoung * other_young / (myYoung + other_young);
            double calculation_area = 0.0;
            const double equiv_shear = equiv_young / (2.0 * (1 + equiv_poisson));

            if (i < (int)mContinuumInitialNeighborsSize) {
                if (r_process_info[USE_INITIAL_BOND_CONTACT_AREA]) {
                    calculation_area = GetInitialBondContactArea(neighbour_iterator_id); //this is the area that we are going to use for the calculation of the forces. It is not the same as the one used in the contact law.
                } else {
                    mContinuumConstitutiveLawArray[i]->GetContactArea(GetRadius(), other_radius, cont_ini_neigh_area, i, calculation_area); //some Constitutive Laws get a value, some others calculate the value.
                }

                if (mContinuumConstitutiveLawArray[i]->GetTypeOfLaw() == "parallel_bond_CL") {
                    mContinuumConstitutiveLawArray[i]->CalculateElasticConstants(kn_el, kt_el, initial_dist, equiv_young, equiv_poisson, calculation_area, this, neighbour_iterator, indentation_particle);
                } else {
                    mContinuumConstitutiveLawArray[i]->CalculateElasticConstants(kn_el, kt_el, initial_dist, equiv_young, equiv_poisson, calculation_area, this, neighbour_iterator, indentation);
                }
            }

            EvaluateDeltaDisplacement(data_buffer, DeltDisp, RelVel, data_buffer.mLocalCoordSystem, data_buffer.mOldLocalCoordSystem, vel, delta_displ);

            if (this->Is(DEMFlags::HAS_ROTATION)) {
                RelativeDisplacementAndVelocityOfContactPointDueToRotationQuaternion(DeltDisp, RelVel, data_buffer.mOldLocalCoordSystem, other_radius, data_buffer.mDt, ang_vel, neighbour_iterator, data_buffer);
            }

            RelativeDisplacementAndVelocityOfContactPointDueToOtherReasons(r_process_info, DeltDisp, RelVel, data_buffer.mOldLocalCoordSystem, data_buffer.mLocalCoordSystem, neighbour_iterator);

            double LocalDeltDisp[3] = {0.0};
            double LocalElasticContactForce[3] = {0.0}; // 0: first tangential, // 1: second tangential, // 2: normal force
            double LocalElasticExtraContactForce[3] = {0.0};
            double GlobalElasticContactForce[3] = {0.0};
            double GlobalElasticExtraContactForce[3] = {0.0};
            double TotalGlobalElasticContactForce[3] = {0.0};
            double OldLocalElasticContactForce[3] = {0.0};

            //FilterNonSignificantDisplacements(DeltDisp, RelVel, indentation);

            GeometryFunctions::VectorGlobal2Local(data_buffer.mLocalCoordSystem, DeltDisp, LocalDeltDisp);

            RotateOldContactForces(data_buffer.mOldLocalCoordSystem, data_buffer.mLocalCoordSystem, mNeighbourElasticContactForces[i]);
            RotateOldContactForces(data_buffer.mOldLocalCoordSystem, data_buffer.mLocalCoordSystem, mNeighbourElasticExtraContactForces[i]);

            GeometryFunctions::VectorGlobal2Local(data_buffer.mLocalCoordSystem, mNeighbourElasticContactForces[i], OldLocalElasticContactForce);

            GlobalElasticContactForce[0] = mNeighbourElasticContactForces[i][0];
            GlobalElasticContactForce[1] = mNeighbourElasticContactForces[i][1];
            GlobalElasticContactForce[2] = mNeighbourElasticContactForces[i][2];

            GeometryFunctions::VectorGlobal2Local(data_buffer.mLocalCoordSystem, GlobalElasticContactForce, LocalElasticContactForce); //TODO: can we remove this? We should overwrite LocalElasticContactForce afterwards

            double ViscoDampingLocalContactForce[3] = {0.0};
            double equiv_visco_damp_coeff_normal;
            double equiv_visco_damp_coeff_tangential;
            double ElasticLocalRotationalMoment[3] = {0.0};
            double ViscoLocalRotationalMoment[3] = {0.0};
            double cohesive_force =  0.0;
            double LocalRelVel[3] = {0.0};
            double LocalContactForce[3] = {0.0};
            double GlobalContactForce[3] = {0.0};
            double LocalUnbondedContactForce[3] = {0.0};// TODO: only works for parallel_bond_CL at the moment

            GeometryFunctions::VectorGlobal2Local(data_buffer.mLocalCoordSystem, RelVel, LocalRelVel);
            
            //*******************Forces calculation****************
            if (i < (int)mContinuumInitialNeighborsSize) {

                mContinuumConstitutiveLawArray[i]->CalculateForces(r_process_info,
                                                                OldLocalElasticContactForce,
                                                                LocalElasticContactForce,
                                                                LocalElasticExtraContactForce,
                                                                data_buffer.mLocalCoordSystem,
                                                                LocalDeltDisp,
                                                                kn_el,
                                                                kt_el,
                                                                contact_sigma,
                                                                contact_tau,
                                                                failure_criterion_state,
                                                                equiv_young,
                                                                equiv_shear,
                                                                indentation,
                                                                indentation_particle,
                                                                calculation_area,
                                                                acumulated_damage,
                                                                this,
                                                                neighbour_iterator,
                                                                i,
                                                                r_process_info[TIME_STEPS],
                                                                sliding,
                                                                equiv_visco_damp_coeff_normal,
                                                                equiv_visco_damp_coeff_tangential,
                                                                LocalRelVel,
                                                                ViscoDampingLocalContactForce);
                
                if (mContinuumConstitutiveLawArray[i]->GetTypeOfLaw() == "parallel_bond_CL") {
                    mContinuumConstitutiveLawArray[i]->GetLocalUnbondedContactForce(LocalUnbondedContactForce, LocalElasticContactForce);
                }

            } else if (indentation_particle > 0.0) {
                const double previous_indentation = indentation_particle + LocalDeltDisp[2];
                mDiscontinuumConstitutiveLaw = pCloneDiscontinuumConstitutiveLawWithNeighbour(data_buffer.mpOtherParticle);
                mDiscontinuumConstitutiveLaw->CalculateForces(r_process_info, OldLocalElasticContactForce, LocalElasticContactForce,
                                                                LocalDeltDisp, LocalRelVel, indentation_particle, previous_indentation,
                                                                ViscoDampingLocalContactForce, cohesive_force, this, data_buffer.mpOtherParticle, sliding, data_buffer.mLocalCoordSystem);
            } else { //Not bonded and no idata_buffer.mpOtherParticlendentation
                LocalElasticContactForce[0] = 0.0;      LocalElasticContactForce[1] = 0.0;      LocalElasticContactForce[2] = 0.0;
                ViscoDampingLocalContactForce[0] = 0.0; ViscoDampingLocalContactForce[1] = 0.0; ViscoDampingLocalContactForce[2] = 0.0;
                cohesive_force= 0.0;
            }

            array_1d<double, 3> other_ball_to_ball_forces(3,0.0);
            ComputeOtherBallToBallForces(other_ball_to_ball_forces);

            if (i < (int)mContinuumInitialNeighborsSize && mContinuumConstitutiveLawArray[i]->GetTypeOfLaw() == "smooth_joint_CL") {
                //double JointLocalCoordSystem[3][3];
                double GlobalJointNormal[3];
                mContinuumConstitutiveLawArray[i]->GetGlobalJointNormal(GlobalJointNormal);
                //GeometryFunctions::RotateCoordToDirection(data_buffer.mLocalCoordSystem, LocalJointNormal, JointLocalCoordSystem);
                //GeometryFunctions::VectorLocal2Global(JointLocalCoordSystem, LocalElasticContactForce, GlobalElasticContactForce);
                GeometryFunctions::VectorLocal2Global(data_buffer.mLocalCoordSystem, LocalElasticContactForce, GlobalElasticContactForce);
                //GlobalElasticContactForce[0] = LocalElasticContactForce[0]; //LocalElasticContactForce is acturlly a global force for "smooth_joint_CL"
                //GlobalElasticContactForce[1] = LocalElasticContactForce[1];
                //GlobalElasticContactForce[2] = LocalElasticContactForce[2];
                GeometryFunctions::RotateVectorToVector(GlobalElasticContactForce, GlobalJointNormal, GlobalElasticContactForce);

            } else {
                GeometryFunctions::VectorLocal2Global(data_buffer.mLocalCoordSystem, LocalElasticContactForce, GlobalElasticContactForce);
            }

            //*******************Add up forces and moments**************
            AddUpForcesAndProjectContinuum(data_buffer.mOldLocalCoordSystem, 
                                            data_buffer.mLocalCoordSystem, 
                                            LocalContactForce, 
                                            LocalElasticContactForce, 
                                            LocalElasticExtraContactForce, 
                                            GlobalContactForce,
                                            GlobalElasticContactForce, 
                                            GlobalElasticExtraContactForce, 
                                            TotalGlobalElasticContactForce, 
                                            ViscoDampingLocalContactForce,
                                            LocalUnbondedContactForce,
                                            0.0, 
                                            other_ball_to_ball_forces, 
                                            rElasticForce, 
                                            rContactForce, 
                                            i, 
                                            r_process_info); //TODO: replace the 0.0 with an actual cohesive force for discontinuum neighbours


            //******************Moments calculation start****************
            if (this->Is(DEMFlags::HAS_ROTATION)) {
                if (i < (int)mContinuumInitialNeighborsSize) {
                    mContinuumConstitutiveLawArray[i]->CalculateMoments(this,
                                                                        neighbour_iterator,
                                                                        equiv_young,
                                                                        data_buffer.mDistance,
                                                                        calculation_area,
                                                                        data_buffer.mLocalCoordSystem,
                                                                        ElasticLocalRotationalMoment,
                                                                        ViscoLocalRotationalMoment,
                                                                        equiv_poisson,
                                                                        indentation,
                                                                        indentation_particle,
                                                                        LocalContactForce[2],
                                                                        GlobalContactForce,
                                                                        data_buffer.mLocalCoordSystem[2],
                                                                        i);
                                                            
                    double LocalContactForceUsedForRolling[3] = {0.0}; //TODO: this can be improved
                    if (mContinuumConstitutiveLawArray[i]->GetTypeOfLaw() == "parallel_bond_CL") {
                        LocalContactForceUsedForRolling[0] = LocalUnbondedContactForce[0];
                        LocalContactForceUsedForRolling[1] = LocalUnbondedContactForce[1];
                        LocalContactForceUsedForRolling[2] = LocalUnbondedContactForce[2];
                    } else {
                        LocalContactForceUsedForRolling[0] = LocalContactForce[0];
                        LocalContactForceUsedForRolling[1] = LocalContactForce[1];
                        LocalContactForceUsedForRolling[2] = LocalContactForce[2];
                    }

                    if (this->Is(DEMFlags::HAS_ROLLING_FRICTION) && !data_buffer.mMultiStageRHS) {
                        if (mRollingFrictionModel->CheckIfThisModelRequiresRecloningForEachNeighbour()){
                            mRollingFrictionModel = pCloneRollingFrictionModelWithNeighbour(data_buffer.mpOtherParticle);
                            mRollingFrictionModel->ComputeRollingFriction(this, 
                                                                        data_buffer.mpOtherParticle, 
                                                                        r_process_info, 
                                                                        LocalContactForceUsedForRolling, 
                                                                        indentation_particle, 
                                                                        mContactMoment, 
                                                                        data_buffer.mLocalCoordSystem[2],
                                                                        mNeighbourRollingFrictionMoments[i]);
                        }
                        else {
                            mRollingFrictionModel->ComputeRollingResistance(this, data_buffer.mpOtherParticle, LocalContactForceUsedForRolling);
                        }
                    }

                } else { //for unbonded particles

                    ComputeMoments(LocalContactForce[2], GlobalContactForce, data_buffer.mLocalCoordSystem[2], data_buffer.mpOtherParticle, indentation_particle, i);

                    if (this->Is(DEMFlags::HAS_ROLLING_FRICTION) && !data_buffer.mMultiStageRHS) {
                        if (mRollingFrictionModel->CheckIfThisModelRequiresRecloningForEachNeighbour()){
                            mRollingFrictionModel = pCloneRollingFrictionModelWithNeighbour(data_buffer.mpOtherParticle);
                            mRollingFrictionModel->ComputeRollingFriction(this, 
                                                                        data_buffer.mpOtherParticle, 
                                                                        r_process_info, 
                                                                        LocalContactForce, 
                                                                        indentation_particle, 
                                                                        mContactMoment, 
                                                                        data_buffer.mLocalCoordSystem[2],
                                                                        mNeighbourRollingFrictionMoments[i]);
                        }
                        else {
                            mRollingFrictionModel->ComputeRollingResistance(this, data_buffer.mpOtherParticle, LocalContactForce);
                        }
                    }
                }
            }
            //*****************Moments calculation end******************

            if (this->Is(DEMFlags::HAS_STRESS_TENSOR) && (i < (int)mContinuumInitialNeighborsSize)) { // We leave apart the discontinuum neighbors (the same for the walls). The neighbor would not be able to do the same if we activate it.
                mContinuumConstitutiveLawArray[i]->AddPoissonContribution(equiv_poisson, data_buffer.mLocalCoordSystem, LocalElasticContactForce[2],
                                                                          calculation_area, mSymmStressTensor, this, neighbour_iterator, r_process_info, i, indentation);
            }

            //*******************Bond failure check*********************
            if (i < (int)mContinuumInitialNeighborsSize) {
                mContinuumConstitutiveLawArray[i]->CheckFailure(i, this, neighbour_iterator, contact_sigma, contact_tau, LocalElasticContactForce,
                                                                    ViscoDampingLocalContactForce, ElasticLocalRotationalMoment, ViscoLocalRotationalMoment);
            }

            if (this->Is(DEMFlags::HAS_ROTATION)) {
                AddUpMomentsAndProject(data_buffer.mLocalCoordSystem, ElasticLocalRotationalMoment, ViscoLocalRotationalMoment);
            }

            
            int contact_mesh_option = r_process_info[CONTACT_MESH_OPTION];
            if (contact_mesh_option == 1 && this->Id() < neighbour_iterator_id && (r_process_info[IS_TIME_TO_PRINT] || r_process_info[IS_TIME_TO_UPDATE_CONTACT_ELEMENT])) {
                if (i < (int)mContinuumInitialNeighborsSize) {
                    double total_local_elastic_contact_force[3] = {0.0};
                    total_local_elastic_contact_force[0] = LocalElasticContactForce[0] + LocalElasticExtraContactForce[0];
                    total_local_elastic_contact_force[1] = LocalElasticContactForce[1] + LocalElasticExtraContactForce[1];
                    total_local_elastic_contact_force[2] = LocalElasticContactForce[2] + LocalElasticExtraContactForce[2];
                    double bond_volume = GetInitialBondVolume(neighbour_iterator_id);
                    CalculateOnContinuumContactElements(i, total_local_elastic_contact_force, contact_sigma, contact_tau, failure_criterion_state, acumulated_damage, time_steps, calculation_area, GlobalContactForce, bond_volume);
                } else {
                    CalculateOnContactElements(i, LocalContactForce, GlobalContactForce);
                }
            } 

            if (this->Is(DEMFlags::HAS_STRESS_TENSOR) /*&& (i < mContinuumInitialNeighborsSize)*/) {
                AddNeighbourContributionToStressTensor(r_process_info, GlobalContactForce, data_buffer.mLocalCoordSystem[2], data_buffer.mDistance, radius_sum, this);
            }

            AddContributionToRepresentativeVolume(data_buffer.mDistance, radius_sum, calculation_area);

            ComputeForceWithNeighbourFinalOperations();

            /*if (i < mContinuumInitialNeighborsSize) {
                DEM_COPY_SECOND_TO_FIRST_3(mArrayOfDeltaDisplacements[i], DeltDisp);
            }*/
        } // for each neighbor

        ComputeBrokenBondsRatio();

        KRATOS_CATCH("")
    } //  ComputeBallToBallContactForceAndMoment

    void SphericContinuumParticle::AddUpForcesAndProjectContinuum(double OldCoordSystem[3][3],
                                                                double LocalCoordSystem[3][3],
                                                                double LocalContactForce[3],
                                                                double LocalElasticContactForce[3],
                                                                double LocalElasticExtraContactForce[3],
                                                                double GlobalContactForce[3],
                                                                double GlobalElasticContactForce[3],
                                                                double GlobalElasticExtraContactForce[3],
                                                                double TotalGlobalElasticContactForce[3],
                                                                double ViscoDampingLocalContactForce[3],
                                                                double LocalUnbondedContactForce[3],
                                                                const double cohesive_force,
                                                                array_1d<double, 3>& other_ball_to_ball_forces,
                                                                array_1d<double, 3>& r_elastic_force,
                                                                array_1d<double, 3>& r_contact_force,
                                                                const unsigned int i_neighbour_count,
                                                                const ProcessInfo& r_process_info)
    {

        for (unsigned int index = 0; index < 3; index++) {
            LocalContactForce[index] = LocalElasticContactForce[index] + ViscoDampingLocalContactForce[index] + other_ball_to_ball_forces[index];
            LocalUnbondedContactForce[index] += other_ball_to_ball_forces[index];
        }
        LocalContactForce[2] -= cohesive_force;
        LocalUnbondedContactForce[2] -= cohesive_force;

        DEM_ADD_SECOND_TO_FIRST(LocalElasticContactForce, other_ball_to_ball_forces);

        GeometryFunctions::VectorLocal2Global(LocalCoordSystem, LocalElasticContactForce, GlobalElasticContactForce);
        GeometryFunctions::VectorLocal2Global(LocalCoordSystem, LocalContactForce, GlobalContactForce);
        GeometryFunctions::VectorLocal2Global(LocalCoordSystem, LocalElasticExtraContactForce, GlobalElasticExtraContactForce);

        // Saving contact forces (We need to, since tangential elastic force is history-dependent)
        DEM_COPY_SECOND_TO_FIRST_3(mNeighbourElasticContactForces[i_neighbour_count], GlobalElasticContactForce)
        DEM_COPY_SECOND_TO_FIRST_3(mNeighbourElasticExtraContactForces[i_neighbour_count], GlobalElasticExtraContactForce)

        TotalGlobalElasticContactForce[0] = GlobalElasticContactForce[0] + GlobalElasticExtraContactForce[0];
        TotalGlobalElasticContactForce[1] = GlobalElasticContactForce[1] + GlobalElasticExtraContactForce[1];
        TotalGlobalElasticContactForce[2] = GlobalElasticContactForce[2] + GlobalElasticExtraContactForce[2];
        DEM_ADD_SECOND_TO_FIRST(r_elastic_force, TotalGlobalElasticContactForce)

        double TotalGlobalContactForce[3];
        TotalGlobalContactForce[0] = GlobalContactForce[0] + GlobalElasticExtraContactForce[0];
        TotalGlobalContactForce[1] = GlobalContactForce[1] + GlobalElasticExtraContactForce[1];
        TotalGlobalContactForce[2] = GlobalContactForce[2] + GlobalElasticExtraContactForce[2];
        DEM_ADD_SECOND_TO_FIRST(r_contact_force, TotalGlobalContactForce )
    }

    void SphericContinuumParticle::ComputeForceWithNeighbourFinalOperations(){}

    void SphericContinuumParticle::ComputeBrokenBondsRatio() {

        int BrokenBondsCounter = 0;

        for (unsigned int i = 0; i < mContinuumInitialNeighborsSize; i++) {
            /*if (i >= mNeighbourElements.size()) {
                BrokenBondsCounter++;
            } else*/
            if (mBondElements[i] == NULL || mIniNeighbourFailureId[i] > 0) {
                BrokenBondsCounter++; 
            }
        }

        if(mContinuumInitialNeighborsSize) {
            GetGeometry()[0].FastGetSolutionStepValue(DAMAGE_RATIO) = double(BrokenBondsCounter) / mContinuumInitialNeighborsSize;
        } else {
            GetGeometry()[0].FastGetSolutionStepValue(DAMAGE_RATIO) = 0.0; // No neighbors present, no damage. Or DAMAGE_RATIO = 1.0.
        }

    }

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

    void SphericContinuumParticle::FinalizeSolutionStep(const ProcessInfo& r_process_info) {
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

            if (mNeighbourElements[i] == NULL || mIniNeighbourFailureId[i]) {
                *mSkinSphere = 1.0;
                break;
            }
        }

        KRATOS_CATCH("")
    }

    void SphericContinuumParticle::ReorderAndRecoverInitialPositionsAndFilter(std::vector<SphericParticle*>& temp_neighbour_elements, const ProcessInfo& r_process_info) {

        KRATOS_TRY

        unsigned int current_neighbors_size = mNeighbourElements.size();
        unsigned int initial_neighbors_size = mIniNeighbourIds.size();
        temp_neighbour_elements.resize(initial_neighbors_size);

        for (unsigned int i = 0; i < initial_neighbors_size; i++) {
            temp_neighbour_elements[i] = NULL;
        }

        // Loop over current neighbors
        for (unsigned int i = 0; i < current_neighbors_size; i++) {
            SphericParticle* i_neighbour = mNeighbourElements[i];
            bool found = false;
            // Loop over initial neighbors
            for (unsigned int kk = 0; kk < initial_neighbors_size; kk++) {

                if (static_cast<int>(i_neighbour->Id()) == mIniNeighbourIds[kk]) {
                    temp_neighbour_elements[kk] = i_neighbour;
                    found = true;
                    break;
                }
            }

            if (!found) {
                double other_radius = i_neighbour->GetInteractionRadius();
                double radius_sum = GetInteractionRadius() + other_radius;
                array_1d<double, 3> other_to_me_vect;
                if (!r_process_info[DOMAIN_IS_PERIODIC]){ // default infinite-domain case
                    noalias(other_to_me_vect) = this->GetGeometry()[0].Coordinates() - i_neighbour->GetGeometry()[0].Coordinates();
                } else { // periodic domain
                    auto& central_node = GetGeometry()[0];
                    auto& neighbour_node = i_neighbour->GetGeometry()[0];
                    double my_coors[3] = {central_node[0], central_node[1], central_node[2]};
                    double other_coors[3] = {neighbour_node[0], neighbour_node[1], neighbour_node[2]};
                    TransformNeighbourCoorsToClosestInPeriodicDomain(r_process_info, my_coors, other_coors);
                    other_to_me_vect[0] = my_coors[0] - other_coors[0];
                    other_to_me_vect[1] = my_coors[1] - other_coors[1];
                    other_to_me_vect[2] = my_coors[2] - other_coors[2];
                }
                double distance = DEM_MODULUS_3(other_to_me_vect);
                double indentation = radius_sum - distance;

                if (indentation > 0.0) {
                    temp_neighbour_elements.push_back(i_neighbour);
                }
            }
        }

        mNeighbourElements.swap(temp_neighbour_elements);

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

    void SphericContinuumParticle::UpdateContinuumNeighboursVector(const ProcessInfo& r_process_info) {}

    void SphericContinuumParticle::ReorderFEMneighbours() {

        KRATOS_TRY

        unsigned int current_neighbors_size = mNeighbourRigidFaces.size();
        unsigned int initial_neighbors_size = mFemIniNeighbourIds.size();

        std::vector<DEMWall*> temp_neighbour_elements(initial_neighbors_size, nullptr);
        std::vector<array_1d<double, 4> > temp_neighbours_weights(initial_neighbors_size, ZeroVector(4));
        std::vector<int> temp_neighbours_contact_types(initial_neighbors_size, 0);

        // Loop over current neighbors
        for (unsigned int i = 0; i < current_neighbors_size; i++) {
            DEMWall* i_neighbour = mNeighbourRigidFaces[i];
            bool found = false;
            // Loop over initial neighbors
            for (unsigned int k = 0; k < initial_neighbors_size; k++) {
                if (static_cast<int>(i_neighbour->Id()) == mFemIniNeighbourIds[k]) {
                    temp_neighbour_elements[k] = i_neighbour;
                    temp_neighbours_weights[k] = mContactConditionWeights[i];
                    temp_neighbours_contact_types[k] = mContactConditionContactTypes[i];
                    found = true;
                    break;
                }
            }

            if (!found) {
                temp_neighbour_elements.push_back(i_neighbour);
                temp_neighbours_weights.push_back(mContactConditionWeights[i]);
                temp_neighbours_contact_types.push_back(mContactConditionContactTypes[i]);
            }
        }

        mNeighbourRigidFaces.swap(temp_neighbour_elements);
        mContactConditionWeights.swap(temp_neighbours_weights);
        mContactConditionContactTypes.swap(temp_neighbours_contact_types);

        KRATOS_CATCH("")
    }

    double SphericContinuumParticle::CalculateMaxSearchDistance(const bool has_mpi, const ProcessInfo& r_process_info) {

        KRATOS_TRY

        double max_local_search = 0.0;

        for (unsigned int i = 0; i < mContinuumInitialNeighborsSize; i++) {
            if (mNeighbourElements[i] == NULL) continue;
            SphericContinuumParticle* r_continuum_ini_neighbour = dynamic_cast<SphericContinuumParticle*>(mNeighbourElements[i]);
            double search_dist = mContinuumConstitutiveLawArray[i]->LocalMaxSearchDistance(i, this, r_continuum_ini_neighbour);
            if (search_dist > max_local_search) max_local_search = search_dist;
        }

        return max_local_search;

        KRATOS_CATCH("")
    }

    bool SphericContinuumParticle::OverlappedParticleRemoval(const ProcessInfo& r_process_info) {

        KRATOS_TRY

        for (unsigned int i = 0; i < mNeighbourElements.size(); i++) {
            if (mNeighbourElements[i] == NULL) continue;
            SphericParticle* ini_cont_neighbour_iterator = mNeighbourElements[i];
            double other_radius = ini_cont_neighbour_iterator->GetRadius();

            array_1d<double, 3> other_to_me_vect;
            auto& central_node = GetGeometry()[0];
            auto& neighbour_node = mNeighbourElements[i]->GetGeometry()[0];
            if (!r_process_info[DOMAIN_IS_PERIODIC]){ // default infinite-domain case
                noalias(other_to_me_vect) = central_node.Coordinates() - neighbour_node.Coordinates();
            } else { // periodic domain
                double my_coors[3] = {central_node[0], central_node[1], central_node[2]};
                double other_coors[3] = {neighbour_node[0], neighbour_node[1], neighbour_node[2]};
                TransformNeighbourCoorsToClosestInPeriodicDomain(r_process_info, my_coors, other_coors);
                other_to_me_vect[0] = my_coors[0] - other_coors[0];
                other_to_me_vect[1] = my_coors[1] - other_coors[1];
                other_to_me_vect[2] = my_coors[2] - other_coors[2];
            }
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
                DenseVector<int>& vector_of_ids_of_neighbours_of_this_neighbour = r_continuum_ini_neighbour->GetValue(NEIGHBOUR_IDS);
                if ((int) this->Id() == vector_of_ids_of_neighbours_of_this_neighbour[j]) index_of_the_neighbour_that_is_me = (int)j;
            }

            if (index_of_the_neighbour_that_is_me == -1) {
                std::string message = "An element (Id " + std::to_string(this->Id()) + ") found a neighbor (had contact area) but the neighbor (Id " \
                                                        + std::to_string(r_continuum_ini_neighbour->Id()) + ") did not have area for that element  ";

                KRATOS_ERROR << message << std::endl;
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
            KRATOS_ERROR_IF(coeff > 1.0) << "The coefficient assigned for virtual mass is larger than one, virtual_mass_coeff= "<< coeff << std::endl;
            else if ((coeff == 1.0) && (r_process_info[VIRTUAL_MASS_OPTION])) { Output = 9.0E09; }
            else {
                if (r_process_info[VIRTUAL_MASS_OPTION]) { mass /= 1 - coeff; }
                double K = GetYoung() * Globals::Pi * GetRadius();
                Output = 0.34 * sqrt(mass / K);
                if (r_process_info[ROTATION_OPTION] == 1) { Output = Output * 0.5; } //factor for critical time step when rotation is allowed.
            }
            return;
        }//CRITICAL DELTA CALCULATION

        SphericParticle::Calculate(rVariable, Output, r_process_info);

        KRATOS_CATCH("")
    }//Calculate

    void SphericContinuumParticle::Initialize(const ProcessInfo& r_process_info){

        SphericParticle::Initialize(r_process_info);

        SetValue(NEIGHBOURS_CONTACT_AREAS, Vector());

        mSkinSphere     = &(this->GetGeometry()[0].FastGetSolutionStepValue(SKIN_SPHERE));
        mContinuumGroup = this->GetGeometry()[0].FastGetSolutionStepValue(COHESIVE_GROUP);
    }

    double SphericContinuumParticle::GetInitialDelta(int index) {
        //if (index < (int) mIniNeighbourDelta.size()) return mIniNeighbourDelta[index];
        if (mIniNeighbourDelta.find(static_cast<int>(index)) != mIniNeighbourDelta.end()){
            return mIniNeighbourDelta[static_cast<int>(index)];
        } else {
            return 0.0;
        }
    }

    /*
    double SphericContinuumParticle::GetInitialBondContactArea(int index) {
        if (mIniBondContactArea.find(static_cast<int>(index)) != mIniBondContactArea.end()){
            return mIniBondContactArea[static_cast<int>(index)];
        } else {
            return 0.0;
        }
    }*/

    double SphericContinuumParticle::GetInitialBondContactArea(int index) {
        std::pair<int, int> key = makeKey(static_cast<int>(this->Id()), static_cast<int>(index));
        if (mCementedContactAreasMapPtr->find(key) != mCementedContactAreasMapPtr->end()){
            return (*mCementedContactAreasMapPtr)[key];
        } else {
            return 0.0;
        }
    }

    double SphericContinuumParticle::GetInitialBondVolume(int index) {
        if (mBondVolume.find(static_cast<int>(index)) != mBondVolume.end()){
            return mBondVolume[static_cast<int>(index)];
        } else {
            return 0.0;
        }
    }

    double SphericContinuumParticle::GetInitialDeltaWithFEM(int index) {
        if(index< (int) mFemIniNeighbourDelta.size()) return mFemIniNeighbourDelta[index];
        else return 0.0;
    }

    void SphericContinuumParticle::CalculateOnContinuumContactElements(size_t i, double LocalElasticContactForce[3],
                                                                        double contact_sigma, 
                                                                        double contact_tau, 
                                                                        double failure_criterion_state, 
                                                                        double acumulated_damage, 
                                                                        int time_steps, 
                                                                        double calculation_area, 
                                                                        double GlobalContactForce[3],
                                                                        double bond_volume)
    {

        KRATOS_TRY
        if (!mBondElements.size()) return; // we skip this function if the vector of bonds hasn't been filled yet.
        //KRATOS_ERROR_IF(i >= mBondElements.size()) << "The index of the bond is larger than the size of the vector of bonds. i = " << i << " and size = " << mBondElements.size() << std::endl;
        ParticleContactElement* bond = mBondElements[i];
        if (bond == NULL) return; //This bond was never created (happens in some MPI cases, see CreateContactElements() in explicit_solve_continumm.h)

        bond->mLocalContactForce[0] = LocalElasticContactForce[0];
        bond->mLocalContactForce[1] = LocalElasticContactForce[1];
        bond->mLocalContactForce[2] = LocalElasticContactForce[2];
        bond->mContactSigma = contact_sigma;
        bond->mContactTau = contact_tau;
        bond->mContactFailure = mIniNeighbourFailureId[i];
        bond->mFailureCriterionState = failure_criterion_state;
        bond->mContactRadius = sqrt(calculation_area / Globals::Pi);
        bond->mBondVolume = bond_volume;

        bond->mGlobalContactForce[0] = GlobalContactForce[0];
        bond->mGlobalContactForce[1] = GlobalContactForce[1];
        bond->mGlobalContactForce[2] = GlobalContactForce[2];

        if ((time_steps == 0) || (acumulated_damage > bond->mUnidimendionalDamage)) {
            bond->mUnidimendionalDamage = acumulated_damage;
        }
        KRATOS_CATCH("")
    }

    void SphericContinuumParticle::GetCementedContactAreasMap(std::map<std::pair<int, int>, double>* CementedContactAreasMap){
        mCementedContactAreasMapPtr = CementedContactAreasMap;
    }

    void SphericContinuumParticle::GetCementedContactPairsSet(std::set<std::pair<int, int>>* CementedContactPairsSet){
        mCementedContactPairsSetPtr = CementedContactPairsSet;
    }

} // namespace Kratos
