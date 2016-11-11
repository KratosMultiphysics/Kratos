// Author: 
// Miguel Angel Celigueta maceli@cimne.upc.edu
//
//README::::look to the key word "VERSION" if you want to find all the points where you have to change something so that you can pass from a kdtree to a bin data search structure;

#include "create_and_destroy.h"
#include "../custom_elements/spheric_continuum_particle.h"
#include "../custom_elements/cluster3D.h"
#include "../custom_utilities/GeometryFunctions.h"
#include "../custom_utilities/AuxiliaryFunctions.h"

namespace Kratos {
    
    ParticleCreatorDestructor::ParticleCreatorDestructor() : mGreatestParticleId(0) {
        mScaleFactor = 1.0;
        mHighPoint[0] = 10e18;
        mHighPoint[1] = 10e18;
        mHighPoint[2] = 10e18;
        mLowPoint[0] = -10e18;
        mLowPoint[1] = -10e18;
        mLowPoint[2] = -10e18;
    }

    //Particle_Creator_Destructor() {};

    /// Destructor.

    ParticleCreatorDestructor::~ParticleCreatorDestructor() {
        mDoSearchNeighbourElements = true; // true by default. It should be set to false by the strategy (friend class) if appropriate
    }

    int ParticleCreatorDestructor::FindMaxNodeIdInModelPart(ModelPart& r_modelpart) {
        KRATOS_TRY
        int max_Id = 1; //GID accepts Id's >= 1
        std::vector<int> thread_maximums(OpenMPUtils::GetNumThreads(),1);

        //#pragma omp parallel for
        for(int i=0; i<(int)r_modelpart.GetCommunicator().LocalMesh().Nodes().size(); i++){
            ModelPart::NodesContainerType::iterator node_it = r_modelpart.GetCommunicator().LocalMesh().NodesBegin() + i;        
            if ((int) (node_it->Id()) > thread_maximums[OpenMPUtils::ThisThread()]) thread_maximums[OpenMPUtils::ThisThread()] = node_it->Id();
        }
        
        for(int i=0; i<OpenMPUtils::GetNumThreads(); i++){
            if(thread_maximums[i] > max_Id) max_Id = thread_maximums[i];
        }
        
        r_modelpart.GetCommunicator().MaxAll(max_Id);
        return max_Id;
        KRATOS_CATCH("")
    }
    
    void ParticleCreatorDestructor::FindAndSaveMaxNodeIdInModelPart(ModelPart& r_modelpart) {
        KRATOS_TRY
        mMaxNodeId = FindMaxNodeIdInModelPart(r_modelpart);
        KRATOS_CATCH("")
    }

    int ParticleCreatorDestructor::FindMaxElementIdInModelPart(ModelPart& r_modelpart) {
        KRATOS_TRY
        int max_Id = 1;
        
        for (ModelPart::ElementsContainerType::iterator element_it = r_modelpart.GetCommunicator().LocalMesh().ElementsBegin();
            element_it != r_modelpart.GetCommunicator().LocalMesh().ElementsEnd();
            element_it++) {

            if ((int) (element_it->Id()) > max_Id) max_Id = element_it->Id();
        }

        r_modelpart.GetCommunicator().MaxAll(max_Id);
        return max_Id;
        KRATOS_CATCH("")
    }
    
    int ParticleCreatorDestructor::FindMaxConditionIdInModelPart(ModelPart& r_modelpart) {
        KRATOS_TRY
        int max_Id = 1;
        
        for (ModelPart::ConditionsContainerType::iterator condition_it = r_modelpart.GetCommunicator().LocalMesh().ConditionsBegin();
            condition_it != r_modelpart.GetCommunicator().LocalMesh().ConditionsEnd();
            condition_it++) {

            if ((int) (condition_it->Id()) > max_Id) max_Id = condition_it->Id();
        }

        r_modelpart.GetCommunicator().MaxAll(max_Id);
        return max_Id;
        KRATOS_CATCH("")
    }
    
    void ParticleCreatorDestructor::RenumberElementIdsFromGivenValue(ModelPart& r_modelpart, const int initial_id) {
        KRATOS_TRY
        const int number_of_elements = r_modelpart.GetCommunicator().LocalMesh().NumberOfElements();
        int total_accumulated_elements = 0;
        
        r_modelpart.GetCommunicator().ScanSum(number_of_elements, total_accumulated_elements);

        int id = total_accumulated_elements - number_of_elements + initial_id;
        auto& local_mesh = r_modelpart.GetCommunicator().LocalMesh();
        
        for (auto element_it = local_mesh.ElementsBegin(); element_it != local_mesh.ElementsEnd(); element_it++) {
            element_it->SetId(id);
            id++;
        }                
        KRATOS_CATCH("")
    }
        
    
    double ParticleCreatorDestructor::rand_normal(const double mean, const double stddev, const double max_radius, const double min_radius) {
        KRATOS_TRY
        if (!stddev) return mean;

        double return_value;

        do {
            double x, y, r;

            do {
                x = 2.0 * rand() / RAND_MAX - 1;
                y = 2.0 * rand() / RAND_MAX - 1;
                r = x * x + y*y;
            } while (r == 0.0 || r > 1.0);

            double d = sqrt(-2.0 * log(r) / r);
            return_value = x * d * stddev + mean;

        } while (return_value < min_radius || return_value > max_radius);

        return return_value;
        KRATOS_CATCH("")
    }
    
    double ParticleCreatorDestructor::rand_lognormal(const double mean, const double stddev, const double max_radius, const double min_radius) {
        KRATOS_TRY
        const double normal_mean = log(mean * mean / sqrt(stddev * stddev + mean * mean));
        const double normal_stddev = sqrt(log(1 + stddev * stddev / (mean * mean)));
        double normally_distributed_value = rand_normal(normal_mean, normal_stddev, log(max_radius), log(min_radius));

        return exp(normally_distributed_value);
        KRATOS_CATCH("")
    }
    
    void ParticleCreatorDestructor::AddRandomPerpendicularVelocityToGivenVelocity(array_1d<double, 3 >& velocity, const double angle_in_degrees) {
        KRATOS_TRY
        double velocity_modulus = sqrt(velocity[0]*velocity[0] + velocity[1]*velocity[1] + velocity[2]*velocity[2]);
        array_1d<double, 3 > unitary_velocity;
        noalias(unitary_velocity) = velocity / velocity_modulus;
        array_1d<double, 3 > normal_1;
        array_1d<double, 3 > normal_2;
        
        if (fabs(unitary_velocity[0])>=0.577) {
            normal_1[0]= - unitary_velocity[1];
            normal_1[1]= unitary_velocity[0];
            normal_1[2]= 0.0;
        }
        else if (fabs(unitary_velocity[1])>=0.577) {
            normal_1[0]= 0.0;
            normal_1[1]= - unitary_velocity[2];
            normal_1[2]= unitary_velocity[1];
        }        
        else {                   
            normal_1[0]= unitary_velocity[2];
            normal_1[1]= 0.0;
            normal_1[2]= - unitary_velocity[0];
        }        
        
        //normalize(normal_1);
        double distance0 = DEM_MODULUS_3(normal_1);
        double inv_distance0 = (distance0 != 0.0) ? 1.0 / distance0 : 0.00;
        normal_1[0] = normal_1[0] * inv_distance0;
        normal_1[1] = normal_1[1] * inv_distance0;
        normal_1[2] = normal_1[2] * inv_distance0;
        
        //CrossProduct(NormalDirection,Vector0,Vector1);
        normal_2[0] = unitary_velocity[1]*normal_1[2] - unitary_velocity[2]*normal_1[1];
        normal_2[1] = unitary_velocity[2]*normal_1[0] - unitary_velocity[0]*normal_1[2];
        normal_2[2] = unitary_velocity[0]*normal_1[1] - unitary_velocity[1]*normal_1[0];     
        
        double angle_in_radians = angle_in_degrees * KRATOS_M_PI / 180;
        double radius = tan(angle_in_radians) * velocity_modulus;
        double radius_square = radius * radius;
        double local_added_velocity_modulus_square = radius_square + 1.0; //just greater than the radius, to get at least one iteration of the while
        array_1d<double, 3> local_added_velocity; local_added_velocity[0] = local_added_velocity[1] = local_added_velocity[2] = 0.0;
        
        while (local_added_velocity_modulus_square > radius_square) {
            //Random in a range: (max - min) * ( (double)rand() / (double)RAND_MAX ) + min
            local_added_velocity[0] = 2*radius * (double)rand() / (double)RAND_MAX - radius;
            local_added_velocity[1] = 2*radius * (double)rand() / (double)RAND_MAX - radius;
            local_added_velocity_modulus_square = local_added_velocity[0]*local_added_velocity[0] + local_added_velocity[1]*local_added_velocity[1];
        }
        
        noalias(velocity) += local_added_velocity[0] * normal_1 + local_added_velocity[1] * normal_2;  
        KRATOS_CATCH("")
    }
    
    void ParticleCreatorDestructor::NodeCreatorWithPhysicalParameters(ModelPart& r_modelpart,
                                                                    Node < 3 > ::Pointer& pnew_node,
                                                                    int aId,
                                                                    Node < 3 > ::Pointer& reference_node,
                                                                    double radius,
                                                                    Properties& params,
                                                                    ModelPart& r_sub_model_part_with_parameters,
                                                                    bool has_sphericity,
                                                                    bool has_rotation,
                                                                    bool initial) {
        KRATOS_TRY
        array_1d<double, 3 > null_vector(3, 0.0);

        double bx = reference_node->X();
        double cy = reference_node->Y();
        double dz = reference_node->Z();

        if (initial) {
            pnew_node = reference_node;
            pnew_node->SetId(aId);
            #pragma omp critical
            {
                r_modelpart.AddNode(pnew_node); // The same node is added to r_modelpart (the calculation model part) 
            }
            pnew_node->FastGetSolutionStepValue(VELOCITY) = null_vector;            
            //(actually it should be the velocity of the inlet layer, which is different from the particles being inserted)
            pnew_node->FastGetSolutionStepValue(PARTICLE_MATERIAL) = params[PARTICLE_MATERIAL] + 100; //So the inlet ghost spheres are not in the same
        }                                                                                             //layer of the inlet newly created spheres
        else {                        
            pnew_node = boost::make_shared< Node<3> >(aId, bx, cy, dz);
            pnew_node->SetSolutionStepVariablesList(&r_modelpart.GetNodalSolutionStepVariablesList());
            pnew_node->SetBufferSize(r_modelpart.GetBufferSize());
            #pragma omp critical
            {
                //pnew_node = r_modelpart.CreateNewNode(aId, bx, cy, dz); //ACTUAL node creation and addition to model part
                r_modelpart.Nodes().push_back(pnew_node);
            }
                        
            array_1d<double, 3 > velocity = r_sub_model_part_with_parameters[VELOCITY];
            double max_rand_deviation_angle = r_sub_model_part_with_parameters[MAX_RAND_DEVIATION_ANGLE];
            AddRandomPerpendicularVelocityToGivenVelocity(velocity, max_rand_deviation_angle);
            pnew_node->FastGetSolutionStepValue(VELOCITY) = velocity;
            pnew_node->FastGetSolutionStepValue(PARTICLE_MATERIAL) = params[PARTICLE_MATERIAL];
        }

        if (has_rotation && pnew_node->SolutionStepsDataHas(PARTICLE_ROTATION_DAMP_RATIO) ) {
            pnew_node->FastGetSolutionStepValue(PARTICLE_ROTATION_DAMP_RATIO) = params[PARTICLE_ROTATION_DAMP_RATIO];
        }

        if (has_sphericity) {
            pnew_node->FastGetSolutionStepValue(PARTICLE_SPHERICITY) = params[PARTICLE_SPHERICITY];
        }

        pnew_node->FastGetSolutionStepValue(RADIUS) = radius;        
        pnew_node->FastGetSolutionStepValue(ANGULAR_VELOCITY) = null_vector;

        pnew_node->AddDof(VELOCITY_X, REACTION_X);
        pnew_node->AddDof(VELOCITY_Y, REACTION_Y);
        pnew_node->AddDof(VELOCITY_Z, REACTION_Z);
        pnew_node->AddDof(ANGULAR_VELOCITY_X, REACTION_X);
        pnew_node->AddDof(ANGULAR_VELOCITY_Y, REACTION_Y);
        pnew_node->AddDof(ANGULAR_VELOCITY_Z, REACTION_Z);

        pnew_node->pGetDof(VELOCITY_X)->FixDof();
        pnew_node->pGetDof(VELOCITY_Y)->FixDof();
        pnew_node->pGetDof(VELOCITY_Z)->FixDof();
        pnew_node->pGetDof(ANGULAR_VELOCITY_X)->FixDof();
        pnew_node->pGetDof(ANGULAR_VELOCITY_Y)->FixDof();
        pnew_node->pGetDof(ANGULAR_VELOCITY_Z)->FixDof();

        pnew_node->Set(DEMFlags::FIXED_VEL_X, true);
        pnew_node->Set(DEMFlags::FIXED_VEL_Y, true);
        pnew_node->Set(DEMFlags::FIXED_VEL_Z, true);
        pnew_node->Set(DEMFlags::FIXED_ANG_VEL_X, true);
        pnew_node->Set(DEMFlags::FIXED_ANG_VEL_Y, true);
        pnew_node->Set(DEMFlags::FIXED_ANG_VEL_Z, true);
        KRATOS_CATCH("")
    }
    
    void ParticleCreatorDestructor::NodeForClustersCreatorWithPhysicalParameters(ModelPart& r_modelpart,
                                                                    Node < 3 > ::Pointer& pnew_node,
                                                                    int aId,
                                                                    Node < 3 > ::Pointer& reference_node,
                                                                    Properties& params,
                                                                    ModelPart& r_sub_model_part_with_parameters,
                                                                    bool has_sphericity,
                                                                    bool has_rotation,
                                                                    bool initial) {
        KRATOS_TRY
        array_1d<double, 3 > null_vector(3, 0.0);

        double bx = reference_node->X();
        double cy = reference_node->Y();
        double dz = reference_node->Z();

        if (initial) {
            pnew_node = reference_node;
            pnew_node->SetId(aId);
            #pragma omp critical
            {
                r_modelpart.AddNode(pnew_node); // The same node is added to r_modelpart (the calculation model part) 
            }
            pnew_node->FastGetSolutionStepValue(VELOCITY) = null_vector;            
            //(actually it should be the velocity of the inlet layer, which is different from the particles being inserted)
            pnew_node->FastGetSolutionStepValue(PARTICLE_MATERIAL) = params[PARTICLE_MATERIAL] + 100; //So the inlet ghost spheres are not in the same
        }                                                                                             //layer of the inlet newly created spheres
        else {                        
            pnew_node = boost::make_shared< Node<3> >(aId, bx, cy, dz);
            pnew_node->SetSolutionStepVariablesList(&r_modelpart.GetNodalSolutionStepVariablesList());
            pnew_node->SetBufferSize(r_modelpart.GetBufferSize());
            #pragma omp critical
            {
                //pnew_node = r_modelpart.CreateNewNode(aId, bx, cy, dz); //ACTUAL node creation and addition to model part
                r_modelpart.Nodes().push_back(pnew_node);
            }
                        
            array_1d<double, 3 > velocity = r_sub_model_part_with_parameters[VELOCITY];
            double max_rand_deviation_angle = r_sub_model_part_with_parameters[MAX_RAND_DEVIATION_ANGLE];
            AddRandomPerpendicularVelocityToGivenVelocity(velocity, max_rand_deviation_angle);
            pnew_node->FastGetSolutionStepValue(VELOCITY) = velocity;
            pnew_node->FastGetSolutionStepValue(PARTICLE_MATERIAL) = params[PARTICLE_MATERIAL];
        }

        if (has_rotation && pnew_node->SolutionStepsDataHas(PARTICLE_ROTATION_DAMP_RATIO)) {
            pnew_node->FastGetSolutionStepValue(PARTICLE_ROTATION_DAMP_RATIO) = params[PARTICLE_ROTATION_DAMP_RATIO];
        }

        if (has_sphericity) {
            pnew_node->FastGetSolutionStepValue(PARTICLE_SPHERICITY) = params[PARTICLE_SPHERICITY];
        }
  
        pnew_node->FastGetSolutionStepValue(ANGULAR_VELOCITY) = null_vector;

        pnew_node->AddDof(VELOCITY_X, REACTION_X);
        pnew_node->AddDof(VELOCITY_Y, REACTION_Y);
        pnew_node->AddDof(VELOCITY_Z, REACTION_Z);
        pnew_node->AddDof(ANGULAR_VELOCITY_X, REACTION_X);
        pnew_node->AddDof(ANGULAR_VELOCITY_Y, REACTION_Y);
        pnew_node->AddDof(ANGULAR_VELOCITY_Z, REACTION_Z);

        pnew_node->pGetDof(VELOCITY_X)->FixDof();
        pnew_node->pGetDof(VELOCITY_Y)->FixDof();
        pnew_node->pGetDof(VELOCITY_Z)->FixDof();
        pnew_node->pGetDof(ANGULAR_VELOCITY_X)->FixDof();
        pnew_node->pGetDof(ANGULAR_VELOCITY_Y)->FixDof();
        pnew_node->pGetDof(ANGULAR_VELOCITY_Z)->FixDof();

        pnew_node->Set(DEMFlags::FIXED_VEL_X, true);
        pnew_node->Set(DEMFlags::FIXED_VEL_Y, true);
        pnew_node->Set(DEMFlags::FIXED_VEL_Z, true);
        pnew_node->Set(DEMFlags::FIXED_ANG_VEL_X, true);
        pnew_node->Set(DEMFlags::FIXED_ANG_VEL_Y, true);
        pnew_node->Set(DEMFlags::FIXED_ANG_VEL_Z, true);
        KRATOS_CATCH("")
    }

    Kratos::Element* ParticleCreatorDestructor::ElementCreatorWithPhysicalParameters(ModelPart& r_modelpart,
                                                                        int r_Elem_Id,
                                                                        Node < 3 > ::Pointer reference_node,
                                                                        Element::Pointer injector_element,
                                                                        Properties::Pointer r_params,
                                                                        ModelPart& r_sub_model_part_with_parameters,
                                                                        const Element& r_reference_element,
                                                                        PropertiesProxy* p_fast_properties,
                                                                        bool has_sphericity,
                                                                        bool has_rotation,
                                                                        bool initial,
                                                                        ElementsContainerType& array_of_injector_elements) {
        KRATOS_TRY
        Node<3>::Pointer pnew_node;

        double radius = r_sub_model_part_with_parameters[RADIUS];
        double max_radius = 1.5 * radius;

        if (initial) {
            radius = max_radius;
        } else {
            double std_deviation = r_sub_model_part_with_parameters[STANDARD_DEVIATION];
            std::string distribution_type = r_sub_model_part_with_parameters[PROBABILITY_DISTRIBUTION];
            double min_radius = 0.5 * radius;
            
            if (distribution_type == "normal") radius = rand_normal(radius, std_deviation, max_radius, min_radius);
            else if (distribution_type == "lognormal") radius = rand_lognormal(radius, std_deviation, max_radius, min_radius);
            else KRATOS_THROW_ERROR(std::runtime_error, "Unknown probability distribution in submodelpart ", r_sub_model_part_with_parameters.Name())
        }

        NodeCreatorWithPhysicalParameters(r_modelpart, pnew_node, r_Elem_Id, reference_node, radius, *r_params, r_sub_model_part_with_parameters, has_sphericity, has_rotation, initial);

        Geometry<Node<3> >::PointsArrayType nodelist;
        
        nodelist.push_back(pnew_node);
        Element::Pointer p_particle = r_reference_element.Create(r_Elem_Id, nodelist, r_params);
        SphericParticle* spheric_p_particle = dynamic_cast<SphericParticle*> (p_particle.get());

        if (initial) {
            array_of_injector_elements.push_back(p_particle);
            p_particle->Set(BLOCKED);
            pnew_node->Set(BLOCKED);
        }
        else {
            array_1d<double, 3 > zero_vector(3, 0.0);
            SphericParticle* injector_spheric_particle = dynamic_cast<SphericParticle*> (injector_element.get());

            if (mDoSearchNeighbourElements) { // there is no risk of contact so there is no need to track overlap
                injector_spheric_particle->mNeighbourElements.push_back(spheric_p_particle);
                injector_spheric_particle->mNeighbourElasticContactForces.push_back(zero_vector);
                spheric_p_particle->mNeighbourElements.push_back(injector_spheric_particle);
                spheric_p_particle->mNeighbourElasticContactForces.push_back(zero_vector);
            }
        }

        p_particle->Set(NEW_ENTITY);
        pnew_node->Set(NEW_ENTITY);
        //p_particle->InitializeSolutionStep(r_modelpart.GetProcessInfo());       

        spheric_p_particle->SetFastProperties(p_fast_properties);

        double density = spheric_p_particle->GetDensity();
        spheric_p_particle->SetRadius(radius);
        spheric_p_particle->SetSearchRadius(radius);
        spheric_p_particle->SetSearchRadiusWithFem(radius);
        double mass = 4.0 / 3.0 * KRATOS_M_PI * density * radius * radius * radius;
        spheric_p_particle->SetMass(mass);

        if (has_rotation) spheric_p_particle->Set(DEMFlags::HAS_ROTATION, true);
        else spheric_p_particle->Set(DEMFlags::HAS_ROTATION, false);

        spheric_p_particle->Initialize(r_modelpart.GetProcessInfo()); 

        #pragma omp critical
        {
            r_modelpart.Elements().push_back(p_particle);
        }
        
        return spheric_p_particle;
        KRATOS_CATCH("")
    }

    void ParticleCreatorDestructor::NodeCreatorForClusters(ModelPart& r_modelpart,
                                                            Node<3>::Pointer& pnew_node,
                                                            int aId,
                                                            array_1d<double, 3>& reference_coordinates,
                                                            double radius,
                                                            Properties& params) {
        KRATOS_TRY
        pnew_node = boost::make_shared< Node<3> >( aId, reference_coordinates[0], reference_coordinates[1], reference_coordinates[2] );
        pnew_node->SetSolutionStepVariablesList(&r_modelpart.GetNodalSolutionStepVariablesList());
        pnew_node->SetBufferSize(r_modelpart.GetBufferSize());

        #pragma omp critical
        {
        //pnew_node = r_modelpart.CreateNewNode(aId, reference_coordinates[0], reference_coordinates[1], reference_coordinates[2]); //ACTUAL node creation and addition to model part
            r_modelpart.Nodes().push_back(pnew_node);
        }
        
        pnew_node->FastGetSolutionStepValue(RADIUS) = radius;
        array_1d<double, 3 > null_vector(3, 0.0);
        pnew_node->FastGetSolutionStepValue(VELOCITY) = null_vector;
        pnew_node->FastGetSolutionStepValue(ANGULAR_VELOCITY) = null_vector;
        pnew_node->FastGetSolutionStepValue(PARTICLE_MATERIAL) = params[PARTICLE_MATERIAL];

        pnew_node->AddDof(VELOCITY_X, REACTION_X);
        pnew_node->AddDof(VELOCITY_Y, REACTION_Y);
        pnew_node->AddDof(VELOCITY_Z, REACTION_Z);
        pnew_node->AddDof(ANGULAR_VELOCITY_X, REACTION_X);
        pnew_node->AddDof(ANGULAR_VELOCITY_Y, REACTION_Y);
        pnew_node->AddDof(ANGULAR_VELOCITY_Z, REACTION_Z);

        pnew_node->pGetDof(VELOCITY_X)->FixDof();
        pnew_node->pGetDof(VELOCITY_Y)->FixDof();
        pnew_node->pGetDof(VELOCITY_Z)->FixDof();
        pnew_node->pGetDof(ANGULAR_VELOCITY_X)->FixDof();
        pnew_node->pGetDof(ANGULAR_VELOCITY_Y)->FixDof();
        pnew_node->pGetDof(ANGULAR_VELOCITY_Z)->FixDof();

        pnew_node->Set(DEMFlags::FIXED_VEL_X, true);
        pnew_node->Set(DEMFlags::FIXED_VEL_Y, true);
        pnew_node->Set(DEMFlags::FIXED_VEL_Z, true);
        pnew_node->Set(DEMFlags::FIXED_ANG_VEL_X, true);
        pnew_node->Set(DEMFlags::FIXED_ANG_VEL_Y, true);
        pnew_node->Set(DEMFlags::FIXED_ANG_VEL_Z, true);
        pnew_node->Set(DEMFlags::BELONGS_TO_A_CLUSTER, true);
        
        //pnew_node->FastGetSolutionStepValue(SPRAYED_MATERIAL) = 0.0;
        KRATOS_CATCH("")
    }

    Kratos::SphericParticle* ParticleCreatorDestructor::SphereCreatorForClusters(ModelPart& r_modelpart,
                                                                                  int r_Elem_Id,
                                                                                  double radius,
                                                                                  array_1d<double, 3>& reference_coordinates,
                                                                                  double cluster_mass,
                                                                                  Properties::Pointer r_params,
                                                                                  const Element& r_reference_element,
                                                                                  const int cluster_id, 
                                                                                  PropertiesProxy* p_fast_properties) {
        KRATOS_TRY
        Node<3>::Pointer pnew_node;
        
        NodeCreatorForClusters(r_modelpart, pnew_node, r_Elem_Id, reference_coordinates, radius, *r_params);

        Geometry<Node <3> >::PointsArrayType nodelist;
        nodelist.push_back(pnew_node);

        Element::Pointer p_particle = r_reference_element.Create(r_Elem_Id, nodelist, r_params);

        Kratos::SphericParticle* spheric_p_particle = dynamic_cast<Kratos::SphericParticle*> (p_particle.get());

        spheric_p_particle->SetFastProperties(p_fast_properties);
        spheric_p_particle->Initialize(r_modelpart.GetProcessInfo());        
        spheric_p_particle->SetRadius(radius);
        spheric_p_particle->SetSearchRadius(radius);
        spheric_p_particle->SetSearchRadiusWithFem(radius);        
        spheric_p_particle->SetMass(cluster_mass);        
        spheric_p_particle->Set(DEMFlags::HAS_ROLLING_FRICTION, false);
        spheric_p_particle->Set(DEMFlags::BELONGS_TO_A_CLUSTER, true);
        spheric_p_particle->SetClusterId(cluster_id);        
        spheric_p_particle->CreateDiscontinuumConstitutiveLaws(r_modelpart.GetProcessInfo());                      

        #pragma omp critical
        {        
            r_modelpart.Elements().push_back(p_particle);
        }
        
        return spheric_p_particle;
        KRATOS_CATCH("")
    }

Kratos::SphericParticle* ParticleCreatorDestructor::SphereCreatorForBreakableClusters(ModelPart& r_modelpart,
                                                                                  int r_Elem_Id,
                                                                                  double radius,
                                                                                  array_1d<double, 3>& reference_coordinates,
                                                                                  Properties::Pointer r_params,
                                                                                  const Element& r_reference_element,
                                                                                  const int cluster_id, 
                                                                                  PropertiesProxy* p_fast_properties) {
        KRATOS_TRY
        Node<3>::Pointer pnew_node;
        
        NodeCreatorForClusters(r_modelpart, pnew_node, r_Elem_Id, reference_coordinates, radius, *r_params);

        Geometry<Node <3> >::PointsArrayType nodelist;
        nodelist.push_back(pnew_node);

        Element::Pointer p_particle = r_reference_element.Create(r_Elem_Id, nodelist, r_params);

        Kratos::SphericParticle* spheric_p_particle = dynamic_cast<Kratos::SphericParticle*> (p_particle.get());

        spheric_p_particle->SetFastProperties(p_fast_properties);
        spheric_p_particle->Initialize(r_modelpart.GetProcessInfo()); 
        spheric_p_particle->SetRadius(radius);
        spheric_p_particle->SetSearchRadius(radius);
        spheric_p_particle->SetSearchRadiusWithFem(radius);        
        spheric_p_particle->SetMass(spheric_p_particle->GetDensity() * spheric_p_particle->CalculateVolume());
        if (spheric_p_particle->Is(DEMFlags::HAS_ROTATION)) {
            spheric_p_particle->GetGeometry()[0].FastGetSolutionStepValue(PARTICLE_MOMENT_OF_INERTIA) = spheric_p_particle->CalculateMomentOfInertia();
        }               
        spheric_p_particle->Set(DEMFlags::HAS_ROLLING_FRICTION, false);
        spheric_p_particle->Set(DEMFlags::BELONGS_TO_A_CLUSTER, true);
        spheric_p_particle->SetClusterId(-1);        
        spheric_p_particle->CreateDiscontinuumConstitutiveLaws(r_modelpart.GetProcessInfo());

        #pragma omp critical
        {        
            r_modelpart.Elements().push_back(p_particle);
        }
        return spheric_p_particle;
        KRATOS_CATCH("")
    }       
    
    void ParticleCreatorDestructor::ClusterCreatorWithPhysicalParameters(ModelPart& r_spheres_modelpart,
                                                                         ModelPart& r_clusters_modelpart,                                                                         
                                                                         int r_Elem_Id,
                                                                         Node<3>::Pointer reference_node,
                                                                         Element::Pointer injector_element,
                                                                         Properties::Pointer r_params,
                                                                         ModelPart& r_sub_model_part_with_parameters,
                                                                         const Element& r_reference_element,
                                                                         PropertiesProxy* p_fast_properties,
                                                                         bool has_sphericity,
                                                                         bool has_rotation,
                                                                         ElementsContainerType& array_of_injector_elements,
                                                                         int& number_of_added_spheres) {
        KRATOS_TRY
        
        ProcessInfo& r_process_info = r_spheres_modelpart.GetProcessInfo();
        
        Node<3>::Pointer pnew_node;

        double radius = r_sub_model_part_with_parameters[RADIUS];
        double min_radius = 0.5 * radius; //TODO: this is a little bit arbitrary
        double max_radius = 1.5 * radius; //TODO: this is a little bit arbitrary

        double OrientationReal;
        array_1d<double, 3> OrientationImag;
        
        double std_deviation = r_sub_model_part_with_parameters[STANDARD_DEVIATION];
        std::string distribution_type = r_sub_model_part_with_parameters[PROBABILITY_DISTRIBUTION];

        if (distribution_type == "normal") radius = rand_normal(radius, std_deviation, max_radius, min_radius);
        else if (distribution_type == "lognormal") radius = rand_lognormal(radius, std_deviation, max_radius, min_radius);
        else KRATOS_THROW_ERROR(std::runtime_error, "Unknown probability distribution in submodelpart ", r_sub_model_part_with_parameters.Name())

        NodeForClustersCreatorWithPhysicalParameters(r_clusters_modelpart, pnew_node, r_Elem_Id, reference_node, *r_params, r_sub_model_part_with_parameters, has_sphericity, has_rotation, false);
        
        pnew_node->FastGetSolutionStepValue(CHARACTERISTIC_LENGTH) = radius * 2.0; //Cluster specific. Can be removed
        
        Geometry<Node<3> >::PointsArrayType nodelist;

        nodelist.push_back(pnew_node);

        Element::Pointer p_new_cluster = r_reference_element.Create(r_Elem_Id, nodelist, r_params);
        Kratos::Cluster3D* p_cluster = dynamic_cast<Kratos::Cluster3D*> (p_new_cluster.get());
                    
        if (r_sub_model_part_with_parameters[RANDOM_ORIENTATION]) {
            OrientationReal    = ((double) rand() / (RAND_MAX));
            OrientationImag[0] = ((double) rand() / (RAND_MAX));
            OrientationImag[1] = ((double) rand() / (RAND_MAX));
            OrientationImag[2] = ((double) rand() / (RAND_MAX));
        }
        else {
            OrientationReal = r_sub_model_part_with_parameters[ORIENTATION_REAL];
            OrientationImag[0] = r_sub_model_part_with_parameters[ORIENTATION_IMAG][0];
            OrientationImag[1] = r_sub_model_part_with_parameters[ORIENTATION_IMAG][1];
            OrientationImag[2] = r_sub_model_part_with_parameters[ORIENTATION_IMAG][2];
        }        

        p_cluster->SetOrientation(OrientationReal, OrientationImag);
        p_cluster->Initialize(r_process_info); 
        
        const bool is_breakable = (*r_params)[BREAKABLE_CLUSTER]; //THIS IS NOT THREAD SAFE!!!
        
        if (!is_breakable) {  
            mMaxNodeId++; //This must be done before CreateParticles because the creation of particles accesses mMaxNodeId to choose what Id is assigned to the new nodes/spheres
        }
        
        ParticleCreatorDestructor* creator_destructor_ptr = this;
        p_cluster->CreateParticles(creator_destructor_ptr, r_spheres_modelpart, p_fast_properties);                        
        
        number_of_added_spheres = p_cluster->GetNumberOfSpheres();
        
        SphericParticle* injector_spheric_particle = dynamic_cast<SphericParticle*> (injector_element.get());

	if (has_rotation) p_cluster->Set(DEMFlags::HAS_ROTATION, true);
        else p_cluster->Set(DEMFlags::HAS_ROTATION, false);

        if (!is_breakable) {  
            #pragma omp critical
            {
            r_clusters_modelpart.Elements().push_back(p_new_cluster);  //We add the cluster element to the clusters modelpart               
            }
        }
        else { 
            p_cluster->GetGeometry()[0].Set(TO_ERASE); //We do not add the cluster to the modelpart (will be erased at the end of this function) and we mark the central node for erasing (none are needed)
            p_cluster->SetContinuumGroupToBreakableClusterSpheres(r_Elem_Id);
            double search_tolerance = 0.02 * radius;
            p_cluster->SetInitialNeighbours(search_tolerance);
            p_cluster->CreateContinuumConstitutiveLaws();
            p_cluster->SetInitialConditionsToSpheres(r_sub_model_part_with_parameters[VELOCITY]);            
        }
        
        array_1d<double, 3 > zero_vector(3, 0.0);
        
        for (unsigned int i=0; i<p_cluster->GetNumberOfSpheres(); i++) { 
            SphericParticle* spheric_p_particle = p_cluster->GetSpheres()[i];

            if (mDoSearchNeighbourElements) { // there is no risk of contact so there is no need to track overlap
                injector_spheric_particle->mNeighbourElements.push_back(spheric_p_particle);
                injector_spheric_particle->mNeighbourElasticContactForces.push_back(zero_vector);
                spheric_p_particle->mNeighbourElements.push_back(injector_spheric_particle);
                spheric_p_particle->mNeighbourElasticContactForces.push_back(zero_vector);
            }

            spheric_p_particle->Set(NEW_ENTITY);
            spheric_p_particle->GetGeometry()[0].Set(NEW_ENTITY);
            spheric_p_particle->SetFastProperties(p_fast_properties);
        }

        p_cluster->Set(NEW_ENTITY);
        pnew_node->Set(NEW_ENTITY);
        KRATOS_CATCH("")
    }


    void ParticleCreatorDestructor::CalculateSurroundingBoundingBox(ModelPart& r_balls_model_part,
                                                                    ModelPart& r_clusters_model_part,
                                                                    ModelPart& r_rigid_faces_model_part,
                                                                    double scale_factor,
                                                                    bool automatic) {
        KRATOS_TRY

        if (automatic) {

            double ref_radius = 0.0;

            if (r_balls_model_part.NumberOfElements(0) == 0 && r_clusters_model_part.NumberOfElements(0) == 0 && r_rigid_faces_model_part.NumberOfElements(0) == 0) {
                KRATOS_THROW_ERROR(std::logic_error, "The Bounding Box cannot be calculated automatically when there are no elements. Kratos stops.", "");
            }

            if (scale_factor < 0.0) {
                KRATOS_THROW_ERROR(std::logic_error, "The enlargement factor for the automatic calculation of the bounding box must be a positive value.", "");
            }

            if (scale_factor < 1.0) {
                std::cout << "\n WARNING" + std::string(2, '\n');
                std::cout << "The enlargement factor for the automatic calculation of the bounding box is less than 1!. \n";
                std::cout << "It is not guaranteed that the model fits inside of it." + std::string(2, '\n');
            }

            if (r_balls_model_part.NumberOfElements(0)) { // loop over spheric elements (balls)
                Configure::ElementsContainerType Elements = r_balls_model_part.GetCommunicator().LocalMesh().Elements();

                ref_radius = (*(Elements.begin().base()))->GetGeometry()[0].FastGetSolutionStepValue(RADIUS);
                const array_1d<double, 3 >& ini_coor = (*(Elements.begin().base()))->GetGeometry()[0].Coordinates();
                mStrictLowPoint = ini_coor;
                mStrictHighPoint = ini_coor;

                for (Configure::ElementsContainerType::iterator particle_pointer_it = Elements.begin(); particle_pointer_it != Elements.end(); ++particle_pointer_it) {
                    const array_1d<double, 3 >& coor = (*(particle_pointer_it.base()))->GetGeometry()[0].Coordinates();
                    double radius = particle_pointer_it->GetGeometry()[0].FastGetSolutionStepValue(RADIUS);
                    ref_radius = (ref_radius < radius) ? radius : ref_radius;

                    for (std::size_t i = 0; i < 3; i++) {
                        mStrictLowPoint[i] = (mStrictLowPoint[i] > coor[i]) ? coor[i] : mStrictLowPoint[i];
                        mStrictHighPoint[i] = (mStrictHighPoint[i] < coor[i]) ? coor[i] : mStrictHighPoint[i];
                    }
                }
            }
                        
                    
            if (r_rigid_faces_model_part.NumberOfConditions(0)) { // loop over rigid faces
                ModelPart::ConditionsContainerType Conditions = r_rigid_faces_model_part.GetCommunicator().LocalMesh().Conditions();

                std::vector <array_1d<double, 3 > > face_coor;
                std::size_t n_nodes = (*(Conditions.begin().base()))->GetGeometry().size();
                face_coor.resize(n_nodes);

                for (std::size_t i = 0; i < n_nodes; ++i) {
                    face_coor[i] = (*(Conditions.begin().base()))->GetGeometry()[i].Coordinates();
                }

                if (r_balls_model_part.NumberOfElements(0) == 0) { // initialize if not initialized already
                    mStrictLowPoint = face_coor[0];
                    mStrictHighPoint = face_coor[0];
                }

                for (ModelPart::ConditionsContainerType::iterator particle_pointer_it = Conditions.begin(); particle_pointer_it != Conditions.end(); ++particle_pointer_it) {
                    std::size_t n_cond_nodes = (*(particle_pointer_it.base()))->GetGeometry().size();
                    face_coor.resize(n_cond_nodes);

                    for (std::size_t i = 0; i < n_cond_nodes; ++i) {
                        face_coor[i] = (*(particle_pointer_it.base()))->GetGeometry()[i].Coordinates();

                        for (std::size_t j = 0; j < 3; j++) {
                            mStrictLowPoint[j] = (mStrictLowPoint[j] > face_coor[i][j]) ? face_coor[i][j] : mStrictLowPoint[j];
                            mStrictHighPoint[j] = (mStrictHighPoint[j] < face_coor[i][j]) ? face_coor[i][j] : mStrictHighPoint[j];
                        }
                    }
                }
            }
                        
            r_balls_model_part.GetCommunicator().MinAll(mStrictLowPoint[0]);
            r_balls_model_part.GetCommunicator().MinAll(mStrictLowPoint[1]);
            r_balls_model_part.GetCommunicator().MinAll(mStrictLowPoint[2]);
            r_balls_model_part.GetCommunicator().MaxAll(mStrictHighPoint[0]);
            r_balls_model_part.GetCommunicator().MaxAll(mStrictHighPoint[1]);
            r_balls_model_part.GetCommunicator().MaxAll(mStrictHighPoint[2]);
            r_balls_model_part.GetCommunicator().MaxAll(ref_radius);
            
             
            array_1d<double, 3 > midpoint;
            noalias(midpoint) = 0.5 * (mStrictHighPoint + mStrictLowPoint);
            mHighPoint = midpoint * (1 - scale_factor) + scale_factor * mStrictHighPoint;
            mLowPoint = midpoint * (1 - scale_factor) + scale_factor * mStrictLowPoint;
                    
            for (std::size_t i = 0; i < 3; i++) {
                mLowPoint[i] -= 2 * ref_radius;
                mHighPoint[i] += 2 * ref_radius;
            }
            
            mStrictDiameter = norm_2(mStrictHighPoint - mStrictLowPoint);
            mDiameter = norm_2(mHighPoint - mLowPoint);
        }// if (automatic)

        else {
            for (int i = 0; i < 3; ++i) {

                if (mHighPoint[i] < mLowPoint[i]) {
                    KRATOS_THROW_ERROR(std::logic_error, "Check limits of the Bounding Box, minimum coordinates exceed maximum coordinates.", "");
                }
            }

            mStrictHighPoint = mHighPoint; // mHighPoint and mLowPoint have been set as an input value
            mStrictLowPoint = mLowPoint;
            mStrictDiameter = norm_2(mStrictHighPoint - mStrictLowPoint);
            mDiameter = norm_2(mHighPoint - mLowPoint);
        }

        KRATOS_CATCH("")
    }    
    
    void ParticleCreatorDestructor::DestroyParticles(ModelPart& r_model_part) {

        KRATOS_TRY
        

        DestroyParticles(r_model_part.GetCommunicator().LocalMesh());
        DestroyParticles(r_model_part.GetCommunicator().GhostMesh());
          
        KRATOS_CATCH("")
    }    
    
    void ParticleCreatorDestructor::DestroyParticles(ModelPart::MeshType& rMesh) {

        KRATOS_TRY
        
        ElementsArrayType& rElements = rMesh.Elements();
        ModelPart::NodesContainerType& rNodes = rMesh.Nodes();       

        if (rElements.size() != rNodes.size()) {
            KRATOS_THROW_ERROR(std::runtime_error, "While removing elements and nodes, the number of elements and the number of nodes are not the same in the ModelPart!", 0);
        }
                                
        int good_elems_counter = 0;
        
        for (int k = 0; k < (int)rElements.size(); k++) {
            Configure::ElementsContainerType::ptr_iterator element_pointer_it = rElements.ptr_begin() + k;
            ModelPart::NodeType& node = (*element_pointer_it)->GetGeometry()[0];
                        
            if (node.IsNot(TO_ERASE) && (*element_pointer_it)->IsNot(TO_ERASE)) {
  	        if (k != good_elems_counter) {
                    *(rElements.ptr_begin() + good_elems_counter) = boost::move(*element_pointer_it);
                }
                good_elems_counter++;                                
            }
            else {
                (*element_pointer_it).reset();
                node.Set(TO_ERASE, true);
            }            
        }   
        int good_nodes_counter = 0;
        
        for (int k = 0; k < (int)rNodes.size(); k++) {
            ModelPart::NodesContainerType::ptr_iterator node_pointer_it = rNodes.ptr_begin() + k;
            if ((*node_pointer_it)->IsNot(TO_ERASE)) {
  	        if (k != good_nodes_counter) {
                    *(rNodes.ptr_begin() + good_nodes_counter) = boost::move(*node_pointer_it);
                }
                good_nodes_counter++;                                
            }

            else (*node_pointer_it).reset();
        }
        
        if (good_elems_counter != good_nodes_counter) {
            KRATOS_THROW_ERROR(std::runtime_error, "While removing elements and nodes, the number of removed elements and the number of removed nodes were not the same!", 0);
        }               
          
        if ((int)rElements.size() != good_elems_counter) {
            rElements.erase(rElements.ptr_begin() + good_elems_counter, rElements.ptr_end());
        }
        
        if ((int)rNodes.size() != good_nodes_counter) {
            rNodes.erase(rNodes.ptr_begin() + good_nodes_counter, rNodes.ptr_end());  
        }    
        KRATOS_CATCH("")
    }  


    void ParticleCreatorDestructor::DestroyContactElements(ModelPart& r_model_part) {
        KRATOS_TRY                        
        
        ElementsArrayType& rElements = r_model_part.GetCommunicator().LocalMesh().Elements();
                            
        int good_elems_counter = 0;
        
        for (int k = 0; k < (int)rElements.size(); k++) {
            Configure::ElementsContainerType::ptr_iterator element_pointer_it = rElements.ptr_begin() + k;
                        
            if ((*element_pointer_it)->IsNot(TO_ERASE)) {
  	        if (k != good_elems_counter) {
                    *(rElements.ptr_begin() + good_elems_counter) = boost::move(*element_pointer_it);
                }
                good_elems_counter++;                                
            }
            else (*element_pointer_it).reset();
        }   
        
        if ((int)rElements.size() != good_elems_counter) rElements.erase(rElements.ptr_begin() + good_elems_counter, rElements.ptr_end());
                
        KRATOS_CATCH("")
    }
    
    void ParticleCreatorDestructor::RemoveUnusedNodesOfTheClustersModelPart(ModelPart& r_clusters_modelpart) {
        KRATOS_TRY
        ModelPart::NodesContainerType& rNodes = r_clusters_modelpart.GetCommunicator().LocalMesh().Nodes();   
        int good_nodes_counter = 0;
        
        for(int k=0; k < (int)rNodes.size(); k++) {
            ModelPart::NodesContainerType::ptr_iterator node_pointer_it = rNodes.ptr_begin() + k;
            if ((*node_pointer_it)->IsNot(TO_ERASE)) {
  	        if(k != good_nodes_counter){
                    *(rNodes.ptr_begin() + good_nodes_counter) = boost::move(*node_pointer_it);
                }
                good_nodes_counter++;                                
            }
            else{ (*node_pointer_it).reset(); }                
        }
        
        if((int)rNodes.size() != good_nodes_counter){
            rNodes.erase(rNodes.ptr_begin() + good_nodes_counter, rNodes.ptr_end());  
        }  
        KRATOS_CATCH("")
    }

    void ParticleCreatorDestructor::MarkDistantParticlesForErasing(ModelPart& r_model_part) {
        MarkParticlesForErasingGivenBoundingBox(r_model_part, mLowPoint, mHighPoint);
    }

    void ParticleCreatorDestructor::MarkParticlesForErasingGivenScalarVariableValue(ModelPart& r_model_part, const Variable<double>& rVariable, double value, double tol) {

        KRATOS_TRY

        Configure::ElementsContainerType& rElements = r_model_part.GetCommunicator().LocalMesh().Elements();
        
        #pragma omp parallel for
        for(int k=0; k<(int)rElements.size(); k++){
            Configure::ElementsContainerType::ptr_iterator particle_pointer_it = rElements.ptr_begin() + k;                

            const double& i_value = (*particle_pointer_it)->GetGeometry()[0].FastGetSolutionStepValue(rVariable);
            bool include = true; // = (erase_flag < 0.5);

            include = include && ((i_value <= value - fabs(tol)) || (i_value >= value + fabs(tol)));

            if (include)
                (*particle_pointer_it)->GetGeometry()[0].Set(TO_ERASE);

        }
        KRATOS_CATCH("")
    }

    void ParticleCreatorDestructor::MarkParticlesForErasingGivenVectorVariableModulus(ModelPart& r_model_part, const Variable<array_1d<double, 3 > >& rVariable, double value, double tol) {

        KRATOS_TRY

        Configure::ElementsContainerType& rElements = r_model_part.GetCommunicator().LocalMesh().Elements();

        #pragma omp parallel for
        for (int k = 0; k < (int)rElements.size(); k++) {
            Configure::ElementsContainerType::ptr_iterator particle_pointer_it = rElements.ptr_begin() + k;  

            array_1d<double, 3 > & i_var = (*particle_pointer_it)->GetGeometry()[0].FastGetSolutionStepValue(rVariable);
            double i_value = sqrt(i_var[0] * i_var[0] + i_var[1] * i_var[1] + i_var[2] * i_var[2]);
            bool include = true; //  = (erase_flag < 0.5);

            include = include && ((i_value <= value - fabs(tol)) || (i_value >= value + fabs(tol)));

            if (include)
                (*particle_pointer_it)->GetGeometry()[0].Set(TO_ERASE);

        }

        KRATOS_CATCH("")

    }

    void ParticleCreatorDestructor::MarkParticlesForErasingGivenBoundingBox(ModelPart& r_model_part, array_1d<double, 3 > low_point, array_1d<double, 3 > high_point) {

        KRATOS_TRY

        ModelPart::NodesContainerType& rNodes = r_model_part.GetCommunicator().LocalMesh().Nodes();
        Configure::ElementsContainerType& rElements = r_model_part.GetCommunicator().LocalMesh().Elements();

        #pragma omp parallel
        {
            #pragma omp for
            for (int k = 0; k < (int)rElements.size(); k++){
                Configure::ElementsContainerType::ptr_iterator particle_pointer_it = rElements.ptr_begin() + k; 

                if ((*particle_pointer_it)->Is(DEMFlags::BELONGS_TO_A_CLUSTER)) continue;
                if ((*particle_pointer_it)->Is(BLOCKED)) continue;

                const array_1d<double, 3 >& coor = (*particle_pointer_it)->GetGeometry()[0].Coordinates();
                bool include = true;

                for (unsigned int i = 0; i < 3; i++) {
                    include = include && (coor[i] >= low_point[i]) && (coor[i] <= high_point[i]);
                }

                if (!include) {
                    (*particle_pointer_it)->GetGeometry()[0].Set(TO_ERASE);
                    (*particle_pointer_it)->Set(TO_ERASE);
                }

            }

            #pragma omp for
            for (int k = 0; k < (int)rNodes.size(); k++) { //TODO: Can we remove this loop? For DEM it is useless. The job was done in the previous loop. Try it.
                ModelPart::NodesContainerType::ptr_iterator node_pointer_it = rNodes.ptr_begin() + k;                

                if ((*node_pointer_it)->Is(DEMFlags::BELONGS_TO_A_CLUSTER)) continue;
                if ((*node_pointer_it)->Is(BLOCKED)) continue;
                const array_1d<double, 3 >& coor = (*node_pointer_it)->Coordinates();
                bool include = true;

                for (unsigned int i = 0; i < 3; i++) {
                    include = include && (coor[i] >= low_point[i]) && (coor[i] <= high_point[i]);
                }

                if (!include) {
                    (*node_pointer_it)->Set(TO_ERASE);
                }

            }
        }

        KRATOS_CATCH("")
    }

    void ParticleCreatorDestructor::MarkContactElementsForErasing(ModelPart& r_model_part, ModelPart& mcontacts_model_part) {
        KRATOS_TRY

        Configure::ElementsContainerType& rElements = r_model_part.GetCommunicator().LocalMesh().Elements();

        #pragma omp parallel for
        for (int k = 0; k < (int)rElements.size(); k++) {
            Configure::ElementsContainerType::ptr_iterator particle_pointer_it = rElements.ptr_begin() + k; 

            if ((*particle_pointer_it)->GetGeometry()[0].Is(TO_ERASE)) {
                Element* p_element = particle_pointer_it->get();
                SphericContinuumParticle* p_continuum_spheric_particle = dynamic_cast<SphericContinuumParticle*> (p_element);
                std::vector<ParticleContactElement*>& array_of_bonds = p_continuum_spheric_particle->mBondElements;
                for (unsigned int i = 0; i < array_of_bonds.size(); i++) {
                    if (array_of_bonds[i] != NULL) { //NULL happens when the initial neighbor was a ghost and had a lower Id than the others
                        array_of_bonds[i]->Set(TO_ERASE);
                    }
                }
            }
        }

        KRATOS_CATCH("")
    }

    void ParticleCreatorDestructor::DestroyParticlesOutsideBoundingBox(ModelPart& r_model_part) {
        KRATOS_TRY
        MarkDistantParticlesForErasing(r_model_part);
        DestroyParticles(r_model_part);
        KRATOS_CATCH("")
    }

    void ParticleCreatorDestructor::MoveParticlesOutsideBoundingBoxBackInside(ModelPart& r_model_part) {
        KRATOS_TRY

        ModelPart::NodesContainerType& rNodes = r_model_part.GetCommunicator().LocalMesh().Nodes();

        #pragma omp parallel for
        for (int k = 0; k < (int)rNodes.size(); k++) { 
            ModelPart::NodesContainerType::ptr_iterator node_pointer_it = rNodes.ptr_begin() + k;

            array_1d<double, 3 >& coor = (*node_pointer_it)->Coordinates();
            array_1d<double, 3 >& displ = (*node_pointer_it)->FastGetSolutionStepValue(DISPLACEMENT);
            double period_0 = mHighPoint[0] - mLowPoint[0];
            double period_1 = mHighPoint[1] - mLowPoint[1];
            double period_2 = mHighPoint[2] - mLowPoint[2];

            if (coor[0] > mHighPoint[0]) {
                displ[0] -= period_0;
                coor[0] -= period_0;
            }
            else if (coor[0] < mLowPoint[0]) {
                displ[0] += period_0;
                coor[0] += period_0;
            }
            if (coor[1] > mHighPoint[1]) {
                displ[1] -= period_1;
                coor[1] -= period_1;
            }
            else if (coor[1] < mLowPoint[1]) {
                displ[1] += period_1;
                coor[1] += period_1;
            }
            if (coor[2] > mHighPoint[2]) {
                displ[2] -= period_2;
                coor[2] -= period_2;
            }
            else if (coor[2] < mLowPoint[2]) {
                displ[2] += period_2;
                coor[2] += period_2;
            }
        }
        KRATOS_CATCH("")
    }

    void ParticleCreatorDestructor::DestroyContactElementsOutsideBoundingBox(ModelPart& r_model_part, ModelPart& mcontacts_model_part) {
        KRATOS_TRY
        MarkContactElementsForErasing(r_model_part, mcontacts_model_part);
        DestroyContactElements(mcontacts_model_part);
        KRATOS_CATCH("")
    }

    array_1d<double, 3 > ParticleCreatorDestructor::GetHighNode() { return (mHighPoint);}
    array_1d<double, 3 > ParticleCreatorDestructor::GetLowNode() { return (mLowPoint); }
    array_1d<double, 3 > ParticleCreatorDestructor::GetStrictHighNode() { return (mStrictHighPoint); }
    array_1d<double, 3 > ParticleCreatorDestructor::GetStrictLowNode() { return (mStrictLowPoint); }
    double ParticleCreatorDestructor::GetDiameter() { return (mDiameter); }
    double ParticleCreatorDestructor::GetStrictDiameter() { return (mStrictDiameter); }
    void ParticleCreatorDestructor::SetHighNode(array_1d<double, 3 > node) { mHighPoint = node; }
    void ParticleCreatorDestructor::SetLowNode(array_1d<double, 3 > node) { mLowPoint = node; }
    unsigned int ParticleCreatorDestructor::GetCurrentMaxNodeId() { return mMaxNodeId; } 
    unsigned int* ParticleCreatorDestructor::pGetCurrentMaxNodeId() { return &mMaxNodeId; }
    void ParticleCreatorDestructor::SetMaxNodeId(unsigned int id) { mMaxNodeId = id; }

    std::string ParticleCreatorDestructor::Info() const { return ""; }
    void ParticleCreatorDestructor::PrintInfo(std::ostream& rOStream) const {}
    void ParticleCreatorDestructor::PrintData(std::ostream& rOStream) const {}

    void ParticleCreatorDestructor::Clear(ModelPart::NodesContainerType::iterator node_it, int step_data_size) {
        KRATOS_TRY
        unsigned int buffer_size = node_it->GetBufferSize();
        for (unsigned int step = 0; step < buffer_size; step++) {
            //getting the data of the solution step
            double* step_data = (node_it)->SolutionStepData().Data(step);
            //copying this data in the position of the vector we are interested in
            for (int j = 0; j < step_data_size; j++) {
                step_data[j] = 0.0;
            }
        }
        KRATOS_CATCH("")
    }

    void ParticleCreatorDestructor::ClearVariables(ModelPart::NodesContainerType::iterator node_it, Variable<array_1d<double, 3 > >& rVariable) {
        KRATOS_TRY
        array_1d<double, 3 > & Aux_var = node_it->FastGetSolutionStepValue(rVariable);
        noalias(Aux_var) = ZeroVector(3);
        KRATOS_CATCH("")
    }

    void ParticleCreatorDestructor::ClearVariables(ParticleIterator particle_it, Variable<double>& rVariable) {}

} //Namespace Kratos
