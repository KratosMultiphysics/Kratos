// Author: 
// Miguel Angel Celigueta maceli@cimne.upc.edu
//
//README::::look to the key word "VERSION" if you want to find all the points where you have to change something so that you can pass from a kdtree to a bin data search structure;

// System includes
#include <string>
#include <iostream>
#include "DEM_application.h"
#include "create_and_destroy.h"

namespace Kratos {

    static double rand_normal(const double mean, const double stddev, const double max_radius, const double min_radius) {

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
    }

    static double rand_lognormal(const double mean, const double stddev, const double max_radius, const double min_radius){
        const double normal_mean = log(mean * mean / sqrt(stddev * stddev + mean * mean));
        const double normal_stddev = sqrt(log(1 + stddev * stddev / (mean * mean)));
        double normally_distributed_value = rand_normal(normal_mean, normal_stddev, log(max_radius), log(min_radius));

        return exp(normally_distributed_value);
    }

    static void AddRandomPerpendicularVelocityToGivenVelocity(array_1d<double, 3 >& velocity, const double angle_in_degrees){
        
        double velocity_modulus = sqrt(velocity[0]*velocity[0] + velocity[1]*velocity[1] + velocity[2]*velocity[2]);
        array_1d<double, 3 > unitary_velocity;
        noalias(unitary_velocity) = velocity / velocity_modulus;
        array_1d<double, 3 > normal_1;
        array_1d<double, 3 > normal_2;
        
        if(fabs(unitary_velocity[0])>=0.577) {
            normal_1[0]= - unitary_velocity[1];
            normal_1[1]= unitary_velocity[0];
            normal_1[2]= 0.0;
        }
        else if(fabs(unitary_velocity[1])>=0.577) {
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
        double inv_distance0 = (distance0 != 0.0) ?  1.0 / distance0 : 0.00;
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
        array_1d<double, 3 > local_added_velocity; local_added_velocity[0] = local_added_velocity[1] = local_added_velocity[2] = 0.0;
        
        while (local_added_velocity_modulus_square > radius_square) {
            //Random in a range: (max - min) * ( (double)rand() / (double)RAND_MAX ) + min
            local_added_velocity[0] = 2*radius * (double)rand() / (double)RAND_MAX - radius;
            local_added_velocity[1] = 2*radius * (double)rand() / (double)RAND_MAX - radius;
            local_added_velocity_modulus_square = local_added_velocity[0]*local_added_velocity[0] + local_added_velocity[1]*local_added_velocity[1];
        }
        
        noalias(velocity) += local_added_velocity[0] * normal_1 + local_added_velocity[1] * normal_2;
        
        
    }

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
    }

    int ParticleCreatorDestructor::FindMaxNodeIdInModelPart(ModelPart& r_modelpart) {

        int max_Id = 1; //GID accepts Id's >= 1
        std::vector<int> thread_maximums(OpenMPUtils::GetNumThreads(),1);

        #pragma omp parallel for
        for(int i=0; i<(int)r_modelpart.GetCommunicator().LocalMesh().Nodes().size(); i++){
            ModelPart::NodesContainerType::iterator node_it = r_modelpart.GetCommunicator().LocalMesh().NodesBegin() + i;        
            if ((int) (node_it->Id()) > thread_maximums[OpenMPUtils::ThisThread()]) thread_maximums[OpenMPUtils::ThisThread()] = node_it->Id();
        }
        
        for(int i=0; i<OpenMPUtils::GetNumThreads(); i++){
            if(thread_maximums[i] > max_Id) max_Id = thread_maximums[i];
        }
        
        r_modelpart.GetCommunicator().MaxAll(max_Id);
        return max_Id;
    }
    
    void ParticleCreatorDestructor::FindAndSaveMaxNodeIdInModelPart(ModelPart& r_modelpart) {
        mMaxNodeId = FindMaxNodeIdInModelPart(r_modelpart);
    }

    int ParticleCreatorDestructor::FindMaxElementIdInModelPart(ModelPart& r_modelpart) {

        int max_Id = 1;
        
        for (ModelPart::ElementsContainerType::iterator element_it = r_modelpart.GetCommunicator().LocalMesh().ElementsBegin();
            element_it != r_modelpart.GetCommunicator().LocalMesh().ElementsEnd();
            element_it++) {

            if ((int) (element_it->Id()) > max_Id) max_Id = element_it->Id();
        }

        r_modelpart.GetCommunicator().MaxAll(max_Id);
        return max_Id;
    }
    
    void ParticleCreatorDestructor::NodeCreatorWithPhysicalParameters(ModelPart& r_modelpart,
                                                                    Node < 3 > ::Pointer& pnew_node,
                                                                    int aId,
                                                                    Node < 3 > ::Pointer& reference_node,
                                                                    double radius,
                                                                    Properties& params,
                                                                    bool has_sphericity,
                                                                    bool has_rotation,
                                                                    bool initial) {
        
        array_1d<double, 3 > null_vector(3, 0.0);

        double bx = reference_node->X();
        double cy = reference_node->Y();
        double dz = reference_node->Z();

        if (initial) {
            pnew_node = reference_node;
            pnew_node->SetId(aId);
            r_modelpart.AddNode(pnew_node); // The same node is added to r_modelpart (the calculation model part)          
            pnew_node->FastGetSolutionStepValue(VELOCITY) = null_vector;            
            //(actually it should be the velocity of the inlet layer, which is different from the particles being inserted)
            pnew_node->FastGetSolutionStepValue(PARTICLE_MATERIAL) = params[PARTICLE_MATERIAL] + 100; //So the inlet ghost spheres are not in the same
        }                                                                                             //layer of the inlet newly created spheres
        else {
            pnew_node = r_modelpart.CreateNewNode(aId, bx, cy, dz); //ACTUAL node creation and addition to model part
            array_1d<double, 3 > velocity = params[VELOCITY];
            double max_rand_deviation_angle = params[MAX_RAND_DEVIATION_ANGLE];
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

        //The next section is commented because initialization to 0 is done automatically of all variables in GetSolutionStepvariable
        /*if (pnew_node->SolutionStepsDataHas(SOLID_FRACTION_PROJECTED)) {
            pnew_node->FastGetSolutionStepValue(SOLID_FRACTION_PROJECTED) = 0.0;
        }
        
        if (pnew_node->SolutionStepsDataHas(MESH_VELOCITY1)) {
            pnew_node->SolutionStepsDataHas(MESH_VELOCITY1) = null_vector
            //pnew_node->FastGetSolutionStepValue(MESH_VELOCITY1_X) = 0.0;
            //pnew_node->FastGetSolutionStepValue(MESH_VELOCITY1_Y) = 0.0;
            //pnew_node->FastGetSolutionStepValue(MESH_VELOCITY1_Z) = 0.0;
        }

        if (pnew_node->SolutionStepsDataHas(DRAG_REACTION)) {
            pnew_node->SolutionStepsDataHas(DRAG_REACTION) = null_vector
            //pnew_node->FastGetSolutionStepValue(DRAG_REACTION_X) = 0.0;
            //pnew_node->FastGetSolutionStepValue(DRAG_REACTION_Y) = 0.0;
            //pnew_node->FastGetSolutionStepValue(DRAG_REACTION_Z) = 0.0;
        }*/

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
    }

    void ParticleCreatorDestructor::ElementCreatorWithPhysicalParameters(ModelPart& r_modelpart,
                                                                        int r_Elem_Id,
                                                                        Node < 3 > ::Pointer reference_node,
                                                                        Element::Pointer injector_element,
                                                                        Properties::Pointer r_params,
                                                                        const Element& r_reference_element,
                                                                        PropertiesProxy* p_fast_properties,
                                                                        bool has_sphericity,
                                                                        bool has_rotation,
                                                                        bool initial,
                                                                        ElementsContainerType& array_of_injector_elements) {

        Node < 3 > ::Pointer pnew_node;

        double radius = (*r_params)[RADIUS];
        double max_radius = 1.5 * radius;

        if (initial) {
            radius = max_radius;
        } else {
            double std_deviation = (*r_params)[STANDARD_DEVIATION];
            std::string distribution_type = (*r_params)[PROBABILITY_DISTRIBUTION];
            double min_radius = 0.5 * radius;
            
            if (distribution_type == "normal") radius = rand_normal(radius, std_deviation, max_radius, min_radius);
            else if (distribution_type == "lognormal") radius = rand_lognormal(radius, std_deviation, max_radius, min_radius);
            else KRATOS_THROW_ERROR(std::runtime_error, "Unknown probability distribution.", "")
        }

        NodeCreatorWithPhysicalParameters(r_modelpart, pnew_node, r_Elem_Id, reference_node, radius, *r_params, has_sphericity, has_rotation, initial);

        Geometry< Node < 3 > >::PointsArrayType nodelist;
        
        nodelist.push_back(pnew_node);
        Element::Pointer p_particle = r_reference_element.Create(r_Elem_Id, nodelist, r_params);
        Kratos::SphericParticle* spheric_p_particle = dynamic_cast<Kratos::SphericParticle*> (p_particle.get());

        if (initial) {
            array_of_injector_elements.push_back(p_particle);
            p_particle->Set(BLOCKED);
            pnew_node->Set(BLOCKED);
        }
        else {
            array_1d<double, 3 > zero_vector(3, 0.0);
            Kratos::SphericParticle* injector_spheric_particle = dynamic_cast<Kratos::SphericParticle*> (injector_element.get());
            injector_spheric_particle->mNeighbourElements.push_back(spheric_p_particle);
            injector_spheric_particle->mNeighbourElasticContactForces.push_back(zero_vector);
            injector_spheric_particle->mNeighbourTotalContactForces.push_back(zero_vector);
            spheric_p_particle->mNeighbourElements.push_back(injector_spheric_particle);
            spheric_p_particle->mNeighbourElasticContactForces.push_back(zero_vector);
            spheric_p_particle->mNeighbourTotalContactForces.push_back(zero_vector);
        }

        p_particle->Set(NEW_ENTITY);
        pnew_node->Set(NEW_ENTITY);
        //p_particle->InitializeSolutionStep(r_modelpart.GetProcessInfo());       

        spheric_p_particle->SetFastProperties(p_fast_properties);

        double density = spheric_p_particle->GetDensity();
        spheric_p_particle->SetRadius(radius);
        double mass = 4.0 / 3.0 * KRATOS_M_PI * density * radius * radius * radius;
        spheric_p_particle->SetMass(mass);

        if (has_rotation) spheric_p_particle->Set(DEMFlags::HAS_ROTATION, true);
        else spheric_p_particle->Set(DEMFlags::HAS_ROTATION, false);

        spheric_p_particle->FullInitialize(r_modelpart.GetProcessInfo()); 

        r_modelpart.Elements().push_back(p_particle);
    }

    void ParticleCreatorDestructor::NodeCreatorForClusters(ModelPart& r_modelpart,
        Node < 3 > ::Pointer& pnew_node,
        int aId,
        array_1d<double, 3>& reference_coordinates,
        double radius,
        Properties& params) {
        double bx = reference_coordinates[0];
        double cy = reference_coordinates[1];
        double dz = reference_coordinates[2];

        #pragma omp critical
        {
        pnew_node = r_modelpart.CreateNewNode(aId, bx, cy, dz); //ACTUAL node creation and addition to model part
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
    }

    Kratos::SphericParticle* ParticleCreatorDestructor::ElementCreatorForClusters(ModelPart& r_modelpart,
                                                                                  int r_Elem_Id,
                                                                                  double radius,
                                                                                  array_1d<double, 3>& reference_coordinates,
                                                                                  double cluster_mass,
                                                                                  Properties::Pointer r_params,
                                                                                  const Element& r_reference_element,
                                                                                  const int cluster_id) {
        Node<3>::Pointer pnew_node;

        NodeCreatorForClusters(r_modelpart, pnew_node, r_Elem_Id, reference_coordinates, radius, *r_params);

        Geometry<Node <3> >::PointsArrayType nodelist;
        nodelist.push_back(pnew_node);

        Element::Pointer p_particle = r_reference_element.Create(r_Elem_Id, nodelist, r_params);
        //p_particle->InitializeSolutionStep(r_modelpart.GetProcessInfo());

        Kratos::SphericParticle* spheric_p_particle = dynamic_cast<Kratos::SphericParticle*> (p_particle.get());

        spheric_p_particle->SetRadius(radius);
        spheric_p_particle->SetMass(cluster_mass);
        spheric_p_particle->MemberDeclarationFirstStep(r_modelpart.GetProcessInfo());
        spheric_p_particle->CreateDiscontinuumConstitutiveLaws(r_modelpart.GetProcessInfo());

        spheric_p_particle->Set(DEMFlags::HAS_ROLLING_FRICTION, false);
        spheric_p_particle->Set(DEMFlags::BELONGS_TO_A_CLUSTER, true);
        spheric_p_particle->SetClusterId(cluster_id);

        #pragma omp critical
        {        
            r_modelpart.Elements().push_back(p_particle);
        }
        
        return spheric_p_particle;
    }    
    
    void ParticleCreatorDestructor::ClusterCreatorWithPhysicalParameters(ModelPart& r_modelpart,
                                                                         ModelPart& r_clusters_modelpart,
                                                                         int r_Elem_Id,
                                                                         Node < 3 > ::Pointer reference_node,
                                                                         Element::Pointer injector_element,
                                                                         Properties::Pointer r_params,
                                                                         const Element& r_reference_element,
                                                                         PropertiesProxy* p_fast_properties,
                                                                         bool has_sphericity,
                                                                         bool has_rotation,
                                                                         bool initial,
                                                                         ElementsContainerType& array_of_injector_elements,
                                                                         int& number_of_added_spheres) {
        
        Node < 3 > ::Pointer pnew_node;

        double radius = (*r_params)[RADIUS];
        double max_radius = 1.5 * radius;
        array_1d<double, 3> euler_angles;
        
        if (initial) {
            radius = max_radius;
        } else {
            double std_deviation = (*r_params)[STANDARD_DEVIATION];
            double min_radius = 0.5 * radius;
            radius = rand_normal(radius, std_deviation, max_radius, min_radius);
        }

        NodeCreatorWithPhysicalParameters(r_clusters_modelpart, pnew_node, r_Elem_Id, reference_node, radius, *r_params, has_sphericity, has_rotation, initial);
        
        pnew_node->FastGetSolutionStepValue(CHARACTERISTIC_LENGTH) = radius * 2.0; //Cluster specific. Can be removed
        
        Geometry< Node < 3 > >::PointsArrayType nodelist;

        nodelist.push_back(pnew_node);

        Element::Pointer p_particle = r_reference_element.Create(r_Elem_Id, nodelist, r_params);

        if (initial) {
            array_of_injector_elements.push_back(p_particle);
            p_particle->Set(BLOCKED);
            pnew_node->Set(BLOCKED);
        }
        else {
            Kratos::Cluster3D* p_cluster = dynamic_cast<Kratos::Cluster3D*> (p_particle.get());
            p_cluster->Initialize();                        
                    
            if ((*r_params)[RANDOM_EULER_ANGLES]) {
                
                double random_factor = 2.0 * KRATOS_M_PI / RAND_MAX;
                euler_angles[0] = random_factor * rand();
                euler_angles[1] = random_factor * rand();
                euler_angles[2] = random_factor * rand();                
            }
            
            else {
                euler_angles[0] = (*r_params)[EULER_ANGLES][0];
                euler_angles[1] = (*r_params)[EULER_ANGLES][1];
                euler_angles[2] = (*r_params)[EULER_ANGLES][2];
            }
            
            p_cluster->SetOrientation(euler_angles);            
            
            ParticleCreatorDestructor* creator_destructor_ptr = this;
            p_cluster->CreateParticles(creator_destructor_ptr, r_modelpart);
            number_of_added_spheres = p_cluster->GetNumberOfSpheres();
            
            Kratos::SphericParticle* injector_spheric_particle = dynamic_cast<Kratos::SphericParticle*> (injector_element.get());
            
            array_1d<double, 3 > zero_vector(3, 0.0);
            
            for (unsigned int i=0; i<p_cluster->GetNumberOfSpheres(); i++) { 
                Kratos::SphericParticle* spheric_p_particle = p_cluster->GetSpheres()[i];
                injector_spheric_particle->mNeighbourElements.push_back(spheric_p_particle);
                injector_spheric_particle->mNeighbourElasticContactForces.push_back(zero_vector);
                injector_spheric_particle->mNeighbourTotalContactForces.push_back(zero_vector);
                spheric_p_particle->mNeighbourElements.push_back(injector_spheric_particle);
                spheric_p_particle->mNeighbourElasticContactForces.push_back(zero_vector);
                spheric_p_particle->mNeighbourTotalContactForces.push_back(zero_vector);
                spheric_p_particle->Set(NEW_ENTITY);
                spheric_p_particle->GetGeometry()[0].Set(NEW_ENTITY);
                //spheric_p_particle->InitializeSolutionStep(r_modelpart.GetProcessInfo());  //protected!!
                spheric_p_particle->SetFastProperties(p_fast_properties);
            }
        }

        p_particle->Set(NEW_ENTITY);
        pnew_node->Set(NEW_ENTITY);            
        
        if (has_rotation) p_particle->Set(DEMFlags::HAS_ROTATION, true);
        else p_particle->Set(DEMFlags::HAS_ROTATION, false);

        r_clusters_modelpart.Elements().push_back(p_particle);                
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
                array_1d<double, 3 > coor = (*(Elements.begin().base()))->GetGeometry()[0].Coordinates();
                mStrictLowPoint = coor;
                mStrictHighPoint = coor;

                for (Configure::ElementsContainerType::iterator particle_pointer_it = Elements.begin(); particle_pointer_it != Elements.end(); ++particle_pointer_it) {
                    coor = (*(particle_pointer_it.base()))->GetGeometry()[0].Coordinates();
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
            
             
            array_1d<double, 3 > midpoint = 0.5 * (mStrictHighPoint + mStrictLowPoint);
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

                //Type definitions
        typedef ModelPart::ElementsContainerType ElementsArrayType;
//        typedef ElementsArrayType::iterator ElementsIterator;

        ElementsArrayType& rElements = r_model_part.GetCommunicator().LocalMesh().Elements();
        ModelPart::NodesContainerType& rNodes = r_model_part.GetCommunicator().LocalMesh().Nodes();

        //ElementsArrayType& rElements                          = r_model_part.Elements();
        //ModelPart::NodesContainerType& rNodes               = r_model_part.Nodes();

        ElementsArrayType temp_particles_container;
        ModelPart::NodesContainerType temp_nodes_container;

        temp_particles_container.reserve(rElements.size());
        temp_nodes_container.reserve(rNodes.size());

        temp_particles_container.swap(rElements);
        temp_nodes_container.swap(rNodes);

        //Add the ones not marked with TO_ERASE
        for (Configure::ElementsContainerType::ptr_iterator particle_pointer_it = temp_particles_container.ptr_begin(); particle_pointer_it != temp_particles_container.ptr_end(); ++particle_pointer_it) {

            if (!(*particle_pointer_it)->GetGeometry()[0].Is(TO_ERASE) && !(*particle_pointer_it)->Is(TO_ERASE)) {
                (rElements).push_back(*particle_pointer_it); //adding the elements

                for (unsigned int i = 0; i < (*particle_pointer_it)->GetGeometry().PointsNumber(); i++) { //GENERAL FOR ELEMENTS OF MORE THAN ONE NODE
                    ModelPart::NodeType::Pointer pNode = (*particle_pointer_it)->GetGeometry().pGetPoint(i);
                    (rNodes).push_back(pNode);
                }
            }

        }
        KRATOS_CATCH("")
    }

    void ParticleCreatorDestructor::DestroyContactElements(ModelPart& r_model_part) {
        KRATOS_TRY

                //Type definitions
                typedef ModelPart::ElementsContainerType ElementsArrayType;
//        typedef ElementsArrayType::iterator ElementsIterator;

        ElementsArrayType& rElements = r_model_part.GetCommunicator().LocalMesh().Elements();

        ElementsArrayType temp_elements_container;

        //Copy the elements and clear the element container
        //temp_elements_container.reserve(pElements->size());
        temp_elements_container.reserve(rElements.size());
        temp_elements_container.swap(rElements);

        //Add the ones not marked with TO_ERASE
        for (Configure::ElementsContainerType::ptr_iterator element_pointer_it = temp_elements_container.ptr_begin(); element_pointer_it != temp_elements_container.ptr_end(); ++element_pointer_it) {

            if (!(*element_pointer_it)->Is(TO_ERASE)) {
                (rElements).push_back(*element_pointer_it); //adding the elements               
            }
        }
        KRATOS_CATCH("")
    }
    /*
    void ParticleCreatorDestructor::MarkInitialNeighboursThatAreBeingRemoved(ModelPart& r_model_part) { //TODO To be removed

        KRATOS_TRY

                typedef ModelPart::ElementsContainerType ElementsArrayType;
        //        typedef ElementsArrayType::iterator ElementsIterator;

        ElementsArrayType& rElements = r_model_part.GetCommunicator().LocalMesh().Elements();

        for (Configure::ElementsContainerType::ptr_iterator particle_pointer_it = rElements.ptr_begin(); particle_pointer_it != rElements.ptr_end(); ++particle_pointer_it) {

            if (!(*particle_pointer_it)->Is(TO_ERASE)) continue;

            Element* p_element = particle_pointer_it->get();
            SphericContinuumParticle* p_continuum_spheric_particle = dynamic_cast<SphericContinuumParticle*> (p_element);

            for (unsigned int i = 0; i < p_continuum_spheric_particle->mContinuumIniNeighbourElements.size(); i++) {
                SphericContinuumParticle* neighbour_i = p_continuum_spheric_particle->mContinuumIniNeighbourElements[i];
                if (neighbour_i == NULL) continue;
                for (unsigned int j = 0; j < neighbour_i->mContinuumIniNeighbourElements.size(); j++) {
                    if (neighbour_i->mContinuumIniNeighbourElements[j] == p_continuum_spheric_particle) {
                        neighbour_i->mContinuumIniNeighbourElements[j] = NULL;
                        break;
                    }
                }
            }

        }

        KRATOS_CATCH("")
    }
    */
    
    void ParticleCreatorDestructor::MarkDistantParticlesForErasing(ModelPart& r_model_part) {
        MarkParticlesForErasingGivenBoundingBox(r_model_part, mLowPoint, mHighPoint);
    }

    void ParticleCreatorDestructor::MarkParticlesForErasingGivenScalarVariableValue(ModelPart& r_model_part, const Variable<double>& rVariable, double value, double tol) {

        KRATOS_TRY

        Configure::ElementsContainerType& rElements = r_model_part.GetCommunicator().LocalMesh().Elements();

        for (Configure::ElementsContainerType::ptr_iterator particle_pointer_it = rElements.ptr_begin();
                particle_pointer_it != rElements.ptr_end(); ++particle_pointer_it) {

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

        for (Configure::ElementsContainerType::ptr_iterator particle_pointer_it = rElements.ptr_begin();
                particle_pointer_it != rElements.ptr_end(); ++particle_pointer_it) {

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

        for (Configure::ElementsContainerType::ptr_iterator particle_pointer_it = rElements.ptr_begin();
                particle_pointer_it != rElements.ptr_end(); ++particle_pointer_it) {

            if ((*particle_pointer_it)->Is(DEMFlags::BELONGS_TO_A_CLUSTER)) continue;

            array_1d<double, 3 > coor = (*particle_pointer_it)->GetGeometry()[0].Coordinates();
            bool include = true;

            for (unsigned int i = 0; i < 3; i++) {
                include = include && (coor[i] >= low_point[i]) && (coor[i] <= high_point[i]);
            }

            if (!include) {
                (*particle_pointer_it)->GetGeometry()[0].Set(TO_ERASE);
                (*particle_pointer_it)->Set(TO_ERASE);
            }

        }

        for (ModelPart::NodesContainerType::ptr_iterator node_pointer_it = rNodes.ptr_begin();
                node_pointer_it != rNodes.ptr_end(); ++node_pointer_it) {

            if ((*node_pointer_it)->Is(DEMFlags::BELONGS_TO_A_CLUSTER)) continue;
            array_1d<double, 3 > coor = (*node_pointer_it)->Coordinates();
            bool include = true;

            for (unsigned int i = 0; i < 3; i++) {
                include = include && (coor[i] >= low_point[i]) && (coor[i] <= high_point[i]);
            }

            if (!include) {
                (*node_pointer_it)->Set(TO_ERASE);
            }

        }

        KRATOS_CATCH("")
    }

    void ParticleCreatorDestructor::MarkContactElementsForErasing(ModelPart& r_model_part, ModelPart& mcontacts_model_part) {
        ///ONLY FOR CONTINUUM!!!

        KRATOS_TRY

        Configure::ElementsContainerType& rElements = r_model_part.GetCommunicator().LocalMesh().Elements();

        for (Configure::ElementsContainerType::ptr_iterator particle_pointer_it = rElements.ptr_begin();
                particle_pointer_it != rElements.ptr_end(); ++particle_pointer_it) {

            if ((*particle_pointer_it)->GetGeometry()[0].Is(TO_ERASE)) {

                Element* p_element = particle_pointer_it->get();
                SphericContinuumParticle* p_continuum_spheric_particle = dynamic_cast<SphericContinuumParticle*> (p_element);
                std::vector<Particle_Contact_Element*>& node_to_neigh_element_pointer = p_continuum_spheric_particle->mBondElements;
                unsigned int number_of_contact_elements = node_to_neigh_element_pointer.size();
                for (unsigned int i = 0; i < number_of_contact_elements; i++) {
                    Particle_Contact_Element* p_to_contact_element = node_to_neigh_element_pointer[i];
                    if (p_to_contact_element != NULL) { //NULL happens when the initial neighbor was a ghost and had a lower Id than the others
                        p_to_contact_element->Set(TO_ERASE);
                    }
                }

            }

        }

        KRATOS_CATCH("")
    }

    void ParticleCreatorDestructor::DestroyParticlesOutsideBoundingBox(ModelPart& r_model_part) {
        MarkDistantParticlesForErasing(r_model_part);
        DestroyParticles(r_model_part);
    }

    void ParticleCreatorDestructor::DestroyContinuumParticlesOutsideBoundingBox(ModelPart& r_model_part) {
        MarkDistantParticlesForErasing(r_model_part);
        //MarkInitialNeighboursThatAreBeingRemoved(r_model_part);
        DestroyParticles(r_model_part);
    }

    void ParticleCreatorDestructor::DestroyContactElementsOutsideBoundingBox(ModelPart& r_model_part, ModelPart& mcontacts_model_part) {
        MarkContactElementsForErasing(r_model_part, mcontacts_model_part);
        DestroyContactElements(mcontacts_model_part);
    }

    array_1d<double, 3 > ParticleCreatorDestructor::GetHighNode() {
        return (mHighPoint);
    }

    array_1d<double, 3 > ParticleCreatorDestructor::GetLowNode() {
        return (mLowPoint);
    }

    array_1d<double, 3 > ParticleCreatorDestructor::GetStrictHighNode() {
        return (mStrictHighPoint);
    }

    array_1d<double, 3 > ParticleCreatorDestructor::GetStrictLowNode() {
        return (mStrictLowPoint);
    }

    double ParticleCreatorDestructor::GetDiameter() {
        return (mDiameter);
    }

    double ParticleCreatorDestructor::GetStrictDiameter() {
        return (mStrictDiameter);
    }

    void ParticleCreatorDestructor::SetHighNode(array_1d<double, 3 > node) {
        mHighPoint = node;
    }

    void ParticleCreatorDestructor::SetLowNode(array_1d<double, 3 > node) {
        mLowPoint = node;
    }

    unsigned int ParticleCreatorDestructor::GetCurrentMaxNodeId() {
        return mMaxNodeId;
    }
    
    unsigned int* ParticleCreatorDestructor::pGetCurrentMaxNodeId() {
        return &mMaxNodeId;
    }

    void ParticleCreatorDestructor::SetMaxNodeId(unsigned int id) {
        mMaxNodeId = id;
    }

    /// Turn back information as a stemplate<class T, std::size_t dim> tring.

    std::string ParticleCreatorDestructor::Info() const {
        return "";
    }

    /// Print information about this object

    void ParticleCreatorDestructor::PrintInfo(std::ostream& rOStream) const {
    }

    /// Print object's data

    void ParticleCreatorDestructor::PrintData(std::ostream& rOStream) const {
    }

    void ParticleCreatorDestructor::Clear(ModelPart::NodesContainerType::iterator node_it, int step_data_size) {
        unsigned int buffer_size = node_it->GetBufferSize();
        for (unsigned int step = 0; step < buffer_size; step++) {
            //getting the data of the solution step
            double* step_data = (node_it)->SolutionStepData().Data(step);
            //copying this data in the position of the vector we are interested in
            for (int j = 0; j < step_data_size; j++) {
                step_data[j] = 0.0;
            }
        }
    }

    void ParticleCreatorDestructor::ClearVariables(ModelPart::NodesContainerType::iterator node_it, Variable<array_1d<double, 3 > >& rVariable) {
        array_1d<double, 3 > & Aux_var = node_it->FastGetSolutionStepValue(rVariable);
        noalias(Aux_var) = ZeroVector(3);
    }

    void ParticleCreatorDestructor::ClearVariables(ParticleIterator particle_it, Variable<double>& rVariable) {
    }

} //Namespace Kratos
