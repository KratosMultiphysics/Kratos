//
// Author: 
// Miguel Angel Celigueta maceli@cimne.upc.edu
//

#include <string>
#include <iostream>

#include "inlet.h"
#include "create_and_destroy.h"
#include "custom_elements/spheric_particle.h"
#include "custom_elements/spheric_continuum_particle.h"
#include "custom_elements/cluster3D.h"
#include "custom_constitutive/DEM_discontinuum_constitutive_law.h"
#include "custom_constitutive/DEM_continuum_constitutive_law.h"
#include "dem_fem_utilities.h"
#include "GeometryFunctions.h"


namespace Kratos {
    
    inline double CalculateNormalizedIndentation(SphericParticle& elem_it_1, SphericParticle& elem_it_2) {
        const array_1d<double,3>& coordinates_1 = elem_it_1.GetGeometry()[0].Coordinates();
        const array_1d<double,3>& coordinates_2 = elem_it_2.GetGeometry()[0].Coordinates();
        
        const double distance = sqrt((coordinates_1[0]- coordinates_2[0]) * (coordinates_1[0] - coordinates_2[0]) +
                                     (coordinates_1[1]- coordinates_2[1]) * (coordinates_1[1] - coordinates_2[1]) +
                                     (coordinates_1[2]- coordinates_2[2]) * (coordinates_1[2] - coordinates_2[2]));

        const double radius_sum = elem_it_1.GetInteractionRadius() + elem_it_2.GetInteractionRadius();
        double indentation = radius_sum - distance;

        indentation /= radius_sum;

        return indentation;
    }

    /// Constructor            
    
    DEM_Inlet::DEM_Inlet(ModelPart& inlet_modelpart): mInletModelPart(inlet_modelpart)
    {                
        mPartialParticleToInsert.resize(inlet_modelpart.NumberOfSubModelParts(), false);
        mLastInjectionTimes.resize(inlet_modelpart.NumberOfSubModelParts(), false);
        mTotalNumberOfDetachedParticles.resize(inlet_modelpart.NumberOfSubModelParts(), false);
        mLayerRemoved.resize(inlet_modelpart.NumberOfSubModelParts(), false);
        
        int smp_iterator_number = 0;   
        for (ModelPart::SubModelPartsContainerType::iterator sub_model_part = inlet_modelpart.SubModelPartsBegin(); sub_model_part != inlet_modelpart.SubModelPartsEnd(); ++sub_model_part) {                        
            mPartialParticleToInsert[smp_iterator_number] = 0.0;
            mLastInjectionTimes[smp_iterator_number] = 0.0;
            mLayerRemoved[smp_iterator_number] = false;
            smp_iterator_number++;
        }                
        
        mFirstTime = true;
        mBallsModelPartHasSphericity = false;  
        mBallsModelPartHasRotation   = false;
    }
            
    /// Destructor
    
    DEM_Inlet::~DEM_Inlet() {}
     
    void DEM_Inlet::InitializeDEM_Inlet(ModelPart& r_modelpart, ParticleCreatorDestructor& creator, const bool using_strategy_for_continuum) {
        
        mStrategyForContinuum = using_strategy_for_continuum;
        
        unsigned int& max_Id=creator.mMaxNodeId;       
        CreatePropertiesProxies(mFastProperties, mInletModelPart);      
        VariablesList r_modelpart_nodal_variables_list = r_modelpart.GetNodalSolutionStepVariablesList();
        
        if (r_modelpart_nodal_variables_list.Has(PARTICLE_SPHERICITY)) mBallsModelPartHasSphericity = true;        
        
        if (r_modelpart.GetProcessInfo()[ROTATION_OPTION]) {
            mBallsModelPartHasRotation = true;
            mInletModelPart.GetProcessInfo()[ROTATION_OPTION] = true;            
        }
        else {
            mInletModelPart.GetProcessInfo()[ROTATION_OPTION] = false;             
        }
               
        int smp_number = 0;
                
        for (ModelPart::SubModelPartsContainerType::iterator smp_it = mInletModelPart.SubModelPartsBegin(); smp_it != mInletModelPart.SubModelPartsEnd(); ++smp_it) {
            ModelPart& mp = *smp_it;

            int mesh_size = smp_it->NumberOfNodes();
            if (!mesh_size) continue;
            ModelPart::NodesContainerType::ContainerType all_nodes = smp_it->NodesArray();
            std::string& identifier = mp[IDENTIFIER];
            mp[INLET_INITIAL_VELOCITY] = mp[LINEAR_VELOCITY]; //This is the velocity of the moving injector particles
            
            array_1d<double, 3>& inlet_velocity = mp[VELOCITY]; //This is the velocity of the injected particles
            
            if ((inlet_velocity[0] == 0.0) &&
                (inlet_velocity[1] == 0.0) &&
                (inlet_velocity[2] == 0.0)) {
                
                KRATOS_THROW_ERROR(std::runtime_error, "The inlet velocity cannot be zero for group ", identifier);
            }
            
            double max_rand_dev_angle = mp[MAX_RAND_DEVIATION_ANGLE];
            if (max_rand_dev_angle < 0.0 || max_rand_dev_angle > 89.5) {
                
                KRATOS_THROW_ERROR(std::runtime_error, "The velocity deviation angle must be between 0 and 90 degrees for group ", identifier);
            }
            
            int general_properties_id = mInletModelPart.GetProperties(mp[PROPERTIES_ID]).Id();  
            PropertiesProxy* p_fast_properties = NULL;
            
            for (unsigned int i = 0; i < mFastProperties.size(); i++) {
                int fast_properties_id = mFastProperties[i].GetId(); 
                if (fast_properties_id == general_properties_id) {  
                    p_fast_properties = &(mFastProperties[i]);
                    break;
                }
                mLastInjectionTimes[smp_number] = mp[INLET_START_TIME];
            }
            
            Element::Pointer dummy_element_pointer;
            std::string ElementNameString;
            if (using_strategy_for_continuum) ElementNameString = "SphericContinuumParticle3D";
            else ElementNameString = "SphericParticle3D";
            const Element& r_reference_element = KratosComponents<Element>::Get(ElementNameString);
            
            Properties::Pointer p_properties = mInletModelPart.pGetProperties(mp[PROPERTIES_ID]);
            
            for (int i = 0; i < mesh_size; i++) {                
                Element* p_element = creator.ElementCreatorWithPhysicalParameters(r_modelpart,
                                                             max_Id+1,
                                                             all_nodes[i],
                                                             dummy_element_pointer,
                                                             p_properties,
                                                             mp,
                                                             r_reference_element,
                                                             p_fast_properties,
                                                             mBallsModelPartHasSphericity,
                                                             mBallsModelPartHasRotation,
                                                             true,
                                                             smp_it->Elements());
		max_Id++;
                if(using_strategy_for_continuum){
                    SphericContinuumParticle* p_continuum_spheric_particle = dynamic_cast<SphericContinuumParticle*>(p_element);
                    p_continuum_spheric_particle->mContinuumInitialNeighborsSize=0;
                    p_continuum_spheric_particle->mInitialNeighborsSize=0;    
                }
            } 
            smp_number++;
        } //for smp_it                                               
    } //InitializeDEM_Inlet

    void DEM_Inlet::DettachElements(ModelPart& r_modelpart, unsigned int& max_Id) {

        vector<unsigned int> ElementPartition;
        OpenMPUtils::CreatePartition(OpenMPUtils::GetNumThreads(), r_modelpart.GetCommunicator().LocalMesh().Elements().size(), ElementPartition);
        typedef ElementsArrayType::iterator ElementIterator;

        #pragma omp parallel for
        for (int k = 0; k < (int)r_modelpart.GetCommunicator().LocalMesh().Elements().size(); k++) {
            ElementIterator elem_it = r_modelpart.GetCommunicator().LocalMesh().Elements().ptr_begin() + k;                                        
            if (elem_it->IsNot(NEW_ENTITY)) continue;

            SphericParticle& spheric_particle = dynamic_cast<SphericParticle&>(*elem_it);
            Node<3>& node_it = elem_it->GetGeometry()[0];

            bool have_just_stopped_touching = true;

            for (unsigned int i = 0; i < spheric_particle.mNeighbourElements.size(); i++) {
                SphericParticle* neighbour_particle = spheric_particle.mNeighbourElements[i];
                if(neighbour_particle == NULL) continue;
                
                Node<3>& neighbour_node = neighbour_particle->GetGeometry()[0];

                const double admissible_indentation_ratio = 0.0;
                const double indentation = CalculateNormalizedIndentation(spheric_particle, *neighbour_particle);
                const bool indentation_is_significant = indentation > admissible_indentation_ratio;
                const bool i_am_injected_he_is_injector = node_it.IsNot(BLOCKED) && neighbour_node.Is(BLOCKED);
                const bool i_am_injector_he_is_injected = node_it.Is(BLOCKED) && neighbour_node.IsNot(BLOCKED);

                if ((i_am_injected_he_is_injector || i_am_injector_he_is_injected) && indentation_is_significant) {
                    have_just_stopped_touching = false;
                    break;
                }
            }

            if (have_just_stopped_touching) {
                if (node_it.IsNot(BLOCKED)) {//The ball must be freed
                    node_it.Set(DEMFlags::FIXED_VEL_X, false);
                    node_it.Set(DEMFlags::FIXED_VEL_Y, false);
                    node_it.Set(DEMFlags::FIXED_VEL_Z, false);
                    node_it.Set(DEMFlags::FIXED_ANG_VEL_X, false);
                    node_it.Set(DEMFlags::FIXED_ANG_VEL_Y, false);
                    node_it.Set(DEMFlags::FIXED_ANG_VEL_Z, false);
                    elem_it->Set(NEW_ENTITY, 0);
                    node_it.Set(NEW_ENTITY, 0);
                    node_it.pGetDof(VELOCITY_X)->FreeDof();
                    node_it.pGetDof(VELOCITY_Y)->FreeDof();
                    node_it.pGetDof(VELOCITY_Z)->FreeDof();
                    node_it.pGetDof(ANGULAR_VELOCITY_X)->FreeDof();
                    node_it.pGetDof(ANGULAR_VELOCITY_Y)->FreeDof();
                    node_it.pGetDof(ANGULAR_VELOCITY_Z)->FreeDof();
                }
                else {
                    //Inlet BLOCKED nodes are ACTIVE when injecting, so when they cease to be in contact with other balls, ACTIVE is set to 'false', as they become available for injecting new elements.
                    node_it.Set(ACTIVE, false);
                    elem_it->Set(ACTIVE, false);
                }
            }            
        }         
        mFirstTime = false;        
    } //Dettach


    void DEM_Inlet::DettachClusters(ModelPart& r_clusters_modelpart, unsigned int& max_Id) {    
                            
        vector<unsigned int> ElementPartition;
        typedef ElementsArrayType::iterator ElementIterator;       

        #pragma omp parallel for
        for (int k = 0; k < (int)r_clusters_modelpart.GetCommunicator().LocalMesh().Elements().size(); k++) {
            ElementIterator elem_it = r_clusters_modelpart.GetCommunicator().LocalMesh().Elements().ptr_begin() + k;                            
            if (elem_it->IsNot(NEW_ENTITY)) continue;

            Kratos::Cluster3D& r_cluster = dynamic_cast<Kratos::Cluster3D&>(*elem_it); 
            Node<3>& cluster_central_node = r_cluster.GetGeometry()[0];
            bool still_touching=false;

            for (unsigned int j = 0; j < r_cluster.GetSpheres().size(); j++) { //loop over the spheres of the cluster
                SphericParticle* spheric_particle = r_cluster.GetSpheres()[j];
                for (unsigned int i = 0; i < spheric_particle->mNeighbourElements.size(); i++) { //loop over the neighbor spheres of each sphere of the cluster
                    SphericParticle* neighbour_iterator = spheric_particle->mNeighbourElements[i];
                    Node<3>& neighbour_node = neighbour_iterator->GetGeometry()[0]; 
                    if (neighbour_node.Is(BLOCKED)) {
                        still_touching = true;
                        break;
                    }              
                }
                if (still_touching) break;
            }

            if (!still_touching) { //The ball must be freed
                cluster_central_node.Set(DEMFlags::FIXED_VEL_X, false);
                cluster_central_node.Set(DEMFlags::FIXED_VEL_Y, false);
                cluster_central_node.Set(DEMFlags::FIXED_VEL_Z, false);
                cluster_central_node.Set(DEMFlags::FIXED_ANG_VEL_X, false);
                cluster_central_node.Set(DEMFlags::FIXED_ANG_VEL_Y, false);
                cluster_central_node.Set(DEMFlags::FIXED_ANG_VEL_Z, false);
                elem_it->Set(NEW_ENTITY, 0);
                cluster_central_node.Set(NEW_ENTITY, 0);
                cluster_central_node.pGetDof(VELOCITY_X)->FreeDof();
                cluster_central_node.pGetDof(VELOCITY_Y)->FreeDof();
                cluster_central_node.pGetDof(VELOCITY_Z)->FreeDof();
                cluster_central_node.pGetDof(ANGULAR_VELOCITY_X)->FreeDof();
                cluster_central_node.pGetDof(ANGULAR_VELOCITY_Y)->FreeDof();
                cluster_central_node.pGetDof(ANGULAR_VELOCITY_Z)->FreeDof();

                for (unsigned int j = 0; j < r_cluster.GetSpheres().size(); j++) { //loop over the spheres of the cluster
                    SphericParticle* spheric_particle = r_cluster.GetSpheres()[j];
                    Node<3>& node_it = spheric_particle->GetGeometry()[0];  
                    spheric_particle->Set(NEW_ENTITY, 0);
                    node_it.Set(NEW_ENTITY, 0);                    
                }                    
            }                
        } 
        mFirstTime=false;
    } //DettachClusters
    
    void DEM_Inlet::CreateElementsFromInletMesh(ModelPart& r_modelpart, ModelPart& r_clusters_modelpart, ParticleCreatorDestructor& creator) {                    
        unsigned int& max_Id=creator.mMaxNodeId; 
        const double current_time = r_modelpart.GetProcessInfo()[TIME];
        DettachElements(r_modelpart, max_Id);
        DettachClusters(r_clusters_modelpart, max_Id);
                
        int smp_number = 0;
        for (ModelPart::SubModelPartsContainerType::iterator smp_it = mInletModelPart.SubModelPartsBegin(); smp_it != mInletModelPart.SubModelPartsEnd(); ++smp_it) {            
            ModelPart& mp = *smp_it;

            if (r_modelpart.GetProcessInfo()[TIME] < mp[INLET_START_TIME]) continue;
            
            const int mesh_size_elements = smp_it->NumberOfElements();
            
            ModelPart::ElementsContainerType::ContainerType all_elements = smp_it->ElementsArray();
                        
            if (r_modelpart.GetProcessInfo()[TIME] > mp[INLET_STOP_TIME]) {
                if (mLayerRemoved[smp_number]) continue;
                for (int i = 0; i < mesh_size_elements; i++) {                   
                    all_elements[i]->Set(TO_ERASE);
                    all_elements[i]->GetGeometry()[0].Set(TO_ERASE);
                }
                mLayerRemoved[smp_number] = true;
                continue;
            }                        

            int total_mesh_size_accross_mpi_processes = mesh_size_elements; //temporary value until reduction is done
            r_modelpart.GetCommunicator().SumAll(total_mesh_size_accross_mpi_processes);
            const double this_mpi_process_portion_of_inlet_mesh = (double) mesh_size_elements / (double) total_mesh_size_accross_mpi_processes;
            double num_part_surface_time = mp[INLET_NUMBER_OF_PARTICLES];
            num_part_surface_time *= this_mpi_process_portion_of_inlet_mesh;
            const double delta_t = current_time - mLastInjectionTimes[smp_number]; // FLUID DELTA_T CAN BE USED ALSO, it will depend on how often we call this function
            double surface = 1.0; //inlet_surface, this should probably be projected to velocity vector
            
            //calculate number of particles to insert from input data
            const double double_number_of_particles_to_insert = num_part_surface_time * delta_t * surface + mPartialParticleToInsert[smp_number];
            int number_of_particles_to_insert = floor(double_number_of_particles_to_insert);
            mPartialParticleToInsert[smp_number] = double_number_of_particles_to_insert - number_of_particles_to_insert;
            
            if (number_of_particles_to_insert) {
                //randomizing mesh
                srand(/*time(NULL)* */r_modelpart.GetProcessInfo()[TIME_STEPS]);
                
                ModelPart::ElementsContainerType::ContainerType inserting_elements(number_of_particles_to_insert);               
                ModelPart::ElementsContainerType::ContainerType valid_elements = smp_it->ElementsArray();
                int valid_elements_length = 0;
               
                for (int i = 0; i < mesh_size_elements; i++) {                    
                    if (all_elements[i]->IsNot(ACTIVE)) { 
		        valid_elements[valid_elements_length] = all_elements[i];   
		        valid_elements_length++; 
                    } // (push_back) //Inlet BLOCKED nodes are ACTIVE when injecting, but once they are not in contact with other balls, ACTIVE can be reseted.
                }

                if (valid_elements_length < number_of_particles_to_insert) {
                    number_of_particles_to_insert = valid_elements_length;
                    //std::cout<<"The number of DEM particles has been reduced to match the available number of nodes of the DEM Inlet mesh"<<std::endl<<std::flush;
                }
               
                for (int i = 0; i < number_of_particles_to_insert; i++) {
                    int pos = rand() % valid_elements_length;
                    //int pos = i;
                    inserting_elements[i] = valid_elements[pos]; //This only works for pos as real position in the vector if 
                    //we use ModelPart::NodesContainerType::ContainerType instead of ModelPart::NodesContainerType
                    valid_elements[pos] = valid_elements[valid_elements_length - 1];
                    valid_elements_length = valid_elements_length - 1;
                }
               
                PropertiesProxy* p_fast_properties = NULL;
                int general_properties_id = mInletModelPart.GetProperties(mp[PROPERTIES_ID]).Id();  
                for (unsigned int i = 0; i < mFastProperties.size(); i++) {
                    int fast_properties_id = mFastProperties[i].GetId(); 
                    if (fast_properties_id == general_properties_id) {  
                        p_fast_properties = &(mFastProperties[i]);
                        break;
                    }
                }
                
                const array_1d<double, 3> angular_velocity = mp[ANGULAR_VELOCITY];
                const double mod_angular_velocity = MathUtils<double>::Norm3(angular_velocity);
                const double angular_velocity_start_time = mp[ANGULAR_VELOCITY_START_TIME];
                const double angular_velocity_stop_time = mp[ANGULAR_VELOCITY_STOP_TIME];
                const double angular_period = mp[ANGULAR_VELOCITY_PERIOD];
                array_1d<double, 3> angular_velocity_changed;
                const double time = r_modelpart.GetProcessInfo()[TIME];
                array_1d<double, 3> new_axes1;
                array_1d<double, 3> new_axes2;
                array_1d<double, 3> new_axes3;
                
                // The objective of the function call that follows is to rotate the inlet velocity, thus preserving perpendicularity with respect to the inlet plane.
                // Here we compute the three axis new_axes1, new_axes2, new_axes3 where we will have to project the initial inlet velocity to obtain the actual inlet velocity.
                GeometryFunctions::RotateGridOfNodes(time, angular_velocity_start_time, angular_velocity_stop_time, angular_velocity_changed,
                                                     angular_period, mod_angular_velocity, angular_velocity, new_axes1, new_axes2, new_axes3);
                 
                array_1d<double, 3> inlet_velocity = mp[INLET_INITIAL_VELOCITY];
                // Dot product to compute the updated inlet velocity from the initial one:
                mp[LINEAR_VELOCITY] = new_axes1 * inlet_velocity[0] + new_axes2 * inlet_velocity[1] + new_axes3 * inlet_velocity[2];
                
                std::string& ElementNameString = mp[ELEMENT_TYPE];
                const Element& r_reference_element = KratosComponents<Element>::Get(ElementNameString);
                
                Properties::Pointer p_properties = mInletModelPart.pGetProperties(mp[PROPERTIES_ID]);
               
                if (mp[CONTAINS_CLUSTERS] == false) {
               
                    for (int i = 0; i < number_of_particles_to_insert; i++) {
                        Element* p_element = creator.ElementCreatorWithPhysicalParameters(r_modelpart,                                                                     
                                                                                        max_Id+1, 
                                                                                        inserting_elements[i]->GetGeometry()(0), 
                                                                                        inserting_elements[i],
                                                                                        p_properties, 
                                                                                        mp,
                                                                                        r_reference_element, 
                                                                                        p_fast_properties, 
                                                                                        mBallsModelPartHasSphericity, 
                                                                                        mBallsModelPartHasRotation, 
                                                                                        false, 
                                                                                        smp_it->Elements());
                        if(mStrategyForContinuum) {
                            SphericContinuumParticle& spheric_cont_particle = dynamic_cast<SphericContinuumParticle&>(*p_element);
                            spheric_cont_particle.mContinuumInitialNeighborsSize = 0;
                        }
                        inserting_elements[i]->Set(ACTIVE); //Inlet BLOCKED nodes are ACTIVE when injecting, but once they are not in contact with other balls, ACTIVE can be reseted. 
                        inserting_elements[i]->GetGeometry()[0].Set(ACTIVE);
                        max_Id++;
                    }        
                }
               
                else {
               
                    for (int i = 0; i < number_of_particles_to_insert; i++) {
                        int number_of_added_spheres = 0;
                        creator.ClusterCreatorWithPhysicalParameters(r_modelpart, 
                                                                     r_clusters_modelpart,
                                                                     max_Id+1, 
                                                                     inserting_elements[i]->GetGeometry()(0), 
                                                                     inserting_elements[i],
                                                                     p_properties,
                                                                     mp,
                                                                     r_reference_element, 
                                                                     p_fast_properties, 
                                                                     mBallsModelPartHasSphericity, 
                                                                     mBallsModelPartHasRotation, 
                                                                     smp_it->Elements(),
                                                                     number_of_added_spheres);
                        inserting_elements[i]->Set(ACTIVE); //Inlet BLOCKED nodes are ACTIVE when injecting, but once they are not in contact with other balls, ACTIVE can be reseted. 
                        inserting_elements[i]->GetGeometry()[0].Set(ACTIVE);
                        max_Id += number_of_added_spheres;
                    }        
                }               
            } //if (number_of_particles_to_insert)
            mLastInjectionTimes[smp_number] = current_time;
            smp_number++;
        } // for smp_it
        
        creator.RemoveUnusedNodesOfTheClustersModelPart(r_clusters_modelpart);
        
    }    //CreateElementsFromInletMesh       

} // namespace Kratos
