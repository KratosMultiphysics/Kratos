//
// Author: 
// Miguel Angel Celigueta maceli@cimne.upc.edu
//

#include <string>
#include <iostream>

#include "inlet.h"
#include "create_and_destroy.h"
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
        const int number_of_submodelparts = inlet_modelpart.NumberOfSubModelParts();
        mPartialParticleToInsert.resize(number_of_submodelparts, false);
        mLastInjectionTimes.resize(number_of_submodelparts, false);
        //mTotalNumberOfDetachedParticles.resize(number_of_submodelparts, false);
        mLayerRemoved.resize(number_of_submodelparts, false);
        mNumberOfParticlesInjected.resize(number_of_submodelparts, false);
        mMassInjected.resize(number_of_submodelparts, false);
        
        int smp_iterator_number = 0;   
        for (ModelPart::SubModelPartsContainerType::iterator sub_model_part = inlet_modelpart.SubModelPartsBegin(); sub_model_part != inlet_modelpart.SubModelPartsEnd(); ++sub_model_part) {                        
            mPartialParticleToInsert[smp_iterator_number] = 0.0;
            mLastInjectionTimes[smp_iterator_number] = 0.0;
            mLayerRemoved[smp_iterator_number] = false;
            //mTotalNumberOfDetachedParticles[smp_iterator_number] = 0.0;
            mNumberOfParticlesInjected[smp_iterator_number] = 0;
            mMassInjected[smp_iterator_number] = 0.0;
            smp_iterator_number++;
        }                
        
        mFirstInjectionIsDone = false;
        mBallsModelPartHasSphericity = false;  
        mBallsModelPartHasRotation   = false;
        mTotalNumberOfParticlesInjected = 0;
        mTotalMassInjected = 0.0;
        SetNormalizedMaxIndentationForRelease(0.0);
        SetNormalizedMaxIndentationForNewParticleCreation(0.0);
        
        mWarningTooSmallInlet = false;
        mWarningTooSmallInletForMassFlow = false;
    }
    
    void DEM_Inlet::CheckSubModelPart(ModelPart& smp) {
        CheckIfSubModelPartHasVariable(smp, IDENTIFIER);        
        CheckIfSubModelPartHasVariable(smp, VELOCITY);
        CheckIfSubModelPartHasVariable(smp, MAX_RAND_DEVIATION_ANGLE);
        CheckIfSubModelPartHasVariable(smp, PROPERTIES_ID);
        CheckIfSubModelPartHasVariable(smp, INLET_START_TIME);
        CheckIfSubModelPartHasVariable(smp, INLET_STOP_TIME);
        CheckIfSubModelPartHasVariable(smp, ELEMENT_TYPE);
        CheckIfSubModelPartHasVariable(smp, INJECTOR_ELEMENT_TYPE);
        CheckIfSubModelPartHasVariable(smp, INLET_NUMBER_OF_PARTICLES);        
        CheckIfSubModelPartHasVariable(smp, CONTAINS_CLUSTERS);                        
        CheckIfSubModelPartHasVariable(smp, RIGID_BODY_MOTION);
        
        if(smp[RIGID_BODY_MOTION]){
            CheckIfSubModelPartHasVariable(smp, LINEAR_VELOCITY);
            CheckIfSubModelPartHasVariable(smp, ANGULAR_VELOCITY);
            CheckIfSubModelPartHasVariable(smp, ANGULAR_VELOCITY_START_TIME);
            CheckIfSubModelPartHasVariable(smp, ANGULAR_VELOCITY_STOP_TIME);
            CheckIfSubModelPartHasVariable(smp, ANGULAR_VELOCITY_PERIOD);
        }
    }
                 
    void DEM_Inlet::InitializeDEM_Inlet(ModelPart& r_modelpart, ParticleCreatorDestructor& creator, const bool using_strategy_for_continuum) {
        
        mStrategyForContinuum = using_strategy_for_continuum;
        unsigned int& max_Id=creator.mMaxNodeId;
        //CreatePropertiesProxies(mFastProperties, mInletModelPart); 
        mFastProperties = PropertiesProxiesManager().GetPropertiesProxies(r_modelpart);
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
            
            CheckSubModelPart(mp);

            int mesh_size = smp_it->NumberOfNodes();
            if (!mesh_size) continue;
            ModelPart::NodesContainerType::ContainerType& all_nodes = smp_it->NodesArray();
            std::string& identifier = mp[IDENTIFIER];
            mp[INLET_INITIAL_VELOCITY] = mp[LINEAR_VELOCITY];    //This is the velocity of the moving injector of particles
            mp[INLET_INITIAL_PARTICLES_VELOCITY] = mp[VELOCITY]; //This is the initial velocity vector of the injected particles
            
            array_1d<double, 3>& inlet_velocity = mp[VELOCITY];
            
            if ((inlet_velocity[0] == 0.0) &&
                (inlet_velocity[1] == 0.0) &&
                (inlet_velocity[2] == 0.0)) {
                
                KRATOS_THROW_ERROR(std::runtime_error, "The inlet velocity cannot be zero for group ", identifier);
            }
            
            double max_rand_dev_angle = mp[MAX_RAND_DEVIATION_ANGLE];
            if (max_rand_dev_angle < 0.0 || max_rand_dev_angle > 89.5) {
                
                KRATOS_THROW_ERROR(std::runtime_error, "The velocity deviation angle must be between 0 and 89.5 degrees for group ", identifier);
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
            std::string& ElementNameString = mp[INJECTOR_ELEMENT_TYPE];
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

                FixInjectorConditions(p_element);
                max_Id++;
                /*if(mStrategyForContinuum){
                    SphericContinuumParticle* p_continuum_spheric_particle = dynamic_cast<SphericContinuumParticle*>(p_element);
                    p_continuum_spheric_particle->mContinuumInitialNeighborsSize=0;
                    p_continuum_spheric_particle->mInitialNeighborsSize=0;    
                }*/
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
            if (elem_it->Is(DEMFlags::BELONGS_TO_A_CLUSTER)) continue;

            SphericParticle& spheric_particle = dynamic_cast<SphericParticle&>(*elem_it);
            Node<3>& node_it = elem_it->GetGeometry()[0];

            bool have_just_stopped_touching = true;

            for (unsigned int i = 0; i < spheric_particle.mNeighbourElements.size(); i++) {
                SphericParticle* neighbour_particle = spheric_particle.mNeighbourElements[i];
                if(neighbour_particle == NULL) continue;
                
                Node<3>& neighbour_node = neighbour_particle->GetGeometry()[0];

                const double indentation = CalculateNormalizedIndentation(spheric_particle, *neighbour_particle);
                const bool indentation_is_significant_for_release = indentation > mNormalizedMaxIndentationForRelease;
                const bool indentation_is_significant_for_injection = indentation > mNormalizedMaxIndentationForNewParticleCreation;
                const bool i_am_injected_he_is_injector = node_it.IsNot(BLOCKED) && neighbour_node.Is(BLOCKED);
                const bool i_am_injector_he_is_injected = node_it.Is(BLOCKED) && neighbour_node.IsNot(BLOCKED);

                if (i_am_injected_he_is_injector && indentation_is_significant_for_release) {
                    have_just_stopped_touching = false;
                    break;
                }

                if (i_am_injector_he_is_injected && indentation_is_significant_for_injection) {
                    have_just_stopped_touching = false;
                    break;
                }
            }

            if (have_just_stopped_touching) {
                if (node_it.IsNot(BLOCKED)) {//The ball must be freed
                    RemoveInjectionConditions(spheric_particle);
                    UpdateTotalThroughput(spheric_particle);
                }
                else {
                    //Inlet BLOCKED nodes are ACTIVE when injecting, so when they cease to be in contact with other balls, ACTIVE is set to 'false', as they become available for injecting new elements.
                    node_it.Set(ACTIVE, false);
                    elem_it->Set(ACTIVE, false);                    
                }
            }            
        }         
    } //Dettach

    void DEM_Inlet::FixInjectorConditions(Element* p_element){}

    void DEM_Inlet::FixInjectionConditions(Element* p_element)
    {
        Node<3>& node = p_element->GetGeometry()[0];
        node.pGetDof(VELOCITY_X)->FixDof();
        node.pGetDof(VELOCITY_Y)->FixDof();
        node.pGetDof(VELOCITY_Z)->FixDof();
        node.pGetDof(ANGULAR_VELOCITY_X)->FixDof();
        node.pGetDof(ANGULAR_VELOCITY_Y)->FixDof();
        node.pGetDof(ANGULAR_VELOCITY_Z)->FixDof();

        node.Set(DEMFlags::FIXED_VEL_X, true);
        node.Set(DEMFlags::FIXED_VEL_Y, true);
        node.Set(DEMFlags::FIXED_VEL_Z, true);
        node.Set(DEMFlags::FIXED_ANG_VEL_X, true);
        node.Set(DEMFlags::FIXED_ANG_VEL_Y, true);
        node.Set(DEMFlags::FIXED_ANG_VEL_Z, true);
    }

    void DEM_Inlet::RemoveInjectionConditions(Element& element)
    {
        Node<3>& node = element.GetGeometry()[0];
        node.Set(DEMFlags::FIXED_VEL_X, false);
        node.Set(DEMFlags::FIXED_VEL_Y, false);
        node.Set(DEMFlags::FIXED_VEL_Z, false);
        node.Set(DEMFlags::FIXED_ANG_VEL_X, false);
        node.Set(DEMFlags::FIXED_ANG_VEL_Y, false);
        node.Set(DEMFlags::FIXED_ANG_VEL_Z, false);
        element.Set(NEW_ENTITY, 0);
        node.Set(NEW_ENTITY, 0);
        node.pGetDof(VELOCITY_X)->FreeDof();
        node.pGetDof(VELOCITY_Y)->FreeDof();
        node.pGetDof(VELOCITY_Z)->FreeDof();
        node.pGetDof(ANGULAR_VELOCITY_X)->FreeDof();
        node.pGetDof(ANGULAR_VELOCITY_Y)->FreeDof();
        node.pGetDof(ANGULAR_VELOCITY_Z)->FreeDof();
    }

    void DEM_Inlet::DettachClusters(ModelPart& r_clusters_modelpart, unsigned int& max_Id) {    
                            
        vector<unsigned int> ElementPartition;
        typedef ElementsArrayType::iterator ElementIterator;       

        #pragma omp parallel for
        for (int k = 0; k < (int)r_clusters_modelpart.GetCommunicator().LocalMesh().Elements().size(); k++) {
            ElementIterator elem_it = r_clusters_modelpart.GetCommunicator().LocalMesh().Elements().ptr_begin() + k;                            
            if (elem_it->IsNot(NEW_ENTITY)) continue;

            Kratos::Cluster3D& r_cluster = dynamic_cast<Kratos::Cluster3D&>(*elem_it); 
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
                RemoveInjectionConditions(r_cluster);
                UpdateTotalThroughput(r_cluster);

                for (unsigned int j = 0; j < r_cluster.GetSpheres().size(); j++) { //loop over the spheres of the cluster
                    SphericParticle* spheric_particle = r_cluster.GetSpheres()[j];
                    Node<3>& node_it = spheric_particle->GetGeometry()[0];  
                    spheric_particle->Set(NEW_ENTITY, 0);
                    node_it.Set(NEW_ENTITY, 0);                    
                }                    
            }                
        } 
    } //DettachClusters
    
    void DEM_Inlet::CreateElementsFromInletMesh(ModelPart& r_modelpart, ModelPart& r_clusters_modelpart, ParticleCreatorDestructor& creator) {                    
        InitializeStep(r_modelpart);
        unsigned int& max_Id=creator.mMaxNodeId;
        const double current_time = r_modelpart.GetProcessInfo()[TIME];
        DettachElements(r_modelpart, max_Id);
        DettachClusters(r_clusters_modelpart, max_Id);
                
        int smp_number = 0;
        for (ModelPart::SubModelPartsContainerType::iterator smp_it = mInletModelPart.SubModelPartsBegin(); smp_it != mInletModelPart.SubModelPartsEnd(); ++smp_it) {            
            ModelPart& mp = *smp_it;

            const double inlet_start_time = mp[INLET_START_TIME];
            if (current_time < inlet_start_time) continue;
            
            const int mesh_size_elements = smp_it->NumberOfElements();
            
            ModelPart::ElementsContainerType::ContainerType& all_elements = smp_it->ElementsArray();
                        
            if (current_time > mp[INLET_STOP_TIME]) {
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
            double num_part_surface_time = GetInputNumberOfParticles(mp);
            num_part_surface_time *= this_mpi_process_portion_of_inlet_mesh;
            const double delta_t = current_time - mLastInjectionTimes[smp_number]; // FLUID DELTA_T CAN BE USED ALSO, it will depend on how often we call this function
            double surface = 1.0; //inlet_surface, this should probably be projected to velocity vector
            
            int number_of_particles_to_insert = 0;
            const double mass_flow = mp[MASS_FLOW];
            const bool imposed_mass_flow_option = mp.Has(IMPOSED_MASS_FLOW_OPTION) && mp[IMPOSED_MASS_FLOW_OPTION];
            if(imposed_mass_flow_option){                
                number_of_particles_to_insert = mesh_size_elements; // The maximum possible, to increase random.                
                
                if(mass_flow) {
                    const double mean_radius = mp[RADIUS];
                    const double density = mInletModelPart.GetProperties(mp[PROPERTIES_ID])[PARTICLE_DENSITY];
                    const double estimated_mass_of_a_particle = density * 4.0/3.0 * KRATOS_M_PI * mean_radius * mean_radius * mean_radius;                               
                    const double maximum_time_until_release = estimated_mass_of_a_particle * mesh_size_elements / mass_flow;
                    const double minimum_velocity = mean_radius * 3.0 / maximum_time_until_release; //The distance necessary to get out of the injector, over the time.
                    array_1d<double, 3> & proposed_velocity = mp[VELOCITY];
                    const double modulus_of_proposed_velocity = DEM_MODULUS_3(proposed_velocity);
                    const double factor = 2.0;
                    DEM_MULTIPLY_BY_SCALAR_3(proposed_velocity, factor * minimum_velocity / modulus_of_proposed_velocity);
                }
            }
            else {           
                //calculate number of particles to insert from input data
                const double double_number_of_particles_to_insert = num_part_surface_time * delta_t * surface + mPartialParticleToInsert[smp_number];

                if (double_number_of_particles_to_insert < INT_MAX){ // otherwise the precision is not enough to see the residuals
                    number_of_particles_to_insert = std::trunc(double_number_of_particles_to_insert);
                    mPartialParticleToInsert[smp_number] = double_number_of_particles_to_insert - number_of_particles_to_insert;
                }

                else {
                    number_of_particles_to_insert = INT_MAX;
                }
            }
            
            if (number_of_particles_to_insert) {                                
                
                //randomizing mesh
                srand(/*time(NULL)* */r_modelpart.GetProcessInfo()[TIME_STEPS]);
                
                ModelPart::ElementsContainerType::ContainerType valid_elements(mesh_size_elements); //This is a new vector we are going to work on 
                int valid_elements_length = 0;
               
                for (int i = 0; i < mesh_size_elements; i++) {                    
                    if (all_elements[i]->IsNot(ACTIVE)) { 
                        valid_elements[valid_elements_length] = all_elements[i];
                        valid_elements_length++;
                    } // (push_back) //Inlet BLOCKED nodes are ACTIVE when injecting, but once they are not in contact with other balls, ACTIVE can be reseted.
                }
                
                 if (valid_elements_length < number_of_particles_to_insert) {                                                
                    number_of_particles_to_insert = valid_elements_length;
                    if(!imposed_mass_flow_option){
                        ThrowWarningTooSmallInlet(mp);
                    }
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
                array_1d<double, 3> angular_velocity_changed, new_axes1, new_axes2, new_axes3;
                
                // The objective of the function call that follows is to rotate the inlet velocity, thus preserving perpendicularity with respect to the inlet plane.
                // Here we compute the three axis new_axes1, new_axes2, new_axes3 where we will have to project the initial inlet velocity to obtain the actual inlet velocity.
                GeometryFunctions::RotateGridOfNodes(current_time, angular_velocity_start_time, angular_velocity_stop_time, angular_velocity_changed,
                                                     angular_period, mod_angular_velocity, angular_velocity, new_axes1, new_axes2, new_axes3);
                 
                array_1d<double, 3> inlet_initial_velocity = mp[INLET_INITIAL_VELOCITY];
                array_1d<double, 3> inlet_initial_particles_velocity = mp[INLET_INITIAL_PARTICLES_VELOCITY];
                // Dot product to compute the updated inlet velocity from the initial one:
                mp[LINEAR_VELOCITY] = new_axes1 * inlet_initial_velocity[0] + new_axes2 * inlet_initial_velocity[1] + new_axes3 * inlet_initial_velocity[2];
                mp[VELOCITY] = new_axes1 * inlet_initial_particles_velocity[0] + new_axes2 * inlet_initial_particles_velocity[1] + new_axes3 * inlet_initial_particles_velocity[2];
                
                std::string& ElementNameString = mp[ELEMENT_TYPE];
                const Element& r_reference_element = KratosComponents<Element>::Get(ElementNameString);
                
                Properties::Pointer p_properties = mInletModelPart.pGetProperties(mp[PROPERTIES_ID]);
                
                const double mass_that_should_have_been_inserted_so_far = mass_flow * (current_time - inlet_start_time);                               

                int i=0;
                for (i = 0; i < number_of_particles_to_insert; i++) {

                    if (imposed_mass_flow_option) {
                        if(GetPartialMassInjectedSoFar(smp_number) >= mass_that_should_have_been_inserted_so_far ) {
                            break;
                        }
                    }

                    int random_pos = rand() % valid_elements_length;
                    
                    if (mp[CONTAINS_CLUSTERS] == false) {
                        SphericParticle* r_spheric_particle = creator.ElementCreatorWithPhysicalParameters(r_modelpart,
                                                                                    max_Id+1, 
                                                                                    valid_elements[random_pos]->GetGeometry()(0),
                                                                                    valid_elements[random_pos],
                                                                                    //This only works for random_pos as real position in the vector if 
                                                                                    //we use ModelPart::NodesContainerType::ContainerType instead of ModelPart::NodesContainerType
                                                                                    p_properties, 
                                                                                    mp,
                                                                                    r_reference_element, 
                                                                                    p_fast_properties, 
                                                                                    mBallsModelPartHasSphericity, 
                                                                                    mBallsModelPartHasRotation, 
                                                                                    false, 
                                                                                    smp_it->Elements());
                    /*if(mStrategyForContinuum) {
                        SphericContinuumParticle& spheric_cont_particle = dynamic_cast<SphericContinuumParticle&>(*p_element);
                        spheric_cont_particle.mContinuumInitialNeighborsSize = 0;
                    }*/

                        FixInjectionConditions(r_spheric_particle);
                        UpdatePartialThroughput(*r_spheric_particle, smp_number);
                        max_Id++;
                    }
                    else {
                        
                        int number_of_added_spheres = 0;
                        creator.ClusterCreatorWithPhysicalParameters(r_modelpart, 
                                                                    r_clusters_modelpart,
                                                                    max_Id+1, 
                                                                    valid_elements[random_pos]->GetGeometry()(0),
                                                                    valid_elements[random_pos],
                                                                   //This only works for random_pos as real position in the vector if 
                                                                   //we use ModelPart::NodesContainerType::ContainerType instead of ModelPart::NodesContainerType
                                                                    p_properties,
                                                                    mp,
                                                                    r_reference_element, 
                                                                    p_fast_properties, 
                                                                    mBallsModelPartHasSphericity, 
                                                                    mBallsModelPartHasRotation, 
                                                                    smp_it->Elements(),
                                                                    number_of_added_spheres,
                                                                    mStrategyForContinuum);
                                               
                        //UpdatePartialThroughput(*p_cluster, smp_number);
                        max_Id += number_of_added_spheres;                        
                    }

                    valid_elements[random_pos]->Set(ACTIVE); //Inlet BLOCKED nodes are ACTIVE when injecting, but once they are not in contact with other balls, ACTIVE can be reseted. 
                    valid_elements[random_pos]->GetGeometry()[0].Set(ACTIVE);                    
                    valid_elements[random_pos] = valid_elements[valid_elements_length - 1]; //Last position overwrites random_pos
                    valid_elements_length--; //we remove last position and next random_pos has one less option
                }  
                
                if(imposed_mass_flow_option) {
                    if (i == number_of_particles_to_insert && (GetPartialMassInjectedSoFar(smp_number) < mass_that_should_have_been_inserted_so_far ) && mFirstInjectionIsDone == true) {
                        ThrowWarningTooSmallInletForMassFlow(mp);                                                
                    }
                }
                
            } //if (number_of_particles_to_insert)
            mLastInjectionTimes[smp_number] = current_time;
            smp_number++;
        } // for smp_it
        
        creator.RemoveUnusedNodesOfTheClustersModelPart(r_clusters_modelpart);
        mFirstInjectionIsDone = true;
        
    }    //CreateElementsFromInletMesh   
        
    void DEM_Inlet::AddRandomPerpendicularComponentToGivenVector(array_1d<double, 3 >& vector, const double angle_in_degrees)
    {
        KRATOS_TRY
        const double vector_modulus = DEM_MODULUS_3(vector);
        array_1d<double, 3 > unitary_vector;
        noalias(unitary_vector) = vector / vector_modulus;
        array_1d<double, 3 > normal_1;
        array_1d<double, 3 > normal_2;

        if (fabs(unitary_vector[0])>=0.577) {
            normal_1[0]= - unitary_vector[1];
            normal_1[1]= unitary_vector[0];
            normal_1[2]= 0.0;
        }
        else if (fabs(unitary_vector[1])>=0.577) {
            normal_1[0]= 0.0;
            normal_1[1]= - unitary_vector[2];
            normal_1[2]= unitary_vector[1];
        }
        else {
            normal_1[0]= unitary_vector[2];
            normal_1[1]= 0.0;
            normal_1[2]= - unitary_vector[0];
        }

        //normalize(normal_1);
        const double distance0 = DEM_MODULUS_3(normal_1);
        const double inv_distance0 = (distance0 != 0.0) ? 1.0 / distance0 : 0.00;
        normal_1[0] *= inv_distance0;
        normal_1[1] *= inv_distance0;
        normal_1[2] *= inv_distance0;

        //CrossProduct(NormalDirection,Vector0,Vector1);
        DEM_SET_TO_CROSS_OF_FIRST_TWO_3(unitary_vector, normal_1, normal_2)

        const double angle_in_radians = angle_in_degrees * KRATOS_M_PI / 180;
        const double radius = tan(angle_in_radians) * vector_modulus;
        const double radius_square = radius * radius;
        double local_added_vector_modulus_square = radius_square + 1.0; //just greater than the radius, to get at least one iteration of the while
        array_1d<double, 3> local_added_vector; local_added_vector[0] = local_added_vector[1] = local_added_vector[2] = 0.0;

        while (local_added_vector_modulus_square > radius_square) {
            //Random in a range: (max - min) * ( (double)rand() / (double)RAND_MAX ) + min
            local_added_vector[0] = 2*radius * (double)rand() / (double)RAND_MAX - radius;
            local_added_vector[1] = 2*radius * (double)rand() / (double)RAND_MAX - radius;
            local_added_vector_modulus_square = local_added_vector[0]*local_added_vector[0] + local_added_vector[1]*local_added_vector[1];
        }

        noalias(vector) += local_added_vector[0] * normal_1 + local_added_vector[1] * normal_2;
        KRATOS_CATCH("")
    }
    
    void DEM_Inlet::ThrowWarningTooSmallInlet(const ModelPart& mp) {
        if(!mWarningTooSmallInlet) {

            std::cout<<std::endl;
            std::cout<<std::endl;
            std::cout<<"WARNING: At Inlet, the number of injected DEM particles has been reduced to match the available number of nodes for injecting, which was too small. Increase the size of inlet called '"<<mp.Name()<<"'."<<std::endl;
            std::cout<<std::endl;
            std::cout<<std::endl<<std::flush;

            mWarningTooSmallInlet = true;
        }              
    }
    
    void DEM_Inlet::ThrowWarningTooSmallInletForMassFlow(const ModelPart& mp) {
        if(!mWarningTooSmallInletForMassFlow) {
            
            std::cout<<std::endl;
            std::cout<<std::endl;
            std::cout<<"WARNING: At Inlet, the mass flow can not be fulfilled because the number of nodes for injecting was too small. Increase the size of inlet called '"<<mp.Name()<<"'."<<std::endl;
            std::cout<<std::endl;
            std::cout<<std::endl<<std::flush;

            mWarningTooSmallInletForMassFlow = true;
        }           
    }

    ModelPart& DEM_Inlet::GetInletModelPart()
    {
        return mInletModelPart;
    }

    void DEM_Inlet::SetNormalizedMaxIndentationForRelease(const double value)
    {
        mNormalizedMaxIndentationForRelease = value;
    }

    void DEM_Inlet::SetNormalizedMaxIndentationForNewParticleCreation(const double value)
    {
        mNormalizedMaxIndentationForNewParticleCreation = value;
    }
    
    int DEM_Inlet::GetPartialNumberOfParticlesInjectedSoFar(const int i)
    {
        return mNumberOfParticlesInjected[i];
    }

    int DEM_Inlet::GetTotalNumberOfParticlesInjectedSoFar()
    {
        return mTotalNumberOfParticlesInjected;
    }
    
    double DEM_Inlet::GetPartialMassInjectedSoFar(const int i)
    {
        return mMassInjected[i];
    }

    double DEM_Inlet::GetTotalMassInjectedSoFar()
    {
        return mTotalMassInjected;
    }

    void DEM_Inlet::UpdateTotalThroughput(SphericParticle& r_spheric_particle)
    {
        ++mTotalNumberOfParticlesInjected;
        mTotalMassInjected += r_spheric_particle.GetMass();
    }
    
    void DEM_Inlet::UpdateTotalThroughput(Cluster3D& r_cluster)
    {
        ++mTotalNumberOfParticlesInjected;
        mTotalMassInjected += r_cluster.GetMass();
    }
    
    void DEM_Inlet::UpdatePartialThroughput(SphericParticle& r_spheric_particle, const int i)
    {
        ++mNumberOfParticlesInjected[i];
        
        mMassInjected[i] += r_spheric_particle.GetMass();
    }
    
    void DEM_Inlet::UpdatePartialThroughput(Cluster3D& r_cluster, const int i)
    {
        ++mNumberOfParticlesInjected[i];
        
        mMassInjected[i] += r_cluster.GetMass();
    }

    double DEM_Inlet::GetInputNumberOfParticles(const ModelPart& mp)
    {
        double num_part_surface_time = mp[INLET_NUMBER_OF_PARTICLES];

        if (num_part_surface_time >= 0){
           return num_part_surface_time;
        }

        else {
            KRATOS_ERROR << "The value of the Model Part variable INLET_NUMBER_OF_PARTICLES is not a positive int: " << num_part_surface_time;
        }
    }


} // namespace Kratos
