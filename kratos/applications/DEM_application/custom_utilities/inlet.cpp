//
//   Project Name:        Kratos
//   Last Modified by:    $Author: Salva $
//   Date:                $Date: 2014-07-28 16:07:33 $
//   Revision:            $Revision: 1.1.1.1 $
//

// System includes
#include <string>
#include <iostream>

// External includes

// Project includes

#include "DEM_application.h"
#include "inlet.h"

namespace Kratos {
    
    /// Constructor            
    
    DEM_Inlet::DEM_Inlet(ModelPart& inlet_modelpart): mInletModelPart(inlet_modelpart)
    {                
        
        PartialParticleToInsert.resize(inlet_modelpart.NumberOfMeshes(), false);
        mLayerRemoved.resize(inlet_modelpart.NumberOfMeshes(), false);
        
        int mesh_iterator_number = 0;   
        
        for (ModelPart::MeshesContainerType::iterator mesh_it  = inlet_modelpart.GetMeshes().begin();
                                                      mesh_it != inlet_modelpart.GetMeshes().end()  ;    ++mesh_it)
        {
         PartialParticleToInsert[mesh_iterator_number] = 0.0;
         mLayerRemoved[mesh_iterator_number] = false;
         mesh_iterator_number++;
        }                
        
        mFirstTime = true;
        mBallsModelPartHasSphericity = false;  
        mBallsModelPartHasRotation   = false;    
        
    }
            
    /// Destructor
    
    DEM_Inlet::~DEM_Inlet() {}
        
    void DEM_Inlet::InitializeDEM_Inlet(ModelPart& r_modelpart, ParticleCreatorDestructor& creator) {
        
        unsigned int& max_Id=creator.mMaxNodeId;       
        
        CreatePropertiesProxies(mFastProperties, mInletModelPart);      
               
        VariablesList r_modelpart_nodal_variables_list = r_modelpart.GetNodalSolutionStepVariablesList();
        if(r_modelpart_nodal_variables_list.Has(PARTICLE_SPHERICITY) )  mBallsModelPartHasSphericity = true;        
        if(r_modelpart.GetProcessInfo()[ROTATION_OPTION] ) {
            mBallsModelPartHasRotation   = true;
            mInletModelPart.GetProcessInfo()[ROTATION_OPTION] = true;            
        }
        else {
            mInletModelPart.GetProcessInfo()[ROTATION_OPTION] = false;             
        }
               
        int mesh_number=0;
        
        //For the next loop, take into account that the meshes of the modelpart is a full array with no gaps. If mesh i is not defined in mdpa, then it is created empty!
        for (ModelPart::MeshesContainerType::iterator mesh_it = mInletModelPart.GetMeshes().begin()+1; mesh_it != mInletModelPart.GetMeshes().end(); ++mesh_it) {      
            
            mesh_number++;
            int mesh_size=mesh_it->NumberOfNodes();
            if (!mesh_size) continue;
            ModelPart::NodesContainerType::ContainerType all_nodes = mesh_it->NodesArray();
            std::string& identifier = mInletModelPart.GetProperties(mesh_number)[IDENTIFIER];
            
            array_1d<double, 3 >& inlet_velocity = mInletModelPart.GetProperties(mesh_number)[VELOCITY];
            
            if ((inlet_velocity[0] == 0.0) &&
                (inlet_velocity[1] == 0.0) &&
                (inlet_velocity[2] == 0.0)) {
                
                KRATOS_ERROR(std::runtime_error, "The inlet velocity cannot be zero for group ", identifier);
            }
            
            double max_rand_dev_angle = mInletModelPart.GetProperties(mesh_number)[MAX_RAND_DEVIATION_ANGLE];
            if ( max_rand_dev_angle < 0.0 || max_rand_dev_angle > 89.5) {
                
                KRATOS_ERROR(std::runtime_error, "The velocity deviation angle must be between 0 and 90 degrees for group ",identifier);
            }
            
            
            int general_properties_id = mInletModelPart.GetProperties(mesh_number).Id();  
            PropertiesProxy* p_fast_properties = NULL;
            
            for (unsigned int i = 0; i < mFastProperties.size(); i++) {
                int fast_properties_id = mFastProperties[i].GetId(); 
                if (fast_properties_id == general_properties_id) {  
                    p_fast_properties = &(mFastProperties[i]);
                    break;
                }
            }   
            
            Element::Pointer dummy_element_pointer;
            std::string ElementNameString = "SphericParticle3D";
            const Element& r_reference_element = KratosComponents<Element>::Get(ElementNameString);
            
            for (int i = 0; i < mesh_size; i++) {                
                creator.ElementCreatorWithPhysicalParameters(r_modelpart,
                                                             max_Id+1,
                                                             all_nodes[i],
                                                             dummy_element_pointer,
                                                             mInletModelPart.pGetProperties(mesh_number),
                                                             r_reference_element,
                                                             p_fast_properties,
                                                             mBallsModelPartHasSphericity,
                                                             mBallsModelPartHasRotation,
                                                             true,
                                                             mesh_it->Elements()); 
		max_Id++;
            }                 
            
        } //for mesh_it                                               
    } //InitializeDEM_Inlet
        
    void DEM_Inlet::DettachElements(ModelPart& r_modelpart, unsigned int& max_Id) {    
                            
        vector<unsigned int> ElementPartition;
        OpenMPUtils::CreatePartition(OpenMPUtils::GetNumThreads(), r_modelpart.GetCommunicator().LocalMesh().Elements().size(), ElementPartition);
        typedef ElementsArrayType::iterator ElementIterator;       

        #pragma omp parallel for
        for (int k = 0; k < OpenMPUtils::GetNumThreads(); k++) {
            ElementIterator it_begin = r_modelpart.GetCommunicator().LocalMesh().Elements().ptr_begin() + ElementPartition[k];
            ElementIterator it_end   = r_modelpart.GetCommunicator().LocalMesh().Elements().ptr_begin() + ElementPartition[k + 1];

            for (ElementsArrayType::iterator elem_it = it_begin; elem_it != it_end; ++elem_it) {
                            
                if(elem_it->IsNot(NEW_ENTITY)) continue;

                Kratos::SphericParticle& spheric_particle = dynamic_cast<Kratos::SphericParticle&>(*elem_it);

                Node<3>& node_it = elem_it->GetGeometry()[0];                  

                bool still_touching=false;

                for (unsigned int i = 0; i < spheric_particle.mNeighbourElements.size(); i++) {
                    SphericParticle* neighbour_iterator = spheric_particle.mNeighbourElements[i];
                    Node<3>& neighbour_node = neighbour_iterator->GetGeometry()[0]; 
                    if ((node_it.IsNot(BLOCKED) && neighbour_node.Is(BLOCKED)) || (neighbour_node.IsNot(BLOCKED) && node_it.Is(BLOCKED))) {
                        still_touching = true;
                        break;
                    }              
                }

                if (!still_touching) { 
                    if (node_it.IsNot(BLOCKED)) {//The ball must be freed
                        node_it.Set(DEMFlags::FIXED_VEL_X, false);
                        node_it.Set(DEMFlags::FIXED_VEL_Y, false);
                        node_it.Set(DEMFlags::FIXED_VEL_Z, false);
                        node_it.Set(DEMFlags::FIXED_ANG_VEL_X, false);
                        node_it.Set(DEMFlags::FIXED_ANG_VEL_Y, false);
                        node_it.Set(DEMFlags::FIXED_ANG_VEL_Z, false);
                        elem_it->Set(NEW_ENTITY, 0);
                        node_it.Set(NEW_ENTITY, 0);
                    }
                    else {
                    //Inlet BLOCKED nodes are ACTIVE when injecting, but once they are not in contact with other balls, ACTIVE can be reseted.             
                    node_it.Set(ACTIVE, false);
                    elem_it->Set(ACTIVE, false);
                    }
                }
            } //loop elements
        } //loop threads
        
        mFirstTime=false;
        
    } //Dettach
    
    void DEM_Inlet::DettachClusters(ModelPart& r_clusters_modelpart, unsigned int& max_Id) {    
                            
        vector<unsigned int> ElementPartition;
        OpenMPUtils::CreatePartition(OpenMPUtils::GetNumThreads(), r_clusters_modelpart.GetCommunicator().LocalMesh().Elements().size(), ElementPartition);
        typedef ElementsArrayType::iterator ElementIterator;       

        #pragma omp parallel for
        for (int k = 0; k < OpenMPUtils::GetNumThreads(); k++) {
            ElementIterator it_begin = r_clusters_modelpart.GetCommunicator().LocalMesh().Elements().ptr_begin() + ElementPartition[k];
            ElementIterator it_end   = r_clusters_modelpart.GetCommunicator().LocalMesh().Elements().ptr_begin() + ElementPartition[k + 1];

            for (ElementsArrayType::iterator elem_it = it_begin; elem_it != it_end; ++elem_it) {
                            
                if(elem_it->IsNot(NEW_ENTITY)) continue;

                Kratos::Cluster3D& r_cluster = dynamic_cast<Kratos::Cluster3D&>(*elem_it); 
                Node<3>& cluster_central_node = r_cluster.GetGeometry()[0];
                bool still_touching=false;
                
                for (unsigned int j = 0; j < r_cluster.GetSpheres().size(); j++) { //loop over the spheres of the cluster
                    Kratos::SphericParticle* spheric_particle = r_cluster.GetSpheres()[j];
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
                    for (unsigned int j = 0; j < r_cluster.GetSpheres().size(); j++) { //loop over the spheres of the cluster
                        Kratos::SphericParticle* spheric_particle = r_cluster.GetSpheres()[j];
                        Node<3>& node_it = spheric_particle->GetGeometry()[0];  
                        spheric_particle->Set(NEW_ENTITY, 0);
                        node_it.Set(NEW_ENTITY, 0);                    
                    }                    
                }                
                
            } //loop clusters
        } //loop threads
        
        mFirstTime=false;
        
    } //DettachClusters
    
    void DEM_Inlet::CreateElementsFromInletMesh(ModelPart& r_modelpart, ModelPart& r_clusters_modelpart, ParticleCreatorDestructor& creator) {                    
        
        unsigned int& max_Id=creator.mMaxNodeId; 
        
        DettachElements(r_modelpart, max_Id);
        DettachClusters(r_clusters_modelpart, max_Id);
                
        int mesh_number = 0;
        for (ModelPart::MeshesContainerType::iterator mesh_it  = mInletModelPart.GetMeshes().begin()+1;
                                                      mesh_it != mInletModelPart.GetMeshes().end()    ; ++mesh_it)
        {            
            mesh_number++;

            if (r_modelpart.GetProcessInfo()[TIME] < mInletModelPart.GetProperties(mesh_number)[INLET_START_TIME]) continue;
            
            int mesh_size_elements = mesh_it->NumberOfElements();
            
            ModelPart::ElementsContainerType::ContainerType all_elements = mesh_it->ElementsArray();
                        
            if (r_modelpart.GetProcessInfo()[TIME] > mInletModelPart.GetProperties(mesh_number)[INLET_STOP_TIME]) {
                if (mLayerRemoved[mesh_number]) continue;
                for (int i = 0; i < mesh_size_elements; i++) {                   
                     all_elements[i]->Set(TO_ERASE);
                     all_elements[i]->GetGeometry()[0].Set(TO_ERASE);
                }
                mLayerRemoved[mesh_number] = true;
                continue;
            }                        
                        
            double num_part_surface_time = mInletModelPart.GetProperties(mesh_number)[INLET_NUMBER_OF_PARTICLES]; 
            double delta_t = r_modelpart.GetProcessInfo()[DELTA_TIME]; // FLUID DELTA_T CAN BE USED ALSO, it will depend on how often we call this function
            double surface = 1.0;//inlet_surface; // this should probably be projected to velocity vector
            
            //calculate number of particles to insert from input data
            double double_number_of_particles_to_insert = num_part_surface_time * delta_t * surface + PartialParticleToInsert[mesh_number - 1];            
            int number_of_particles_to_insert = floor(double_number_of_particles_to_insert);
            PartialParticleToInsert[mesh_number - 1] = double_number_of_particles_to_insert - number_of_particles_to_insert;
            
            if (number_of_particles_to_insert) {
                //randomizing mesh
                srand( time(NULL)*r_modelpart.GetProcessInfo()[TIME_STEPS] );
                
                ModelPart::ElementsContainerType::ContainerType inserting_elements(number_of_particles_to_insert);               
                ModelPart::ElementsContainerType::ContainerType valid_elements = mesh_it->ElementsArray();
                int valid_elements_length = 0;
               
                for (int i = 0; i < mesh_size_elements; i++) {                    
                    if( all_elements[i]->IsNot(ACTIVE) ) { 
		        valid_elements[valid_elements_length]=all_elements[i];   
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
                int general_properties_id = mInletModelPart.GetProperties(mesh_number).Id();  
                for (unsigned int i = 0; i < mFastProperties.size(); i++) {
                    int fast_properties_id = mFastProperties[i].GetId(); 
                    if (fast_properties_id == general_properties_id) {  
                        p_fast_properties = &(mFastProperties[i]);
                        break;
                    }
                }     
                                  
               std::string& ElementNameString = mInletModelPart.GetProperties(mesh_number)[ELEMENT_TYPE];
               const Element& r_reference_element = KratosComponents<Element>::Get(ElementNameString);
               
               if( mInletModelPart.GetProperties(mesh_number)[CONTAINS_CLUSTERS] == false ) {
               
                    for (int i = 0; i < number_of_particles_to_insert; i++) {
                        creator.ElementCreatorWithPhysicalParameters(r_modelpart,                                                                     
                                                                     max_Id+1, 
                                                                     inserting_elements[i]->GetGeometry()(0), 
                                                                     inserting_elements[i],
                                                                     mInletModelPart.pGetProperties(mesh_number), 
                                                                     r_reference_element, 
                                                                     p_fast_properties, 
                                                                     mBallsModelPartHasSphericity, 
                                                                     mBallsModelPartHasRotation, 
                                                                     false, 
                                                                     mesh_it->Elements());
                        inserting_elements[i]->Set(ACTIVE); //Inlet BLOCKED nodes are ACTIVE when injecting, but once they are not in contact with other balls, ACTIVE can be reseted. 
                        inserting_elements[i]->GetGeometry()[0].Set(ACTIVE);
                        max_Id++;
                    }        
               }
               
               else {
               
                    for (int i = 0; i < number_of_particles_to_insert; i++) {
                        int number_of_added_spheres=0;
                        creator.ClusterCreatorWithPhysicalParameters(r_modelpart, 
                                                                     r_clusters_modelpart,
                                                                     max_Id+1, 
                                                                     inserting_elements[i]->GetGeometry()(0), 
                                                                     inserting_elements[i],
                                                                     mInletModelPart.pGetProperties(mesh_number), 
                                                                     r_reference_element, 
                                                                     p_fast_properties, 
                                                                     mBallsModelPartHasSphericity, 
                                                                     mBallsModelPartHasRotation, 
                                                                     false, 
                                                                     mesh_it->Elements(),
                                                                     number_of_added_spheres);
                        inserting_elements[i]->Set(ACTIVE); //Inlet BLOCKED nodes are ACTIVE when injecting, but once they are not in contact with other balls, ACTIVE can be reseted. 
                        inserting_elements[i]->GetGeometry()[0].Set(ACTIVE);
                        max_Id += number_of_added_spheres;
                    }        
               }
               
           } //if (number_of_particles_to_insert)
            
        } // for mesh_it
        
    }    //CreateElementsFromInletMesh       

} // namespace Kratos.