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
    
    DEM_Inlet::DEM_Inlet(ModelPart& inlet_modelpart): InletModelPart(inlet_modelpart)
    {                
        
        PartialParticleToInsert.resize(inlet_modelpart.NumberOfMeshes(),false);
        mLayerRemoved.resize(inlet_modelpart.NumberOfMeshes(),false);
        
        int mesh_iterator_number=0;   
        
        for (ModelPart::MeshesContainerType::iterator mesh_it = inlet_modelpart.GetMeshes().begin();
                                                      mesh_it != inlet_modelpart.GetMeshes().end();    ++mesh_it)
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
    
    DEM_Inlet::~DEM_Inlet() {};
        
    void DEM_Inlet::InitializeDEM_Inlet(ModelPart& r_modelpart, ParticleCreatorDestructor& creator, const std::string& ElementNameString) {
        
        unsigned int& max_Id=creator.mMaxNodeId;       
        
        CreatePropertiesProxies(mFastProperties, InletModelPart);      
               
        VariablesList r_modelpart_nodal_variables_list = r_modelpart.GetNodalSolutionStepVariablesList();
        if(r_modelpart_nodal_variables_list.Has(PARTICLE_SPHERICITY) )  mBallsModelPartHasSphericity = true;        
        if(r_modelpart.GetProcessInfo()[ROTATION_OPTION] ) {
            mBallsModelPartHasRotation   = true;
            InletModelPart.GetProcessInfo()[ROTATION_OPTION] = true;            
        }
        else {
            InletModelPart.GetProcessInfo()[ROTATION_OPTION] = false;             
        }
        
        const Element& r_reference_element = KratosComponents<Element>::Get(ElementNameString);
        
        int mesh_number=0;
        for (ModelPart::MeshesContainerType::iterator mesh_it = InletModelPart.GetMeshes().begin()+1;
                                               mesh_it != InletModelPart.GetMeshes().end();    ++mesh_it)
        {            
            mesh_number++;
            int mesh_size=mesh_it->NumberOfNodes();
            ModelPart::NodesContainerType::ContainerType all_nodes = mesh_it->NodesArray();
            
            int general_properties_id = InletModelPart.GetProperties(mesh_number).Id();  
            PropertiesProxy* p_fast_properties = NULL;
            
            for (unsigned int i = 0; i < mFastProperties.size(); i++) {
                int fast_properties_id = mFastProperties[i].GetId(); 
                if (fast_properties_id == general_properties_id) {  
                    p_fast_properties = &(mFastProperties[i]);
                    break;
                }
            }   
            
            for (int i = 0; i < mesh_size; i++) {                
                creator.ElementCreatorWithPhysicalParameters(r_modelpart,
                                                             max_Id+1,
                                                             all_nodes[i],
                                                             InletModelPart.pGetProperties(mesh_number),
                                                             r_reference_element,
                                                             p_fast_properties,
                                                             mBallsModelPartHasSphericity,
                                                             mBallsModelPartHasRotation,
                                                             true); 
		max_Id++;
            }      
            
        } //for mesh_it                                               
    } //InitializeDEM_Inlet
        
    void DEM_Inlet::DettachElements(ModelPart& r_modelpart, unsigned int& max_Id) {    
                            
        vector<unsigned int> ElementPartition;
        OpenMPUtils::CreatePartition(OpenMPUtils::GetNumThreads(), r_modelpart.GetCommunicator().LocalMesh().Elements().size(), ElementPartition);
        typedef ElementsArrayType::iterator ElementIterator;       

        //#pragma omp parallel for
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
            } //loop nodes
        } //loop threads
        
        mFirstTime=false;
        
    } //DettachElements
    
    void DEM_Inlet::CreateElementsFromInletMesh(ModelPart& r_modelpart,
                                                ModelPart& inlet_modelpart,
                                                ParticleCreatorDestructor& creator,
                                                const std::string& ElementNameString)
    {                    
        //unsigned int max_Id=0; 
        unsigned int& max_Id=creator.mMaxNodeId; 
        
        DettachElements(r_modelpart, max_Id);                
                
        int mesh_number = 0;
        for (ModelPart::MeshesContainerType::iterator mesh_it = inlet_modelpart.GetMeshes().begin()+1;
                                               mesh_it != inlet_modelpart.GetMeshes().end();    ++mesh_it)
        {            
            mesh_number++;

            if(r_modelpart.GetProcessInfo()[TIME] < InletModelPart.GetProperties(mesh_number)[INLET_START_TIME]) continue;
            
            int mesh_size=mesh_it->NumberOfNodes();
            ModelPart::NodesContainerType::ContainerType all_nodes = mesh_it->NodesArray();
                        
            if (r_modelpart.GetProcessInfo()[TIME] > InletModelPart.GetProperties(mesh_number)[INLET_STOP_TIME]) {
                if (mLayerRemoved[mesh_number]) continue;
                for (int i = 0; i < mesh_size; i++) {                   
		     all_nodes[i]->Set(TO_ERASE);		   
                }
                mLayerRemoved[mesh_number] = true;
                continue;
            }                        
                        
            double num_part_surface_time = InletModelPart.GetProperties(mesh_number)[INLET_NUMBER_OF_PARTICLES]; 
            double delta_t = r_modelpart.GetProcessInfo()[DELTA_TIME]; // FLUID DELTA_T CAN BE USED ALSO, it will depend on how often we call this function
            double surface = 1.0;//inlet_surface; // this should probably be projected to velocity vector
            
            //calculate number of particles to insert from input data
            double double_number_of_particles_to_insert = num_part_surface_time * delta_t * surface + PartialParticleToInsert[mesh_number - 1];            
            int number_of_particles_to_insert = floor(double_number_of_particles_to_insert);
            PartialParticleToInsert[mesh_number - 1] = double_number_of_particles_to_insert - number_of_particles_to_insert;
            
            if (number_of_particles_to_insert) {
                //randomizing mesh
                srand( time(NULL)*r_modelpart.GetProcessInfo()[TIME_STEPS] );
                ModelPart::NodesContainerType::ContainerType inserting_nodes(number_of_particles_to_insert);               
                ModelPart::NodesContainerType::ContainerType valid_nodes = mesh_it->NodesArray();
                int valid_nodes_length = 0;
               
                for (int i = 0; i < mesh_size; i++){
                    if( all_nodes[i]->IsNot(ACTIVE) ) { 
		        valid_nodes[valid_nodes_length]=all_nodes[i];   
		        valid_nodes_length++; 
		    } // (push_back) //Inlet BLOCKED nodes are ACTIVE when injecting, but once they are not in contact with other balls, ACTIVE can be reseted. 
                }

                if (valid_nodes_length < number_of_particles_to_insert) {
                    number_of_particles_to_insert = valid_nodes_length;
                    //std::cout<<"The number of DEM particles has been reduced to match the available number of nodes of the DEM Inlet mesh"<<std::endl<<std::flush;
                }
               
                for (int i = 0; i < number_of_particles_to_insert; i++) {
		    int pos = rand() % valid_nodes_length;
		    //int pos = i;
                    inserting_nodes[i] = valid_nodes[pos]; //This only works for pos as real position in the vector if 
                    //we use ModelPart::NodesContainerType::ContainerType 
                    //instead of ModelPart::NodesContainerType
                    valid_nodes[pos] = valid_nodes[valid_nodes_length - 1];
                    valid_nodes_length = valid_nodes_length - 1;
                }
               
                PropertiesProxy* p_fast_properties = NULL;
                int general_properties_id = InletModelPart.GetProperties(mesh_number).Id();  
                for (unsigned int i = 0; i < mFastProperties.size(); i++) {
                    int fast_properties_id = mFastProperties[i].GetId(); 
                    if (fast_properties_id == general_properties_id) {  
                        p_fast_properties = &(mFastProperties[i]);
                        break;
                    }
                }     
                                             
               const Element& r_reference_element = KratosComponents<Element>::Get(ElementNameString);
               
               for (int i = 0; i < number_of_particles_to_insert; i++) {
                   creator.ElementCreatorWithPhysicalParameters(r_modelpart, max_Id+1, inserting_nodes[i], InletModelPart.pGetProperties(mesh_number), r_reference_element, p_fast_properties, mBallsModelPartHasSphericity, mBallsModelPartHasRotation, false);
                   inserting_nodes[i]->Set(ACTIVE); //Inlet BLOCKED nodes are ACTIVE when injecting, but once they are not in contact with other balls, ACTIVE can be reseted. 
                   max_Id++;
               }                                                                   
               
           } //if (number_of_particles_to_insert)
            
        } // for mesh_it
        
    }    //CreateElementsFromInletMesh       

} // namespace Kratos.