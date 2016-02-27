//
//   Project Name:        Kratos
//   Last Modified by:    $Author: Salva $
//   Date:                $Date: 2006-11-27 16:07:33 $
//   Revision:            $Revision: 1.1.1.1 $
//
//


// System includes
#include <string>
#include <iostream>
#include <stdlib.h>

// External includes

// Project includes
#include "cluster3D.h"
#include "custom_utilities/GeometryFunctions.h"
#include "custom_utilities/AuxiliaryFunctions.h"
#include "DEM_application_variables.h"
#include "includes/kratos_flags.h"
#include "includes/variables.h"

namespace Kratos {
        
    // using namespace GeometryFunctions;

    Cluster3D::Cluster3D() : Element() {}
            
      
    Cluster3D::Cluster3D(IndexType NewId, GeometryType::Pointer pGeometry)
    : Element(NewId, pGeometry) {}
      
      
    Cluster3D::Cluster3D(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
    : Element(NewId, pGeometry, pProperties) {}

      
    Cluster3D::Cluster3D(IndexType NewId, NodesArrayType const& ThisNodes)
    : Element(NewId, ThisNodes) {}

      
    Element::Pointer Cluster3D::Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const {
          
        return Element::Pointer(new Cluster3D(NewId, GetGeometry().Create(ThisNodes), pProperties));
        
    }      

    // Destructor
    Cluster3D::~Cluster3D() {
    
        for (unsigned int i=0; i<mListOfCoordinates.size(); i++) {
            mListOfSphericParticles[i]->Set(DEMFlags::BELONGS_TO_A_CLUSTER, false);
            mListOfSphericParticles[i]->Set(TO_ERASE, true);                        
        }    
        
        mListOfSphericParticles.clear();
        mListOfCoordinates.clear();  
        mListOfRadii.clear();  
    }

      
    void Cluster3D::Initialize() {
        
        if (GetGeometry()[0].GetDof(VELOCITY_X).IsFixed())         GetGeometry()[0].Set(DEMFlags::FIXED_VEL_X,true);
        else                                                          GetGeometry()[0].Set(DEMFlags::FIXED_VEL_X,false);
        if (GetGeometry()[0].GetDof(VELOCITY_Y).IsFixed())         GetGeometry()[0].Set(DEMFlags::FIXED_VEL_Y,true);
        else                                                          GetGeometry()[0].Set(DEMFlags::FIXED_VEL_Y,false);
        if (GetGeometry()[0].GetDof(VELOCITY_Z).IsFixed())         GetGeometry()[0].Set(DEMFlags::FIXED_VEL_Z,true);
        else                                                          GetGeometry()[0].Set(DEMFlags::FIXED_VEL_Z,false);
        if (GetGeometry()[0].GetDof(ANGULAR_VELOCITY_X).IsFixed()) GetGeometry()[0].Set(DEMFlags::FIXED_ANG_VEL_X,true);
        else                                                          GetGeometry()[0].Set(DEMFlags::FIXED_ANG_VEL_X,false);
        if (GetGeometry()[0].GetDof(ANGULAR_VELOCITY_Y).IsFixed()) GetGeometry()[0].Set(DEMFlags::FIXED_ANG_VEL_Y,true);
        else                                                          GetGeometry()[0].Set(DEMFlags::FIXED_ANG_VEL_Y,false);
        if (GetGeometry()[0].GetDof(ANGULAR_VELOCITY_Z).IsFixed()) GetGeometry()[0].Set(DEMFlags::FIXED_ANG_VEL_Z,true);
        else                                                          GetGeometry()[0].Set(DEMFlags::FIXED_ANG_VEL_Z,false);
        
        CustomInitialize();
    }
    
    void Cluster3D::CustomInitialize() {}
    
    void Cluster3D::SetOrientation(const array_1d<double, 3>& euler_angles) {        
        noalias( this->GetGeometry()[0].FastGetSolutionStepValue(EULER_ANGLES) ) = euler_angles;        
    }
      
    void Cluster3D::CreateParticles(ParticleCreatorDestructor* p_creator_destructor, ModelPart& dem_model_part){
        
        KRATOS_TRY 
        
        int cluster_id = (int) this->Id();
        
        unsigned int max_Id = 0;        
        unsigned int* p_max_Id=p_creator_destructor->pGetCurrentMaxNodeId();  //must have been found
          
        std::string ElementNameString = "SphericParticle3D";
            
        //We now create a spheric particle and keep it as a reference to an Element
        const Element& r_reference_element = KratosComponents<Element>::Get(ElementNameString);
        
        Node<3>& central_node = GetGeometry()[0]; //CENTRAL NODE OF THE CLUSTER
          
        const array_1d<double, 3>& euler_angles = central_node.FastGetSolutionStepValue(EULER_ANGLES);
        
        const double mass = central_node.FastGetSolutionStepValue(NODAL_MASS);
        array_1d<double, 3> coordinates_of_sphere;
        array_1d<double, 3> global_relative_coordinates;
        double radius_of_sphere;
        
        double rotation_matrix[3][3];
        GeometryFunctions::GetRotationMatrix(euler_angles, rotation_matrix);
                
        for (unsigned int i=0; i<mListOfCoordinates.size(); i++) {
            
            GeometryFunctions::VectorLocal2Global(rotation_matrix, mListOfCoordinates[i], global_relative_coordinates);
            coordinates_of_sphere[0]= central_node.Coordinates()[0] + global_relative_coordinates[0];
            coordinates_of_sphere[1]= central_node.Coordinates()[1] + global_relative_coordinates[1];
            coordinates_of_sphere[2]= central_node.Coordinates()[2] + global_relative_coordinates[2];    
            
            radius_of_sphere = mListOfRadii[i];   
            #pragma omp critical
            {
            (*p_max_Id)++;
            max_Id = *p_max_Id;
            }
             
            Kratos::SphericParticle* new_sphere = p_creator_destructor->ElementCreatorForClusters(dem_model_part, 
                                                                                                  max_Id, 
                                                                                                  radius_of_sphere, 
                                                                                                  coordinates_of_sphere, 
                                                                                                  mass,
                                                                                                  this->pGetProperties(), 
                                                                                                  r_reference_element,
                                                                                                  cluster_id);
            
            p_creator_destructor->SetMaxNodeId(max_Id);       
            mListOfSphericParticles[i] = new_sphere;
        }
                
        KRATOS_CATCH("")
    }
    
    
    void Cluster3D::UpdateLinearDisplacementAndVelocityOfSpheres(){
        Node<3>& central_node = GetGeometry()[0]; //CENTRAL NODE OF THE CLUSTER
        array_1d<double, 3>& cluster_velocity = central_node.FastGetSolutionStepValue(VELOCITY);
        array_1d<double, 3> global_relative_coordinates;
        array_1d<double, 3 > & EulerAngles = central_node.FastGetSolutionStepValue(EULER_ANGLES);
        double rotation_matrix[3][3];
        
        GeometryFunctions::GetRotationMatrix(EulerAngles, rotation_matrix);
        
        for (unsigned int i=0; i<mListOfCoordinates.size(); i++) {
            
            GeometryFunctions::VectorLocal2Global(rotation_matrix, mListOfCoordinates[i], global_relative_coordinates);
            Node<3>& sphere_node = mListOfSphericParticles[i]->GetGeometry()[0]; 
            array_1d<double, 3>& sphere_position = sphere_node.Coordinates();
            array_1d<double, 3>& delta_displacement = sphere_node.FastGetSolutionStepValue(DELTA_DISPLACEMENT);
            array_1d<double, 3> previous_position; 
            noalias(previous_position) = sphere_position;
            noalias(sphere_position)= central_node.Coordinates() + global_relative_coordinates;
            noalias(delta_displacement) = sphere_position - previous_position;
            noalias(sphere_node.FastGetSolutionStepValue(VELOCITY)) = cluster_velocity;
        }        
    }
    
    void Cluster3D::UpdatePositionOfSpheres(const double RotationMatrix[3][3]) {
        
        Node<3>& central_node = GetGeometry()[0]; //CENTRAL NODE OF THE CLUSTER
        array_1d<double, 3> global_relative_coordinates;      
        array_1d<double, 3> linear_vel_due_to_rotation;
        array_1d<double, 3>& cluster_velocity = central_node.FastGetSolutionStepValue(VELOCITY);
        array_1d<double, 3>& cluster_angular_velocity = central_node.FastGetSolutionStepValue(ANGULAR_VELOCITY);
        array_1d<double, 3>& cluster_delta_rotation = central_node.FastGetSolutionStepValue(DELTA_ROTATION);
        
        array_1d<double, 3> previous_position;        
        
        for (unsigned int i=0; i<mListOfCoordinates.size(); i++) {
            
            GeometryFunctions::VectorLocal2Global(RotationMatrix, mListOfCoordinates[i], global_relative_coordinates);
            Node<3>& sphere_node = mListOfSphericParticles[i]->GetGeometry()[0]; 
            array_1d<double, 3>& sphere_position = sphere_node.Coordinates();
            array_1d<double, 3>& delta_displacement = sphere_node.FastGetSolutionStepValue(DELTA_DISPLACEMENT);
            noalias(previous_position) = sphere_position;
            noalias(sphere_position)= central_node.Coordinates() + global_relative_coordinates;
            noalias(delta_displacement) = sphere_position - previous_position;
            
            GeometryFunctions::CrossProduct( cluster_angular_velocity, global_relative_coordinates, linear_vel_due_to_rotation );
            
            array_1d<double, 3>& velocity = sphere_node.FastGetSolutionStepValue(VELOCITY);
            
            noalias(velocity) = cluster_velocity + linear_vel_due_to_rotation;                                    
            noalias(sphere_node.FastGetSolutionStepValue(ANGULAR_VELOCITY)) = cluster_angular_velocity;
            noalias(sphere_node.FastGetSolutionStepValue(DELTA_ROTATION)) = cluster_delta_rotation;
        }                        
    }
    
    void Cluster3D::CollectForcesAndTorquesFromSpheres() {
        
        Node<3>& central_node = GetGeometry()[0]; //CENTRAL NODE OF THE CLUSTER
        array_1d<double, 3>& center_forces = central_node.FastGetSolutionStepValue(TOTAL_FORCES);        
        array_1d<double, 3>& center_torque = central_node.FastGetSolutionStepValue(PARTICLE_MOMENT);
        center_forces[0] = center_forces[1]= center_forces[2]= center_torque[0]= center_torque[1]= center_torque[2]= 0.0;
        
        array_1d<double, 3> center_to_sphere_vector;
        array_1d<double, 3> additional_torque;
        
        for (unsigned int i=0; i<mListOfSphericParticles.size(); i++) {
            
            if (mListOfSphericParticles[i]->mNeighbourElements.size()==0 && mListOfSphericParticles[i]->mNeighbourRigidFaces.size()==0) continue; //Assuming the sphere only adds contact forces to the cluster
            
            Node<3>& sphere_node = mListOfSphericParticles[i]->GetGeometry()[0]; 
            array_1d<double, 3>& particle_forces = sphere_node.FastGetSolutionStepValue(TOTAL_FORCES);  
            center_forces[0] += particle_forces[0];
            center_forces[1] += particle_forces[1];
            center_forces[2] += particle_forces[2];
            
            array_1d<double, 3>& particle_torque = sphere_node.FastGetSolutionStepValue(PARTICLE_MOMENT); 
            center_torque[0] += particle_torque[0];
            center_torque[1] += particle_torque[1];
            center_torque[2] += particle_torque[2];
                        
            //Now adding the torque due to the eccentric forces (spheres are not on the center of the cluster)
            array_1d<double, 3>& sphere_position = sphere_node.Coordinates();
            center_to_sphere_vector[0] = sphere_position[0] - central_node.Coordinates()[0];
            center_to_sphere_vector[1] = sphere_position[1] - central_node.Coordinates()[1];
            center_to_sphere_vector[2] = sphere_position[2] - central_node.Coordinates()[2];
            GeometryFunctions::CrossProduct( center_to_sphere_vector, particle_forces, additional_torque );
            center_torque[0] += additional_torque[0];
            center_torque[1] += additional_torque[1];
            center_torque[2] += additional_torque[2];
        }    
    }
    
    void Cluster3D::GetClustersForce(const array_1d<double,3>& gravity) {
        
        CollectForcesAndTorquesFromSpheres();
        ComputeAdditionalForces(gravity);
    }
    
    void Cluster3D::ComputeAdditionalForces(const array_1d<double,3>& gravity){
        const double mass = GetGeometry()[0].FastGetSolutionStepValue(NODAL_MASS);
        noalias(GetGeometry()[0].FastGetSolutionStepValue(TOTAL_FORCES)) += mass * gravity;                        
    }    
}  // namespace Kratos.

