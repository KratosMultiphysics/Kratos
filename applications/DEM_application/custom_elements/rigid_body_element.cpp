// Created by: Salva Latorre, latorre@cimne.upc.edu

// System includes
#include <stdlib.h>

// Project includes
#include "rigid_body_element.h"
#include "custom_utilities/GeometryFunctions.h"
#include "custom_utilities/AuxiliaryFunctions.h"
#include "DEM_application_variables.h"
#include "includes/kratos_flags.h"
#include "includes/variables.h"
#include "geometries/point_3d.h"

namespace Kratos {

    RigidBodyElement3D::RigidBodyElement3D() : Element() {
    }
            
    RigidBodyElement3D::RigidBodyElement3D(IndexType NewId, GeometryType::Pointer pGeometry)
    : Element(NewId, pGeometry) {
        mpIntegrationScheme = NULL;
    }
      
    RigidBodyElement3D::RigidBodyElement3D(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
    : Element(NewId, pGeometry, pProperties) {
        mpIntegrationScheme = NULL;
    }
      
    RigidBodyElement3D::RigidBodyElement3D(IndexType NewId, NodesArrayType const& ThisNodes)
    : Element(NewId, ThisNodes) {
        mpIntegrationScheme = NULL;
    }
    
    Element::Pointer RigidBodyElement3D::Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const {
        return Element::Pointer(new RigidBodyElement3D(NewId, GetGeometry().Create(ThisNodes), pProperties));
    }      
    
    // Destructor
    RigidBodyElement3D::~RigidBodyElement3D() {
    
        mListOfCoordinates.clear();  
        if (mpIntegrationScheme != NULL) delete mpIntegrationScheme;
        // Destroy triangles?
    }
      
    void RigidBodyElement3D::Initialize(ProcessInfo& r_process_info, ModelPart& rigid_body_element_sub_model_part) {
        
        GetGeometry()[0].Set(DEMFlags::FIXED_VEL_X, false);
        GetGeometry()[0].Set(DEMFlags::FIXED_VEL_Y, false);
        GetGeometry()[0].Set(DEMFlags::FIXED_VEL_Z, false);
        GetGeometry()[0].Set(DEMFlags::FIXED_ANG_VEL_X, false);
        GetGeometry()[0].Set(DEMFlags::FIXED_ANG_VEL_Y, false);
        GetGeometry()[0].Set(DEMFlags::FIXED_ANG_VEL_Z, false);
              
        mInertias = ZeroVector(3);
        
        mInertias[0] = rigid_body_element_sub_model_part[RIGID_BODY_INERTIAS][0];
        mInertias[1] = rigid_body_element_sub_model_part[RIGID_BODY_INERTIAS][1];
        mInertias[2] = rigid_body_element_sub_model_part[RIGID_BODY_INERTIAS][2];
        mMass = rigid_body_element_sub_model_part[RIGID_BODY_MASS]; 
        
        GetGeometry()[0].FastGetSolutionStepValue(ORIENTATION) = Quaternion<double>(1.0, 0.0, 0.0, 0.0);
        GetGeometry()[0].FastGetSolutionStepValue(NODAL_MASS) = mMass;    
        GetGeometry()[0].FastGetSolutionStepValue(EXTERNAL_APPLIED_FORCE)[0] = rigid_body_element_sub_model_part[EXTERNAL_APPLIED_FORCE][0];
        GetGeometry()[0].FastGetSolutionStepValue(EXTERNAL_APPLIED_FORCE)[1] = rigid_body_element_sub_model_part[EXTERNAL_APPLIED_FORCE][1];
        GetGeometry()[0].FastGetSolutionStepValue(EXTERNAL_APPLIED_FORCE)[2] = rigid_body_element_sub_model_part[EXTERNAL_APPLIED_FORCE][2];
        GetGeometry()[0].FastGetSolutionStepValue(EXTERNAL_APPLIED_MOMENT)[0] = rigid_body_element_sub_model_part[EXTERNAL_APPLIED_MOMENT][0];
        GetGeometry()[0].FastGetSolutionStepValue(EXTERNAL_APPLIED_MOMENT)[1] = rigid_body_element_sub_model_part[EXTERNAL_APPLIED_MOMENT][1];
        GetGeometry()[0].FastGetSolutionStepValue(EXTERNAL_APPLIED_MOMENT)[2] = rigid_body_element_sub_model_part[EXTERNAL_APPLIED_MOMENT][2];
        
        mFloatingFlag = rigid_body_element_sub_model_part[FLOATING_OPTION];
                
        //DEMIntegrationScheme::Pointer& integration_scheme = GetProperties()[DEM_INTEGRATION_SCHEME_POINTER];
        //SetIntegrationScheme(integration_scheme);
    }   
    
    void RigidBodyElement3D::SetIntegrationScheme(DEMIntegrationScheme::Pointer& integration_scheme){
        mpIntegrationScheme = integration_scheme->CloneRaw();
    }
    
    void RigidBodyElement3D::CustomInitialize() {

        const array_1d<double,3>& reference_inertias = mInertias;                                
        GetGeometry()[0].FastGetSolutionStepValue(NODAL_MASS) = mMass;
        GetGeometry()[0].FastGetSolutionStepValue(PRINCIPAL_MOMENTS_OF_INERTIA)[0] = reference_inertias[0];
        GetGeometry()[0].FastGetSolutionStepValue(PRINCIPAL_MOMENTS_OF_INERTIA)[1] = reference_inertias[1];
        GetGeometry()[0].FastGetSolutionStepValue(PRINCIPAL_MOMENTS_OF_INERTIA)[2] = reference_inertias[2];
        array_1d<double, 3> base_principal_moments_of_inertia = GetGeometry()[0].FastGetSolutionStepValue(PRINCIPAL_MOMENTS_OF_INERTIA);
        
        Quaternion<double>& Orientation = GetGeometry()[0].FastGetSolutionStepValue(ORIENTATION);                
        Orientation.normalize();

        array_1d<double, 3> angular_velocity = GetGeometry()[0].FastGetSolutionStepValue(ANGULAR_VELOCITY);       
        array_1d<double, 3> angular_momentum;
        double LocalTensor[3][3];
        double GlobalTensor[3][3];
        GeometryFunctions::ConstructLocalTensor(base_principal_moments_of_inertia, LocalTensor);
        GeometryFunctions::QuaternionTensorLocal2Global(Orientation, LocalTensor, GlobalTensor);                   
        GeometryFunctions::ProductMatrix3X3Vector3X1(GlobalTensor, angular_velocity, angular_momentum);
        noalias(this->GetGeometry()[0].FastGetSolutionStepValue(ANGULAR_MOMENTUM)) = angular_momentum;

        array_1d<double, 3> local_angular_velocity;
        GeometryFunctions::QuaternionVectorGlobal2Local(Orientation, angular_velocity, local_angular_velocity);
        noalias(this->GetGeometry()[0].FastGetSolutionStepValue(LOCAL_ANGULAR_VELOCITY)) = local_angular_velocity;
    }
    
    void RigidBodyElement3D::SetOrientation(const Quaternion<double> Orientation) {
        this->GetGeometry()[0].FastGetSolutionStepValue(ORIENTATION) = Orientation;
    }

    void RigidBodyElement3D::Calculate(const Variable<double>& rVariable, double& Output, const ProcessInfo& r_process_info) {}

    void RigidBodyElement3D::UpdateLinearDisplacementAndVelocityOfNodes() {
        
        Node<3>& central_node = GetGeometry()[0]; //CENTRAL NODE OF THE RBE
        array_1d<double, 3>& rigid_element_velocity = central_node.FastGetSolutionStepValue(VELOCITY);
        array_1d<double, 3> global_relative_coordinates;
        Quaternion<double>& Orientation = central_node.FastGetSolutionStepValue(ORIENTATION);
        NodesArrayType::iterator i_begin = mListOfNodes.begin();
        NodesArrayType::iterator i_end = mListOfNodes.end();
        int iter = 0;
                
        for (ModelPart::NodeIterator i = i_begin; i != i_end; ++i) {
            GeometryFunctions::QuaternionVectorLocal2Global(Orientation, mListOfCoordinates[iter], global_relative_coordinates);
            array_1d<double, 3>& node_position = i->Coordinates();

            array_1d<double, 3>& delta_displacement = i->FastGetSolutionStepValue(DELTA_DISPLACEMENT);
            array_1d<double, 3> previous_position; 
            noalias(previous_position) = node_position;
            noalias(node_position)= central_node.Coordinates() + global_relative_coordinates;
            noalias(delta_displacement) = node_position - previous_position;
            noalias(i->FastGetSolutionStepValue(VELOCITY)) = rigid_element_velocity;

            iter++;
        }
    }    
    
    void RigidBodyElement3D::UpdatePositionOfNodes() {
        
        Node<3>& central_node = GetGeometry()[0]; //CENTRAL NODE OF THE RBE
        array_1d<double, 3> global_relative_coordinates;      
        array_1d<double, 3> linear_vel_due_to_rotation;
        array_1d<double, 3>& rigid_body_velocity = central_node.FastGetSolutionStepValue(VELOCITY);
        array_1d<double, 3>& rigid_body_angular_velocity = central_node.FastGetSolutionStepValue(ANGULAR_VELOCITY);
        Quaternion<double>& Orientation = central_node.FastGetSolutionStepValue(ORIENTATION);
              
        array_1d<double, 3> previous_position;       
        NodesArrayType::iterator i_begin = mListOfNodes.begin();
        NodesArrayType::iterator i_end = mListOfNodes.end();
        int iter = 0;
                
        for (ModelPart::NodeIterator i = i_begin; i != i_end; ++i) {
            GeometryFunctions::QuaternionVectorLocal2Global(Orientation, mListOfCoordinates[iter], global_relative_coordinates);
            array_1d<double, 3>& node_position = i->Coordinates();
            array_1d<double, 3>& delta_displacement = i->FastGetSolutionStepValue(DELTA_DISPLACEMENT);
            noalias(previous_position) = node_position;
            noalias(node_position)= central_node.Coordinates() + global_relative_coordinates;
            noalias(delta_displacement) = node_position - previous_position;
            GeometryFunctions::CrossProduct(rigid_body_angular_velocity, global_relative_coordinates, linear_vel_due_to_rotation);
            array_1d<double, 3>& velocity = i->FastGetSolutionStepValue(VELOCITY);
            noalias(velocity) = rigid_body_velocity + linear_vel_due_to_rotation;                                    
            iter++;
        }
    }   

    void RigidBodyElement3D::CollectForcesAndTorquesFromTheNodesOfARigidBodyElement() {
        
        KRATOS_TRY
                
        Node<3>& central_node = GetGeometry()[0]; //CENTRAL NODE OF THE RBE
        array_1d<double, 3>& center_forces = central_node.FastGetSolutionStepValue(TOTAL_FORCES);
        array_1d<double, 3>& center_torque = central_node.FastGetSolutionStepValue(PARTICLE_MOMENT);
        center_forces[0] = center_forces[1] = center_forces[2] = center_torque[0] = center_torque[1] = center_torque[2] = 0.0;

        array_1d<double, 3> center_to_node_vector;
        array_1d<double, 3> additional_torque;

        NodesArrayType::iterator i_begin = mListOfNodes.begin();
        NodesArrayType::iterator i_end = mListOfNodes.end();
                
        for (ModelPart::NodeIterator i = i_begin; i != i_end; ++i) {
                        
            array_1d<double, 3>& node_forces = i->FastGetSolutionStepValue(ELASTIC_FORCES); //CONTACT_FORCES, ELASTIC_FORCES???
            center_forces[0] += node_forces[0];
            center_forces[1] += node_forces[1];
            center_forces[2] += node_forces[2];
                        
            array_1d<double, 3>& node_position = i->Coordinates();
            
            center_to_node_vector[0] = node_position[0] - central_node.Coordinates()[0];
            center_to_node_vector[1] = node_position[1] - central_node.Coordinates()[1];
            center_to_node_vector[2] = node_position[2] - central_node.Coordinates()[2];
            GeometryFunctions::CrossProduct(center_to_node_vector, node_forces, additional_torque);
            center_torque[0] += additional_torque[0];
            center_torque[1] += additional_torque[1];
            center_torque[2] += additional_torque[2];
        }
        
        KRATOS_CATCH("")
    }
    
    void RigidBodyElement3D::ComputeBuoyancyEffects() {
    
        KRATOS_TRY
        
        if (!mFloatingFlag) return;
        const double water_density = 1000;
        const double gravity = 9.81;
        const double water_level = 0.0;

        for (unsigned int i = 0; i != mListOfRigidFaces.size(); ++i) {
            double mean_pressure = 0.0;
            double node_Z_coordinate = 0.0;
            double rigid_face_area = 0.0;
            Point rigid_face_centroid;
            array_1d<double, 3> normal_to_rigid_face = ZeroVector(3);
            array_1d<double, 3> normal_rigid_face_force = ZeroVector(3);
            array_1d<double, 3> rigid_body_centroid_to_rigid_face_controid_vector = ZeroVector(3);
            array_1d<double, 3> buoyancy_moment = ZeroVector(3);
            unsigned int rigid_face_size = mListOfRigidFaces[i]->GetGeometry().size();
                    
            for (unsigned int j = 0; j < rigid_face_size; j++) {
                node_Z_coordinate = mListOfRigidFaces[i]->GetGeometry()[j].Coordinates()[2];
                mean_pressure += ((node_Z_coordinate >= water_level) ? 0.0 : -node_Z_coordinate * water_density * gravity);
            }
            
            rigid_face_centroid = mListOfRigidFaces[i]->GetGeometry().Center();
            if (rigid_face_size) mean_pressure /= rigid_face_size;
            else std::cout << "A rigid face with no nodes was found!";
            mListOfRigidFaces[i]->CalculateNormal(normal_to_rigid_face);
            rigid_face_area = mListOfRigidFaces[i]->GetGeometry().Area();
            normal_rigid_face_force[0] = mean_pressure * rigid_face_area * normal_to_rigid_face[0];
            normal_rigid_face_force[1] = mean_pressure * rigid_face_area * normal_to_rigid_face[1];
            normal_rigid_face_force[2] = mean_pressure * rigid_face_area * normal_to_rigid_face[2];
            
            for (unsigned int i = 0; i < rigid_face_size; i++) {
                rigid_body_centroid_to_rigid_face_controid_vector[0] = rigid_face_centroid.Coordinates()[0] - GetGeometry()[0].Coordinates()[0];
                rigid_body_centroid_to_rigid_face_controid_vector[1] = rigid_face_centroid.Coordinates()[1] - GetGeometry()[0].Coordinates()[1];
                rigid_body_centroid_to_rigid_face_controid_vector[2] = rigid_face_centroid.Coordinates()[2] - GetGeometry()[0].Coordinates()[2];
                if (GeometryFunctions::DotProduct(rigid_body_centroid_to_rigid_face_controid_vector, normal_to_rigid_face) > 0.0) {
                    DEM_MULTIPLY_BY_SCALAR_3(normal_rigid_face_force, -1.0)
                }
            }
            
            array_1d<double, 3>& total_forces = GetGeometry()[0].FastGetSolutionStepValue(TOTAL_FORCES);
            DEM_ADD_SECOND_TO_FIRST(total_forces, normal_rigid_face_force)
            array_1d<double, 3>& total_moments = GetGeometry()[0].FastGetSolutionStepValue(PARTICLE_MOMENT);
            GeometryFunctions::CrossProduct(rigid_body_centroid_to_rigid_face_controid_vector, normal_rigid_face_force, buoyancy_moment);
            DEM_ADD_SECOND_TO_FIRST(total_moments, buoyancy_moment)   
        }

        KRATOS_CATCH("")
    }
    
    void RigidBodyElement3D::GetRigidBodyElementsForce(const array_1d<double,3>& gravity) {
        
        KRATOS_TRY
        
        CollectForcesAndTorquesFromTheNodesOfARigidBodyElement();
        ComputeAdditionalForces(gravity);
        AddUpAllForcesAndMoments();
        
        KRATOS_CATCH("")
    }
    
    void RigidBodyElement3D::ComputeAdditionalForces(const array_1d<double,3>& gravity) {

        KRATOS_TRY
        // Gravity
        noalias(GetGeometry()[0].FastGetSolutionStepValue(TOTAL_FORCES)) += mMass * gravity;
        
        // Other forces and moments
        if (mFloatingFlag) { // This is only when having a ship
            ComputeBuoyancyEffects();
            if (0) { // This is only when having a ship
                ComputeEngineForce();
                ComputeWaterDragForce();
            }
        }
        
        KRATOS_CATCH("")
    }
    
    void RigidBodyElement3D::ComputeEngineForce() {
    
        KRATOS_TRY
        
        // We are assuming the ship is moving in the X direction
        const double engine_power = 60000000; // 60MW Arktika-class icebreaker
        const double max_engine_force = 60000000; // with 20MN the ship almost couldn't make it through the ice
        const double threshold_velocity = 1.0; // It was set to 3.0 m/s before, which corresponded to a maximum force of 20MN
        const double engine_performance = 1.0;

        array_1d<double, 3>& external_applied_force = GetGeometry()[0].FastGetSolutionStepValue(EXTERNAL_APPLIED_FORCE);
        const array_1d<double, 3> velocity = GetGeometry()[0].FastGetSolutionStepValue(VELOCITY);

        if ((GetGeometry()[0].FastGetSolutionStepValue(VELOCITY)[0]) < threshold_velocity) external_applied_force[0] = engine_performance * max_engine_force;
        else if (velocity[0]) external_applied_force[0] = engine_performance * engine_power / velocity[0];
                
        KRATOS_CATCH("")
    }
    
    void RigidBodyElement3D::ComputeWaterDragForce() {
        
        KRATOS_TRY
        
        const double drag_constant_X = 500000; // Such that the X maximum velocity is 11 m/s, which corresponds to Arktika-class icebreakers
        const double drag_constant_Y = 240000000; // Such that the Y maximum velocity is 0.5 m/s
        const double drag_constant_Z = 240000000; // Such that the Z maximum velocity is 0.5 m/s
        array_1d<double, 3>& external_applied_force  = GetGeometry()[0].FastGetSolutionStepValue(EXTERNAL_APPLIED_FORCE);
        const array_1d<double, 3> velocity = GetGeometry()[0].FastGetSolutionStepValue(VELOCITY);
        
        // Drag forces due to water. We are assuming the ship is moving in the X direction
        // Quadratic laws were chosen. They may be linear
        external_applied_force[0] += ((velocity[0] >= 0.0) ? -drag_constant_X * velocity[0] * velocity[0] : drag_constant_X * velocity[0] * velocity[0]);
        external_applied_force[1] += ((velocity[1] >= 0.0) ? -drag_constant_Y * velocity[1] * velocity[1] : drag_constant_Y * velocity[1] * velocity[1]);
        external_applied_force[2] += ((velocity[2] >= 0.0) ? -drag_constant_Z * velocity[2] * velocity[2] : drag_constant_Z * velocity[2] * velocity[2]);

        KRATOS_CATCH("")
    }
    
    void RigidBodyElement3D::AddUpAllForcesAndMoments() {
    
        KRATOS_TRY
        
        const array_1d<double, 3> external_applied_force  = GetGeometry()[0].FastGetSolutionStepValue(EXTERNAL_APPLIED_FORCE);
        const array_1d<double, 3> external_applied_torque = GetGeometry()[0].FastGetSolutionStepValue(EXTERNAL_APPLIED_MOMENT);
        noalias(GetGeometry()[0].FastGetSolutionStepValue(TOTAL_FORCES)) += external_applied_force;
        noalias(GetGeometry()[0].FastGetSolutionStepValue(PARTICLE_MOMENT)) += external_applied_torque;
        
        KRATOS_CATCH("")
    }
    
    void RigidBodyElement3D::SetInitialConditionsToNodes(const array_1d<double,3>& velocity) {
        for (unsigned int i=0; i<mListOfCoordinates.size(); i++) {
            this->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY) = velocity;
        } 
    }
    
    void RigidBodyElement3D::Move(const double delta_t, const bool rotation_option, const double force_reduction_factor, const int StepFlag ) {
       
        SymplecticEulerScheme A;
        A.MoveRigidBodyElement(this, GetGeometry()[0], delta_t, rotation_option, force_reduction_factor, StepFlag);
        //GetIntegrationScheme().MoveRigidBodyElement(this, GetGeometry()[0], delta_t, rotation_option, force_reduction_factor, StepFlag);            
    }   
} // namespace Kratos
