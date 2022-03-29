// System includes
#include <string>
#include <iostream>
#include <stdlib.h>

// Project includes
#include "inlet_element.h"
#include "custom_utilities/GeometryFunctions.h"
#include "DEM_application_variables.h"
#include "includes/variables.h"

namespace Kratos {

    InletElement3D::InletElement3D() : RigidBodyElement3D() {}

    InletElement3D::InletElement3D(IndexType NewId, GeometryType::Pointer pGeometry)
    : RigidBodyElement3D(NewId, pGeometry) {}

    InletElement3D::InletElement3D(IndexType NewId, GeometryType::Pointer pGeometry, PropertiesType::Pointer pProperties)
    : RigidBodyElement3D(NewId, pGeometry, pProperties) {}

    InletElement3D::InletElement3D(IndexType NewId, NodesArrayType const& ThisNodes)
    : RigidBodyElement3D(NewId, ThisNodes) {}

    Element::Pointer InletElement3D::Create(IndexType NewId, NodesArrayType const& ThisNodes, PropertiesType::Pointer pProperties) const {
        return Element::Pointer(new InletElement3D(NewId, GetGeometry().Create(ThisNodes), pProperties));
    }

    InletElement3D::~InletElement3D() {}




    void InletElement3D::Initialize(const ProcessInfo& r_process_info) {
        /*
        This will be called once from the solver_Strategy -> rigid_body_element->Initialize
        */

        auto& central_node = GetGeometry()[0];

        central_node.Set(DEMFlags::FIXED_VEL_X, true);
        central_node.Set(DEMFlags::FIXED_VEL_Y, true);
        central_node.Set(DEMFlags::FIXED_VEL_Z, true);
        central_node.Set(DEMFlags::FIXED_ANG_VEL_X, true);
        central_node.Set(DEMFlags::FIXED_ANG_VEL_Y, true);
        central_node.Set(DEMFlags::FIXED_ANG_VEL_Z, true);

        DEMIntegrationScheme::Pointer& translational_integration_scheme = GetProperties()[DEM_TRANSLATIONAL_INTEGRATION_SCHEME_POINTER];
        DEMIntegrationScheme::Pointer& rotational_integration_scheme = GetProperties()[DEM_ROTATIONAL_INTEGRATION_SCHEME_POINTER];
        SetIntegrationScheme(translational_integration_scheme, rotational_integration_scheme);
    }

    //*************************************************
    /// from inlet.cpp
    inline double CalculateNormalizedIndentation(SphericParticle& elem_it_1, SphericParticle& elem_it_2) {
        const array_1d<double,3>& coordinates_1 = elem_it_1.GetGeometry()[0].Coordinates();
        const array_1d<double,3>& coordinates_2 = elem_it_2.GetGeometry()[0].Coordinates();

        const double distance = std::sqrt((coordinates_1[0]- coordinates_2[0]) * (coordinates_1[0] - coordinates_2[0]) +
                                          (coordinates_1[1]- coordinates_2[1]) * (coordinates_1[1] - coordinates_2[1]) +
                                          (coordinates_1[2]- coordinates_2[2]) * (coordinates_1[2] - coordinates_2[2]));

        const double radius_sum = elem_it_1.GetInteractionRadius() + elem_it_2.GetInteractionRadius();
        double indentation = radius_sum - distance;

        indentation /= radius_sum;

        return indentation;
    }

    bool SortSubModelPartsByName(ModelPart* A, ModelPart* B) {
        return (A->Name() < B->Name());
    }






    //*************************************************
    ///reference - maybe used
    void InletElement3D::UpdateLinearDisplacementAndVelocityOfNodes() {
        // and update the vector of the injected particles in the same function?

        Node<3>& central_node = GetGeometry()[0];
        array_1d<double, 3>& rigid_body_velocity = central_node.FastGetSolutionStepValue(VELOCITY);
        array_1d<double, 3> global_relative_coordinates;
        Quaternion<double>& Orientation = central_node.FastGetSolutionStepValue(ORIENTATION);

        NodesArrayType::iterator i_begin = mListOfNodes.begin();
        NodesArrayType::iterator i_end = mListOfNodes.end();
        int iter = 0;

        for (ModelPart::NodeIterator i = i_begin; i != i_end; ++i) {

            GeometryFunctions::QuaternionVectorLocal2Global(Orientation, mListOfCoordinates[iter], global_relative_coordinates);
            array_1d<double, 3>& node_position = i->Coordinates();
            array_1d<double, 3>& delta_displacement = i->FastGetSolutionStepValue(DELTA_DISPLACEMENT);
            array_1d<double, 3>& displacement = i->FastGetSolutionStepValue(DISPLACEMENT);
            array_1d<double, 3> previous_position;
            noalias(previous_position) = node_position;
            noalias(node_position)= central_node.Coordinates() + global_relative_coordinates;
            noalias(delta_displacement) = node_position - previous_position;
            noalias(displacement) += delta_displacement;
            noalias(i->FastGetSolutionStepValue(VELOCITY)) = rigid_body_velocity;
            iter++;
        }
    }

    void InletElement3D::UpdateAngularDisplacementAndVelocityOfNodes() {
        // and update the vector of the injected particles in the same function?

        Node<3>& central_node = GetGeometry()[0];
        array_1d<double, 3> global_relative_coordinates;
        array_1d<double, 3> linear_vel_due_to_rotation;
        array_1d<double, 3>& rigid_body_velocity = central_node.FastGetSolutionStepValue(VELOCITY);
        array_1d<double, 3>& rigid_body_angular_velocity = central_node.FastGetSolutionStepValue(ANGULAR_VELOCITY);
        array_1d<double, 3>& rigid_body_delta_rotation = central_node.FastGetSolutionStepValue(DELTA_ROTATION);
        Quaternion<double>& Orientation = central_node.FastGetSolutionStepValue(ORIENTATION);

        array_1d<double, 3> previous_position;

        NodesArrayType::iterator i_begin = mListOfNodes.begin();
        NodesArrayType::iterator i_end = mListOfNodes.end();
        int iter = 0;

        for (ModelPart::NodeIterator i = i_begin; i != i_end; ++i) {

            GeometryFunctions::QuaternionVectorLocal2Global(Orientation, mListOfCoordinates[iter], global_relative_coordinates);
            GeometryFunctions::CrossProduct( rigid_body_angular_velocity, global_relative_coordinates, linear_vel_due_to_rotation );
            array_1d<double, 3>& velocity = i->FastGetSolutionStepValue(VELOCITY);
            noalias(velocity) = rigid_body_velocity + linear_vel_due_to_rotation;
            noalias(i->FastGetSolutionStepValue(ANGULAR_VELOCITY)) = rigid_body_angular_velocity;
            noalias(i->FastGetSolutionStepValue(DELTA_ROTATION)) = rigid_body_delta_rotation;
            iter++;
        }
    }

} // namespace Kratos
