//
// Author:
//
// Project includes
#include "glued_to_wall_scheme.h"

namespace Kratos {

    GluedToWallScheme::GluedToWallScheme(Condition* p_wall, SphericParticle* p_sphere) {

        mCondition = p_wall;

        array_1d<double, 3> wall_to_sphere_node;
        noalias(wall_to_sphere_node) = p_sphere->GetGeometry()[0].Coordinates() - p_wall->GetGeometry()[0].Coordinates();

        noalias(mInitialNormalToWall) = mCondition->GetGeometry().UnitNormal(p_wall->GetGeometry()[0]);

        mDistanceSignedWithNormal = DEM_INNER_PRODUCT_3(wall_to_sphere_node, mInitialNormalToWall);

        DEM_MULTIPLY_BY_SCALAR_3(mInitialNormalToWall, mDistanceSignedWithNormal);
        array_1d<double, 3> inner_point;
        noalias(inner_point) = p_wall->GetGeometry()[0].Coordinates() + wall_to_sphere_node - mInitialNormalToWall;

        array_1d<double, 3> point_local_coordinates;
        mCondition->GetGeometry().PointLocalCoordinates(point_local_coordinates, inner_point);
        mShapeFunctionsValues.resize(3);
        mCondition->GetGeometry().ShapeFunctionsValues(mShapeFunctionsValues, point_local_coordinates);
    }

    void GluedToWallScheme::SetTranslationalIntegrationSchemeInProperties(Properties::Pointer pProp, bool verbose) const {
        pProp->SetValue(DEM_TRANSLATIONAL_INTEGRATION_SCHEME_POINTER, this->CloneShared());
    }

    void GluedToWallScheme::SetRotationalIntegrationSchemeInProperties(Properties::Pointer pProp, bool verbose) const {
        pProp->SetValue(DEM_ROTATIONAL_INTEGRATION_SCHEME_POINTER, this->CloneShared());
    }

    void GluedToWallScheme::Move(Node<3> & i, const double delta_t, const double force_reduction_factor, const int StepFlag) {
        if (i.Is(DEMFlags::BELONGS_TO_A_CLUSTER)) return;

        array_1d<double, 3> inner_point;
        noalias(inner_point) = mShapeFunctionsValues[0] * mCondition->GetGeometry()[0] + mShapeFunctionsValues[1] * mCondition->GetGeometry()[1] + mShapeFunctionsValues[2] * mCondition->GetGeometry()[2];

        noalias(mCurrentNormalToWall) = mCondition->GetGeometry().UnitNormal(mCondition->GetGeometry()[0]);
        DEM_MULTIPLY_BY_SCALAR_3(mCurrentNormalToWall, mDistanceSignedWithNormal);

        array_1d<double, 3> & coordinates = i.Coordinates();
        array_1d<double, 3> old_coordinates;
        noalias(old_coordinates) = coordinates;
        noalias(coordinates) = inner_point + mCurrentNormalToWall;
        array_1d<double, 3> & delta_displacement = i.FastGetSolutionStepValue(DELTA_DISPLACEMENT);
        noalias(delta_displacement) = coordinates - old_coordinates;
        noalias(i.FastGetSolutionStepValue(DISPLACEMENT)) += delta_displacement;

        //Finding velocity of inner_point:
        array_1d<double, 3> velocity_of_inner_point = ZeroVector(3);
        noalias(velocity_of_inner_point) += mShapeFunctionsValues[0] * mCondition->GetGeometry()[0].FastGetSolutionStepValue(VELOCITY);
        noalias(velocity_of_inner_point) += mShapeFunctionsValues[1] * mCondition->GetGeometry()[1].FastGetSolutionStepValue(VELOCITY);
        noalias(velocity_of_inner_point) += mShapeFunctionsValues[2] * mCondition->GetGeometry()[2].FastGetSolutionStepValue(VELOCITY);

        //Finding angular velocity of inner_point:
        BoundedMatrix<double,9,3>  matrix_a;
        BoundedVector<double,9> vector_b;
        std::vector<array_1d<double, 3> > r;
        r.resize(3);

        for (int i=0; i<3; i++) {
            noalias(r[i]) = mCondition->GetGeometry()[i] - inner_point;
            matrix_a(3*i,0)   =  0.0;
            matrix_a(3*i,1)   =  r[i][2];
            matrix_a(3*i,2)   = -r[i][1];
            matrix_a(3*i+1,0) = -r[i][2];
            matrix_a(3*i+1,1) =  0.0;
            matrix_a(3*i+1,2) =  r[i][0];
            matrix_a(3*i+2,0) =  r[i][1];
            matrix_a(3*i+2,1) = -r[i][0];
            matrix_a(3*i+2,2) =  0.0;
	    array_1d<double, 3>& velocity_of_node_i = mCondition->GetGeometry()[i].FastGetSolutionStepValue(VELOCITY);
            vector_b[3*i]   = velocity_of_node_i[0] - velocity_of_inner_point[0];
            vector_b[3*i+1] = velocity_of_node_i[1] - velocity_of_inner_point[1];
            vector_b[3*i+2] = velocity_of_node_i[2] - velocity_of_inner_point[2];
        }

        BoundedMatrix<double,3,9>  trans_matrix_a = trans(matrix_a);

        BoundedMatrix<double,3,3> new_lhs = prod(trans_matrix_a, matrix_a);
        BoundedVector<double,3> new_rhs = prod(trans_matrix_a, vector_b);
        double det = 0.0;
        Matrix inverse_new_lhs(3,3);
        MathUtils<double>::InvertMatrix3(new_lhs, inverse_new_lhs, det);
        array_1d<double, 3>&  angular_velocity = i.FastGetSolutionStepValue(ANGULAR_VELOCITY);
        noalias(angular_velocity) = prod(inverse_new_lhs, new_rhs);

        //Once the angular velocity is found, we use it to calculate the total velocity of the sphere:
        array_1d<double, 3> linear_vel_of_sphere_due_to_rotation;
        MathUtils<double>::CrossProduct(linear_vel_of_sphere_due_to_rotation, angular_velocity, mCurrentNormalToWall);
        noalias(i.FastGetSolutionStepValue(VELOCITY)) = velocity_of_inner_point + linear_vel_of_sphere_due_to_rotation;
    }

    void GluedToWallScheme::Rotate(Node<3> & i, const double delta_t, const double moment_reduction_factor, const int StepFlag) {
        if (i.Is(DEMFlags::BELONGS_TO_A_CLUSTER)) return;
        //NOTE: ANGULAR_VELOCITY was calculated in the 'Move' function.
        array_1d<double, 3> rotation_vector;
        MathUtils<double>::CrossProduct(rotation_vector, mInitialNormalToWall, mCurrentNormalToWall);
        const double norm = MathUtils<double>::Norm3(rotation_vector);
        const double angle_between_normals = asin( norm / ( MathUtils<double>::Norm3(mInitialNormalToWall) * MathUtils<double>::Norm3(mCurrentNormalToWall) ) );
        DEM_MULTIPLY_BY_SCALAR_3(rotation_vector, angle_between_normals);
        array_1d<double, 3>& particle_rotation_angle =  i.FastGetSolutionStepValue(PARTICLE_ROTATION_ANGLE);
        noalias(i.FastGetSolutionStepValue(DELTA_ROTATION)) = rotation_vector - particle_rotation_angle;
        noalias(particle_rotation_angle) = rotation_vector;
    }

} //namespace Kratos
