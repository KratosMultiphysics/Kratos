// Author: Miguel Angel Celigueta, Salva Latorre

#include "excavator_utility.h"

namespace Kratos {

    ExcavatorUtility::ExcavatorUtility(ModelPart& rModelPart,
                    const double angular_velocity_of_arm_x,
                    const double coordinates_of_arm_articulation_y,
                    const double coordinates_of_arm_articulation_z,
                    const double arm_start_time,
                    const double arm_stop_time,
                    const double angular_velocity_of_bucket_x,
                    const double initial_coordinates_of_bucket_articulation_y,
                    const double initial_coordinates_of_bucket_articulation_z,
                    const double bucket_start_time,
                    const double bucket_stop_time,
                    const double time_to_lift_the_bucket,
                    const double time_to_stop_lifting_the_bucket,
                    const double bucket_lifting_velocity_z): mrModelPart(rModelPart) {
        
        mW1[0] = angular_velocity_of_arm_x;
        mW1[1] = 0.0;
        mW1[2] = 0.0;
        const double X_coord_of_YZ_of_symmetry = -2.252; // Hardcoded at the moment
        mCoordinatesOfStatorCenter[0] = X_coord_of_YZ_of_symmetry;
        mCoordinatesOfStatorCenter[1] = coordinates_of_arm_articulation_y;
        mCoordinatesOfStatorCenter[2] = coordinates_of_arm_articulation_z;
        mInitialCoordinatesOfRotorCenter[0] = X_coord_of_YZ_of_symmetry;
        mInitialCoordinatesOfRotorCenter[1] = initial_coordinates_of_bucket_articulation_y;
        mInitialCoordinatesOfRotorCenter[2] = initial_coordinates_of_bucket_articulation_z;
        const array_1d<double,3> centers_distance = mInitialCoordinatesOfRotorCenter - mCoordinatesOfStatorCenter;
        mEccentricity = sqrt(centers_distance[1] * centers_distance[1] + centers_distance[2] * centers_distance[2]);
        mW2[0] = angular_velocity_of_bucket_x;
        mW2[1] = 0.0;
        mW2[2] = 0.0;
        mArmStartTime = arm_start_time;
        mArmStopTime = arm_stop_time;
        mBucketStartTime = bucket_start_time;
        mBucketStopTime = bucket_stop_time;
        mTimeLiftBucket = time_to_lift_the_bucket;
        mTimeStopLiftBucket = time_to_stop_lifting_the_bucket;
        mBucketLiftingVelocity[0] = 0.0;
        mBucketLiftingVelocity[1] = 0.0;
        mBucketLiftingVelocity[2] = bucket_lifting_velocity_z;
    }

    /// Destructor
    ExcavatorUtility::~ExcavatorUtility() {}

    void ExcavatorUtility::ExecuteBeforeSolutionLoop() {
        
        KRATOS_TRY;
        KRATOS_CATCH("");
    }

    void ExcavatorUtility::ExecuteInitializeSolutionStep() {
        
        KRATOS_TRY;
        ProcessInfo& rCurrentProcessInfo = mrModelPart.GetProcessInfo();
        const double& rCurrentTime = rCurrentProcessInfo[TIME];

        if (!mrModelPart.NodesBegin()->SolutionStepsDataHas(DISPLACEMENT))
            KRATOS_ERROR << "DISPLACEMENT variable is not in move rotor submodelpart.";

        if (!mrModelPart.NodesBegin()->SolutionStepsDataHas(VELOCITY))
            KRATOS_ERROR << "VELOCITY variable is not in move rotor submodelpart.";

        //UPDATE POSITION OF ROTOR AXIS:
        const double initial_angle = atan2((mInitialCoordinatesOfRotorCenter[2] - mCoordinatesOfStatorCenter[2]), (mInitialCoordinatesOfRotorCenter[1] - mCoordinatesOfStatorCenter[1]));
        double rotated_angle1;
        static double final_rotated_angle1;
        if (rCurrentTime < mArmStopTime) {
            rotated_angle1 = mW1[0] * (rCurrentTime - mArmStartTime);
            final_rotated_angle1 = rotated_angle1;
        } else {
            rotated_angle1 = final_rotated_angle1;
            mW1[0] = 0.0;
        }
        
        const double current_angle1 = initial_angle + rotated_angle1;
        array_1d<double,3> vector_to_rotor_center;
        vector_to_rotor_center[0] = 0.0;
        vector_to_rotor_center[1] = mEccentricity * std::cos(current_angle1);
        vector_to_rotor_center[2] = mEccentricity * std::sin(current_angle1);
        const array_1d<double,3> current_rotor_position = mCoordinatesOfStatorCenter + vector_to_rotor_center;
        noalias(mrModelPart[ROTATION_CENTER]) = current_rotor_position;
        //UPDATE VELOCITY OF ROTOR AXIS:
        array_1d<double,3> rotor_velocity;
        MathUtils<double>::CrossProduct(rotor_velocity, mW1, vector_to_rotor_center);

        //UPDATE LOCAL AXES (ROTATE THEM AROUND (0,0,1) THE ROTATED ANGLE )
        double rotated_angle2;
        static double final_rotated_angle2;
        
        if (rCurrentTime < mBucketStartTime) {
            rotated_angle2 = 0.0;
        } else if (rCurrentTime < mBucketStopTime) {
            rotated_angle2 = mW2[0] * (rCurrentTime - mBucketStartTime);
            final_rotated_angle2 = rotated_angle2;
        } else {
            rotated_angle2 = final_rotated_angle2;
            mW2[0] = 0.0;
        }
   
        array_1d<double,3> current_local_axis_1;
        array_1d<double,3> current_local_axis_2;
        array_1d<double,3> current_local_axis_3;

        array_1d<double,3> vertical_vector; vertical_vector[0] = 1.0; vertical_vector[1] = 0.0; vertical_vector[2] = 0.0;
        array_1d<double,3> initial_local_axis_1; initial_local_axis_1[0] = 0.0; initial_local_axis_1[1] = 1.0; initial_local_axis_1[2] = 0.0; //(local axes are assumed oriented as global axes at the beginning)
        array_1d<double,3> initial_local_axis_2; initial_local_axis_2[0] = 0.0; initial_local_axis_2[1] = 0.0; initial_local_axis_2[2] = 1.0; //(local axes are assumed oriented as global axes at the beginning)

        GeometryFunctions::RotateAVectorAGivenAngleAroundAUnitaryVector(initial_local_axis_1, vertical_vector, rotated_angle1 + rotated_angle2, current_local_axis_1);
        GeometryFunctions::RotateAVectorAGivenAngleAroundAUnitaryVector(initial_local_axis_2, vertical_vector, rotated_angle1 + rotated_angle2, current_local_axis_2);

        //UPDATE POSITION AND VELOCITY OF ALL NODES
        array_1d<double,3> from_rotor_center_to_node;
        array_1d<double,3> backup_coordinates;
        
        for (ModelPart::NodesContainerType::iterator node_i = mrModelPart.NodesBegin(); node_i != mrModelPart.NodesEnd(); node_i++) {
            //Get local coordinates at the beginning (local axes are assumed oriented as global axes at the beginning)
            array_1d<double,3> local_coordinates;
            local_coordinates[1] = node_i->Y0() - mInitialCoordinatesOfRotorCenter[1];
            local_coordinates[2] = node_i->Z0() - mInitialCoordinatesOfRotorCenter[2];

            //Use local coordinates with the updated local axes
            noalias(from_rotor_center_to_node) = local_coordinates[1] * current_local_axis_1 + local_coordinates[2] * current_local_axis_2;

            array_1d<double,3>& current_node_position = node_i->Coordinates();
            backup_coordinates = current_node_position;
            current_node_position[1] = current_rotor_position[1] + local_coordinates[1] * current_local_axis_1[1] + local_coordinates[2] * current_local_axis_2[1];
            current_node_position[2] = current_rotor_position[2] + local_coordinates[1] * current_local_axis_1[2] + local_coordinates[2] * current_local_axis_2[2];
            
            if ((rCurrentTime > mTimeLiftBucket) && (rCurrentTime <= mTimeStopLiftBucket)) current_node_position[2] = current_node_position[2] + mBucketLiftingVelocity[2] * (rCurrentTime - mTimeLiftBucket);
            if (rCurrentTime > mTimeStopLiftBucket) current_node_position[2] += mBucketLiftingVelocity[2] * (mTimeStopLiftBucket - mTimeLiftBucket);

            noalias(node_i->FastGetSolutionStepValue(DISPLACEMENT)) = node_i->Coordinates() - node_i->GetInitialPosition();            
            noalias(node_i->FastGetSolutionStepValue(DELTA_DISPLACEMENT)) = node_i->Coordinates() - backup_coordinates;

            array_1d<double,3>& current_node_velocity = node_i->FastGetSolutionStepValue(VELOCITY);
            array_1d<double,3> aux;
            MathUtils<double>::CrossProduct(aux, mW2, from_rotor_center_to_node);
            noalias(current_node_velocity) = rotor_velocity + aux;
            if ((rCurrentTime > mTimeLiftBucket) && (rCurrentTime <= mTimeStopLiftBucket)) current_node_velocity[2] += mBucketLiftingVelocity[2];

        }//end of loop over nodes
        KRATOS_CATCH("");
    }

    /// Turn back information as a string.
    std::string ExcavatorUtility::Info() const {
        
        std::stringstream buffer;
        buffer << "ExcavatorUtility" ;
        return buffer.str();
    }

    /// Print information about this object.
    void ExcavatorUtility::PrintInfo(std::ostream& rOStream) const {rOStream << "ExcavatorUtility";}

    /// Print object's data.
    void ExcavatorUtility::PrintData(std::ostream& rOStream) const {}
} // namespace Kratos
