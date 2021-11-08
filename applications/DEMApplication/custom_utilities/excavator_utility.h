#ifndef EXCAVATOR_UTILITY_H
#define EXCAVATOR_UTILITY_H

#include "GeometryFunctions.h"

namespace Kratos {

    class KRATOS_API(DEM_APPLICATION) ExcavatorUtility {

        public:

        KRATOS_CLASS_POINTER_DEFINITION(ExcavatorUtility);
        
        ExcavatorUtility(ModelPart& rModelPart,
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
                    const double bucket_lifting_velocity_z);
        
        /// Destructor
        virtual ~ExcavatorUtility();

        void ExecuteBeforeSolutionLoop();

        void ExecuteInitializeSolutionStep();

        /// Turn back information as a string
        virtual std::string Info() const;

        /// Print information about this object
        virtual void PrintInfo(std::ostream& rOStream) const;

        /// Print object's data
        virtual void PrintData(std::ostream& rOStream) const;

        protected:

        ModelPart&                                 mrModelPart; //EL SUBMODELPART QUE TOQUE
        array_1d<double,3>                                 mW1;
        array_1d<double,3>                                 mW2;
        double                                   mEccentricity;
        array_1d<double,3>                         mLocalAxis1;
        array_1d<double,3>                         mLocalAxis2;
        array_1d<double,3>                         mLocalAxis3;
        array_1d<double,3>    mInitialCoordinatesOfRotorCenter;
        array_1d<double,3>          mCoordinatesOfStatorCenter;
        double mArmStartTime;
        double mBucketStartTime;
        double mArmStopTime;
        double mBucketStopTime;
        double mTimeLiftBucket;
        double mTimeStopLiftBucket;
        array_1d<double,3> mBucketLiftingVelocity;

        private:

        /// Assignment operator
        ExcavatorUtility & operator=(ExcavatorUtility const& rOther);

        }; // Class ExcavatorUtility
} // namespace Kratos

#endif // EXCAVATOR_UTILITY_H