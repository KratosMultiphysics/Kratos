//    |  /           |
//    ' /   __| _` | __|  _ \   __|
//    . \  |   (   | |   (   |\__ `
//   _|\_\_|  \__,_|\__|\___/ ____/
//                   Multi-Physics
//
//  License:		 BSD License
//					 Kratos default license: kratos/license.txt
//
//  Main authors:    Ruben Zorrilla
//
//

// System includes

// External includes

// Project includes

// Application includes
#include "acceleration_limitation_utilities.h"

namespace Kratos
{

    /* Public functions *******************************************************/

    /**
     * @brief Iteration over all nodes in ModelPart to find and reduce excessive accelerations.
     * Model part is specified in the constructor
     */

    void AccelerationLimitationUtilities::Execute() {

        const double dt = mrModelPart.GetProcessInfo().GetValue(DELTA_TIME);

        ModelPart::NodesContainerType rNodes = mrModelPart.Nodes();

        /* double min_x_delta_v = 0.0;
        double min_y_delta_v = 0.0;
        double min_z_delta_v = 0.0;
        double max_x_delta_v = 0.0;
        double max_y_delta_v = 0.0;
        double max_z_delta_v = 0.0; */

        double max_delta_v = 0.0;

        #pragma omp parallel for reduction(max:max_delta_v) //reduction(max:max_x_delta_v) reduction(min:min_x_delta_v) reduction(max:max_y_delta_v) reduction(min:min_y_delta_v) reduction(max:max_z_delta_v) reduction(min:min_z_delta_v)
        for(int count = 0; count < static_cast<int>(rNodes.size()); count++)
        {
            ModelPart::NodesContainerType::iterator i = rNodes.begin() + count;

            if (i->Is(INLET)){
                const array_1d<double, 3> &v  = i->FastGetSolutionStepValue( VELOCITY, 0 );
                const array_1d<double, 3> &vn = i->FastGetSolutionStepValue( VELOCITY, 1 );
                const array_1d<double, 3> delta_v = v - vn;

                const double norm_delta_v = std::sqrt( delta_v[0]*delta_v[0] + delta_v[1]*delta_v[1] + delta_v[2]*delta_v[2] );

                max_delta_v = max_delta_v > norm_delta_v ? max_delta_v : norm_delta_v;

                //if (i->IsFixed(VELOCITY_X)){
                    //max_x_delta_v = max_x_delta_v > delta_v[0] ? max_x_delta_v : delta_v[0];
                    //min_x_delta_v = min_x_delta_v < delta_v[0] ? min_x_delta_v : delta_v[0];
                //}

                //if (i->IsFixed(VELOCITY_Y)){
                    //max_y_delta_v = max_y_delta_v > delta_v[1] ? max_y_delta_v : delta_v[1];
                    //min_y_delta_v = min_y_delta_v < delta_v[1] ? min_y_delta_v : delta_v[1];
                //}

                //if (i->IsFixed(VELOCITY_Z)){
                    //max_z_delta_v = max_z_delta_v > delta_v[2] ? max_z_delta_v : delta_v[2];
                    //min_z_delta_v = min_z_delta_v < delta_v[2] ? min_z_delta_v : delta_v[2];
                //}
            }
        }

        /* KRATOS_INFO("max_x_delta_v") << max_x_delta_v;
        KRATOS_INFO("min_x_delta_v") << min_x_delta_v;

        max_x_delta_v *= mMaximalAccelaration;
        min_x_delta_v *= mMaximalAccelaration;
        max_y_delta_v *= mMaximalAccelaration;
        min_y_delta_v *= mMaximalAccelaration;
        max_z_delta_v *= mMaximalAccelaration;
        min_z_delta_v *= mMaximalAccelaration; */

        KRATOS_INFO("max_delta_v") << max_delta_v << std::endl;

        #pragma omp parallel for
        for(int count = 0; count < static_cast<int>(rNodes.size()); count++)
        {
            ModelPart::NodesContainerType::iterator i = rNodes.begin() + count;

            // retrieving velocities of current and last time step
            array_1d<double, 3> &v  = i->FastGetSolutionStepValue( VELOCITY, 0 );
            const array_1d<double, 3> &vn = i->FastGetSolutionStepValue( VELOCITY, 1 );

            array_1d<double, 3> delta_v = v - vn;
            const double norm_delta_v = std::sqrt( delta_v[0]*delta_v[0] + delta_v[1]*delta_v[1] + delta_v[2]*delta_v[2] );

            const double alpha = norm_delta_v / ( dt * mMaximalAccelaration * 9.81 );

            if ( alpha > 1.0 && !i->Is(INLET) && norm_delta_v > mMaximalAccelaration*max_delta_v){//!i->Is(INLET) ){
                // setting a new and "reasonable" velocity by scaling
                v = vn + ( 1.0 / alpha ) * delta_v;
                /*if ( !i->IsFixed(VELOCITY_X) &&  (delta_v[0] > max_x_delta_v || delta_v[0] < min_x_delta_v) ){
                    i->FastGetSolutionStepValue( VELOCITY_X, 0 ) = v[0];
                }
                if ( !i->IsFixed(VELOCITY_Y) && (delta_v[1] > max_y_delta_v || delta_v[1] < min_y_delta_v) ){
                    i->FastGetSolutionStepValue( VELOCITY_Y, 0 ) = v[1];
                }
                if ( !i->IsFixed(VELOCITY_Z) && (delta_v[2] > max_z_delta_v || delta_v[2] < min_z_delta_v) ){
                    i->FastGetSolutionStepValue( VELOCITY_Z, 0 ) = v[2];
                }*/
            }
        }

    }


    /**
     * @brief Set a new value for the maximal acceleration
     * The value is given as a multiple of the gravitational acceleration (1g = 9.81 m/s²)
     *
     * @param newMaxAcc new maximal acceleration
     */

    void AccelerationLimitationUtilities::SetLimitAsMultipleOfGravitionalAcceleration( double& newMaxAcc ){

        mMaximalAccelaration = newMaxAcc;
    }


    /* External functions *****************************************************/

    /// output stream function
    inline std::ostream& operator << (
        std::ostream& rOStream,
        const AccelerationLimitationUtilities& rThis ) {

        rThis.PrintData(rOStream);
        return rOStream;
    }

}
