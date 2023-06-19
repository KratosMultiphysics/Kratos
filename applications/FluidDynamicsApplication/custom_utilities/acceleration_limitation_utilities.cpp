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
#include "utilities/parallel_utilities.h"
#include "utilities/reduction_utilities.h"

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

        array_1d<double,3> delta_v;
        double max_delta_v = block_for_each<MaxReduction<double>>(mrModelPart.Nodes(), delta_v, [&](Node& rNode, array_1d<double,3>& rDeltaVTLS){
            double norm_delta_v_fixed = 0.0;
            if ( rNode.Is(INLET) || rNode.IsFixed(VELOCITY_X) || rNode.IsFixed(VELOCITY_Y) || rNode.IsFixed(VELOCITY_Z) ){
                const array_1d<double, 3> &r_v  = rNode.FastGetSolutionStepValue( VELOCITY, 0 );
                const array_1d<double, 3> &r_vn = rNode.FastGetSolutionStepValue( VELOCITY, 1 );
                rDeltaVTLS = r_v - r_vn;

                norm_delta_v_fixed = norm_2(rDeltaVTLS);
            }

            return norm_delta_v_fixed;
        });

        // Synchronize maximum fixed acceleration between processes
        max_delta_v = mrModelPart.GetCommunicator().GetDataCommunicator().MaxAll(max_delta_v);

        block_for_each(mrModelPart.Nodes(), [&](Node& rNode){
            // retrieving velocities of current and last time step
            array_1d<double, 3> &v  = rNode.FastGetSolutionStepValue( VELOCITY, 0 );
            const array_1d<double, 3> &vn = rNode.FastGetSolutionStepValue( VELOCITY, 1 );

            array_1d<double, 3> delta_v = v - vn;
            const double norm_delta_v = std::sqrt( delta_v[0]*delta_v[0] + delta_v[1]*delta_v[1] + delta_v[2]*delta_v[2] );

            const double alpha = norm_delta_v / ( dt * mMaximalAccelaration * 9.81 );

            if ( alpha > 1.0 && !rNode.Is(INLET) && norm_delta_v > mMaximalAccelaration*max_delta_v ){
                // setting a new and "reasonable" velocity by scaling
                v = vn + ( 1.0 / alpha ) * delta_v;
                if (!rNode.IsFixed(VELOCITY_X) ){
                    rNode.FastGetSolutionStepValue( VELOCITY_X, 0 ) = v[0];
                }
                if (!rNode.IsFixed(VELOCITY_Y) ){
                    rNode.FastGetSolutionStepValue( VELOCITY_Y, 0 ) = v[1];
                }
                if (!rNode.IsFixed(VELOCITY_Z) ){
                    rNode.FastGetSolutionStepValue( VELOCITY_Z, 0 ) = v[2];
                }
            }
        });
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
