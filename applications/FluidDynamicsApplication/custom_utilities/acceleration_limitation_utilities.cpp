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

        double dt = mrModelPart.GetProcessInfo().GetValue(DELTA_TIME);

        ModelPart::NodesContainerType rNodes = mrModelPart.Nodes();
        #pragma omp parallel for
        for(int count = 0; count < static_cast<int>(rNodes.size()); count++)
        {
            ModelPart::NodesContainerType::iterator i = rNodes.begin() + count;

            // retrieving velocities of current and last time step
            array_1d<double, 3> &v  = i->FastGetSolutionStepValue( VELOCITY, 0 );
            const array_1d<double, 3> vn = i->FastGetSolutionStepValue( VELOCITY, 1 );

            array_1d<double, 3> delta_v = v - vn;
            double norm_delta_v = sqrt( delta_v[0]*delta_v[0] + delta_v[1]*delta_v[1] + delta_v[2]*delta_v[2] );

            const double alpha = norm_delta_v / ( dt * mMaximalAccelaration * 9.81 );

            if ( alpha > 1.0){
                // setting a new and "reasonable" velocity by scaling
                v = vn + ( 1.0 / alpha ) * delta_v;
                i->FastGetSolutionStepValue( VELOCITY, 0 ) = v;
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
